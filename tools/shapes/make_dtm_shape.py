"""Convert a global radius DTM (simple cylindrical, planetocentric, east longitude
from -180 deg at the left edge, 32-bit float radii in meters) into a body mesh for
assets/shapes/, e.g. the DLR Dawn HAMO DTM of Vesta (USGS Astrogeology).

Sources and conventions are listed in assets/shapes/SOURCES.md. Download the
DTM there first, then (from the repository root):

    python tools/shapes/make_dtm_shape.py <dtm.tif> <output.mesh> [grid_deg]

e.g. Vesta_Dawn_HAMO_DTM_DLR_Global_48ppd.tif -> assets/shapes/vesta.mesh (1.5 deg grid)

The DTM is averaged over cells of grid_deg (pixels below 100 km count as no data);
each grid vertex takes the mean of the four cells around it, a pole the mean of
its whole cell row. Vertices without data are filled by diffusion over the grid
from their neighbours. The mesh keeps the grid (see axmesh.grid_mesh), so it
carries map coordinates for a texture in the same planetocentric system. The
albedo is uniform.

Output format: see axmesh.py. Requires Pillow (only to read the TIFF layout).
"""

import math
import sys
from array import array

from PIL import Image

from axmesh import grid_mesh, volume, write_mesh

Image.MAX_IMAGE_PIXELS = None

MIN_RADIUS_M = 1.0e5
FILL_ITERATIONS = 2000


def cell_means(path, cell_px):
    """Mean radius (km) of the valid pixels in each cell, or None; rows from the north."""
    with Image.open(path) as img:
        w, h = img.size
        tiles = img.tile
        if img.mode != 'F' or any(t.codec_name != 'raw' or t.extents[2] - t.extents[0] != w or
                                  t.extents[3] - t.extents[1] != 1 for t in tiles):
            sys.exit('expected an uncompressed float32 TIFF stored one row per strip')
        offsets = {t.extents[1]: t.offset for t in tiles}
    if w % cell_px or h % cell_px:
        sys.exit(f'{w} x {h} pixels are not whole {cell_px}-pixel cells')
    cols, rows = w // cell_px, h // cell_px
    sums = [[0.0] * cols for _ in range(rows)]
    counts = [[0] * cols for _ in range(rows)]
    with open(path, 'rb') as f:
        for y in range(h):
            f.seek(offsets[y])
            row = array('f')
            row.frombytes(f.read(4 * w))
            if sys.byteorder != 'little':
                row.byteswap()
            s, n = sums[y // cell_px], counts[y // cell_px]
            for c in range(cols):
                valid = [v for v in row[c * cell_px:(c + 1) * cell_px] if v > MIN_RADIUS_M]
                s[c] += sum(valid)
                n[c] += len(valid)
    return [[sums[r][c] / counts[r][c] / 1000.0 if counts[r][c] else None for c in range(cols)]
            for r in range(rows)], cols, rows


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    grid_deg = float(sys.argv[3]) if len(sys.argv) == 4 else 1.5
    with Image.open(sys.argv[1]) as img:
        px_per_deg = img.size[0] / 360.0
    cell_px = round(grid_deg * px_per_deg)
    cells, cols, rows = cell_means(sys.argv[1], cell_px)

    # Vertex grid: latitude rows from -90 (south) to 90; columns at east longitudes
    # 0, grid_deg, ... The cells' columns start at -180 deg (the map's left edge).
    lats = [-90.0 + k * grid_deg for k in range(rows + 1)]
    shift = cols // 2  # cell column of east longitude 0

    def mean(values):
        values = [v for v in values if v is not None]
        return sum(values) / len(values) if values else None

    radius = []
    for j in range(rows + 1):
        r_south = rows - j  # cell row just south of this vertex row (cell rows count from the north)
        if j in (0, rows):
            row = cells[rows - 1] if j == 0 else cells[0]
            radius.append([mean(row)] * cols)
            continue
        line = []
        for i in range(cols):
            c_east = (i + shift) % cols
            c_west = (c_east - 1) % cols
            line.append(mean([cells[r_south - 1][c_west], cells[r_south - 1][c_east],
                              cells[r_south][c_west], cells[r_south][c_east]]))
        radius.append(line)

    missing = [(j, i) for j in range(rows + 1) for i in range(cols) if radius[j][i] is None]
    known = [radius[j][i] for j in range(rows + 1) for i in range(cols) if radius[j][i] is not None]
    for j, i in missing:
        radius[j][i] = sum(known) / len(known)
    for _ in range(FILL_ITERATIONS if missing else 0):
        for j, i in missing:
            nbrs = [radius[j][(i - 1) % cols], radius[j][(i + 1) % cols]]
            if j > 0:
                nbrs.append(radius[j - 1][i])
            if j < rows:
                nbrs.append(radius[j + 1][i])
            radius[j][i] = sum(nbrs) / len(nbrs)

    positions, map_u, triangles = grid_mesh(lats, cols, lambda j, i: radius[j][i])
    v = volume(positions, triangles)
    extent = [max(p[k] for p in positions) - min(p[k] for p in positions) for k in range(3)]
    print(f'{rows + 1} x {cols + 1} grid ({grid_deg} deg), {len(positions)} vertices, {len(triangles)} triangles; '
          f'{len(missing)} vertices filled ({100.0 * len(missing) / (len(radius) * cols):.1f}%)')
    print(f'volume {v:.4g} km^3 (equal-volume diameter {2.0 * (3.0 * v / (4.0 * math.pi)) ** (1.0 / 3.0):.2f} km); '
          f'extent {extent[0]:.2f} x {extent[1]:.2f} x {extent[2]:.2f} km')
    write_mesh(sys.argv[2], positions, [1.0] * len(positions), triangles, map_u)


if __name__ == '__main__':
    main()
