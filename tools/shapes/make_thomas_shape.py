"""Convert a Peter Thomas shape model (PDS SBN "Small Body Optical Shape Models",
ast-sat.thomas.shape-models) into a body mesh for assets/shapes/.

Sources and conventions are listed in assets/shapes/SOURCES.md. Download the
table there first, then (from the repository root):

    python tools/shapes/make_thomas_shape.py <table.tab> <output.mesh>

e.g. m2deimos.tab -> assets/shapes/deimos.mesh

The table is a grid of planetocentric latitude, longitude and radius (km), with
both 0 and 360 deg longitude rows. Its longitudes are WEST-positive: on Thomas'
Phobos table the Stickney crater (49 deg W) shows as a depression at 49, not at
311 (checked against the IAU coordinates). The mesh keeps the grid: one vertex per
grid point, east longitude increasing with the column, so the seam column is
duplicated (map u = 0 and 1) and the poles are rows of coincident vertices; the
renderer averages normals over coincident vertices. Degenerate triangles at the
poles are dropped. The albedo is uniform (the surface map is a separate texture,
in the same planetocentric coordinates).

Output format: see axmesh.py (with map coordinates). Only the Python standard
library is used.
"""

import math
import sys

from axmesh import write_mesh


def read_table(path):
    radius = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3:
                lat, lon_west, r = (float(x) for x in parts)
                radius[(round(lat, 6), round(lon_west % 360.0, 6))] = r
    lats = sorted({lat for lat, _ in radius})
    lons = sorted({lon for _, lon in radius})
    return radius, lats, lons


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    radius, lats, lons = read_table(sys.argv[1])
    step = lons[1] - lons[0]
    columns = round(360.0 / step)
    if any(abs(lons[k] - k * step) > 1e-6 for k in range(len(lons))) or len(lons) != columns:
        sys.exit('expected evenly spaced longitudes 0 .. 360 - step')

    positions, map_u, triangles = [], [], []
    for lat in lats:
        for i in range(columns + 1):
            lon_east = (i % columns) * step  # the seam column repeats column 0 exactly
            r = radius[(lat, round((360.0 - lon_east) % 360.0, 6))]
            if abs(lat) == 90.0:
                p = (0.0, 0.0, math.copysign(r, lat))  # one exact point for the whole pole row
            else:
                la, lo = math.radians(lat), math.radians(lon_east)
                p = (r * math.cos(la) * math.cos(lo), r * math.cos(la) * math.sin(lo), r * math.sin(la))
            positions.append(p)
            map_u.append(i / columns)

    # Counter-clockwise seen from outside (latitude rows south to north, east longitude increasing).
    width = columns + 1
    for j in range(len(lats) - 1):
        for i in range(columns):
            a = j * width + i
            b, c = a + 1, a + width
            d = c + 1
            if lats[j] != -90.0:
                triangles.append((a, b, d))
            if lats[j + 1] != 90.0:
                triangles.append((a, d, c))

    volume = sum(
        (positions[a][0] * (positions[b][1] * positions[c][2] - positions[b][2] * positions[c][1]) +
         positions[a][1] * (positions[b][2] * positions[c][0] - positions[b][0] * positions[c][2]) +
         positions[a][2] * (positions[b][0] * positions[c][1] - positions[b][1] * positions[c][0])) / 6.0
        for a, b, c in triangles)
    print(f'{len(lats)} x {columns + 1} grid, {len(positions)} vertices, {len(triangles)} triangles; '
          f'volume {volume:.1f} km^3 (equal-volume radius {(3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0):.3f} km)')
    write_mesh(sys.argv[2], positions, [1.0] * len(positions), triangles, map_u)


if __name__ == '__main__':
    main()
