"""Bake the relief maps of the Moon, Mars and the Earth (assets/textures/*_relief.jpg).

    python -E -P make_relief.py moon  <ldem_16.img>      <out.jpg>
    python -E -P make_relief.py mars  <megt90n000eb.img> <out.jpg>
    python -E -P make_relief.py earth <ETOPO5.DAT> <ne_50m_land.zip> <ne_50m_lakes.zip>
                                      <ne_50m_antarctic_ice_shelves_polys.zip> <earth.jpg> <out.jpg>

The height maps are averaged down to the 2048 x 1024 of the albedo maps (left edge at
180 W) and turned into slopes and a water mask by relief_map.py. Downloads, licenses
and what is checked are in assets/textures/SOURCES.md.

    Moon   LOLA LDEM_16: 16 pixels per degree, 16-bit little-endian, 0.5 m per count
           above a 1737.4 km sphere, from 0 E eastward.
    Mars   MOLA MEGDR megt90n000eb: 16 pixels per degree, 16-bit big-endian, meters
           above the areoid, planetocentric, from 0 E eastward.
    Earth  ETOPO5: 5 arc minutes, 16-bit big-endian, meters, rows from 90 N (none for
           the south pole), from 0 E eastward. The water is not "below sea level"
           (the Qattara depression is dry, the Great Lakes are not at sea level) but
           Natural Earth's 1:50m land without its lakes, with the Antarctic ice
           shelves as land; where the albedo map is bright there (sea ice) it is not
           open water either. Water is level: its slopes are zero.

Prints a few known heights from each source before baking.

Requires numpy and Pillow.
"""

import io
import os
import struct
import sys
import zipfile

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import relief_map  # noqa: E402

WIDTH = 2048


def read_grid(path, dtype, rows, cols):
    data = np.fromfile(path, dtype=dtype)
    if data.size != rows * cols:
        sys.exit("%s: expected %d x %d samples" % (path, cols, rows))
    return data.reshape(rows, cols).astype(np.float64)


def at(grid, lat_deg, east_lon_deg):
    """The sample of a grid that starts at 90 N, 0 E."""
    rows, cols = grid.shape
    return grid[min(rows - 1, int((90.0 - lat_deg) / 180.0 * rows)), int((east_lon_deg % 360.0) / 360.0 * cols) % cols]


def to_map(grid):
    """From a grid starting at 0 E to the albedo maps' layout (left edge 180 W),
    averaged down to WIDTH."""
    grid = np.roll(grid, grid.shape[1] // 2, axis=1)
    image = Image.fromarray(grid.astype(np.float32), "F").resize((WIDTH, WIDTH // 2), Image.BOX)
    return np.asarray(image, dtype=np.float64)


def shapefile_rings(zip_path):
    """(points, is_hole) of every polygon ring of the .shp in a Natural Earth zip.
    Outer rings run clockwise, holes counter-clockwise (ESRI Shapefile Technical
    Description, 1998)."""
    with zipfile.ZipFile(zip_path) as z:
        name = next(n for n in z.namelist() if n.endswith(".shp"))
        data = z.read(name)
    pos = 100
    while pos < len(data):
        length = struct.unpack(">i", data[pos + 4:pos + 8])[0] * 2
        record = data[pos + 8:pos + 8 + length]
        pos += 8 + length
        if struct.unpack("<i", record[:4])[0] != 5:
            continue
        parts, count = struct.unpack("<ii", record[36:44])
        starts = list(struct.unpack("<%di" % parts, record[44:44 + 4 * parts])) + [count]
        points = np.frombuffer(record, dtype="<f8", count=2 * count, offset=44 + 4 * parts).reshape(count, 2)
        for a, b in zip(starts[:-1], starts[1:]):
            ring = points[a:b]
            area = 0.5 * np.sum(ring[:-1, 0] * ring[1:, 1] - ring[1:, 0] * ring[:-1, 1])
            yield ring, area > 0.0


def earth_water(land_zip, lakes_zip, shelves_zip, albedo_path):
    scale = 4
    w, h = WIDTH * scale, WIDTH // 2 * scale
    image = Image.new("L", (w, h), 0)
    draw = ImageDraw.Draw(image)

    def fill(zip_path, outer, hole):
        rings = list(shapefile_rings(zip_path))
        for want_hole, value in ((False, outer), (True, hole)):
            for ring, is_hole in rings:
                if is_hole == want_hole:
                    xy = [((lon + 180.0) / 360.0 * w, (90.0 - lat) / 180.0 * h) for lon, lat in ring]
                    draw.polygon(xy, fill=value)

    fill(land_zip, 255, 0)
    fill(lakes_zip, 0, 255)
    fill(shelves_zip, 255, 255)
    land = np.asarray(image.resize((WIDTH, WIDTH // 2), Image.BOX), dtype=np.float64) / 255.0
    # Sea ice in the albedo map (July 2004): bright where the open sea is nearly black.
    albedo = np.asarray(Image.open(albedo_path).convert("L").resize((WIDTH, WIDTH // 2), Image.BOX), dtype=np.float64)
    ice = np.clip((albedo - 60.0) / 60.0, 0.0, 1.0)
    return (1.0 - land) * (1.0 - ice)


def main():
    args = sys.argv[1:]
    if len(args) == 3 and args[0] == "moon":
        grid = read_grid(args[1], "<i2", 2880, 5760) * 0.5e-3
        print("Tycho's floor %.2f km, its rim to the east %.2f km; Mare Imbrium %.2f km; far-side highlands %.2f km" %
              (at(grid, -43.3, -11.2), at(grid, -43.3, -9.4), at(grid, 35.0, -17.0), at(grid, 5.4, 201.4)))
        relief_map.write(args[2], to_map(grid), 1737.4)
    elif len(args) == 3 and args[0] == "mars":
        grid = read_grid(args[1], ">i2", 2880, 5760) * 1e-3
        print("Olympus Mons %.2f km (highest sample %.2f km); Hellas %.2f km (lowest %.2f km); Valles Marineris %.2f km" %
              (at(grid, 18.65, 226.2), grid.max(), at(grid, -42.4, 70.5), grid.min(), at(grid, -13.9, 300.8)))
        relief_map.write(args[2], to_map(grid), 3389.5)
    elif len(args) == 7 and args[0] == "earth":
        grid = read_grid(args[1], ">i2", 2160, 4320) * 1e-3
        print("Mount Everest's cell %.2f km; Mariana Trench %.2f km; Tibet %.2f km; Dead Sea %.2f km" %
              (at(grid, 27.99, 86.93), at(grid, 11.35, 142.2), at(grid, 33.0, 88.0), at(grid, 31.5, 35.5)))
        water = earth_water(args[2], args[3], args[4], args[5])
        relief_map.write(args[6], np.maximum(to_map(grid), 0.0), 6371.0, water)
    else:
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
