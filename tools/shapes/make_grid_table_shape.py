"""Convert a shape model given as a table of planetocentric latitude, longitude and
radius on a regular grid into a body mesh for assets/shapes/:
  - Peter Thomas' models (PDS SBN "Small Body Optical Shape Models",
    ast-sat.thomas.shape-models): columns latitude, longitude, radius;
  - Philip Stooke's models (PDS SBN "Stooke Small Body Shape Models",
    small_bodies.stooke.shape-models): columns longitude, latitude, radius
    (pass --lon-first).

Sources and conventions are listed in assets/shapes/SOURCES.md. Download the
table there first, then (from the repository root):

    python tools/shapes/make_grid_table_shape.py [--lon-first] [--east] <table.tab> <output.mesh>

e.g. m1phobos.tab -> phobos.mesh, m2deimos.tab -> deimos.mesh,
     --lon-first j5amalthea.tab -> amalthea.mesh, --lon-first j14thebe.tab -> thebe.mesh,
     --lon-first --east 1682q1halley.tab -> halley.mesh

Both sets list 0 and 360 deg longitude rows (Halley's table only 360), and for
satellites their longitudes are WEST-positive: on Thomas' Phobos table the Stickney crater (49 deg W) shows as
a depression at 49, not at 311 (checked against the IAU coordinates); Stooke's
labels state it (and that asteroids and comets are east-positive, V2.0). The mesh keeps the grid (see axmesh.grid_mesh) with east
longitude increasing with the column. The albedo is uniform (a surface map, if
any, is a separate texture in the same planetocentric coordinates).

Output format: see axmesh.py (with map coordinates). Only the Python standard
library is used.
"""

import math
import sys

from axmesh import grid_mesh, volume, write_mesh


def read_table(path, lon_first, east):
    radius = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3:
                a, b, r = (float(x) for x in parts)
                lat, lon = (b, a) if lon_first else (a, b)
                lon_west = -lon if east else lon
                radius[(round(lat, 6), round(lon_west % 360.0, 6))] = r
    lats = sorted({lat for lat, _ in radius})
    lons = sorted({lon for _, lon in radius})
    return radius, lats, lons


def main():
    flags = {'--lon-first', '--east'}
    args = [a for a in sys.argv[1:] if a not in flags]
    if len(args) != 2:
        sys.exit(__doc__)
    radius, lats, lons = read_table(args[0], '--lon-first' in sys.argv[1:], '--east' in sys.argv[1:])
    step = lons[1] - lons[0]
    columns = round(360.0 / step)
    if any(abs(lons[k] - k * step) > 1e-6 for k in range(len(lons))) or len(lons) != columns:
        sys.exit('expected evenly spaced longitudes 0 .. 360 - step')
    if lats[0] != -90.0 or lats[-1] != 90.0:
        sys.exit('expected latitudes from -90 to 90')

    def radius_at(j, i):
        lon_east = i * step
        return radius[(lats[j], round((360.0 - lon_east) % 360.0, 6))]

    positions, map_u, triangles = grid_mesh(lats, columns, radius_at)
    v = volume(positions, triangles)
    print(f'{len(lats)} x {columns + 1} grid, {len(positions)} vertices, {len(triangles)} triangles; '
          f'volume {v:.1f} km^3 (equal-volume radius {(3.0 * v / (4.0 * math.pi)) ** (1.0 / 3.0):.3f} km)')
    write_mesh(args[1], positions, [1.0] * len(positions), triangles, map_u)


if __name__ == '__main__':
    main()
