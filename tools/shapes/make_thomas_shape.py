"""Convert a Peter Thomas shape model (PDS SBN "Small Body Optical Shape Models",
ast-sat.thomas.shape-models) into a body mesh for assets/shapes/.

Sources and conventions are listed in assets/shapes/SOURCES.md. Download the
table there first, then (from the repository root):

    python tools/shapes/make_thomas_shape.py <table.tab> <output.mesh>

e.g. m1phobos.tab -> assets/shapes/phobos.mesh, m2deimos.tab -> assets/shapes/deimos.mesh

The table is a grid of planetocentric latitude, longitude and radius (km), with
both 0 and 360 deg longitude rows. Its longitudes are WEST-positive: on Thomas'
Phobos table the Stickney crater (49 deg W) shows as a depression at 49, not at
311 (checked against the IAU coordinates). The mesh keeps the grid (see
axmesh.grid_mesh) with east longitude increasing with the column. The albedo is
uniform (the surface map is a separate texture, in the same planetocentric
coordinates).

Output format: see axmesh.py (with map coordinates). Only the Python standard
library is used.
"""

import math
import sys

from axmesh import grid_mesh, volume, write_mesh


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

    def radius_at(j, i):
        lon_east = i * step
        return radius[(lats[j], round((360.0 - lon_east) % 360.0, 6))]

    positions, map_u, triangles = grid_mesh(lats, columns, radius_at)
    v = volume(positions, triangles)
    print(f'{len(lats)} x {columns + 1} grid, {len(positions)} vertices, {len(triangles)} triangles; '
          f'volume {v:.1f} km^3 (equal-volume radius {(3.0 * v / (4.0 * math.pi)) ** (1.0 / 3.0):.3f} km)')
    write_mesh(sys.argv[2], positions, [1.0] * len(positions), triangles, map_u)


if __name__ == '__main__':
    main()
