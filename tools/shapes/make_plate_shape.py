"""Convert a plate model (PDS SBN "Saturn Small Moon Shape Models", Thomas, Cassini
ISS; e.g. hyperion_30k_plt.tab) into a body mesh for assets/shapes/.

Sources and conventions are listed in assets/shapes/SOURCES.md. Download the
table there first, then (from the repository root):

    python tools/shapes/make_plate_shape.py <plates.tab> <output.mesh>

Format of the table: a line with the vertex and plate counts, one line of x, y, z
(km, body-fixed) per vertex, then one line of three vertex indices (from 0) per
plate, counter-clockwise seen from outside (checked here: the signed volume must be
positive). The mesh has a uniform albedo and no map coordinates.

Output format: see axmesh.py. Only the Python standard library is used.
"""

import math
import sys

from axmesh import volume, write_mesh


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    with open(sys.argv[1]) as f:
        lines = [line.split() for line in f if line.strip()]
    vertex_count, plate_count = int(lines[0][0]), int(lines[0][1])
    positions = [tuple(float(x) for x in line[:3]) for line in lines[1:1 + vertex_count]]
    triangles = [tuple(int(x) for x in line[:3]) for line in lines[1 + vertex_count:1 + vertex_count + plate_count]]
    if len(positions) != vertex_count or len(triangles) != plate_count:
        sys.exit('truncated table')
    if any(not 0 <= i < vertex_count for t in triangles for i in t):
        sys.exit('plate vertex index out of range (expected indices from 0)')
    v = volume(positions, triangles)
    if v <= 0.0:
        sys.exit('plates are not counter-clockwise seen from outside')
    extent = [max(p[k] for p in positions) - min(p[k] for p in positions) for k in range(3)]
    print(f'{vertex_count} vertices, {plate_count} plates; volume {v:.4g} km^3 '
          f'(equal-volume radius {(3.0 * v / (4.0 * math.pi)) ** (1.0 / 3.0):.2f} km); '
          f'extent {extent[0]:.1f} x {extent[1]:.1f} x {extent[2]:.1f} km')
    write_mesh(sys.argv[2], positions, [1.0] * vertex_count, triangles)


if __name__ == '__main__':
    main()
