"""Convert a plate model (PDS SBN "Saturn Small Moon Shape Models", Thomas, Cassini
ISS; e.g. hyperion_30k_plt.tab), a triangle mesh in Wavefront OBJ form (e.g. the
VLT/SPHERE MPCD model of Pallas) or a VRML IndexedFaceSet (the Rosetta shape models of
67P, Steins and Lutetia on PDS SBN, *_cart.wrl) into a body mesh for assets/shapes/.

Sources and conventions are listed in assets/shapes/SOURCES.md. Download the
file there first, then (from the repository root):

    python tools/shapes/make_plate_shape.py <plates.tab | shape.obj> <output.mesh>

Format of the table: a line with the vertex and plate counts, one line of x, y, z
(km, body-fixed) per vertex, then one line of three vertex indices (from 0) per
plate. OBJ files: "v x y z" lines (km, body-fixed) and "f a b c" lines (indices
from 1; texture and normal indices are ignored). VRML: the "point [ ... ]" list (km,
body-fixed) and the "coordIndex [ a b c -1 ... ]" list (indices from 0). Either way the plates must be
counter-clockwise seen from outside (checked here: the signed volume must be
positive). The mesh has a uniform albedo and no map coordinates.

Output format: see axmesh.py. Only the Python standard library is used.
"""

import math
import sys

from axmesh import volume, write_mesh


def read_plate_table(path):
    with open(path) as f:
        lines = [line.split() for line in f if line.strip()]
    vertex_count, plate_count = int(lines[0][0]), int(lines[0][1])
    positions = [tuple(float(x) for x in line[:3]) for line in lines[1:1 + vertex_count]]
    triangles = [tuple(int(x) for x in line[:3]) for line in lines[1 + vertex_count:1 + vertex_count + plate_count]]
    if len(positions) != vertex_count or len(triangles) != plate_count:
        sys.exit('truncated table')
    return positions, triangles


def read_obj(path):
    positions, triangles = [], []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == 'v':
                positions.append(tuple(float(x) for x in fields[1:4]))
            elif fields[0] == 'f':
                if len(fields) != 4:
                    sys.exit('only triangular faces are supported')
                triangles.append(tuple(int(x.split('/')[0]) - 1 for x in fields[1:4]))
    return positions, triangles


def read_vrml(path):
    with open(path) as f:
        text = f.read()
    point = text.index('point', text.index('Coordinate'))
    body = text[text.index('[', point) + 1:text.index(']', point)]
    numbers = [float(x) for line in body.splitlines() if not line.lstrip().startswith('#') for x in line.split()]
    positions = [tuple(numbers[k:k + 3]) for k in range(0, len(numbers), 3)]
    index = text.index('coordIndex')
    ids = [int(x) for x in text[text.index('[', index) + 1:text.index(']', index)].split()]
    triangles = []
    face = []
    for i in ids:
        if i >= 0:
            face.append(i)
            continue
        if len(face) != 3:
            sys.exit('only triangular faces are supported')
        triangles.append(tuple(face))
        face = []
    return positions, triangles


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    name = sys.argv[1].lower()
    read = read_obj if name.endswith('.obj') else read_vrml if name.endswith('.wrl') else read_plate_table
    positions, triangles = read(sys.argv[1])
    vertex_count, plate_count = len(positions), len(triangles)
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
