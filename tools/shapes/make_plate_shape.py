"""Convert a plate model (PDS SBN "Saturn Small Moon Shape Models", Thomas, Cassini
ISS; e.g. hyperion_30k_plt.tab), a triangle mesh in Wavefront OBJ form (e.g. the
VLT/SPHERE MPCD model of Pallas), a VRML IndexedFaceSet (the Rosetta shape models of
67P, Steins and Lutetia on PDS SBN, *_cart.wrl) or an implicitly connected quadrilateral
(ICQ) table (Gaskell's shape models on PDS SBN, e.g. phoebe_quad64q.tab; see read_icq)
into a body mesh for assets/shapes/.

Sources and conventions are listed in assets/shapes/SOURCES.md. Download the
file there first, then (from the repository root):

    python tools/shapes/make_plate_shape.py <plates.tab | shape.obj | shape.wrl | icq.tab> <output.mesh>

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


def read_icq(path):
    """Gaskell's implicitly connected quadrilateral format (document icqmodel.asc of the
    PDS SBN shape model bundles): Q, then vertices v(i, j, f) for faces f = 1..6, rows
    j = 0..Q, columns i = 0..Q. Edge and corner vertices are stored once per face; the
    duplicates are merged with Gaskell's index map, and each quadrilateral is split along
    its shorter diagonal, as in his conversion code."""
    with open(path) as f:
        q = int(f.readline().split()[0])
        rows = [line.split() for line in f if line.strip()]
    if len(rows) != 6 * (q + 1) ** 2:
        sys.exit('truncated ICQ table')
    v = {}
    n = {}
    for k, row in enumerate(rows):
        face, rest = divmod(k, (q + 1) ** 2)
        j, i = divmod(rest, q + 1)
        v[i, j, face + 1] = tuple(float(x) for x in row[:3])
        n[i, j, face + 1] = k
    for i in range(1, q):
        n[i, q, 6] = n[q - i, q, 4]
        n[i, 0, 6] = n[i, q, 2]
        n[i, 0, 5] = n[q, q - i, 1]
        n[i, 0, 4] = n[q - i, 0, 1]
        n[i, 0, 3] = n[0, i, 1]
        n[i, 0, 2] = n[i, q, 1]
    for j in range(1, q):
        n[q, j, 6] = n[j, q, 5]
        n[q, j, 5] = n[0, j, 4]
        n[q, j, 4] = n[0, j, 3]
        n[q, j, 3] = n[0, j, 2]
        n[0, j, 6] = n[q - j, q, 3]
        n[0, j, 5] = n[q, j, 2]
    for a, b in [((0, 0, 3), (0, 0, 1)), ((q, 0, 4), (0, 0, 1)), ((0, 0, 2), (0, q, 1)),
                 ((q, 0, 3), (0, q, 1)), ((0, 0, 4), (q, 0, 1)), ((q, 0, 5), (q, 0, 1)),
                 ((0, 0, 5), (q, q, 1)), ((q, 0, 2), (q, q, 1)), ((0, 0, 6), (0, q, 2)),
                 ((q, q, 3), (0, q, 2)), ((0, q, 5), (q, q, 2)), ((q, 0, 6), (q, q, 2)),
                 ((q, q, 4), (0, q, 3)), ((0, q, 6), (0, q, 3)), ((q, q, 5), (0, q, 4)),
                 ((q, q, 6), (0, q, 4))]:
        n[a] = n[b]

    def cross2(a, b):
        return ((a[1] * b[2] - a[2] * b[1]) ** 2 + (a[2] * b[0] - a[0] * b[2]) ** 2
                + (a[0] * b[1] - a[1] * b[0]) ** 2)

    triangles = []
    for face in range(1, 7):
        for i in range(q):
            for j in range(q):
                if cross2(v[i, j, face], v[i + 1, j + 1, face]) <= cross2(v[i + 1, j, face], v[i, j + 1, face]):
                    triangles.append((n[i, j, face], n[i + 1, j + 1, face], n[i + 1, j, face]))
                    triangles.append((n[i, j, face], n[i, j + 1, face], n[i + 1, j + 1, face]))
                else:
                    triangles.append((n[i, j, face], n[i, j + 1, face], n[i + 1, j, face]))
                    triangles.append((n[i + 1, j, face], n[i, j + 1, face], n[i + 1, j + 1, face]))

    # Keep only the vertices the triangles use (6 Q^2 + 2 of them).
    used = sorted({k for t in triangles for k in t})
    index = {k: m for m, k in enumerate(used)}
    positions = [tuple(float(x) for x in rows[k][:3]) for k in used]
    return positions, [tuple(index[k] for k in t) for t in triangles]


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    name = sys.argv[1].lower()
    first = open(sys.argv[1]).readline().split()
    icq = len(first) == 1  # an ICQ table starts with Q alone; a plate table with two counts
    read = (read_obj if name.endswith('.obj') else read_vrml if name.endswith('.wrl')
            else read_icq if icq else read_plate_table)
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
