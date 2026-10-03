"""Convert the New Horizons shape model of Arrokoth (Porter et al. 2024, PDS SBN)
into a body mesh for assets/shapes/.

Sources and conventions are listed in assets/shapes/SOURCES.md. Download the
originals there first, then (from the repository root):

    python tools/shapes/make_arrokoth.py <source_dir> [output]

Expected files in <source_dir>:
    arrokoth_porter_2024_v01.obj          shape (km, body-fixed, 40960 facets)
    albedo_arrokoth4_fp36h2_masked1.png   albedo map the OBJ's texture coordinates refer to

The albedo map covers only what New Horizons saw (roughly the southern
hemisphere); the rest holds one constant mask value. Each vertex takes the
albedo under its texture coordinate; unimaged vertices are filled by diffusion
over the mesh from their imaged neighbours, fading to the mean albedo farther
away, so the unseen side renders as a plain surface. The albedo is stored
relative to the mean of the imaged part (the scene's color sets the tint and
overall brightness).

Output format (little endian), see src/scene/shape_model.cpp:
    char[8]  magic "AXMESH1\\0"
    uint32   vertex count
    uint32   triangle count
    float32  positions (x, y, z km) per vertex
    float32  relative albedo per vertex
    uint32   vertex indices (3 per triangle, counter-clockwise seen from outside)

Requires Pillow.
"""

import os
import struct
import sys

from PIL import Image

OBJ = 'arrokoth_porter_2024_v01.obj'
ALBEDO = 'albedo_arrokoth4_fp36h2_masked1.png'
DEFAULT_OUTPUT = os.path.join('assets', 'shapes', 'arrokoth.mesh')

# PNG pixel -> albedo, from the albedo product label (albedo_arrokoth4_fp36h2_masked1.lblx).
PIXEL_SCALE = 1361975.975
PIXEL_OFFSET = 0.03188

FILL_ITERATIONS = 300
FILL_PULL_TO_MEAN = 0.03  # per iteration: unimaged areas fade to the mean over ~1/0.03 rings of vertices
MIN_RELATIVE, MAX_RELATIVE = 0.25, 4.0


def read_obj(path):
    positions, uvs, triangles = [], [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == 'v':
                positions.append(tuple(float(x) for x in parts[1:4]))
            elif parts[0] == 'vt':
                uvs.append((float(parts[1]), float(parts[2])))
            elif parts[0] == 'f':
                corners = [tuple(int(i) for i in p.split('/')[:2]) for p in parts[1:]]
                if len(corners) != 3 or any(v != t for v, t in corners):
                    sys.exit('expected triangles whose vertex and texture indices agree')
                triangles.append(tuple(v - 1 for v, _ in corners))
    if len(uvs) != len(positions):
        sys.exit('expected one texture coordinate per vertex')
    return positions, uvs, triangles


def sample_albedo(path, uvs):
    """Albedo under each texture coordinate, or None where the map is masked.
    The PNG is stored top row first (it is the FITS array flipped for display),
    so v = 0 is the bottom row, as usual for OBJ texture coordinates."""
    with Image.open(path) as img:
        img = img.convert('I')
        w, h = img.size
        px = img.load()
        counts = {}
        for y in range(h):
            for x in range(w):
                counts[px[x, y]] = counts.get(px[x, y], 0) + 1
        mask = max(counts, key=counts.get)  # the unimaged background (~84% of the map)
        out = []
        for u, v in uvs:
            x = min(w - 1, max(0, int(u * w)))
            y = min(h - 1, max(0, int((1.0 - v) * h)))
            p = px[x, y]
            out.append(None if p == mask else p / PIXEL_SCALE + PIXEL_OFFSET)
    return out


def fill_unimaged(albedo, triangles):
    neighbours = [set() for _ in albedo]
    for a, b, c in triangles:
        neighbours[a].update((b, c))
        neighbours[b].update((a, c))
        neighbours[c].update((a, b))
    known = [a for a in albedo if a is not None]
    mean = sum(known) / len(known)
    unknown = [i for i, a in enumerate(albedo) if a is None]
    values = [mean if a is None else a for a in albedo]
    for _ in range(FILL_ITERATIONS):
        updated = {}
        for i in unknown:
            nbr = sum(values[j] for j in neighbours[i]) / len(neighbours[i])
            updated[i] = (1.0 - FILL_PULL_TO_MEAN) * nbr + FILL_PULL_TO_MEAN * mean
        for i, value in updated.items():
            values[i] = value
    return values, mean, len(known)


def main():
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    src_dir = sys.argv[1]
    output = sys.argv[2] if len(sys.argv) == 3 else DEFAULT_OUTPUT

    positions, uvs, triangles = read_obj(os.path.join(src_dir, OBJ))
    albedo = sample_albedo(os.path.join(src_dir, ALBEDO), uvs)
    values, mean, imaged = fill_unimaged(albedo, triangles)
    relative = [min(MAX_RELATIVE, max(MIN_RELATIVE, a / mean)) for a in values]

    lo = [min(p[k] for p in positions) for k in range(3)]
    hi = [max(p[k] for p in positions) for k in range(3)]
    print(f'{len(positions)} vertices, {len(triangles)} triangles; '
          f'extent {hi[0] - lo[0]:.3f} x {hi[1] - lo[1]:.3f} x {hi[2] - lo[2]:.3f} km')
    print(f'imaged vertices: {imaged} ({100.0 * imaged / len(positions):.1f}%), mean albedo {mean:.5f}; '
          f'relative albedo {min(relative):.2f} .. {max(relative):.2f}')

    os.makedirs(os.path.dirname(output) or '.', exist_ok=True)
    with open(output, 'wb') as f:
        f.write(b'AXMESH1\0')
        f.write(struct.pack('<II', len(positions), len(triangles)))
        f.write(struct.pack(f'<{3 * len(positions)}f', *(c for p in positions for c in p)))
        f.write(struct.pack(f'<{len(relative)}f', *relative))
        f.write(struct.pack(f'<{3 * len(triangles)}I', *(i for t in triangles for i in t)))
    print(f'-> {output} ({os.path.getsize(output) // 1024} KB)')


if __name__ == '__main__':
    main()
