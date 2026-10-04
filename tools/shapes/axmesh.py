"""Writer for the body mesh files in assets/shapes/ (read by src/scene/shape_model.cpp).

Format (little endian):
    char[8]  magic "AXMESH2\\0"
    uint32   vertex count
    uint32   triangle count
    uint32   flags: bit 0 = map coordinates present
    float32  positions (x, y, z km, body-fixed: x toward the prime meridian, z along the pole)
    float32  relative albedo per vertex
    float32  map u per vertex (if flag bit 0): east longitude / 360 deg, continuous
             across each triangle (seam vertices are duplicated with u = 0 and 1);
             the latitude is taken from the vertex direction
    uint32   vertex indices (3 per triangle, counter-clockwise seen from outside)
"""

import os
import struct


def write_mesh(path, positions, albedo, triangles, map_u=None):
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'wb') as f:
        f.write(b'AXMESH2\0')
        f.write(struct.pack('<III', len(positions), len(triangles), 1 if map_u is not None else 0))
        f.write(struct.pack(f'<{3 * len(positions)}f', *(c for p in positions for c in p)))
        f.write(struct.pack(f'<{len(albedo)}f', *albedo))
        if map_u is not None:
            f.write(struct.pack(f'<{len(map_u)}f', *map_u))
        f.write(struct.pack(f'<{3 * len(triangles)}I', *(i for t in triangles for i in t)))
    print(f'-> {path} ({os.path.getsize(path) // 1024} KB)')
