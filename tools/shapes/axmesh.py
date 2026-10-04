"""Writer for the body mesh files in assets/shapes/ (read by src/scene/shape_model.cpp),
and a builder for meshes given as a latitude/longitude grid of radii.

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

import math
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


def grid_mesh(lats, columns, radius):
    """A mesh from planetocentric radii on a grid: `lats` ascending from -90 to 90 deg,
    `columns` evenly spaced east longitudes from 0 deg, radius(j, i) in km for
    latitude row j and column i (0 <= i < columns). One vertex per grid point; the
    seam column is duplicated (map u = 0 and 1) and each pole is a row of coincident
    vertices (the renderer averages normals over coincident vertices). Degenerate
    triangles at the poles are dropped. Returns positions, map_u, triangles."""
    positions, map_u, triangles = [], [], []
    step = 360.0 / columns
    for j, lat in enumerate(lats):
        for i in range(columns + 1):
            r = radius(j, i % columns)  # the seam column repeats column 0 exactly
            if abs(lat) == 90.0:
                p = (0.0, 0.0, math.copysign(r, lat))  # one exact point for the whole pole row
            else:
                la, lo = math.radians(lat), math.radians((i % columns) * step)
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
    return positions, map_u, triangles


def volume(positions, triangles):
    """Signed volume (km^3; positive for a closed, outward-facing mesh)."""
    total = 0.0
    for a, b, c in triangles:
        (ax, ay, az), (bx, by, bz), (cx, cy, cz) = positions[a], positions[b], positions[c]
        total += (ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)) / 6.0
    return total
