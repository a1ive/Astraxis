"""Bake NASA VTAD dwarf-planet GLB base colors to planetocentric maps.

    python tools/textures/extract_dwarf_textures.py ref/dwarf assets/textures

Requires numpy and Pillow. Inputs, hashes and usage terms are recorded in
assets/textures/SOURCES.md. These are illustrative surfaces. The source meshes
and Haumea's normal map are not runtime assets.
"""

import argparse
import hashlib
import io
import json
from pathlib import Path
import struct

import numpy as np
from PIL import Image


MODELS = {
    "Eris_1_2326.glb": "1b97f1b9eea244b047a72c2a3ac0db315715b7d421033d860d669ffd78681d54",
    "Haumea_1_1000.glb": "7b31737593054bd4ef8f61138a6ad1a2577bec03b84d1c3e76e8ac78632a87f8",
    "Makemake_1_1430.glb": "e7d5fa27f50afae10f0e7290cdf7a14ca856c0b9e99b2db02f25b75441fbdbef",
}


def read_model(path):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != MODELS[path.name]:
        raise ValueError(f"{path}: unexpected source SHA-256")
    size = struct.unpack_from("<I", data, 12)[0]
    doc = json.loads(data[20:20 + size])
    blob = data[28 + size:]

    def accessor(index):
        a = doc["accessors"][index]
        view = doc["bufferViews"][a["bufferView"]]
        components = {"VEC3": 3, "VEC2": 2, "SCALAR": 1}[a["type"]]
        dtype = {5126: "<f4", 5123: "<u2"}[a["componentType"]]
        return np.frombuffer(blob, dtype=dtype, count=a["count"] * components,
                             offset=view.get("byteOffset", 0) + a.get("byteOffset", 0)
                             ).reshape(-1, components).copy()

    primitive = doc["meshes"][0]["primitives"][0]
    positions = accessor(primitive["attributes"]["POSITION"]).astype(np.float64)
    # Eris and Makemake have Y as north; Haumea's short/polar axis is Z.
    # Map Y-up to body Z-up with a rotation (x, y, z) -> (x, -z, y).
    if not path.name.startswith("Haumea"):
        positions = positions[:, [0, 2, 1]] * [1, -1, 1]
    uv = accessor(primitive["attributes"]["TEXCOORD_0"])
    indices = accessor(primitive["indices"]).reshape(-1, 3)
    material = doc["materials"][primitive["material"]]["pbrMetallicRoughness"]
    texture = doc["textures"][material["baseColorTexture"]["index"]]
    image = doc["images"][texture["source"]]
    view = doc["bufferViews"][image["bufferView"]]
    offset = view.get("byteOffset", 0)
    picture = Image.open(io.BytesIO(blob[offset:offset + view["byteLength"]]))
    if picture.getchannel("A").getextrema() != (255, 255):
        raise ValueError(f"{path}: base color is not opaque")
    rgb = np.asarray(picture.convert("RGB"), dtype=np.float64) / 255
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    return positions, uv, indices, linear


def sample(image, uv):
    """glTF linear filtering with its default repeat wrap, in linear RGB."""
    h, w = image.shape[:2]
    xy = uv * [w, h] - 0.5
    ij = np.floor(xy).astype(int)
    f = xy - ij
    x, y = ij.T
    a = image[y % h, x % w] * (1 - f[:, 0, None]) + image[y % h, (x + 1) % w] * f[:, 0, None]
    b = image[(y + 1) % h, x % w] * (1 - f[:, 0, None]) + image[(y + 1) % h, (x + 1) % w] * f[:, 0, None]
    return a * (1 - f[:, 1, None]) + b * f[:, 1, None]


def bake(positions, uv, triangles, image, width=4096):
    """Intersect radial rays with source triangles, then interpolate their UVs.

    Work within a spherical cap around each triangle, including caps crossing
    the longitude seam or a pole. Using actual vertex positions also accounts
    for Haumea's ellipsoidal UVs when producing a planetocentric latitude map.
    """
    height = width // 2
    latitudes = np.pi / 2 - (np.arange(height) + 0.5) * np.pi / height
    longitudes = (np.arange(width) + 0.5) * 2 * np.pi / width
    colors = np.zeros((height, width, 3), dtype=np.float64)
    covered = np.zeros((height, width), dtype=bool)
    for tri in triangles:
        vertices = positions[tri]
        unit = vertices / np.linalg.norm(vertices, axis=1)[:, None]
        center = unit.sum(axis=0)
        center /= np.linalg.norm(center)
        radius = np.arccos(np.clip(unit @ center, -1, 1)).max() + 1e-7
        lat = np.arcsin(center[2])
        lon = np.arctan2(center[1], center[0])
        rows = np.flatnonzero(np.abs(latitudes - lat) <= radius)
        if abs(lat) + radius >= np.pi / 2:
            cols = np.arange(width)
        else:
            half_lon = np.arcsin(np.sin(radius) / np.cos(lat))
            delta = (longitudes - lon + np.pi) % (2 * np.pi) - np.pi
            cols = np.flatnonzero(np.abs(delta) <= half_lon)
        rr, cc = np.meshgrid(rows, cols, indexing="ij")
        rr, cc = rr.ravel(), cc.ravel()
        phi, lam = latitudes[rr], longitudes[cc]
        rays = np.column_stack((np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)))
        v0, v1, v2 = vertices
        e1, e2 = v1 - v0, v2 - v0
        normal = np.cross(e1, e2)
        determinant = rays @ normal
        if np.dot(normal, normal) < 1e-16:
            continue
        distance = np.divide(v0 @ normal, determinant, out=np.zeros_like(determinant),
                             where=np.abs(determinant) > 1e-16)
        hit = distance[:, None] * rays - v0
        b1 = hit @ np.cross(e2, normal) / (normal @ normal)
        b2 = hit @ np.cross(normal, e1) / (normal @ normal)
        mask = (distance > 0) & (b1 >= -1e-8) & (b2 >= -1e-8) & (b1 + b2 <= 1 + 1e-8)
        tex_uv = uv[tri[0]] + b1[mask, None] * (uv[tri[1]] - uv[tri[0]]) + b2[mask, None] * (uv[tri[2]] - uv[tri[0]])
        colors[rr[mask], cc[mask]] = sample(image, tex_uv)
        covered[rr[mask], cc[mask]] = True
    if not covered.all():
        raise ValueError(f"Source mesh leaves {np.count_nonzero(~covered)} pixels uncovered")
    # Supersample at twice the final resolution; average radiance before sRGB.
    colors = colors.reshape(height // 2, 2, width // 2, 2, 3).mean(axis=(1, 3))
    srgb = np.where(colors <= 0.0031308, colors * 12.92, 1.055 * colors ** (1 / 2.4) - 0.055)
    return Image.fromarray(np.uint8(np.clip(srgb * 255 + 0.5, 0, 255)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename in MODELS:
        name = filename.split("_")[0].lower()
        result = bake(*read_model(args.source_dir / filename))
        output = args.output_dir / f"{name}.jpg"
        result.save(output, quality=90, subsampling=0)
        print(f"{output}: {result.width}x{result.height}", flush=True)


if __name__ == "__main__":
    main()
