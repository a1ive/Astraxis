"""Bake the fictional Pandora surface, relief and cloud maps (NumPy and Pillow).

    python -E -P tools/textures/make_avatar_pandora.py assets/textures
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
import relief_map
from make_kerbol_textures import (color, continents, fbm, linear_to_srgb8, mix,
                                  octaves_for, ramp, smoothstep, sphere_grid)


def surface(p, lat, width):
    land, elevation, depth = continents(p, lat, 2154, 0.62, width)
    sea = ramp([(0.0, color('#509b98')), (0.08, color('#206b7c')),
                (0.24, color('#124862')), (0.8, color('#102c50'))], depth)
    moisture = fbm(p, 423, 3.0, 5) + 0.17 * np.cos(2 * lat)
    forest = mix(color('#62834b'), color('#235e4d'), smoothstep(-0.38, 0.45, moisture))
    ridge = (1 - np.abs(fbm(p, 800, 5, octaves_for(5, width)))) ** 5
    height_km = elevation ** 1.2 * (2.2 + 4.8 * ridge)
    ground = mix(forest, color('#8e837a'), smoothstep(1.8, 4.2, height_km))
    ground = mix(ground, color('#cdd5cf'), smoothstep(4.5, 6.8, height_km))
    ground = mix(ground, color('#a2ac7e'), 0.55 * smoothstep(0.025, 0.0, elevation))
    return np.where(land[..., None], ground, sea), height_km, ~land


def cloud_cover(p, lat, width):
    # Twist local patches of a continuous spherical field to curl the cloud fronts.
    q = p.copy()
    for lon_deg, lat_deg, radius_deg, twist in [
        (-132, 36, 16, 4.0), (-51, -32, 13, -4.4), (38, 47, 16, 4.5),
        (107, -44, 17, -3.8), (161, 19, 10, 3.0),
    ]:
        lon, phi, radius = np.radians([lon_deg, lat_deg, radius_deg])
        axis = np.array([np.cos(phi) * np.cos(lon), np.cos(phi) * np.sin(lon), np.sin(phi)])
        along = np.einsum('...i,i->...', q, axis)
        angle = twist * np.exp((along - 1) / (1 - np.cos(radius)))
        cosine = np.cos(angle)[..., None]
        q = q * cosine + np.cross(axis, q) * np.sin(angle)[..., None] + along[..., None] * axis * (1 - cosine)

    systems = fbm(q, 826, 2.6, 4)
    stretched = q * np.array([1.0, 1.0, 3.7])
    stretched[..., 2] += 0.55 * systems
    fronts = fbm(stretched, 827, 4, octaves_for(15, width), gain=0.68)
    small_clouds = fbm(q, 828, 18, octaves_for(18, width), gain=0.68)
    tropical = np.exp(-((np.degrees(lat) + 14 * systems) / 15) ** 2)
    return smoothstep(-0.07, 0.55, 0.65 * systems + fronts + 0.35 * small_clouds + 0.10 * tropical)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    out = parser.parse_args().output
    out.mkdir(parents=True, exist_ok=True)
    width = 2048
    p, lat, _ = sphere_grid(width)
    albedo, height_km, water = surface(p, lat, width)
    Image.fromarray(linear_to_srgb8(albedo)).save(out / 'avatar_pandora.jpg', quality=92, optimize=True)
    relief_map.write(out / 'avatar_pandora_relief.jpg', height_km, 5723.5, water)
    coverage = cloud_cover(p, lat, width)
    Image.fromarray(np.rint(255 * coverage).astype(np.uint8)).save(
        out / 'avatar_pandora_clouds.jpg', quality=92, optimize=True)


if __name__ == '__main__':
    main()
