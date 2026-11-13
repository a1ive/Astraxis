"""Generate illustrative cloud maps for Beta Pic b, c and d (NumPy, Pillow).

Run from the repository root: python tools/textures/make_betapic_textures.py
"""

from pathlib import Path

import numpy as np
from PIL import Image

from make_kerbol_textures import color, fbm, linear_to_srgb8, mix, sphere_grid


def cloud_map(width, seed, palette, band_share, shear):
    p, lat, _ = sphere_grid(width)
    # Zonal shear stretches the clouds along longitude. Spherical noise keeps
    # the same pattern across the longitude seam and around the poles.
    stretched = p * [1.0, 1.0, shear]
    billows = fbm(p, seed, 3.5, 6)
    streaks = fbm(stretched, seed + 1, 2.8, 5)
    latitude = lat + 0.035 * billows + 0.025 * streaks
    zonal = np.stack([latitude * 5.0, np.full_like(lat, 0.37), np.full_like(lat, 0.61)], axis=-1)
    bands = fbm(zonal, seed + 2, 1.0, 4, gain=0.45)
    small = fbm(p + 0.11 * billows[..., None], seed + 3, 12.0, 5)
    cover = np.clip(0.5 + band_share * bands + (0.80 - band_share) * streaks +
                    0.20 * billows + 0.13 * small, 0.0, 1.0)
    return mix(color(palette[0]), color(palette[1]), cover)


def main():
    out = Path(__file__).resolve().parents[2] / "assets/textures"
    # Display palettes, not measured visible colors. b: thick warm clouds;
    # c: subtler bands; d: a cooler T-type giant with patchier cloud cover.
    planets = {
        "b": (341, ("#885747", "#e9bd9c"), 0.55, 4.7),
        "c": (578, ("#817064", "#dec8af"), 0.42, 3.8),
        "d": (829, ("#656c7e", "#b4b9c5"), 0.19, 2.5),
    }
    for name, params in planets.items():
        pixels = linear_to_srgb8(cloud_map(1536, *params))
        Image.fromarray(pixels).save(out / f"beta_pic_{name}.jpg", quality=90, optimize=True)


if __name__ == "__main__":
    main()
