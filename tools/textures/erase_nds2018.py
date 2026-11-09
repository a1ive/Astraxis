"""Takes the dark spot NDS-2018 and its bright companion out of neptune.jpg.

    python -E -P erase_nds2018.py <neptune.jpg>      (needs numpy and Pillow)

make_ice_giant_textures.py adds the clouds of one Hubble map (OPAL 2019a) to the zonal
profile, and that map shows the spot. The scene now draws the dark spots itself, each
only while it existed, so the texture must not have one all the time. Run this after
make_ice_giant_textures.py; it rewrites the file.

Within a soft ellipse around the darkest place of the blue channel between 8 and 30 N, the
departures from the zonal median are removed. The fainter secondary feature of 2019 (about
4 percent in F467M, Wong et al. 2022 Table 3) stays: it does not show in true color.
"""

import sys

import numpy as np
from PIL import Image

# Half widths of the erased region (deg): the spot is about 20 x 9 deg in 2019.
HALF_LON_DEG = 26.0
HALF_LAT_DEG = 11.0
SEARCH_LAT_DEG = (8.0, 30.0)
SMOOTH_DEG = 3.0


def to_linear(srgb):
    return np.where(srgb <= 0.04045, srgb / 12.92, ((srgb + 0.055) / 1.055) ** 2.4)


def to_srgb(linear):
    return np.where(linear <= 0.0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - 0.055)


def smooth_lon_lat(field, sigma_px):
    radius = int(3 * sigma_px)
    kernel = np.exp(-0.5 * (np.arange(-radius, radius + 1) / sigma_px) ** 2)
    kernel /= kernel.sum()
    out = np.apply_along_axis(lambda row: np.convolve(np.concatenate((row[-radius:], row, row[:radius])), kernel, 'valid'), 1, field)
    return np.apply_along_axis(lambda col: np.convolve(col, kernel, 'same'), 0, out)


def main():
    path = sys.argv[1]
    linear = to_linear(np.asarray(Image.open(path).convert('RGB'), dtype=np.float64) / 255.0)
    height, width, _ = linear.shape
    lat = 90.0 - (np.arange(height) + 0.5) * 180.0 / height
    lon = (np.arange(width) + 0.5) * 360.0 / width

    zonal = np.median(linear, axis=1, keepdims=True)
    departure = linear / zonal - 1.0

    blue = smooth_lon_lat(departure[:, :, 2], SMOOTH_DEG * width / 360.0)
    rows = (lat > SEARCH_LAT_DEG[0]) & (lat < SEARCH_LAT_DEG[1])
    search = np.where(rows[:, None], blue, 1.0)
    y, x = np.unravel_index(np.argmin(search), search.shape)
    print(f'darkest at {lat[y]:+.1f} deg, {lon[x]:.1f} E: {100.0 * blue[y, x]:.1f} percent in blue')

    dlon = (lon[None, :] - lon[x] + 180.0) % 360.0 - 180.0
    dlat = lat[:, None] - lat[y]
    r = np.sqrt((dlon / HALF_LON_DEG) ** 2 + (dlat / HALF_LAT_DEG) ** 2)
    keep = np.clip((r - 0.7) / 0.3, 0.0, 1.0)  # 0 inside, 1 outside, soft
    keep = keep * keep * (3.0 - 2.0 * keep)

    result = zonal * (1.0 + departure * keep[..., None])
    out = np.round(to_srgb(np.clip(result, 0.0, 1.0)) * 255.0).astype(np.uint8)
    Image.fromarray(out).save(path, quality=92)
    print(f'erased {100.0 * np.mean(r < 1.0):.2f} percent of the map')


if __name__ == '__main__':
    main()
