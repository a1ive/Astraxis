"""Make uranus.jpg and neptune.jpg from the Hubble OPAL global maps.

    python make_ice_giant_textures.py <opal_dir> <output_dir> [width]

<opal_dir> holds the OPAL FITS maps in F657N, F547M and F467M
(hlsp_opal_hst_wfc3-uvis_<planet>-<rotation>_<filter>_v1_globalmap.fits) and the
per-year README files (hlsp_opal_hst_wfc3-uvis_<planet>-<year>_all_v*_readme.txt),
which give the factors that turn the FITS values into I/F. Download them as
described in assets/textures/SOURCES.md.

Each OPAL map only covers the hemisphere that faced Earth, and its edge (high
emission angles) is smeared into colored streaks, so no single map makes a
texture. Instead:

1. Every rotation's F657N/F547M/F467M maps become linear I/F (red, green, blue;
   Minnaert-corrected by OPAL, so limb darkening is already removed). EXCLUDE
   lists rotations whose color is off.
2. Pixels within EDGE_LAT_DEG of the map's edge in latitude (or EDGE_LON_DEG of a
   gap in longitude) are dropped; the zonal median of the rest, per rotation,
   and then the median over rotations give the latitude profile.
3. Latitudes no rotation covers (Uranus south of ~45 S, Neptune north of ~55 N)
   take the profile of the opposite latitude, scaled to meet the covered profile
   at the boundary (the scale fades out over FILL_BLEND_DEG).
4. One rotation's map (DETAIL) adds its departures from its own zonal median,
   times its gain, faded out towards its edge, so the bands carry real clouds.
5. The result is resampled to planetocentric latitude (OPAL maps are
   planetographic; the scene's texture lookup is planetocentric) with east
   longitude increasing to the right from 0 E at the left edge, then encoded as
   sRGB. I/F is reflectance relative to sunlight, so it is used as albedo
   directly.

Requires numpy and Pillow.
"""

import glob
import math
import os
import re
import sys

import numpy as np
from PIL import Image

FILTERS = ('F657N', 'F547M', 'F467M')  # red, green, blue
EDGE_LAT_DEG = 22.0
EDGE_LON_DEG = 22.0
FILL_BLEND_DEG = 20.0
DETAIL_FADE_DEG = 8.0
DETAIL_CLIP = (-0.35, 0.5)
# Uranus 2022b's F657N map is 1.4 times brighter than 2022a's (same I/F factor in the
# README) and every other year's, so the composite turns pink.
EXCLUDE = {'uranus': {'2022b'}}

# planet: (equatorial, polar radius in km as in the OPAL READMEs, detail rotation,
#          east longitude increases to the right in the detail map, detail gain).
# Neptune's detail map has the northern dark spot NDS-2018 (~20 N), at full contrast;
# 2018b shows it too, but next to a mosaic seam whose streaks reach 30 deg in.
PLANETS = {
    'uranus': (25559.0, 24973.0, '2025a', False, 0.75),   # 2021+: right edge = 0 E, increasing left
    'neptune': (24764.0, 24341.0, '2019a', True, 1.0),   # left edge = 360 W, decreasing right
}


def read_fits(path):
    with open(path, 'rb') as f:
        data = f.read()
    header = {}
    offset = 0
    done = False
    while not done:
        block = data[offset:offset + 2880].decode('ascii')
        offset += 2880
        for i in range(36):
            card = block[i * 80:(i + 1) * 80]
            if card[:8].strip() == 'END':
                done = True
                break
            if card[8:10] == '= ':
                header[card[:8].strip()] = card[10:].split('/')[0].strip().strip("'").strip()
    if header.get('BITPIX') != '-32' or header.get('NAXIS') != '2':
        raise ValueError(f'{path}: expected a 2-D float32 image')
    nx, ny = int(header['NAXIS1']), int(header['NAXIS2'])
    image = np.frombuffer(data[offset:offset + nx * ny * 4], dtype='>f4').reshape(ny, nx)
    return image.astype(np.float64)  # row 0 = +90 deg latitude


def read_scales(opal_dir, planet):
    """I/F factor per (year, filter), from the READMEs' Minnaert/I/F tables."""
    scales = {}
    for path in glob.glob(os.path.join(opal_dir, f'*_{planet}-*_readme.txt')):
        year = re.search(planet + r'-(\d{4})', os.path.basename(path)).group(1)
        with open(path, encoding='latin-1') as f:
            for line in f:
                words = line.split()
                if words and words[0] in FILTERS:
                    scales.setdefault(year, {})[words[0]] = float(words[-1])
    return scales


def load_rotations(opal_dir, planet):
    scales = read_scales(opal_dir, planet)
    rotations = {}
    pattern = os.path.join(opal_dir, f'hlsp_opal_hst_wfc3-uvis_{planet}-*_f547m_v1_globalmap.fits')
    for green in sorted(glob.glob(pattern)):
        rotation = re.search(planet + r'-(\d{4}[a-z])_', os.path.basename(green)).group(1)
        year = rotation[:4]
        if rotation in EXCLUDE.get(planet, ()):
            continue
        paths = [green.replace('_f547m_', '_' + f.lower() + '_') for f in FILTERS]
        if not all(os.path.exists(p) for p in paths) or not all(f in scales.get(year, {}) for f in FILTERS):
            continue
        maps = [read_fits(p) * scales[year][f] for p, f in zip(paths, FILTERS)]
        rgb = np.stack(maps, axis=-1)[:, :-1]  # drop the repeated 360 deg column
        rgb[~np.isfinite(rgb)] = 0.0
        rotations[rotation] = rgb
    return rotations


def run_min(valid, radius, axis, wrap):
    """True where every pixel within radius along axis is valid."""
    pad = [(0, 0), (0, 0)]
    pad[axis] = (radius, radius)
    padded = np.pad(valid, pad, mode='wrap' if wrap else 'edge').astype(np.int32)
    c = np.cumsum(padded, axis=axis)
    c = np.concatenate([np.zeros_like(np.take(c, [0], axis=axis)), c], axis=axis)
    n = valid.shape[axis]
    window = np.take(c, np.arange(2 * radius + 1, n + 2 * radius + 1), axis=axis) - np.take(c, np.arange(n), axis=axis)
    return window == 2 * radius + 1


def good_mask(rgb):
    ppd = rgb.shape[1] / 360.0
    valid = np.all(rgb > 0.0, axis=-1)
    # Off the planet's edge rows run past the poles: count them as valid (edge padding).
    good = run_min(valid, int(round(EDGE_LAT_DEG * ppd)), 0, False)
    return good & run_min(valid, int(round(EDGE_LON_DEG * ppd)), 1, True)


def zonal_median(rgb, good):
    rows = rgb.shape[0]
    profile = np.full((rows, 3), np.nan)
    for y in range(rows):
        if good[y].mean() > 0.3:
            profile[y] = np.median(rgb[y, good[y]], axis=0)
    return profile


def smooth_rows(profile, sigma_rows):
    """Gaussian smoothing along latitude of a profile without gaps."""
    k = np.arange(-3 * int(sigma_rows) - 1, 3 * int(sigma_rows) + 2)
    w = np.exp(-0.5 * (k / sigma_rows) ** 2)
    w /= w.sum()
    padded = np.pad(profile, ((len(k) // 2, len(k) // 2), (0, 0)), mode='edge')
    return np.stack([np.convolve(padded[:, c], w, mode='valid') for c in range(3)], axis=-1)


def fill_profile(profile, lat):
    """Fills uncovered latitudes from the opposite ones, scaled to join smoothly."""
    covered = np.all(np.isfinite(profile), axis=-1)
    idx = np.where(covered)[0]
    first, last = idx[0], idx[-1]
    if not covered[first:last + 1].all():
        raise ValueError('the profile has a gap between covered latitudes')
    out = profile.copy()
    rows = len(lat)
    for edge, side in ((first, -1), (last, 1)):
        if (side < 0 and edge == 0) or (side > 0 and edge == rows - 1):
            continue
        mirror_edge = rows - 1 - edge
        if not covered[mirror_edge]:
            raise ValueError('the opposite latitude is not covered either')
        ratio = profile[edge] / profile[mirror_edge]
        y = edge + side
        while 0 <= y < rows:
            fade = max(0.0, 1.0 - abs(lat[y] - lat[edge]) / FILL_BLEND_DEG)
            out[y] = profile[rows - 1 - y] * (1.0 + (ratio - 1.0) * fade)
            y += side
        print(f'    filled beyond {lat[edge]:+.1f} deg from the opposite latitudes (scale {np.round(ratio, 3)})')
    return out


def blur_mask(mask, radius):
    """Box-blurred mask (twice, about a triangle of the given half width)."""
    m = mask.astype(np.float64)
    for _ in range(2):
        for axis, mode in ((0, 'edge'), (1, 'wrap')):
            pad = [(0, 0), (0, 0)]
            pad[axis] = (radius, radius)
            c = np.cumsum(np.pad(m, pad, mode=mode), axis=axis)
            c = np.concatenate([np.zeros_like(np.take(c, [0], axis=axis)), c], axis=axis)
            n = m.shape[axis]
            m = (np.take(c, np.arange(2 * radius + 1, n + 2 * radius + 1), axis=axis) -
                 np.take(c, np.arange(n), axis=axis)) / (2 * radius + 1)
    return m


def make_planet(opal_dir, planet, width):
    a_km, b_km, detail_rotation, east_right, detail_gain = PLANETS[planet]
    rotations = load_rotations(opal_dir, planet)
    print(f'  {planet}: {len(rotations)} rotations ({", ".join(sorted(rotations))})')
    rows = next(iter(rotations.values())).shape[0]
    lat_g = 90.0 - np.arange(rows) * 180.0 / (rows - 1)  # planetographic, row centers

    per_rotation = np.stack([zonal_median(m, good_mask(m)) for m in rotations.values()])
    count = np.sum(np.all(np.isfinite(per_rotation), axis=-1), axis=0)
    profile = np.nanmedian(np.where(count[None, :, None] > 0, per_rotation, 0.0), axis=0)
    profile[count == 0] = np.nan
    covered = np.where(count > 0)[0]
    print(f'    profile covers {lat_g[covered[-1]]:+.1f} .. {lat_g[covered[0]]:+.1f} deg '
          f'(1 to {count.max()} rotations per latitude)')
    profile = fill_profile(profile, lat_g)
    profile = smooth_rows(profile, 2.0)

    detail = rotations[detail_rotation]
    if not east_right:
        detail = detail[:, ::-1]
    good = good_mask(detail)
    own = zonal_median(detail, good)
    with np.errstate(invalid='ignore', divide='ignore'):
        departure = np.clip(detail / own[:, None, :] - 1.0, *DETAIL_CLIP)
    weight = blur_mask(good & np.all(np.isfinite(own), axis=-1)[:, None],
                       int(round(DETAIL_FADE_DEG * detail.shape[1] / 360.0)))
    weight[~good] = 0.0
    departure = np.nan_to_num(departure) * weight[..., None]
    texture_g = profile[:, None, :] * (1.0 + detail_gain * departure)

    # Resample to planetocentric latitude and the output size.
    height = width // 2
    lat_c = 90.0 - (np.arange(height) + 0.5) * 180.0 / height
    lat_out_g = np.degrees(np.arctan((a_km / b_km) ** 2 * np.tan(np.radians(lat_c))))
    fy = (90.0 - lat_out_g) * (rows - 1) / 180.0
    src_w = texture_g.shape[1]
    fx = (np.arange(width) + 0.5) * src_w / width - 0.5
    y0 = np.clip(np.floor(fy).astype(int), 0, rows - 2)
    ty = (fy - y0)[:, None, None]
    x0 = np.floor(fx).astype(int)
    tx = (fx - x0)[None, :, None]
    x0 %= src_w
    x1 = (x0 + 1) % src_w
    top = texture_g[y0][:, x0] * (1 - tx) + texture_g[y0][:, x1] * tx
    bottom = texture_g[y0 + 1][:, x0] * (1 - tx) + texture_g[y0 + 1][:, x1] * tx
    linear = np.clip(top * (1 - ty) + bottom * ty, 0.0, 1.0)
    srgb = np.where(linear <= 0.0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - 0.055)
    eq = linear[height // 2].mean(axis=0)
    print(f'    equator I/F {np.round(eq, 3)}, sRGB #' + ''.join(f'{round(v * 255):02x}' for v in
          np.where(eq <= 0.0031308, 12.92 * eq, 1.055 * eq ** (1 / 2.4) - 0.055)))
    return Image.fromarray(np.round(srgb * 255).astype(np.uint8))


def main():
    if len(sys.argv) not in (3, 4):
        print(__doc__)
        sys.exit(1)
    opal_dir, out_dir = sys.argv[1], sys.argv[2]
    width = int(sys.argv[3]) if len(sys.argv) == 4 else 2048
    os.makedirs(out_dir, exist_ok=True)
    for planet in PLANETS:
        image = make_planet(opal_dir, planet, width)
        path = os.path.join(out_dir, planet + '.jpg')
        image.save(path, quality=90)
        print(f'  -> {path}')


if __name__ == '__main__':
    main()
