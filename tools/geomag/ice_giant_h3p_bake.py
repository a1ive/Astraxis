"""Bake the JWST H3+ maps of Uranus and Neptune into body-fixed brightness maps.

Usage (numpy; run with python -E -P):
    ice_giant_h3p_bake.py uranus <uranus_superim_fits_2.5.fits.zip> <results_lon_*.fits dir> <output.txt>
    ice_giant_h3p_bake.py neptune <figure_data dir> <west longitude offset, deg> <output.txt>

Both print what the scene needs besides the map. Sources, licenses and the
assumptions of the fits: assets/geomag/SOURCES.md.
"""
import math
from pathlib import Path
import sys
import zipfile

import numpy as np

BLOCK = 2880


def fits_images(raw):
    """The image HDUs of a FITS file: (EXTNAME or PRIMARY, array)."""
    pos, images = 0, []
    while pos < len(raw):
        cards, done = {}, False
        while not done and pos < len(raw):
            block = raw[pos:pos + BLOCK]
            pos += BLOCK
            for i in range(0, BLOCK, 80):
                card = block[i:i + 80].decode("ascii", "replace")
                if card.startswith("END"):
                    done = True
                    break
                if card[8:10] == "= ":
                    cards[card[:8].strip()] = card[10:].split("/")[0].strip().strip("'").strip()
        if not done:
            break
        bitpix = int(cards.get("BITPIX", 8))
        dims = [int(cards[f"NAXIS{k}"]) for k in range(1, int(cards.get("NAXIS", 0)) + 1)]
        count = math.prod(dims) if dims else 0
        size = count * abs(bitpix) // 8
        if count:
            dtype = {-64: ">f8", -32: ">f4", 8: "u1", 16: ">i2", 32: ">i4", 64: ">i8"}[bitpix]
            images.append((cards.get("EXTNAME", "PRIMARY"), np.frombuffer(raw[pos:pos + size], dtype=dtype).reshape(dims[::-1]).astype(float)))
        pos += math.ceil(size / BLOCK) * BLOCK
    return images


def write_map(path, values, value_range):
    """South-to-north rows of west longitude cell centers, in the scene's body frame."""
    height, width = values.shape
    text = f"AXAURORAMAP1 {width} {height} {value_range:g}\n"
    text += "".join(" ".join(f"{v:.3f}" for v in row) + "\n" for row in np.clip(values, 0.0, value_range))
    Path(path).write_text(text, encoding="ascii", newline="\n")


def smoothstep(a, b, x):
    t = np.clip((x - a) / (b - a), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def normalized_smooth(values, weights, sigma_cells):
    """Gaussian smoothing over the valid cells of a south-to-north, wrapping lat-lon
    map, by distance: a cell's width in longitude shrinks with cos(latitude)."""
    height, width = values.shape
    lat = (np.arange(height) + 0.5) / height * math.pi - math.pi / 2.0
    total = np.where(weights > 0.0, values, 0.0) * weights
    weight = weights.copy()
    # Along each row, in the frequency domain (the row wraps); wider near the pole.
    frequency = np.fft.rfftfreq(width)
    for j in range(height):
        sigma = min(sigma_cells / max(math.cos(lat[j]), 1e-3), width / 2.0)
        response = np.exp(-2.0 * (math.pi * frequency * sigma) ** 2)
        total[j] = np.fft.irfft(np.fft.rfft(total[j]) * response, width)
        weight[j] = np.fft.irfft(np.fft.rfft(weight[j]) * response, width)
    radius = int(math.ceil(3.0 * sigma_cells))
    kernel = np.exp(-0.5 * (np.arange(-radius, radius + 1) / sigma_cells) ** 2)
    kernel /= kernel.sum()
    along = lambda column: np.convolve(column, kernel, "same")
    total = np.apply_along_axis(along, 0, total)
    weight = np.apply_along_axis(along, 0, weight)
    return np.divide(total, weight, out=np.zeros_like(total), where=weight > 1e-6), weight


def uranus(zip_path, profile_dir, output):
    # Melin 2025 Figure 3a: total emission (slice 5, W m-2 sr-1) on 2.5 deg cells,
    # columns of ULS west longitude; the rows run south to north in the IAU sense
    # (the sunlit, observed pole is ULS south, IAU north).
    with zipfile.ZipFile(zip_path) as z:
        cube = fits_images(z.read("uranus_superim_fits_2.5.fits"))[0][1]
    emission = cube[4] * 1e6  # uW m-2 sr-1
    height, width = emission.shape
    lat = -90.0 + (np.arange(height) + 0.5) * 180.0 / height
    lon_uls = (np.arange(width) + 0.5) * 360.0 / width
    valid = np.isfinite(emission) & (emission != 0.0)

    # Tiranti 2026: emission profiles every 10 deg of ULS west longitude near the
    # limb, within 25 deg of the equator. Altitudes above 5025 km are unreliable.
    grid = 475.0 + 350.0 * np.arange(14)
    profiles = {}
    for f in Path(profile_dir).glob("results_lon_*.fits"):
        hdus = dict(fits_images(f.read_bytes()))
        order = np.argsort(hdus["ALTITUDES"])
        alt = hdus["ALTITUDES"][order]
        intensity = hdus["INTENSITY_PROFILES"][order]
        keep = np.isfinite(intensity)
        profiles[int(f.stem.split("_")[-1])] = np.interp(grid, alt[keep], intensity[keep])
    assert sorted(profiles) == list(range(5, 360, 10)), "expected 36 longitudes"
    stack = np.array([profiles[k] for k in sorted(profiles)])
    mean_profile = stack.mean(axis=0)
    column = stack.sum(axis=1)
    limb_lon = np.arange(5.0, 360.0, 10.0)
    limb = np.interp(lon_uls, limb_lon, np.convolve(np.r_[column[-1], column, column[0]], [0.25, 0.5, 0.25], "valid"), period=360.0)

    # The projection's cells near its edge fall off and are noisy: smooth over
    # ~ the spatial resolution, and trust the map only north of the equator.
    trust = smoothstep(-5.0, 10.0, lat)
    smooth, _ = normalized_smooth(emission, valid.astype(float) * trust[:, None], 1.3)

    # Toward 25 S, the limb's longitude curve, scaled to the map where both are
    # observed; south of it, no data and no emission.
    band = (lat >= 15.0) & (lat <= 35.0)
    scale = np.median(smooth[band]) / np.median(limb)
    pearson = np.corrcoef(smooth[band].mean(axis=0), limb)[0, 1]
    blend = smoothstep(2.0, 18.0, lat)[:, None]
    merged = blend * smooth + (1.0 - blend) * scale * limb[None, :]
    merged, _ = normalized_smooth(merged, np.ones_like(merged), 1.5)
    merged *= smoothstep(-35.0, -25.0, lat)[:, None]

    # To the scene's frame (IAU new, Lamy 2025): east longitude = ULS west - 50.6;
    # rows of west longitude cell centers, as for the other maps.
    east = lon_uls - 50.6
    west_cells = (np.arange(width) + 0.5) * 360.0 / width
    baked = np.zeros_like(merged)
    for j in range(height):
        baked[j] = np.interp(-west_cells, east, merged[j], period=360.0)
    write_map(output, baked, 1.0)

    # Zero one step below and above: AuroraVertical::Table.
    table = np.r_[0.0, mean_profile, 0.0]
    table /= np.trapezoid(table, np.r_[grid[0] - 350.0, grid, grid[-1] + 350.0])
    print("limb profile, per km, from 125 km every 350 km:")
    print(" ".join(f"{v:.4e}" for v in table))
    print(f"limb column max/min {column.max() / column.min():.2f}; correlation with the map at 15-35 N {pearson:.2f}")
    print(f"map median north of 5 N {np.median(baked[lat >= 5.0]):.3f}, max {baked.max():.3f} uW m-2 sr-1, clipped above 1: {(baked > 1.0).mean():.4f}")


def neptune(directory, west_offset, output):
    # Melin 2025 Figures 2 and 4: the H3+ radiance per spaxel, with the corners
    # (west longitude, planetocentric latitude) of the spaxels fully on the disk.
    # The spaxels are too noisy to show (signal-to-noise ~3), so fit what the
    # paper reports: a smooth background per view plus one localized enhancement.
    d = Path(directory)
    views = []
    for name, c, h, cloud, coords, poly in (
        ("1", 1, "figure2c_neptune_longitude1_h3p", "figure2a_neptune_longitude1_clouds", "figure4_coords_longitude1", "figure4_polygons_longitude1"),
        ("2", 2, "figure2d_neptune_longitude2_h3p", "figure2b_neptune_longitude2_clouds", "figure4_coords_longitude2", "figure4_polygons_longitude2"),
    ):
        h3p = fits_images((d / f"{h}.fits").read_bytes())[0][1] * 1e6
        clouds = fits_images((d / f"{cloud}.fits").read_bytes())[0][1]
        xy = fits_images((d / f"{coords}.fits").read_bytes())[0][1].astype(int)
        corners = fits_images((d / f"{poly}.fits").read_bytes())[0][1].copy()
        # Spaxels across the longitude seam were offset by -3600 deg for plotting.
        corners[..., 0] = np.mod(corners[..., 0], 360.0)
        lon_ref = corners[:, 0, 0][:, None]
        lon = lon_ref[:, 0] + np.mean(((corners[..., 0] - lon_ref) + 180.0) % 360.0 - 180.0, axis=1)
        lat = corners[..., 1].mean(axis=1)
        value = h3p[xy[:, 0], xy[:, 1]]
        keep = clouds[xy[:, 0], xy[:, 1]] < 20.0
        views.append((xy[keep].astype(float), np.mod(lon[keep], 360.0), lat[keep], value[keep]))
        print(f"view {name}: {keep.sum()} spaxels without clouds, median {np.median(value[keep]):.1f} uW m-2 sr-1 um-1")

    n = sum(len(v[3]) for v in views)
    def solve(lon0, lat0, sigma_lon, sigma_lat):
        # Linear in the two views' planes and the spot's amplitude.
        rows, values, offset = np.zeros((n, 7)), np.zeros(n), 0
        for k, (xy, lon, lat, value) in enumerate(views):
            m = len(value)
            rows[offset:offset + m, 3 * k] = 1.0
            rows[offset:offset + m, 3 * k + 1] = (xy[:, 0] - 21.0) / 15.0
            rows[offset:offset + m, 3 * k + 2] = (xy[:, 1] - 18.0) / 15.0
            dlon = ((lon - lon0) + 180.0) % 360.0 - 180.0
            rows[offset:offset + m, 6] = np.exp(-0.5 * ((dlon / sigma_lon) ** 2 + ((lat - lat0) / sigma_lat) ** 2))
            values[offset:offset + m] = value
            offset += m
        coeffs, *_ = np.linalg.lstsq(rows, values, rcond=None)
        return float(np.sum((rows @ coeffs - values) ** 2)), coeffs

    best = None
    for lon0 in np.arange(180.0, 301.0, 5.0):
        for lat0 in np.arange(-70.0, -19.0, 5.0):
            # The two views don't bound the width in longitude: at most the paper's 80 deg.
            for sigma_lon in (10.0, 15.0, 20.0, 25.0, 30.0, 34.0):
                for sigma_lat in (5.0, 7.5, 10.0, 12.5, 15.0, 20.0):
                    chi, coeffs = solve(lon0, lat0, sigma_lon, sigma_lat)
                    if coeffs[6] > 0.0 and (best is None or chi < best[0]):
                        best = (chi, lon0, lat0, sigma_lon, sigma_lat, coeffs)
    chi, lon0, lat0, sigma_lon, sigma_lat, coeffs = best
    # The background under the spot, in the view that shows it.
    center = np.mean(views[1][0], axis=0)
    background = coeffs[3] + coeffs[4] * (center[0] - 21.0) / 15.0 + coeffs[5] * (center[1] - 18.0) / 15.0
    ratio = 1.0 + coeffs[6] / background
    print(f"spot: {lon0:.0f} W, {lat0:.0f} deg, sigma {sigma_lon:.0f} x {sigma_lat:.0f} deg, amplitude {coeffs[6]:.1f} over background {background:.1f}: {ratio:.2f} x")
    print(f"residual rms {math.sqrt(chi / n):.1f} (off-disk noise ~5 uW m-2 sr-1 um-1 per smoothed spaxel)")

    # The map: relative to the background, in the scene's body frame.
    height, width = 72, 144
    lat = -90.0 + (np.arange(height) + 0.5) * 180.0 / height
    west = (np.arange(width) + 0.5) * 360.0 / width
    dlon = ((west[None, :] - (lon0 + west_offset)) + 180.0) % 360.0 - 180.0
    spot = np.exp(-0.5 * ((dlon / sigma_lon) ** 2 + ((lat[:, None] - lat0) / sigma_lat) ** 2))
    write_map(output, 1.0 + (ratio - 1.0) * spot, 4.0)


if __name__ == "__main__":
    if sys.argv[1] == "uranus":
        uranus(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "neptune":
        neptune(sys.argv[2], float(sys.argv[3]), sys.argv[4])
    else:
        sys.exit(__doc__)
