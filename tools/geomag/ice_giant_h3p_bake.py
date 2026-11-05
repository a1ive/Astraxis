"""Bake the JWST H3+ maps of Uranus and Neptune into body-fixed brightness maps.

Usage (numpy; run with python -E -P):
    ice_giant_h3p_bake.py uranus <uranus_superim_fits_2.5.fits.zip> <results_lon_*.fits dir> <L_shells.txt> <dip_equator.txt> <output.txt>
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


def write_map(path, values, value_range, dark=None):
    """South-to-north rows of west longitude cell centers, in the scene's body frame;
    a second layer after the first: what is left without sunlight."""
    height, width = values.shape
    text = f"AXAURORAMAP{1 if dark is None else 2} {width} {height} {value_range:g}\n"
    for layer in [values] if dark is None else [values, dark]:
        text += "".join(" ".join(f"{v:.3f}" for v in row) + "\n" for row in np.clip(layer, 0.0, value_range))
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


def read_table(path):
    rows = [line.split() for line in Path(path).read_text(encoding="ascii").splitlines() if line.strip() and not line.startswith("#")]
    return np.array(rows, dtype=float)


def unit_vectors(latitude_deg, west_deg):
    lat, lon = np.radians(latitude_deg), np.radians(west_deg)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


def uranus_magnetic_latitude(shells_path, equator_path, uls_latitude_deg, uls_west_deg):
    """Signed invariant latitude asin(sqrt(1 - 1/L)) in deg at the given points, from the
    contours of Melin 2025's tables for Q3: + around ULS north, 0 on the dip equator. A
    thin-plate spline on the sphere through the contours' points. Q3, not the scene's AH5:
    Melin finds the dark chevron lines up with the dip equator of Q3 and UH5, more poorly
    with AH5's, and the map agrees (darker by 3.9%, 2.9%, 0.9% on the equator than 10-15
    deg off it)."""
    q3 = 1.0
    shells = read_table(shells_path)
    shells = shells[shells[:, 0] == q3]
    equator = read_table(equator_path)
    equator = equator[equator[:, 0] == q3]

    # The file's labels are in the wrong order: their notebook contours Rmax at the
    # ascending levels [3, 5, 10, 20] and names them ['20', '10', '5', '3']. Tracing Q3
    # from the contours gives Rmax ~ 20, 10, 5, 2.9 for the labels 3, 5, 10, 20.
    true_l = {3.0: 20.0, 5.0: 10.0, 10.0: 5.0, 20.0: 3.0}
    spans = [np.ptp(shells[(shells[:, 1] == label) & (shells[:, 2] == 2.0), 4]) for label in sorted(true_l)]
    assert spans == sorted(spans), "the label order changed: 'L = 3' should be the innermost contour"

    points, values = [unit_vectors(equator[:, 3], equator[:, 2])], [np.zeros(len(equator))]
    for label, l_shell in true_l.items():
        for pole, sign in ((1.0, 1.0), (2.0, -1.0)):
            contour = shells[(shells[:, 1] == label) & (shells[:, 2] == pole)]
            points.append(unit_vectors(contour[:, 4], contour[:, 3]))
            values.append(np.full(len(contour), sign * math.sqrt(1.0 - 1.0 / l_shell)))
    points, values = np.concatenate(points), np.concatenate(values)

    # Nodes at least 4 deg apart along each contour: 1 deg steps only make the system
    # ill-conditioned. Greedy, so that the result does not depend on a random draw.
    nodes = []
    limit = math.cos(math.radians(4.0))
    for p, v in zip(points, values):
        if not any(v == w and p @ q > limit for q, w in nodes):
            nodes.append((p, v))
    node_points = np.array([n[0] for n in nodes])
    node_values = np.array([n[1] for n in nodes])

    def kernel(a, b):
        d2 = np.maximum(2.0 - 2.0 * a @ b.T, 1e-12)
        return 0.5 * d2 * np.log(d2)

    n = len(nodes)
    system = np.zeros((n + 4, n + 4))
    system[:n, :n] = kernel(node_points, node_points) + 1e-5 * np.eye(n)
    system[:n, n] = 1.0
    system[:n, n + 1:] = node_points
    system[n:, :n] = system[:n, n:].T
    weights = np.linalg.solve(system, np.r_[node_values, np.zeros(4)])
    query = unit_vectors(uls_latitude_deg, uls_west_deg)
    shape = query.shape[:-1]
    query = query.reshape(-1, 3)
    coordinate = kernel(query, node_points) @ weights[:n] + weights[n] + query @ weights[n + 1:]
    return np.degrees(np.arcsin(np.clip(coordinate.reshape(shape), -1.0, 1.0)))


def interp_over_nan(x):
    index = np.arange(len(x))
    ok = np.isfinite(x)
    return np.interp(index, index[ok], x[ok])


def uranus_unobserved_model(smooth, weights, observed_share, lat, magnetic_lat):
    """log emission = constant + a(|latitude|) + b(magnetic latitude), fitted to the
    observed cells: the extrapolation for where nothing was observed. Held out in
    sectors of the observed map, a(|latitude|) does most of the work and b adds little
    on average; b carries the dip equator's dark band and the auroral shells.

    Also what particle precipitation alone leaves, without sunlight: on the auroral
    shells (Melin 2025: L < 3 is non-auroral), what the observed map holds above the
    non-auroral background, or where it was not observed, the shells' mean excess."""
    log_map = np.log(np.maximum(smooth, 1e-3))
    # Melin 2025: L < 3 is non-auroral.
    l3_deg = math.degrees(math.asin(math.sqrt(1.0 - 1.0 / 3.0)))
    l2_deg = math.degrees(math.asin(math.sqrt(1.0 - 1.0 / 2.0)))
    # 7.5 deg of latitude: the dark caps' width; 5 deg of magnetic latitude: the band's.
    lat_edges = np.arange(0.0, 90.1, 7.5)
    mag_edges = np.arange(-90.0, 90.1, 5.0)
    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    mag_centers = 0.5 * (mag_edges[:-1] + mag_edges[1:])
    lat_cells = np.broadcast_to(lat[:, None], smooth.shape)
    lat_index = np.minimum(np.digitize(np.abs(lat_cells), lat_edges) - 1, len(lat_centers) - 1)
    mag_index = np.minimum(np.digitize(magnetic_lat, mag_edges) - 1, len(mag_centers) - 1)

    # Only cells inside the map's trusted latitudes count toward the profiles.
    use = weights > 0.5
    constant = log_map[use].mean()
    a, b = np.zeros(len(lat_centers)), np.zeros(len(mag_centers))
    min_cells = 5
    for _ in range(8):
        for k in range(len(a)):
            cells = use & (lat_index == k)
            a[k] = (log_map[cells] - constant - b[mag_index[cells]]).mean() if cells.sum() >= min_cells else np.nan
        a = interp_over_nan(a)
        for k in range(len(b)):
            cells = use & (mag_index == k)
            b[k] = (log_map[cells] - constant - a[lat_index[cells]]).mean() if cells.sum() >= min_cells else np.nan
        # The inner shells of the ULS north pole lie entirely in the unobserved part: the
        # same shells around the other pole, for lack of anything better. Not the same
        # field lines: Melin finds the two aurorae not conjugate, and the map's brightness
        # does not correlate between the feet of a line.
        missing = ~np.isfinite(b)
        b[missing] = b[::-1][missing]
        b = interp_over_nan(b)

    # Only the shape of the band: a fitted bin's noise is not a magnetic structure.
    smoothed_b = np.convolve(np.r_[b[0], b, b[-1]], [0.25, 0.5, 0.25], "valid")
    shell_term = np.interp(magnetic_lat, mag_centers, smoothed_b)
    fill = constant + np.interp(np.abs(lat_cells), lat_centers, a) + shell_term
    non_auroral = np.exp(shell_term[use & (np.abs(magnetic_lat) < l3_deg)]).mean()
    on_shells = smoothstep(l2_deg, l3_deg, np.abs(magnetic_lat))
    shell_mean = np.exp(constant) * np.maximum(np.exp(shell_term) - non_auroral, 0.0)
    background = np.exp(constant + np.interp(np.abs(lat_cells), lat_centers, a)) * non_auroral
    observed_excess = np.maximum(smooth - background, 0.0)
    dark = on_shells * (observed_share * observed_excess + (1.0 - observed_share) * shell_mean)
    print(f"unobserved model: a(|lat|) {np.exp(a.min()):.2f} .. {np.exp(a.max()):.2f}, b(magnetic lat) {np.exp(smoothed_b.min()):.2f} .. {np.exp(smoothed_b.max()):.2f}, "
          f"non-auroral b {non_auroral:.3f}, aurora without sunlight: observed up to {(on_shells * observed_share * observed_excess).max():.3f}, "
          f"shell mean up to {(on_shells * shell_mean).max():.3f} uW m-2 sr-1")
    return np.exp(fill), dark


def uranus(zip_path, profile_dir, shells_path, equator_path, output):
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
    # observed; south of it, nothing was observed: a model of the magnetic shells.
    band = (lat >= 15.0) & (lat <= 35.0)
    scale = np.median(smooth[band]) / np.median(limb)
    pearson = np.corrcoef(smooth[band].mean(axis=0), limb)[0, 1]
    blend = smoothstep(2.0, 18.0, lat)[:, None]
    merged = blend * smooth + (1.0 - blend) * scale * limb[None, :]
    # ULS latitude is minus the rows' IAU latitude.
    west_grid, lat_grid = np.meshgrid(lon_uls, lat)
    magnetic_lat = uranus_magnetic_latitude(shells_path, equator_path, -lat_grid, west_grid)
    # The observed cells that the projection's edge does not reach: see trust above.
    observed = valid & (lat_grid >= 8.0)
    fill, dark = uranus_unobserved_model(smooth, observed.astype(float), blend, lat, magnetic_lat)
    observed_part = smoothstep(-35.0, -25.0, lat)[:, None]
    merged = observed_part * merged + (1.0 - observed_part) * fill
    merged, _ = normalized_smooth(merged, np.ones_like(merged), 1.5)
    # All of the above is as if sunlit: the limb's curve gives no latitude dependence and
    # was scaled to the sunlit map, the fill mirrors the sunlit hemisphere. The scene mixes
    # it with the dark layer by the sun of the moment (AuroraPattern::BodyMapSeasons).
    dark, _ = normalized_smooth(dark, np.ones_like(dark), 1.5)

    # To the scene's frame (IAU new, Lamy 2025): east longitude = ULS west - 50.6;
    # rows of west longitude cell centers, as for the other maps.
    east = lon_uls - 50.6
    west_cells = (np.arange(width) + 0.5) * 360.0 / width
    baked = np.zeros_like(merged)
    baked_dark = np.zeros_like(dark)
    for j in range(height):
        baked[j] = np.interp(-west_cells, east, merged[j], period=360.0)
        baked_dark[j] = np.interp(-west_cells, east, dark[j], period=360.0)
    write_map(output, baked, 1.0, baked_dark)

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
    # Melin 2025, Methods: the central meridians, deg W. The sun is within 2 deg of JWST.
    for name, central_west, h, cloud, coords, poly in (
        ("1", 131.0, "figure2c_neptune_longitude1_h3p", "figure2a_neptune_longitude1_clouds", "figure4_coords_longitude1", "figure4_polygons_longitude1"),
        ("2", 303.2, "figure2d_neptune_longitude2_h3p", "figure2b_neptune_longitude2_clouds", "figure4_coords_longitude2", "figure4_polygons_longitude2"),
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
        # From -1 at the dawn limb to +1 at the dusk limb: east of noon is afternoon.
        hour_angle = np.radians(central_west - lon)
        dusk = np.cos(np.radians(lat)) * np.sin(hour_angle)
        views.append((xy[keep].astype(float), np.mod(lon[keep], 360.0), lat[keep], value[keep], dusk[keep]))
        print(f"view {name}: {keep.sum()} spaxels without clouds, median {np.median(value[keep]):.1f} uW m-2 sr-1 um-1")

    n = sum(len(v[3]) for v in views)
    def solve(lon0, lat0, sigma_lon, sigma_lat):
        # Linear in the two views' planes and the spot's amplitude.
        rows, values, offset = np.zeros((n, 7)), np.zeros(n), 0
        for k, (xy, lon, lat, value, _) in enumerate(views):
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

    # The paper sees a rise from dawn to dusk and gives no number: a background per view,
    # one gradient and the spot, now in the sun's frame instead of the views' planes.
    rows, values, offset = np.zeros((n, 4)), np.zeros(n), 0
    for k, (_, lon, lat, value, dusk) in enumerate(views):
        m = len(value)
        rows[offset:offset + m, k] = 1.0
        rows[offset:offset + m, 2] = dusk
        dlon = ((lon - lon0) + 180.0) % 360.0 - 180.0
        rows[offset:offset + m, 3] = np.exp(-0.5 * ((dlon / sigma_lon) ** 2 + ((lat - lat0) / sigma_lat) ** 2))
        values[offset:offset + m] = value
        offset += m
    gradient, *_ = np.linalg.lstsq(rows, values, rcond=None)
    mean_background = 0.5 * (gradient[0] + gradient[1])
    print(f"dawn to dusk: gradient / background = {gradient[2] / mean_background:.3f} (aurora_model.cpp), backgrounds {gradient[0]:.1f}, {gradient[1]:.1f}; "
          f"spot over it {1.0 + gradient[3] / mean_background:.2f} x; correlations with dusk "
          + ", ".join(f"{np.corrcoef(v[3], v[4])[0, 1]:.2f}" for v in views))

    # The map: relative to the background, in the scene's body frame.
    height, width = 72, 144
    lat = -90.0 + (np.arange(height) + 0.5) * 180.0 / height
    west = (np.arange(width) + 0.5) * 360.0 / width
    dlon = ((west[None, :] - (lon0 + west_offset)) + 180.0) % 360.0 - 180.0
    spot = np.exp(-0.5 * ((dlon / sigma_lon) ** 2 + ((lat[:, None] - lat0) / sigma_lat) ** 2))
    write_map(output, 1.0 + (ratio - 1.0) * spot, 4.0)


if __name__ == "__main__":
    if sys.argv[1] == "uranus":
        uranus(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "neptune":
        neptune(sys.argv[2], float(sys.argv[3]), sys.argv[4])
    else:
        sys.exit(__doc__)
