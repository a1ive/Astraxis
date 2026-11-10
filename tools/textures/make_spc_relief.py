"""Bake relief from the Cassini SPC radius GeoTIFFs of Weirich et al. (2025).

    python -E -P tools/textures/make_spc_relief.py <moon> <data_dir> <out.jpg>

Use mimas_radius_g.tif or rhea_radius_g.tif, or the three files
<moon>_{eq,np,sp}radius_g.tif for Tethys and Dione. Pass the bundle's data
directory for Mimas and Rhea, or data/global for Tethys and Dione.
Requires numpy, Pillow and tifffile.
"""

from pathlib import Path
import sys

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parent))
import relief_map


WIDTH = 2048
# Mean radius used in the scene, and the reference ellipsoid (Thomas 2010, table 1).
MOONS = {
    "mimas": (198.2, (207.8, 196.7, 190.6)),
    "tethys": (531.1, (538.4, 528.3, 526.3)),
    "dione": (561.4, (563.4, 561.3, 559.6)),
    "rhea": (763.5, (765.0, 763.1, 762.4)),
}
POLAR_TILES = {"tethys": 531000.0, "dione": 561400.0}


def sample_radius(path, x, y, shape, period_x=None):
    """Bilinear samples in projected metres, returned as weighted radii and coverage."""
    with tifffile.TiffFile(path) as tif:
        page = tif.pages[0]
        if page.shape != shape or page.dtype.kind != "f":
            sys.exit(f"{path}: expected a float radius grid of shape {shape}")
        grid = page.asarray()
        dx, dy, _ = page.tags["ModelPixelScaleTag"].value
        i, j, _, x0, y0, _ = page.tags["ModelTiepointTag"].value

    rows, cols = shape
    u = (x - x0) / dx + i - 0.5
    v = np.clip((y0 - y) / dy + j - 0.5, 0, rows - 1)
    if period_x is None:
        u = np.clip(u, 0, cols - 1)
        left = np.floor(u).astype(int)
        right = np.minimum(left + 1, cols - 1)
        du = u - left
    else:
        # Some grids fall slightly short of 360 degrees. Interpolate across the
        # entire gap between the last and first centres, measured in source pixels.
        period = max(cols, period_x / dx)
        u %= period
        left = np.floor(np.minimum(u, cols - 1)).astype(int)
        right = (left + 1) % cols
        du = (u - left) / np.where(right == 0, period - (cols - 1), 1.0)
    top = np.floor(v).astype(int)
    bottom = np.minimum(top + 1, rows - 1)
    dv = v - top
    radius = np.zeros(np.broadcast_shapes(u.shape, v.shape))
    coverage = np.zeros_like(radius)
    for row, col, weight in ((top, left, (1 - du) * (1 - dv)), (top, right, du * (1 - dv)),
                             (bottom, left, (1 - du) * dv), (bottom, right, du * dv)):
        samples = grid[row, col]
        # The PDS missing constant is negative. Physical radii are positive.
        valid = (samples > 0.0) & np.isfinite(samples)
        radius += np.where(valid, samples, 0.0) * weight
        coverage += weight * valid
    return radius, coverage


def bake(name, source, out):
    radius_km, (a, b, c) = MOONS[name]
    rows = WIDTH // 2
    lat = (0.5 - (np.arange(rows)[:, None] + 0.5) / rows) * np.pi
    lon = ((np.arange(WIDTH)[None, :] + 0.5) / WIDTH - 0.5) * (2.0 * np.pi)

    polar_radius = POLAR_TILES.get(name)
    projection_radius = polar_radius or radius_km * 1000.0
    eq_path = source / (f"{name}_eqradius_g.tif" if polar_radius else f"{name}_radius_g.tif")
    # The equatorial product has its central meridian at 180 E. The albedo starts at 180 W.
    x = projection_radius * (lon % (2.0 * np.pi) - np.pi)
    radius, coverage = sample_radius(eq_path, x, projection_radius * lat,
                                     (679 if polar_radius else 1111, 2222), period_x=2.0 * np.pi * projection_radius)
    if polar_radius:
        # The equatorial map ends at 55 deg and the polar squares reach about 55 deg
        # along their axes. A narrow blend uses the last row or column at each edge.
        lat_deg = np.degrees(lat[:, 0])
        polar_weight = np.clip((np.abs(lat_deg) - 54.8) / 0.7, 0.0, 1.0)
        for pole, prefix in ((1, "np"), (-1, "sp")):
            mask = (pole * lat_deg > 0.0) & (polar_weight > 0.0)
            rho = 2.0 * polar_radius * np.tan((np.pi / 2.0 - pole * lat[mask]) / 2.0)
            cap, cap_coverage = sample_radius(source / f"{name}_{prefix}radius_g.tif",
                                              rho * np.sin(lon), -pole * rho * np.cos(lon), (444, 444))
            weight = polar_weight[mask, None]
            radius[mask] = (1 - weight) * radius[mask] + weight * cap
            coverage[mask] = (1 - weight) * coverage[mask] + weight * cap_coverage

    radius_km_map = np.divide(radius * 1e-3, coverage, out=np.zeros_like(radius), where=coverage > 0.0)
    reference = 1.0 / np.sqrt((np.cos(lat) * np.cos(lon) / a) ** 2 +
                              (np.cos(lat) * np.sin(lon) / b) ** 2 + (np.sin(lat) / c) ** 2)
    height = np.where(coverage > 0.0, radius_km_map - reference, 0.0)
    print(f"{name}: {coverage.mean():.2%} valid radius coverage")
    relief_map.write(out, height, radius_km, coverage=coverage)


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in MOONS:
        sys.exit(__doc__)
    bake(sys.argv[1], Path(sys.argv[2]), sys.argv[3])
