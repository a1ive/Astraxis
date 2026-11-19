"""Model sky seen from the centre of M87 (assets/sky/m87.jpg, assets/sky/m87_stars.csv).

From inside a giant elliptical the sky is its own starlight: M87's measured V-band
surface-brightness profile (Kormendy et al. 2009, Table 3) is deprojected (Abel) into
a luminosity density j(r) (L_sun/pc^3), assumed spherical on the circularized radius.
Seen from the centre that sky is the same in every direction: a uniform glow, plus
point stars drawn from an old stellar population. Everything else is fainter than
the glow by 7 magnitudes or more and is left out (Virgo's other galaxies, M87's
globular clusters, the Milky Way). The brightness is calibrated and compressed as
for the Milky Way skies (make_galaxy_sky.py, galaxy_calibration.json); see
assets/sky/SOURCES.md.

Usage (needs numpy and Pillow):
    python make_m87_sky.py [assets/scenes/m87.toml]
"""

import math
import os
import sys
import tomllib

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_galaxy_sky as g  # noqa: E402

DISTANCE_PC = 1.68e7     # [EHT19VI], as m87.toml (KFCB used 17.14 Mpc; only the angles are taken)
A_V = 0.072              # [KFCB] Table 1: Galactic extinction toward NGC 4486 (their profile is not corrected)
PC_PER_ARCSEC = DISTANCE_PC * math.pi / 648000.0
R_MIN_PC = 1.0           # stars and diffuse light from 1 pc out (~0.012", inside the innermost measured point)

# [KFCB] Kormendy, Fisher, Cornell & Bender 2009, ApJS 182, 216, Table 3: NGC 4486 = M 87
# composite major-axis surface photometry: radius (arcsec), mu_V (mag/arcsec^2, not
# corrected for Galactic extinction), ellipticity (None where the table gives none).
PROFILE = [
    (0.017, 16.266, 0.160), (0.044, 16.358, 0.161), (0.088, 16.511, 0.161), (0.176, 16.589, 0.161),
    (0.220, 16.646, 0.162), (0.264, 16.700, 0.162), (0.308, 16.746, 0.162), (0.352, 16.788, 0.138),
    (0.396, 16.838, 0.131), (0.440, 16.889, 0.118), (0.484, 16.927, 0.109), (0.548, 16.953, 0.097),
    (0.604, 16.966, 0.094), (0.660, 16.996, 0.084), (0.727, 17.031, 0.090), (0.795, 17.062, 0.087),
    (0.867, 17.091, 0.079), (0.950, 17.114, 0.075), (1.038, 17.134, 0.071), (1.147, 17.165, 0.072),
    (1.254, 17.195, 0.049), (1.365, 17.210, 0.030), (1.515, 17.241, 0.023), (1.669, 17.270, 0.015),
    (1.825, 17.290, 0.007), (1.998, 17.318, 0.018), (2.196, 17.346, 0.015), (2.419, 17.371, 0.012),
    (2.640, 17.399, 0.008), (2.835, 17.418, 0.005), (3.218, 17.470, 0.012), (3.823, 17.538, 0.010),
    (4.546, 17.613, 0.017), (5.413, 17.715, 0.021), (6.092, 17.790, 0.028), (7.118, 17.913, 0.023),
    (7.780, 17.991, 0.028), (8.610, 18.086, 0.020), (9.441, 18.183, 0.026), (10.304, 18.277, 0.030),
    (11.552, 18.409, 0.026), (12.322, 18.489, 0.030), (13.715, 18.622, 0.030), (15.109, 18.749, 0.032),
    (16.615, 18.879, 0.036), (18.249, 19.009, 0.036), (19.971, 19.143, 0.040), (21.945, 19.284, 0.043),
    (23.961, 19.415, None),
    (26.318, 19.560, 0.051), (29.040, 19.716, 0.052), (31.750, 19.860, 0.059), (35.015, 20.016, 0.064),
    (38.371, 20.161, 0.072), (42.073, 20.318, 0.076), (45.779, 20.455, 0.078), (50.855, 20.631, 0.082),
    (56.040, 20.794, 0.086), (61.094, 20.936, 0.096), (67.531, 21.097, 0.100), (72.277, 21.217, 0.099),
    (77.179, 21.331, 0.109), (84.918, 21.499, 0.114), (93.972, 21.693, 0.128), (104.954, 21.912, 0.139),
    (116.011, 22.116, 0.153), (127.938, 22.317, 0.157), (139.798, 22.515, 0.171), (154.170, 22.714, 0.185),
    (166.341, 22.870, 0.206), (180.926, 23.019, 0.222), (200.909, 23.220, 0.237), (222.587, 23.420, 0.254),
    (242.103, 23.573, 0.275), (265.053, 23.742, 0.293), (293.990, 23.934, 0.303), (321.366, 24.096, 0.313),
    (346.737, 24.257, 0.329), (381.651, 24.441, 0.337), (419.276, 24.658, 0.348), (462.914, 24.820, 0.370),
    (502.343, 25.011, 0.381), (541.377, 25.090, 0.388), (593.608, 25.288, 0.398), (653.131, 25.486, 0.427),
    (719.449, 25.697, 0.447), (794.328, 25.917, 0.454), (878.348, 26.100, 0.447), (946.237, 26.328, 0.457),
    (1046.325, 26.620, 0.464), (1145.513, 26.848, 0.454), (1230.269, 26.995, 0.443), (1336.595, 27.180, 0.439),
    (1479.109, 27.305, 0.436), (1621.810, 27.535, 0.433), (1778.279, 27.715, 0.429), (1995.262, 27.755, 0.422),
    (2443.700, 28.045, None),
]


def surface_profile():
    """Circularized radius (pc) and surface brightness (L_sun,V / pc^2), extinction-corrected.
    mu = M_V,sun + 21.572 - 2.5 lg I."""
    eps = np.array([e for _, _, e in PROFILE], dtype=float)
    known = ~np.isnan(np.array([np.nan if e is None else e for _, _, e in PROFILE], dtype=float))
    eps = np.interp(np.arange(len(PROFILE)), np.flatnonzero(known), eps[known].astype(float))
    r = np.array([p[0] for p in PROFILE]) * np.sqrt(1.0 - eps) * PC_PER_ARCSEC
    mu = np.array([p[1] for p in PROFILE]) - A_V
    return r, 10 ** (-0.4 * (mu - g.MV_SUN - 21.572))


def log_slopes(R, I):
    """Power-law pieces: I = I_k (R / R_k)^s_k on [R_k, R_k+1); the last slope continues
    outward, the first one inward."""
    s = np.diff(np.log(I)) / np.diff(np.log(R))
    return np.append(s, s[-1])


def surface_at(R, I, s, x):
    k = np.clip(np.searchsorted(R, x, side='right') - 1, 0, len(R) - 1)
    return I[k] * (x / R[k]) ** s[k], s[k]


def deproject(R, I, r):
    """j(r) = -1/pi int_r^inf (dI/dR) dR / sqrt(R^2 - r^2), with dI/dR from the power-law
    pieces, integrated in u = acosh(R / r) (no singularity at R = r)."""
    s = log_slopes(R, I)
    out = np.empty_like(r)
    u = np.linspace(0.0, 14.0, 6000)  # R up to ~6e5 r
    du = u[1] - u[0]
    for i, ri in enumerate(r):
        x = ri * np.cosh(u)
        Ix, sx = surface_at(R, I, s, x)
        # dI/dR dR / sqrt(R^2 - r^2) = (s I / R) (r sinh u du) / (r sinh u) = s I / R du
        out[i] = -float(np.sum(sx * Ix / x) * du) / math.pi
    return out


def old_luminosity_function(cal):
    """The solar neighbourhood's luminosity function (make_galaxy_sky) keeping the stars an
    old (~10 Gyr) population still has: redder than its turnoff (B - V >= 0.6) and fainter
    than the tip of its giant branches (M_V >= -2.5). Stars per kpc^3, templates, and the
    V emissivity (L_sun/pc^3) of that mix."""
    phi, templates, _ = g.luminosity_function(cal)
    phi_old, kept = [], []
    for p, t in zip(phi, templates):
        keep = (t[:, 1] >= 0.6) & (t[:, 0] >= -2.5) if len(t) else np.zeros(0, dtype=bool)
        phi_old.append(p * keep.mean() if len(t) else 0.0)
        kept.append(t[keep])
    phi_old = np.array(phi_old)
    j_old = float((phi_old / 1e9 * 10 ** (-0.4 * (g.MAG_MID - g.MV_SUN))).sum())
    return phi_old, kept, j_old


def light_weighted_bv(phi, templates):
    num = den = 0.0
    for p, t in zip(phi, templates):
        if len(t) == 0 or p <= 0:
            continue
        bv = t[:, 1][~np.isnan(t[:, 1])]
        L = p * 10 ** (-0.4 * t[:, 0][~np.isnan(t[:, 1])])
        num += float((L * bv).sum()) / len(bv)
        den += float(L.sum()) / len(bv)
    return num / den


def blackbody_srgb(kelvin):
    """Linear sRGB of a black body (max component 1): core/color.cpp blackbody_linear_srgb."""
    t = min(max(kelvin, 1667.0), 25000.0)
    if t <= 4000.0:
        x = -0.2661239e9 / t ** 3 - 0.2343589e6 / t ** 2 + 0.8776956e3 / t + 0.179910
    else:
        x = -3.0258469e9 / t ** 3 + 2.1070379e6 / t ** 2 + 0.2226347e3 / t + 0.240390
    if t <= 2222.0:
        y = -1.1063814 * x ** 3 - 1.34811020 * x ** 2 + 2.18555832 * x - 0.20219683
    elif t <= 4000.0:
        y = -0.9549476 * x ** 3 - 1.37418593 * x ** 2 + 2.09137015 * x - 0.16748867
    else:
        y = 3.0817580 * x ** 3 - 5.87338670 * x ** 2 + 3.75112997 * x - 0.37001483
    X, Y, Z = x / y, 1.0, (1.0 - x - y) / y
    rgb = np.array([3.2406 * X - 1.5372 * Y - 0.4986 * Z, -0.9689 * X + 1.8758 * Y + 0.0415 * Z,
                    0.0557 * X - 0.2040 * Y + 1.0570 * Z])
    rgb = np.maximum(rgb, 0.0)
    return rgb / rgb.max()


def temperature_from_bv(bv):
    """core/color.cpp temperature_from_bv (Ballesteros 2012)."""
    return 4600.0 * (1.0 / (0.92 * bv + 1.7) + 1.0 / (0.92 * bv + 0.62))


class Galaxy:
    """j(r) on a log grid from R_MIN_PC out, and its cumulative integrals."""

    def __init__(self):
        R, I = surface_profile()
        self.r = np.geomspace(R_MIN_PC, 3.0 * R[-1], 600)
        self.j = deproject(R, I, self.r)
        self.R, self.I = R, I

    def column(self, weight=None):
        """int j ds (L_sun/pc^2) from R_MIN_PC outward, times weight(s) if given."""
        w = self.j if weight is None else self.j * weight
        return float(np.trapezoid(w, self.r))

    def volume_cdf(self, r_max):
        """Cumulative int j 4 pi r^2 dr (L_sun) on the grid, up to r_max."""
        r = np.minimum(self.r, r_max)
        f = self.j * 4 * math.pi * self.r ** 2 * (self.r <= r_max)
        return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))])


def expected_counts(galaxy, phi, j_old, point_mag):
    """Point stars per magnitude bin brighter than point_mag: phi scaled by j / j_old over
    the volume within which each bin is brighter than the limit."""
    counts = np.zeros(len(phi))
    for b, lo in enumerate(g.MAG_EDGES[:-1]):
        if phi[b] <= 0:
            continue
        D = 10 ** ((point_mag - lo) / 5 + 1)  # pc
        if D <= R_MIN_PC:
            continue
        light = galaxy.volume_cdf(D)[-1]  # L_sun
        counts[b] = phi[b] / 1e9 * light / j_old  # phi per pc^3 at emissivity j_old
    return counts


def sample_stars(galaxy, phi, templates, j_old, point_mag, rng):
    counts = expected_counts(galaxy, phi, j_old, point_mag)
    dirs, V, bv = [], [], []
    for b, lo in enumerate(g.MAG_EDGES[:-1]):
        n = rng.poisson(counts[b])
        if n == 0 or len(templates[b]) == 0:
            continue
        D = 10 ** ((point_mag - lo) / 5 + 1)
        cdf = galaxy.volume_cdf(D)
        r = np.interp(rng.uniform(0, cdf[-1], n), cdf, np.minimum(galaxy.r, D))
        u = rng.normal(size=(n, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        t = templates[b][rng.integers(len(templates[b]), size=n)]
        dirs.append(u)
        V.append(t[:, 0] + 5 * np.log10(r / 10.0))
        bv.append(t[:, 1])
    dirs, V, bv = np.concatenate(dirs), np.concatenate(V), np.concatenate(bv)
    keep = V <= point_mag
    return dirs[keep], V[keep], bv[keep]


def choose_point_mag(galaxy, phi, j_old):
    if expected_counts(galaxy, phi, j_old, g.NAKED_EYE_MAG).sum() <= g.MAX_POINT_STARS:
        return g.NAKED_EYE_MAG
    lo, hi = -15.0, g.NAKED_EYE_MAG
    for _ in range(20):
        mid = 0.5 * (lo + hi)
        if expected_counts(galaxy, phi, j_old, mid).sum() > g.MAX_POINT_STARS:
            hi = mid
        else:
            lo = mid
    return round(lo, 2)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(g.ASSETS, 'scenes', 'm87.toml')
    with open(path, 'rb') as f:
        scene = tomllib.load(f)
    sky = scene['sky']
    viewer = sky['viewer_ra_dec_distance_pc']
    cal = g.load_calibration()
    rng = np.random.default_rng(g.NOISE_SEED)

    galaxy = Galaxy()
    phi, templates, j_old = old_luminosity_function(cal)
    point_mag = choose_point_mag(galaxy, phi, j_old)
    resolved = g.RESOLVED_MAG if point_mag >= g.NAKED_EYE_MAG else point_mag
    # The diffuse light: stars fainter than the point-star limit at each distance.
    f_unres = g.unresolved_fraction(phi, galaxy.r / 1000.0, resolved)
    column = galaxy.column(f_unres)  # L_sun/pc^2 = L_sun/pc^3 x pc
    bv = light_weighted_bv(phi, templates)
    color = blackbody_srgb(temperature_from_bv(bv))
    # Map units as the Milky Way skies: light (L_sun/pc^3 x kpc) x map_per_light, its luminance
    # kept, its color that of the old population.
    lum_per_light = float(np.array(cal['map_per_light']) @ np.array([0.2126, 0.7152, 0.0722]))
    level = column / 1000.0 * lum_per_light / float(color @ np.array([0.2126, 0.7152, 0.0722]))
    width, height = 512, 256
    model = np.broadcast_to(color * level, (height, width, 3)).copy()

    dirs, V, star_bv = sample_stars(galaxy, phi, templates, j_old, point_mag, rng)
    white = float(model.max())
    brightness, mag_offset, ratio = g.brightness_for(model, cal)
    g.write_map(os.path.join(g.ASSETS, sky['milky_way']), model, white)
    header = [
        f'Model sky for {scene.get("name", path)}: stars brighter than V = {point_mag} seen from the centre of M87.',
        'Generated by tools/sky/make_m87_sky.py (see assets/sky/SOURCES.md); statistical, not real stars.',
        f'viewer_ra_dec_distance_pc = {viewer[0]}, {viewer[1]}, {viewer[2]}',
        f'milky_way_brightness = {brightness * white:.4f}',
        f'display_mag_offset = {mag_offset:.3f}',
    ]
    g.write_stars(os.path.join(g.ASSETS, sky['stars']), dirs, V, star_bv, header)

    j1 = float(np.interp(10.0, galaxy.r, galaxy.j))
    j100 = float(np.interp(100.0, galaxy.r, galaxy.j))
    mu_sky = g.MV_SUN + 21.572 - 2.5 * math.log10(galaxy.column())
    print(f'j(10 pc) = {j1:.4g}, j(100 pc) = {j100:.4g} L_sun/pc^3; total column {galaxy.column():.4g} L_sun/pc^2 '
          f'(sky {mu_sky:.2f} V mag/arcsec^2), diffuse {column:.4g}')
    print(f'old population: B-V = {bv:.3f} (light-weighted), j_old = {j_old:.4g} L_sun/pc^3 at the Sun\'s density')
    print(f'{len(V)} stars to V = {point_mag} (brightest V = {V.min():.2f}); mean sky {ratio:.4g} x the Sun\'s; '
          f'milky_way_brightness = {brightness * white:.4f}; stars shown {mag_offset:.2f} mag fainter')


if __name__ == '__main__':
    main()
