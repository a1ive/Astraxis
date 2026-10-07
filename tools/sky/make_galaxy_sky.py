"""Model Milky Way skies for the scenes far from the Sun (assets/sky/).

From the Sun the sky is real (assets/stars/hyg.csv, assets/textures/milky_way.jpg).
From a star 0.5-10 kpc away it is modelled here: a smooth Galaxy (stellar disks,
bulge, dust disk) with statistical dust clumps, seen from the scene's viewer, plus
point stars drawn from the solar neighbourhood's luminosity function. The model,
its sources and the calibration are described in assets/sky/SOURCES.md.

Usage (needs numpy and Pillow):
    python make_galaxy_sky.py calibrate [--out DIR]
        Fits the model's emission and dust to milky_way.jpg as seen from the Sun and
        writes tools/sky/galaxy_calibration.json; DIR gets comparison images.
    python make_galaxy_sky.py sun [--out DIR]
        Validation: the model sky from the Sun (star counts; a map in DIR).
    python make_galaxy_sky.py render <scene.toml> [...]
        For each scene whose [sky] has viewer_ra_dec_distance_pc, milky_way and stars,
        writes the map and the star list (paths relative to assets/) and prints the
        milky_way_brightness to put in the scene file.
"""

import json
import math
import multiprocessing as mp
import os
import sys
import tomllib

import numpy as np
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ASSETS = os.path.join(ROOT, 'assets')
CALIBRATION = os.path.join(ROOT, 'tools', 'sky', 'galaxy_calibration.json')

# ---------------------------------------------------------------- geometry
# Galactic frame (ESA SP-1200, vol. 1, sec. 1.5.3), as core/math.hpp galactic_to_icrf.
NGP_RA, NGP_DEC, L_NCP = 192.85948, 27.12825, 122.93192
R0 = 8.2467       # kpc, Sun - Galactic centre [GR20] (as sgr_a.toml)
DS = R0 / 8.0     # [DS01] lengths are for R0 = 8 kpc: scaled by this
Z_SUN = 0.0146 * DS  # kpc, the Sun above the plane [DS01] Table 2

# Stellar thin disk [DS01] Table 2: exp(-R/r) sech^2(z/h), cut off beyond r_c by
# exp(-(R - r_c) / (r / 5)) (eq. 27). [DS01] fit the disk only beyond 0.35 R0
# (the bulge region was excluded); inside, it is tapered here like their dust
# hole (Gaussian of width 0.25 R0), leaving the centre to the bulge.
THIN_R, THIN_H, THIN_RC = 2.264 * DS, 0.2822 * DS, 10.52 * DS
DISK_HOLE = 0.35 * R0
# Thick disk [BHG16] sec. 5.1.3 / 5.2.2: 4% of the thin disk's local density,
# exponential, z_T = 900 pc, R_T = 2.0 kpc.
THICK_F, THICK_Z, THICK_R = 0.04, 0.900, 2.0
# Bulge/bar [RO12] Table 2, model S+E (the best two-ellipsoid fit to 2MASS star
# counts): a boxy sech^2 ellipsoid and an exponential one, both with the major
# axis 12.9 deg from the Sun - centre line, near end at positive longitude:
# density norm * f(Rs) * cutoff, Rs^c_par = [(|x|/x0)^c_perp + (|y|/y0)^c_perp]^(c_par/c_perp)
# + (|z|/z0)^c_par, cutoff exp(-((R_xy - R_c) / 0.5 kpc)^2) beyond R_c (sec. 2.2).
# (Red-clump studies find a larger angle, ~27 deg [BHG16]; [RO12] is used as a
# consistent set.) Columns: f, x0, y0, z0, norm, R_c, c_par, c_perp.
BAR_ANGLE = math.radians(12.9)
BULGE = [
    ('sech2', 1.46, 0.49, 0.39, 35.45, 3.43, 3.007, 3.329),
    ('exp', 4.44, 1.31, 0.80, 2.27, 6.83, 2.786, 3.917),
]
# Dust disk [DS01] Table 1: rho0 exp(-R/h_r) sech^2(z/h_d), h_d flaring linearly
# beyond r_f; inside 0.5 R0 a Gaussian hole (sec. 3.2). Density in MJy/sr/kpc,
# V opacity 0.0180 (MJy/sr)^-1 (Table 2): tau_V per kpc. The spiral arms and the
# local arm of [DS01] are left out; DUST_SCALE (calibrated) absorbs them on average.
DUST_RHO0, DUST_HR = 1098.0, 2.26 * DS
DUST_H0, DUST_H1, DUST_RF = 0.1344 * DS, 0.0148, 4.40 * DS
DUST_KAPPA_V = 0.0180
# Extinction in the map's blue, green and red channels relative to V: [CCM89]
# with R_V = 3.1 at 440, 550 and 640 nm (assumed channel wavelengths).
R_V = 3.1


def ccm89(wavelength_um):
    y = 1.0 / wavelength_um - 1.82
    a = np.polyval((0.32999, -0.7753, 0.01979, 0.72085, -0.02427, -0.50447, 0.17699, 1.0), y)
    b = np.polyval((-2.09002, 5.3026, -0.62251, -5.38434, 1.07233, 2.28305, 1.41338, 0.0), y)
    return a + b / R_V


EXT_RGB = np.array([ccm89(0.640), ccm89(0.550), ccm89(0.440)])

# Statistical dust clumps (artistic, see SOURCES.md): the dust density times
# exp(SIGMA g) / E[exp(SIGMA g)], g a sum of ridged gradient-noise octaves in
# Galactocentric space (cell size in kpc, amplitude ~ sqrt(cell)). Octaves finer
# than FOOTPRINT (an angle: the final maps' pixel) fade out with distance, so far
# dust is as smooth as a pixel sees it - the same 3D field at any render size.
SIGMA = 2.0
OCTAVES = [(0.32, 1.0), (0.16, 0.71), (0.08, 0.5), (0.04, 0.35), (0.02, 0.25), (0.01, 0.18)]
FOOTPRINT = math.pi / 1024
NOISE_SEED = 20261007

# Light of stars fainter than this (from the viewer) is in the diffuse map: the
# SVS map leaves out the Hipparcos/Tycho stars (assets/textures/SOURCES.md).
RESOLVED_MAG = 11.5
NAKED_EYE_MAG = 6.5         # point stars, as scene/star_catalog.hpp kNakedEyeMag
MAG_GAMMA = 0.56            # star flux compression, render/starfield_pass.cpp kMagnitudeGamma
DEFAULT_BRIGHTNESS = 0.25   # scene/scene.hpp SceneSky::milky_way_brightness (Earth's map)

# Galaxies outside ours, carried over from the real map ([SIMBAD] centre and
# size; distances [PI19], [GR20s], [LI21]): (ra, dec, radius deg, distance kpc).
EXTERNAL = [
    (80.8942, -69.7561, 6.0, 49.59),   # LMC
    (13.1583, -72.8003, 3.0, 62.44),   # SMC
    (10.6847, 41.2688, 2.2, 761.0),    # M31
]

# Ray marching: log-spaced samples from 1 pc to 40 kpc.
STEPS = 192
S_EDGES = np.geomspace(0.001, 40.0, STEPS + 1)
S_MID = np.sqrt(S_EDGES[1:] * S_EDGES[:-1])
S_DS = np.diff(S_EDGES)


def unit_ra_dec(ra_deg, dec_deg):
    ra, dec = np.broadcast_arrays(np.radians(ra_deg), np.radians(dec_deg))
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)


def galactic_to_icrf():
    z = unit_ra_dec(NGP_RA, NGP_DEC)
    p = np.array([0.0, 0.0, 1.0]) - z * z[2]
    p /= np.linalg.norm(p)
    q = np.cross(z, p)
    lr = math.radians(L_NCP)
    x = math.cos(lr) * p - math.sin(lr) * q
    return np.stack([x, np.cross(z, x), z], axis=1)  # columns: Galactic axes in ICRF


G2I = galactic_to_icrf()
I2G = G2I.T

# ---------------------------------------------------------------- model


def inner_taper(R):
    return np.where(R < DISK_HOLE, np.exp(-((R - DISK_HOLE) / (0.25 * R0)) ** 2), 1.0)


def thin_disk(R, Z):
    Rh = np.maximum(R, DISK_HOLE)
    d = np.exp(-Rh / THIN_R) / np.cosh(Z / THIN_H) ** 2 * inner_taper(R)
    return np.where(R > THIN_RC, d * np.exp(-(R - THIN_RC) / (THIN_R / 5.0)), d)


def thick_disk(R, Z):
    Rh = np.maximum(R, DISK_HOLE)
    return (THICK_F * math.exp(-R0 / THIN_R) * np.exp(-(Rh - R0) / THICK_R - np.abs(Z) / THICK_Z) *
            inner_taper(R))


def disk_light(R, Z):
    return thin_disk(R, Z) + thick_disk(R, Z)


def bulge_light(X, Y, Z):
    """The two [RO12] ellipsoids, each at its own normalization."""
    c, s = math.cos(BAR_ANGLE), math.sin(BAR_ANGLE)
    along = np.abs(-c * X + s * Y)
    across = np.abs(s * X + c * Y)
    az = np.abs(Z)
    rxy = np.hypot(along, across)
    parts = []
    for shape, x0, y0, z0, norm, rc, c_par, c_perp in BULGE:
        rs = (((along / x0) ** c_perp + (across / y0) ** c_perp) ** (c_par / c_perp) + (az / z0) ** c_par) ** (1 / c_par)
        f = 1.0 / np.cosh(rs) ** 2 if shape == 'sech2' else np.exp(-rs)
        parts.append(norm * f * np.where(rxy > rc, np.exp(-((rxy - rc) / 0.5) ** 2), 1.0))
    return parts


# Emission components, each with its own fitted RGB amplitude (calibration).
COMPONENTS = ['disk', 'bulge_sech2', 'bulge_exp']


def component_light(X, Y, Z, R):
    return [disk_light(R, Z)] + bulge_light(X, Y, Z)


def dust_tau_v(R, Z):
    """Smooth V optical depth per kpc."""
    Rh = np.maximum(R, 0.5 * R0)
    hd = np.where(Rh > DUST_RF, DUST_H0 + DUST_H1 * (Rh - DUST_RF), DUST_H0)
    rho = DUST_RHO0 * np.exp(-Rh / DUST_HR) / np.cosh(Z / hd) ** 2
    hole = np.exp(-((R / R0) - 0.5) ** 2 / 0.25 ** 2)
    return DUST_KAPPA_V * np.where(R < 0.5 * R0, rho * hole, rho)


def hash3(i, j, k, seed):
    h = (i * 73856093) ^ (j * 19349663) ^ (k * 83492791) ^ seed
    h &= 0xFFFFFFFF
    h ^= h >> 16
    h = (h * 0x45D9F3B) & 0xFFFFFFFF
    h ^= h >> 16
    h = (h * 0x45D9F3B) & 0xFFFFFFFF
    h ^= h >> 16
    return h


def _random_unit_vectors(n, seed):
    v = np.random.default_rng(seed).normal(size=(n, 3))
    return (v / np.linalg.norm(v, axis=1)[:, None]).astype(np.float32)


GRADIENTS = _random_unit_vectors(256, NOISE_SEED)


def gradient_noise(x, y, z, seed):
    """Perlin-style gradient noise (random gradients, quintic interpolation)."""
    xi, yi, zi = np.floor(x), np.floor(y), np.floor(z)
    fx, fy, fz = (x - xi).astype(np.float32), (y - yi).astype(np.float32), (z - zi).astype(np.float32)
    sx, sy, sz = (f * f * f * (f * (f * 6 - 15) + 10) for f in (fx, fy, fz))
    xi, yi, zi = xi.astype(np.int64), yi.astype(np.int64), zi.astype(np.int64)
    out = np.zeros(x.shape, np.float32)
    for dx in (0, 1):
        wx = sx if dx else 1 - sx
        for dy in (0, 1):
            wxy = wx * (sy if dy else 1 - sy)
            for dz in (0, 1):
                g = GRADIENTS[hash3(xi + dx, yi + dy, zi + dz, seed) & 255]
                dot = g[:, 0] * (fx - dx) + g[:, 1] * (fy - dy) + g[:, 2] * (fz - dz)
                out += wxy * (sz if dz else 1 - sz) * dot
    return out


def _random_rotation(seed):
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    return q * np.sign(np.diag(r))


# Each octave turned by its own random rotation, so no lattice direction shows.
ROTATIONS = [_random_rotation(NOISE_SEED + o) for o in range(len(OCTAVES))]


def _octave_stats():
    """Standardization of one ridged octave, v = (mean|n| - |n|) / std|n| (zero mean,
    unit variance, peaking on the surfaces where the noise n crosses zero: sheets,
    which look like filaments in projection), and its cumulant generating function
    log E[exp(t v)] on a grid of t, to keep the clumping mean-preserving."""
    rng = np.random.default_rng(1)
    p = rng.uniform(-50, 50, (3, 400000))
    a = np.abs(gradient_noise(p[0], p[1], p[2], 1)).astype(np.float64)
    mean, std = float(a.mean()), float(a.std())
    v = (mean - a) / std
    t = np.linspace(0.0, 4.0, 81)
    k = np.log(np.mean(np.exp(t[:, None] * v[None, ::4]), axis=1))
    return mean, std, t, k


RIDGE_MEAN, RIDGE_STD, CGF_T, CGF_K = _octave_stats()
OCTAVE_NORM = math.sqrt(sum(a * a for _, a in OCTAVES))


def clump_factor(X, Y, Z, footprint):
    """Mean-preserving clumping of the dust density: exp(SIGMA g) / E[exp(SIGMA g)],
    g the normalized sum of the (independent) octaves still resolved."""
    log_f = np.zeros(X.shape, np.float32)
    for o, (cell, amp) in enumerate(OCTAVES):
        fade = np.clip(1.5 - footprint / cell, 0.0, 1.0).astype(np.float32)
        if not fade.any():
            continue
        live = fade > 0
        rot = ROTATIONS[o] / cell
        x, y, z = X[live], Y[live], Z[live]
        n = gradient_noise(rot[0, 0] * x + rot[0, 1] * y + rot[0, 2] * z,
                           rot[1, 0] * x + rot[1, 1] * y + rot[1, 2] * z,
                           rot[2, 0] * x + rot[2, 1] * y + rot[2, 2] * z, NOISE_SEED + 7919 * o)
        coef = SIGMA * amp * fade[live] / OCTAVE_NORM
        v = (RIDGE_MEAN - np.abs(n)) / RIDGE_STD
        log_f[live] += coef * v - np.interp(coef, CGF_T, CGF_K).astype(np.float32)
    return np.exp(log_f)


def galactocentric(P):
    X = P[..., 0] - R0
    Y = P[..., 1]
    Z = P[..., 2] + Z_SUN
    return X, Y, Z, np.hypot(X, Y)


# ---------------------------------------------------------------- luminosity function


def load_hyg():
    """Stars of assets/stars/hyg.csv with a distance: ra, dec, distance (pc), V, B-V (nan if unknown)."""
    rows = []
    with open(os.path.join(ASSETS, 'stars', 'hyg.csv'), encoding='utf-8') as f:
        header = True
        for line in f:
            if line.startswith('#'):
                continue
            if header:
                header = False
                continue
            p = line.rstrip('\n').split(',')
            if not p[4]:
                continue
            rows.append((float(p[2]), float(p[3]), float(p[4]), float(p[5]), float(p[6]) if p[6] else math.nan))
    return np.array(rows)


MAG_EDGES = np.arange(-9.0, 7.51, 0.5)
LF_MAX_PC = 500.0    # Hipparcos distances are usable to here


def mc_volume(center, radius, rng, n=20000):
    """Volume of a sphere (kpc^3) weighted by the disk density relative to the Sun's."""
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    p = center + u * (radius * rng.uniform(size=n) ** (1 / 3))[:, None]
    X, Y, Z, R = galactocentric(p)
    sun = disk_light(np.array(R0), np.array(Z_SUN))
    return 4 / 3 * math.pi * radius ** 3 * float(np.mean(disk_light(R, Z) / sun))


def luminosity_function(cal=None):
    """Stars per kpc^3 (at the Sun's density) per bin of absolute V, and per-bin
    templates (M, B-V). With a calibration, absolute magnitudes are corrected for
    the model's (smooth) extinction."""
    stars = load_hyg()
    d_kpc = stars[:, 2] / 1000.0
    M = stars[:, 3] - 5 * np.log10(stars[:, 2] / 10.0)
    if cal:
        dirs = I2G @ unit_ra_dec(stars[:, 0], stars[:, 1]).T
        M -= extinction_v(np.zeros(3), dirs.T, d_kpc, cal, smooth=True)
    rng = np.random.default_rng(5)
    phi, templates = [], []
    for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
        dc = min(10 ** ((NAKED_EYE_MAG - hi) / 5 + 1), LF_MAX_PC) / 1000.0
        sel = (M >= lo) & (M < hi) & (d_kpc < dc) & (stars[:, 3] <= NAKED_EYE_MAG)
        phi.append(sel.sum() / mc_volume(np.zeros(3), dc, rng))
        templates.append(np.stack([M[sel], stars[sel, 4]], axis=1))
    return np.array(phi), templates


def unresolved_fraction(phi):
    """Fraction of the stellar light (V) at each S_MID in stars fainter than RESOLVED_MAG."""
    mid = 0.5 * (MAG_EDGES[:-1] + MAG_EDGES[1:])
    light = phi * 10 ** (-0.4 * mid)
    m_lim = RESOLVED_MAG - 5 * np.log10(S_MID * 1000.0 / 10.0)
    frac = np.clip((MAG_EDGES[1:][None, :] - m_lim[:, None]) / 0.5, 0.0, 1.0)  # part of each bin fainter
    return (frac * light[None, :]).sum(1) / light.sum()


# ---------------------------------------------------------------- rendering

_W = {}


def _init(params):
    global SIGMA
    _W.update(params)
    SIGMA = params.get('sigma', SIGMA)  # overrides reach the worker processes


def _render_chunk(args):
    dirs = args
    viewer, cal, f_unres = _W['viewer'], _W['cal'], _W['f_unres']
    P = viewer[None, None, :] + dirs[:, None, :] * S_MID[None, :, None]
    X, Y, Z, R = galactocentric(P)
    lights = component_light(X, Y, Z, R)
    footprint = np.broadcast_to(S_MID[None, :] * FOOTPRINT, R.shape)
    dtau = dust_tau_v(R, Z) * cal['dust_scale'] * S_DS[None, :]
    if _W['clumps']:
        dtau *= clump_factor(X, Y, Z, footprint)
    tau = np.cumsum(dtau, axis=1) - 0.5 * dtau
    w = (f_unres * S_DS)[None, :]
    out = np.empty((len(lights), len(dirs), 3))
    for c in range(3):
        t = np.exp(-tau * EXT_RGB[c]) * w
        for i, j in enumerate(lights):
            out[i, :, c] = (j * t).sum(1)
    return out, dtau.sum(1)


def pixel_dirs(width, height):
    """Galactic unit vectors of an equirectangular ICRF map (RA 0h at the centre,
    increasing to the left; north up), as the renderer's milky_way.frag.hlsl."""
    ra = 360.0 * (0.5 - (np.arange(width) + 0.5) / width)
    dec = 90.0 - 180.0 * (np.arange(height) + 0.5) / height
    d = unit_ra_dec(ra[None, :], dec[:, None])
    return d.reshape(-1, 3) @ G2I  # ICRF -> Galactic (row vectors)


def render(viewer, width, height, cal, f_unres, clumps=True, pool_size=None):
    """Light of each component ((N, H, W, 3), unit amplitudes) and V optical depth to infinity."""
    dirs = pixel_dirs(width, height)
    chunks = [dirs[i:i + 2048] for i in range(0, len(dirs), 2048)]
    params = {'viewer': viewer, 'cal': cal, 'f_unres': f_unres, 'clumps': clumps}
    with mp.Pool(pool_size, _init, (params,)) as pool:
        parts = pool.map(_render_chunk, chunks)
    comps = np.concatenate([p[0] for p in parts], axis=1).reshape(len(COMPONENTS), height, width, 3)
    tau = np.concatenate([p[1] for p in parts]).reshape(height, width)
    return comps, tau


def extinction_v(viewer, dirs, dist, cal, smooth=False, steps=48):
    """A_V (mag) from viewer (heliocentric Galactic, kpc) along dirs to dist (kpc)."""
    t = (np.arange(steps) + 0.5) / steps
    s = dist[:, None] * t[None, :]
    P = viewer[None, None, :] + dirs[:, None, :] * s[:, :, None]
    X, Y, Z, R = galactocentric(P)
    k = dust_tau_v(R, Z) * cal['dust_scale']
    if not smooth:
        k = k * clump_factor(X, Y, Z, np.maximum(s * 0.002, 0.005))
    return 1.0857 * (k * (dist / steps)[:, None]).sum(1)


# ---------------------------------------------------------------- maps


def srgb_decode(v):
    v = v.astype(np.float32) / 255.0
    return np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)


def srgb_encode(x):
    x = np.clip(x, 0.0, 1.0)
    return np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(x, 1 / 2.4) - 0.055)


def load_real_map(factor=1):
    img = np.asarray(Image.open(os.path.join(ASSETS, 'textures', 'milky_way.jpg')).convert('RGB'))
    lin = srgb_decode(img)
    if factor > 1:
        h, w = lin.shape[0] // factor, lin.shape[1] // factor
        lin = lin[:h * factor, :w * factor].reshape(h, factor, w, factor, 3).mean(axis=(1, 3))
    return lin


def solid_angle_weights(height, width):
    dec = np.radians(90.0 - 180.0 * (np.arange(height) + 0.5) / height)
    return np.repeat(np.cos(dec)[:, None], width, axis=1)


def map_icrf_dirs(width, height):
    ra = 360.0 * (0.5 - (np.arange(width) + 0.5) / width)
    dec = 90.0 - 180.0 * (np.arange(height) + 0.5) / height
    return unit_ra_dec(ra[None, :], dec[:, None])


def sample_map(lin, dirs_icrf):
    """Bilinear lookup of an equirectangular ICRF map in the given directions."""
    h, w = lin.shape[:2]
    ra = np.degrees(np.arctan2(dirs_icrf[..., 1], dirs_icrf[..., 0]))
    dec = np.degrees(np.arcsin(np.clip(dirs_icrf[..., 2], -1, 1)))
    u = ((0.5 - ra / 360.0) % 1.0) * w - 0.5
    v = np.clip((90.0 - dec) / 180.0 * h - 0.5, 0, h - 1.001)
    u0, v0 = np.floor(u).astype(int), np.floor(v).astype(int)
    fu, fv = (u - u0)[..., None], (v - v0)[..., None]
    u0 %= w
    u1 = (u0 + 1) % w
    return ((lin[v0, u0] * (1 - fu) + lin[v0, u1] * fu) * (1 - fv) +
            (lin[v0 + 1, u0] * (1 - fu) + lin[v0 + 1, u1] * fu) * fv)


def external_galaxies(real, viewer_icrf, width, height, transmission):
    """The LMC, SMC and M31 from the real map (minus its local background), as
    seen from the viewer: each treated as a sheet across its line of sight from
    the Sun, at its distance; dimmed by the viewer's extinction."""
    dirs = map_icrf_dirs(width, height)
    out = np.zeros((height, width, 3))
    all_dirs = map_icrf_dirs(real.shape[1], real.shape[0])
    for ra, dec, radius, dist in EXTERNAL:
        c = unit_ra_dec(ra, dec)
        cos_r = math.cos(math.radians(radius))
        ang = all_dirs @ c
        ring = (ang < cos_r) & (ang > math.cos(math.radians(radius * 1.3)))
        # A high percentile, so that the Milky Way's faint foreground stars do not survive.
        background = np.percentile(real[ring], 90, axis=0)
        along = dirs @ c
        t = (dist - viewer_icrf @ c) / np.maximum(along, 1e-6)
        P = viewer_icrf[None, None, :] + dirs * t[..., None]
        from_sun = P / np.linalg.norm(P, axis=-1, keepdims=True)
        cos_a = from_sun @ c
        a = np.degrees(np.arccos(np.clip(cos_a, -1, 1))) / radius
        t = np.clip((1.0 - a) / 0.4, 0.0, 1.0)
        weight = t * t * (3 - 2 * t) * (along > 0)
        live = weight > 0
        sample = np.zeros((height, width, 3))
        sample[live] = np.maximum(sample_map(real, from_sun[live]) - background, 0.0)
        out += sample * weight[..., None]
    return out * transmission[..., None]


def combine(cal, comps):
    return np.einsum('nhwc,nc->hwc', comps, np.array([cal['amplitudes'][n] for n in COMPONENTS]))


def luminance(lin):
    return lin @ np.array([0.2126, 0.7152, 0.0722])


def mean_luminance(lin):
    w = solid_angle_weights(*lin.shape[:2])
    return float((luminance(lin) * w).sum() / w.sum())


def save_preview(lin, path, scale):
    Image.fromarray((srgb_encode(lin * scale) * 255 + 0.5).astype(np.uint8)).save(path)


# ---------------------------------------------------------------- calibration


def block_mean(img, n):
    h, w = img.shape[0] // n, img.shape[1] // n
    return img[:h * n, :w * n].reshape(h, n, w, n, *img.shape[2:]).mean(axis=(1, 3))


def fit_log(comps, real, weight):
    """Per channel: amplitudes a_i >= 0 minimizing the weighted squared
    log(sum a_i comp_i) - log(real) (damped Gauss-Newton on log a_i)."""
    n = len(comps)
    amps = np.zeros((n, 3))
    w = weight.ravel()
    for c in range(3):
        C = comps[..., c].reshape(n, -1)
        y = np.log(real[..., c].ravel())
        la = np.full(n, float(np.median(y - np.log(C[0]))))
        for i in range(1, n):  # start each bulge part at ~1% of the disk near its peak
            la[i] = la[0] + float(np.log(C[0][np.argmax(C[i])] / C[i].max())) - 4.0
        for _ in range(300):
            A = np.exp(la)
            mod = A @ C
            r = np.log(mod) - y
            J = (C * A[:, None] / mod).T
            H = J.T @ (J * w[:, None])
            step = np.linalg.solve(H + 1e-6 * np.trace(H) * np.eye(n), J.T @ (r * w))
            step = np.clip(step, -1.0, 1.0)
            la -= step
            la = np.maximum(la, la[0] - 40.0)
            if np.abs(step).max() < 1e-6:
                break
        amps[:, c] = np.exp(la)
    return {name: [float(v) for v in amps[i]] for i, name in enumerate(COMPONENTS)}


def calibrate(out_dir):
    """Fits the model to the real map, both averaged in 2.8 deg cells."""
    w, h = 512, 256
    real = block_mean(load_real_map(4096 // w), 4)
    phi, _ = luminosity_function()
    f_unres = unresolved_fraction(phi)
    weight = block_mean(solid_angle_weights(h, w), 4)
    dirs = map_icrf_dirs(w // 4, h // 4)
    for ra, dec, radius, _ in EXTERNAL:  # not part of the model
        weight[dirs @ unit_ra_dec(ra, dec) > math.cos(math.radians(radius * 1.3))] = 0.0
    real = np.maximum(real, 1e-5)
    best = None
    for k in np.geomspace(0.7, 2.8, 13):
        cal = {'dust_scale': float(k)}
        comps, _ = render(np.zeros(3), w, h, cal, f_unres)
        comps = np.stack([block_mean(x, 4) for x in comps])
        cal['amplitudes'] = fit_log(comps, real, weight)
        model = combine(cal, comps)
        r = np.log(luminance(model)) - np.log(luminance(real))
        cost = float(np.sqrt((r * r * weight).sum() / weight.sum()))
        print(f'dust scale {k:.3f}: rms log residual {cost:.4f}')
        if best is None or cost < best[0]:
            best = (cost, cal, model)
    cost, cal, model = best
    cal['rms_log_residual'] = cost
    tau_sun = float(dust_tau_v(np.array(R0), np.array(Z_SUN)) * cal['dust_scale'])
    cal['a_v_per_kpc_at_sun'] = 1.0857 * tau_sun
    # The mean sky brightness from the Sun, rendered as the scenes are (with clumps).
    comps, tau = render(np.zeros(3), 512, 256, cal, f_unres)
    sun_sky = combine(cal, comps) + external_galaxies(load_real_map(2), np.zeros(3), 512, 256,
                                                             np.exp(-tau * EXT_RGB[1]))
    cal['mean_luminance_sun'] = mean_luminance(sun_sky)
    cal['mean_luminance_sun_real'] = mean_luminance(load_real_map(8))
    print(json.dumps(cal, indent=2))
    with open(CALIBRATION, 'w', encoding='utf-8', newline='\n') as f:
        json.dump(cal, f, indent=2)
        f.write('\n')
    if out_dir:
        np.save(os.path.join(out_dir, 'calibration_real.npy'), real)
        np.save(os.path.join(out_dir, 'calibration_model.npy'), model)
        np.save(os.path.join(out_dir, 'calibration_sun_clumps.npy'), sun_sky)


def load_calibration():
    with open(CALIBRATION, encoding='utf-8') as f:
        return json.load(f)


# ---------------------------------------------------------------- stars


def sample_stars(viewer, cal, rng):
    """Point stars brighter than NAKED_EYE_MAG seen from viewer (heliocentric Galactic, kpc):
    returns ICRF unit vectors, V and B-V."""
    phi, templates = luminosity_function(cal)
    positions, mags, colors = [], [], []
    for b, (lo, _) in enumerate(zip(MAG_EDGES[:-1], MAG_EDGES[1:])):
        if phi[b] <= 0 or len(templates[b]) == 0:
            continue
        D = min(10 ** ((NAKED_EYE_MAG - lo) / 5 + 1), 15000.0) / 1000.0
        n = rng.poisson(phi[b] * mc_volume(viewer, D, rng))
        if n == 0:
            continue
        sun = disk_light(np.array(R0), np.array(Z_SUN))
        got = []
        while sum(len(g) for g in got) < n:
            m = max(4 * n, 1000)
            u = rng.normal(size=(m, 3))
            u /= np.linalg.norm(u, axis=1)[:, None]
            p = viewer + u * (D * rng.uniform(size=m) ** (1 / 3))[:, None]
            X, Y, Z, R = galactocentric(p)
            ratio = disk_light(R, Z) / sun
            got.append(p[rng.uniform(size=m) * max(ratio.max(), 1e-12) < ratio])
        p = np.concatenate(got)[:n]
        t = templates[b][rng.integers(len(templates[b]), size=n)]
        positions.append(p)
        mags.append(t[:, 0])
        colors.append(t[:, 1])
    p = np.concatenate(positions)
    M = np.concatenate(mags)
    bv0 = np.concatenate(colors)
    rel = p - viewer
    dist = np.linalg.norm(rel, axis=1)
    dirs = rel / dist[:, None]
    av = np.concatenate([extinction_v(viewer, dirs[i:i + 4096], dist[i:i + 4096], cal)
                         for i in range(0, len(dirs), 4096)])
    V = M + 5 * np.log10(dist * 1000.0 / 10.0) + av
    bv = bv0 + av / R_V
    keep = V <= NAKED_EYE_MAG
    return dirs[keep] @ I2G, V[keep], bv[keep]


def write_stars(path, dirs_icrf, V, bv, header):
    ra = np.degrees(np.arctan2(dirs_icrf[:, 1], dirs_icrf[:, 0])) % 360.0
    dec = np.degrees(np.arcsin(np.clip(dirs_icrf[:, 2], -1, 1)))
    order = np.argsort(V)
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        for line in header:
            f.write(f'# {line}\n')
        f.write('ra_deg,dec_deg,vmag,bv\n')
        for i in order:
            b = '' if math.isnan(bv[i]) else f'{bv[i]:.3f}'
            f.write(f'{ra[i]:.4f},{dec[i]:.4f},{V[i]:.2f},{b}\n')


# ---------------------------------------------------------------- scenes


def viewer_from(ra, dec, dist_pc):
    icrf = unit_ra_dec(ra, dec) * dist_pc / 1000.0
    return icrf, I2G @ icrf


def sky_for(viewer_ra_dec_pc, width, height, cal, rng):
    viewer_icrf, viewer = viewer_from(*viewer_ra_dec_pc)
    phi, _ = luminosity_function()
    comps, tau = render(viewer, width, height, cal, unresolved_fraction(phi))
    model = combine(cal, comps)
    real = load_real_map(1)
    model += external_galaxies(real, viewer_icrf, width, height, np.exp(-tau * EXT_RGB[1]))
    dirs, V, bv = sample_stars(viewer, cal, rng)
    return model, dirs, V, bv


def brightness_for(model, cal):
    """Displayed brightness: Earth's default, with the overall level compressed like
    the stars' fluxes (relative to the model sky at the Sun)."""
    ratio = mean_luminance(model) / cal['mean_luminance_sun']
    return DEFAULT_BRIGHTNESS * ratio ** (MAG_GAMMA - 1.0), ratio


def write_map(path, model, white):
    rng = np.random.default_rng(0)
    enc = srgb_encode(model / white) * 255 + rng.uniform(-0.5, 0.5, model.shape)
    Image.fromarray(np.clip(enc + 0.5, 0, 255).astype(np.uint8)).save(path, quality=92, subsampling=0)


def render_scene(path, width=2048, height=1024):
    with open(path, 'rb') as f:
        scene = tomllib.load(f)
    sky = scene.get('sky', {})
    viewer = sky.get('viewer_ra_dec_distance_pc')
    if not viewer or 'stars' not in sky or 'milky_way' not in sky:
        sys.exit(f'{path}: [sky] needs viewer_ra_dec_distance_pc, milky_way and stars')
    cal = load_calibration()
    rng = np.random.default_rng(NOISE_SEED)
    model, dirs, V, bv = sky_for(viewer, width, height, cal, rng)
    white = float(np.percentile(model.max(axis=2), 99.99))
    brightness, ratio = brightness_for(model, cal)
    write_map(os.path.join(ASSETS, sky['milky_way']), model, white)
    header = [
        f'Model sky for {scene.get("name", path)}: stars brighter than V = {NAKED_EYE_MAG} seen from the viewer.',
        'Generated by tools/sky/make_galaxy_sky.py (see assets/sky/SOURCES.md); statistical, not real stars.',
        f'viewer_ra_dec_distance_pc = {viewer[0]}, {viewer[1]}, {viewer[2]}',
        f'milky_way_brightness = {brightness * white:.4f}',
    ]
    write_stars(os.path.join(ASSETS, sky['stars']), dirs, V, bv, header)
    print(f'{scene.get("name", path)}: {len(V)} stars (brightest V = {V.min():.2f}); '
          f'mean sky {ratio:.2f} x the Sun\'s; milky_way_brightness = {brightness * white:.4f}')


def sun(out_dir):
    cal = load_calibration()
    rng = np.random.default_rng(NOISE_SEED)
    model, dirs, V, bv = sky_for((0.0, 0.0, 0.0), 1024, 512, cal, rng)
    b = np.degrees(np.arcsin((dirs @ G2I)[:, 2]))
    real = load_hyg()
    real_b = np.degrees(np.arcsin((unit_ra_dec(real[:, 0], real[:, 1]) @ G2I)[:, 2]))
    real_v = real[:, 3] <= NAKED_EYE_MAG
    print(f'model from the Sun: {len(V)} stars to V = {NAKED_EYE_MAG}, {np.mean(np.abs(b) < 10):.3f} within 10 deg '
          f'of the plane, {np.sum(V < 2):.0f} brighter than V = 2')
    print(f'real (hyg.csv with distances): {real_v.sum()} stars, '
          f'{np.mean(np.abs(real_b[real_v]) < 10):.3f} within 10 deg, {np.sum(real[:, 3] < 2):.0f} brighter than V = 2')
    if out_dir:
        scale = 4.0 / np.percentile(luminance(model), 99.5)
        save_preview(model, os.path.join(out_dir, 'model_from_sun.png'), scale)


def main():
    args = sys.argv[1:]
    out_dir = None
    if '--out' in args:
        i = args.index('--out')
        out_dir = args[i + 1]
        del args[i:i + 2]
    if not args:
        sys.exit(__doc__)
    if args[0] == 'calibrate':
        calibrate(out_dir)
    elif args[0] == 'sun':
        sun(out_dir)
    elif args[0] == 'render' and len(args) > 1:
        for path in args[1:]:
            render_scene(path)
    else:
        sys.exit(__doc__)


if __name__ == '__main__':
    main()
