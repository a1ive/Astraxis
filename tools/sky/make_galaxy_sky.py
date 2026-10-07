"""Model Milky Way skies for the scenes far from the Sun (assets/sky/).

From the Sun the sky is real (assets/stars/hyg.csv, assets/textures/milky_way.jpg).
From a star 0.5-10 kpc away, or from the Galactic Centre, it is modelled here: the
Galaxy's stars as a V-band emissivity in L_sun/pc^3 (disks, bulge, nuclear disk
and cluster), dust (a disk, the nuclear molecular disk and the circumnuclear disk)
with statistical clumps, seen from the scene's viewer, plus point stars drawn from
the solar neighbourhood's luminosity function. The model, its sources and the
calibration are described in assets/sky/SOURCES.md.

Usage (needs numpy and Pillow):
    python make_galaxy_sky.py calibrate [--out DIR]
        Fits the conversion from emissivity to map units and the disk's dust to
        milky_way.jpg as seen from the Sun; writes tools/sky/galaxy_calibration.json
        (DIR gets the arrays compared).
    python make_galaxy_sky.py sun [--out DIR]
        Validation: the model sky from the Sun (star counts; a map in DIR).
    python make_galaxy_sky.py render <scene.toml> [...]
        For each scene whose [sky] has viewer_ra_dec_distance_pc, milky_way and stars,
        writes the map and the star list (paths relative to assets/).
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
SGR_A = (266.41683333, -29.00781611)  # [RB04], as sgr_a.toml; at R0

# ---------------------------------------------------------------- stars
# Everything is a V-band emissivity in L_sun / pc^3 (M_V,sun = 4.81 [WI18]).
MV_SUN = 4.81
# Stellar thin disk [DS01] Table 2: exp(-R/r) sech^2(z/h), cut off beyond r_c by
# exp(-(R - r_c) / (r / 5)) (eq. 27). [DS01] fit the disk only beyond 0.35 R0
# (the bulge region was excluded); inside, it is tapered here like their dust
# hole (Gaussian of width 0.25 R0), leaving the centre to the bulge. Normalized
# to the local emissivity measured from hyg.csv (luminosity_function).
THIN_R, THIN_H, THIN_RC = 2.264 * DS, 0.2822 * DS, 10.52 * DS
DISK_HOLE = 0.35 * R0
# Thick disk [BHG16] sec. 5.1.3 / 5.2.2: 4% of the thin disk's local density,
# exponential, z_T = 900 pc, R_T = 2.0 kpc.
THICK_F, THICK_Z, THICK_R = 0.04, 0.900, 2.0
# Mass-to-light ratio of the old populations below (bulge, nuclear disk and
# cluster): the Galaxy's M/L_V = 1.70 [BHG16] Table 2.
ML_V = 1.70
# Bulge/bar: the boxy sech^2 ellipsoid of [RO12] Table 2 model S+E (major axis
# 12.9 deg from the Sun - centre line, near end at positive longitude), density
# f(Rs) with Rs^c_par = [(|x|/x0)^c_perp + (|y|/y0)^c_perp]^(c_par/c_perp) +
# (|z|/z0)^c_par, cut off as exp(-((R_xy - R_c) / 0.5 kpc)^2) beyond R_c
# (sec. 2.2). Its second, exponential ellipsoid is left out: fitted to the real
# map from the Sun its amplitude goes to zero (it would light up the sky well
# above and below the bulge). Mass: 1.55e10 M_sun inside the VVV box
# (+-2.2 x +-1.4 x +-1.2 kpc), the middle of [BHG16] sec. 4.2.4's 1.4-1.7e10.
# (Red-clump studies find a larger bar angle, ~27 deg [BHG16].)
BAR_ANGLE = math.radians(12.9)
BULGE_SHAPE = (1.46, 0.49, 0.39, 3.43, 3.007, 3.329)  # x0, y0, z0, R_c, c_par, c_perp
BULGE_MASS_BOX = 1.55e10
# Nuclear stellar disk [LA02] sec. 5.2, eq. 11: emissivity ~ exp(ln(1/2) (R/R_half)^n),
# two components R_half = 120 and 220 pc (taken equal at the centre), n = 5;
# vertically the same form with half width 45 pc and n = 1.4; 1.4e9 M_sun (Table 7).
NSD_RADII, NSD_N, NSD_HZ, NSD_NZ, NSD_MASS = (0.120, 0.220), 5.0, 0.045, 1.4, 1.4e9
# Nuclear star cluster [SC18] sec. 4.4 / Table 2 (mean): 3D Nuker law
# rho = rho(rb) 2^((beta-gamma)/alpha) (r/rb)^-gamma (1 + (r/rb)^alpha)^((gamma-beta)/alpha),
# rb = 3.1 pc, gamma = 1.13, beta = 3.5, alpha = 10; rho(1 pc) = 1.5e5 M_sun/pc^3
# (Table 3, for 2.5e7 M_sun [SC14]); flattened along the Galactic plane with q = 0.71
# [SC14]. Centred on Sgr A*.
NSC_RB, NSC_GAMMA, NSC_BETA, NSC_ALPHA, NSC_RHO1, NSC_Q = 0.0031, 1.13, 3.5, 10.0, 1.5e5, 0.71

# ---------------------------------------------------------------- dust
# Dust disk [DS01] Table 1: rho0 exp(-R/h_r) sech^2(z/h_d), h_d flaring linearly
# beyond r_f; inside 0.5 R0 a Gaussian hole (sec. 3.2). Density in MJy/sr/kpc,
# V opacity 0.0180 (MJy/sr)^-1 (Table 2): tau_V per kpc. The spiral arms and the
# local arm of [DS01] are left out; dust_scale (calibrated) absorbs them on average.
DUST_RHO0, DUST_HR = 1098.0, 2.26 * DS
DUST_H0, DUST_H1, DUST_RF = 0.1344 * DS, 0.0148, 4.40 * DS
DUST_KAPPA_V = 0.0180
# Gas near the centre, as hydrogen mass; A_V = N_H / 2.21e21 cm^-2 [GO09]. One
# M_sun/pc^3 of hydrogen is 40.5 atoms/cm^3: A_V per kpc = 40.5 * 3.086e21 / 2.21e21
# per M_sun/pc^3.
AV_PER_KPC_PER_MSUN_PC3 = 40.5 * 3.086e21 / 2.21e21
# Nuclear molecular disk [LA02] abstract / sec. 5.5: 2e7 M_sun of hydrogen, a warm
# inner disk of radius 110 pc and a cold outer torus (> 80% of the mass) of the
# nuclear disk's size; taken here as 20% in exp(ln(1/2) (R/110 pc)^5) and 80% in
# the difference of that and the 220 pc profile, with the nuclear disk's height.
NMD_MASS, NMD_INNER = 2e7, 0.2
# Circumnuclear disk [GE10] sec. 3.2: dense clumps (0.2-0.3 pc) at R = 1.5-4 pc,
# tilted 20-30 deg to the Galactic plane; a few 1e4 M_sun of gas from its dust
# emission (virial estimates reach 1e6). Taken: 3e4 M_sun of hydrogen, uniform in
# 1.5-4 pc with soft edges, Gaussian half width 0.25 pc, tilted 25 deg about the
# Galactic y axis (the sense is assumed); clumps finer than the Galaxy's (below).
CND_R, CND_H, CND_MASS, CND_TILT = (0.0015, 0.004), 0.00025, 3e4, math.radians(25.0)

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
CND_OCTAVES = [(0.0016, 1.0), (0.0008, 0.71), (0.0004, 0.5), (0.0002, 0.35), (0.0001, 0.25)]
FOOTPRINT = math.pi / 1024
NOISE_SEED = 20261007

# ---------------------------------------------------------------- stars and maps
# Light of stars fainter than this (from the viewer) is in the diffuse map: the
# SVS map leaves out the Hipparcos/Tycho stars (assets/textures/SOURCES.md).
RESOLVED_MAG = 11.5
NAKED_EYE_MAG = 6.5         # point stars, as scene/star_catalog.hpp kNakedEyeMag
MAX_POINT_STARS = 40000     # beyond this (the Galactic Centre) the point/diffuse limit moves brighter
MAG_GAMMA = 0.56            # star flux compression, render/starfield_pass.cpp kMagnitudeGamma
DEFAULT_BRIGHTNESS = 0.25   # scene/scene.hpp SceneSky::milky_way_brightness (Earth's map)
# The overall sky level (relative to the Sun's) is shown compressed like star fluxes,
# up to this many times the Sun's; beyond it the eye adapts: map and stars together
# get one exposure factor (the Galactic Centre, ~4600 times the Sun's sky).
SKY_LEVEL_CAP = 3.0

# Galaxies outside ours, carried over from the real map ([SIMBAD] centre and
# size; distances [PI19], [GR20s], [LI21]): (ra, dec, radius deg, distance kpc).
EXTERNAL = [
    (80.8942, -69.7561, 6.0, 49.59),   # LMC
    (13.1583, -72.8003, 3.0, 62.44),   # SMC
    (10.6847, 41.2688, 2.2, 761.0),    # M31
]


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
NUCLEUS = I2G @ unit_ra_dec(*SGR_A) * R0  # Sgr A*, heliocentric Galactic (kpc)


def ray_grid(s_min, steps):
    """Log-spaced samples along a ray from s_min to 40 kpc: midpoints and lengths."""
    edges = np.geomspace(s_min, 40.0, steps + 1)
    return np.sqrt(edges[1:] * edges[:-1]), np.diff(edges)


def galactocentric(P):
    X = P[..., 0] - R0
    Y = P[..., 1]
    Z = P[..., 2] + Z_SUN
    return X, Y, Z, np.hypot(X, Y)


# ---------------------------------------------------------------- stellar emissivity


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


def disk_shape(R, Z):
    return thin_disk(R, Z) + thick_disk(R, Z)


DISK_AT_SUN = float(disk_shape(np.array(R0), np.array(Z_SUN)))


def bulge_shape(X, Y, Z):
    c, s = math.cos(BAR_ANGLE), math.sin(BAR_ANGLE)
    along = np.abs(-c * X + s * Y)
    across = np.abs(s * X + c * Y)
    x0, y0, z0, rc, c_par, c_perp = BULGE_SHAPE
    rs = (((along / x0) ** c_perp + (across / y0) ** c_perp) ** (c_par / c_perp) + (np.abs(Z) / z0) ** c_par) ** (1 / c_par)
    rxy = np.hypot(along, across)
    return np.where(rxy > rc, np.exp(-((rxy - rc) / 0.5) ** 2), 1.0) / np.cosh(rs) ** 2


def _bulge_density_scale():
    """M_sun/pc^3 per unit of bulge_shape: BULGE_MASS_BOX inside the VVV box (in the bar frame)."""
    x = np.linspace(-2.2, 2.2, 221)
    y = np.linspace(-1.4, 1.4, 141)
    z = np.linspace(-1.2, 1.2, 121)
    A, B, Z = np.meshgrid(x, y, z, indexing='ij')
    c, s = math.cos(BAR_ANGLE), math.sin(BAR_ANGLE)
    total = bulge_shape(-c * A + s * B, s * A + c * B, Z).sum() * (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]) * 1e9
    return BULGE_MASS_BOX / total


BULGE_RHO = _bulge_density_scale()


def half_profile(x, half, n):
    return np.exp(math.log(0.5) * (x / half) ** n)


def nsd_shape(Rn, Zn):
    radial = 0.5 * (half_profile(Rn, NSD_RADII[0], NSD_N) + half_profile(Rn, NSD_RADII[1], NSD_N))
    return radial * half_profile(np.abs(Zn), NSD_HZ, NSD_NZ)


def _volume_integral(shape_rz, r_max, z_max):
    """Integral over an axisymmetric shape(R, z) in pc^3."""
    r = np.linspace(0.0, r_max, 1201)[1:]
    z = np.linspace(-z_max, z_max, 801)
    Rg, Zg = np.meshgrid(r, z, indexing='ij')
    return float((shape_rz(Rg, Zg) * 2 * math.pi * Rg).sum() * (r[1] - r[0]) * (z[1] - z[0]) * 1e9)


NSD_RHO = NSD_MASS / _volume_integral(nsd_shape, 0.6, 0.4)  # M_sun/pc^3 at the centre


def nsc_density(rn):
    """M_sun/pc^3 at (flattened) radius rn (kpc)."""
    x = np.maximum(rn, 1e-5) / NSC_RB
    g, b, a = NSC_GAMMA, NSC_BETA, NSC_ALPHA
    shape = x ** -g * (1 + x ** a) ** ((g - b) / a)
    x1 = 0.001 / NSC_RB
    return NSC_RHO1 * shape / (x1 ** -g * (1 + x1 ** a) ** ((g - b) / a))


def nuclear_coords(P):
    """Sgr A*-centred Galactic axes (kpc): in-plane radius, height."""
    d = P - NUCLEUS
    return np.hypot(d[..., 0], d[..., 1]), d[..., 2], d


def emissivity(P, j_sun):
    """V-band stellar emissivity (L_sun/pc^3) at heliocentric Galactic positions P (kpc)."""
    X, Y, Z, R = galactocentric(P)
    Rn, Zn, _ = nuclear_coords(P)
    mass = BULGE_RHO * bulge_shape(X, Y, Z) + NSD_RHO * nsd_shape(Rn, Zn)
    mass += nsc_density(np.sqrt(Rn ** 2 + (Zn / NSC_Q) ** 2))
    return j_sun * disk_shape(R, Z) / DISK_AT_SUN + mass / ML_V


# ---------------------------------------------------------------- dust


def dust_disk_tau(R, Z):
    """Smooth [DS01] V optical depth per kpc (before dust_scale)."""
    Rh = np.maximum(R, 0.5 * R0)
    hd = np.where(Rh > DUST_RF, DUST_H0 + DUST_H1 * (Rh - DUST_RF), DUST_H0)
    rho = DUST_RHO0 * np.exp(-Rh / DUST_HR) / np.cosh(Z / hd) ** 2
    hole = np.exp(-((R / R0) - 0.5) ** 2 / 0.25 ** 2)
    return DUST_KAPPA_V * np.where(R < 0.5 * R0, rho * hole, rho)


def nmd_shape(Rn, Zn):
    inner = half_profile(Rn, 0.110, 5.0)
    torus = np.maximum(half_profile(Rn, 0.220, 5.0) - inner, 0.0)
    return (NMD_INNER * inner / NMD_NORM[0] + (1 - NMD_INNER) * torus / NMD_NORM[1]) * \
        half_profile(np.abs(Zn), NSD_HZ, NSD_NZ)


NMD_NORM = (_volume_integral(lambda r, z: half_profile(r, 0.110, 5.0) * half_profile(np.abs(z), NSD_HZ, NSD_NZ), 0.6, 0.4),
            _volume_integral(lambda r, z: np.maximum(half_profile(r, 0.220, 5.0) - half_profile(r, 0.110, 5.0), 0.0) *
                             half_profile(np.abs(z), NSD_HZ, NSD_NZ), 0.6, 0.4))


def cnd_frame(d):
    """Nuclear offsets d (kpc) in the CND's frame: in-plane radius, height."""
    c, s = math.cos(CND_TILT), math.sin(CND_TILT)
    x = c * d[..., 0] - s * d[..., 2]
    z = s * d[..., 0] + c * d[..., 2]
    return np.hypot(x, d[..., 1]), z


def cnd_shape(Rc, Zc):
    lo, hi = CND_R
    edge = 0.0003
    radial = 1 / (1 + np.exp(-(Rc - lo) / edge * 4)) / (1 + np.exp((Rc - hi) / edge * 4))
    return radial * np.exp(-0.5 * (Zc / (CND_H / 1.1774)) ** 2)


CND_RHO = CND_MASS / _volume_integral(cnd_shape, 0.006, 0.002)  # M_sun/pc^3 of hydrogen


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


def clump_factor(X, Y, Z, footprint, octaves=OCTAVES, seed=NOISE_SEED):
    """Mean-preserving clumping of a dust density: exp(SIGMA g) / E[exp(SIGMA g)],
    g the normalized sum of the (independent) octaves still resolved. Each octave
    is turned by its own random rotation, so no lattice direction shows."""
    norm = math.sqrt(sum(a * a for _, a in octaves))
    log_f = np.zeros(X.shape, np.float32)
    for o, (cell, amp) in enumerate(octaves):
        fade = np.clip(1.5 - footprint / cell, 0.0, 1.0).astype(np.float32)
        if not fade.any():
            continue
        live = fade > 0
        rot = _random_rotation(seed + o) / cell
        x, y, z = X[live], Y[live], Z[live]
        n = gradient_noise(rot[0, 0] * x + rot[0, 1] * y + rot[0, 2] * z,
                           rot[1, 0] * x + rot[1, 1] * y + rot[1, 2] * z,
                           rot[2, 0] * x + rot[2, 1] * y + rot[2, 2] * z, seed + 7919 * o)
        coef = SIGMA * amp * fade[live] / norm
        v = (RIDGE_MEAN - np.abs(n)) / RIDGE_STD
        log_f[live] += coef * v - np.interp(coef, CGF_T, CGF_K).astype(np.float32)
    return np.exp(log_f)


def dust_tau(P, cal, footprint=None):
    """V optical depth per kpc at heliocentric positions P (kpc); clumped unless
    footprint (kpc, the resolution at each point) is None."""
    X, Y, Z, R = galactocentric(P)
    Rn, Zn, d = nuclear_coords(P)
    disk = dust_disk_tau(R, Z) * cal['dust_scale']
    # Inside the CND the gas is hot and ionized (the central cavity [GE10] sec. 3.1-3.2).
    cavity = 1 / (1 + np.exp(-(np.hypot(Rn, Zn) - CND_R[0]) / 0.0001))
    nmd = NMD_MASS * nmd_shape(Rn, Zn) * cavity
    near = np.hypot(Rn, Zn) < 0.01  # the CND only matters within 10 pc
    cnd = np.zeros(np.shape(R))
    if near.any():
        cnd[near] = CND_RHO * cnd_shape(*cnd_frame(d[near]))
    if footprint is not None:
        clumps = clump_factor(X, Y, Z, footprint)  # one interstellar medium for both
        disk = disk * clumps
        nmd = nmd * clumps
        if near.any():
            cnd[near] *= clump_factor(X[near], Y[near], Z[near], footprint[near], CND_OCTAVES, NOISE_SEED + 101)
    return disk + (nmd + cnd) * AV_PER_KPC_PER_MSUN_PC3 / 1.0857


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
MAG_MID = 0.5 * (MAG_EDGES[:-1] + MAG_EDGES[1:])
LF_MAX_PC = 500.0    # Hipparcos distances are usable to here


def disk_volume(center, radius, rng, n=20000):
    """Volume of a sphere (kpc^3) weighted by the disk density relative to the Sun's."""
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    p = center + u * (radius * rng.uniform(size=n) ** (1 / 3))[:, None]
    X, Y, Z, R = galactocentric(p)
    return 4 / 3 * math.pi * radius ** 3 * float(np.mean(disk_shape(R, Z) / DISK_AT_SUN))


def luminosity_function(cal):
    """Stars per kpc^3 (at the Sun's density) per bin of absolute V, per-bin templates
    (M, B-V), and the local V emissivity (L_sun/pc^3). The volume of each bin is
    the sphere within which it is complete in hyg.csv (V <= 6.5, at most 500 pc),
    weighted by the disk's density. Absolute magnitudes are corrected for the
    model's mean (unclumped) extinction, which the model then applies again."""
    stars = load_hyg()
    d_kpc = stars[:, 2] / 1000.0
    M = stars[:, 3] - 5 * np.log10(stars[:, 2] / 10.0)
    dirs = unit_ra_dec(stars[:, 0], stars[:, 1]) @ G2I
    M -= np.concatenate([extinction_v(np.zeros(3), dirs[i:i + 4096], d_kpc[i:i + 4096], cal, clumped=False)
                         for i in range(0, len(dirs), 4096)])
    rng = np.random.default_rng(5)
    phi, templates = [], []
    for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
        dc = min(10 ** ((NAKED_EYE_MAG - hi) / 5 + 1), LF_MAX_PC) / 1000.0
        sel = (M >= lo) & (M < hi) & (d_kpc < dc) & (stars[:, 3] <= NAKED_EYE_MAG)
        phi.append(sel.sum() / disk_volume(np.zeros(3), dc, rng))
        templates.append(np.stack([M[sel], stars[sel, 4]], axis=1))
    phi = np.array(phi)
    j_sun = float((phi / 1e9 * 10 ** (-0.4 * (MAG_MID - MV_SUN))).sum())
    return phi, templates, j_sun


def unresolved_fraction(phi, s_mid, resolved_mag):
    """Fraction of the stellar light (V) at distances s_mid (kpc) in stars fainter than resolved_mag."""
    light = phi * 10 ** (-0.4 * MAG_MID)
    m_lim = resolved_mag - 5 * np.log10(s_mid * 1000.0 / 10.0)
    frac = np.clip((MAG_EDGES[1:][None, :] - m_lim[:, None]) / 0.5, 0.0, 1.0)  # part of each bin fainter
    return (frac * light[None, :]).sum(1) / light.sum()


# ---------------------------------------------------------------- rendering

_W = {}


def _init(params):
    _W.update(params)


def _render_chunk(dirs):
    viewer, cal, s_mid, s_ds = _W['viewer'], _W['cal'], _W['s_mid'], _W['s_ds']
    P = viewer[None, None, :] + dirs[:, None, :] * s_mid[None, :, None]
    j = emissivity(P, _W['j_sun'])
    footprint = np.broadcast_to(s_mid[None, :] * FOOTPRINT, j.shape) if _W['clumps'] else None
    dtau = dust_tau(P, cal, footprint) * s_ds[None, :]
    tau = np.cumsum(dtau, axis=1) - 0.5 * dtau
    w = j * (_W['f_unres'] * s_ds)[None, :]
    out = np.empty((len(dirs), 3))
    for c in range(3):
        out[:, c] = (w * np.exp(-tau * EXT_RGB[c])).sum(1)
    return out, dtau.sum(1)


def pixel_dirs(width, height):
    """Galactic unit vectors of an equirectangular ICRF map (RA 0h at the centre,
    increasing to the left; north up), as the renderer's milky_way.frag.hlsl."""
    ra = 360.0 * (0.5 - (np.arange(width) + 0.5) / width)
    dec = 90.0 - 180.0 * (np.arange(height) + 0.5) / height
    d = unit_ra_dec(ra[None, :], dec[:, None])
    return d.reshape(-1, 3) @ G2I  # ICRF -> Galactic (row vectors)


def render(viewer, width, height, cal, lf, resolved_mag=RESOLVED_MAG, s_min=0.001, steps=192, clumps=True):
    """Unresolved starlight (H, W, 3) in L_sun/pc^3 * kpc (times V-band transmission
    per channel) and the V optical depth to infinity (H, W)."""
    phi, _, j_sun = lf
    s_mid, s_ds = ray_grid(s_min, steps)
    dirs = pixel_dirs(width, height)
    chunks = [dirs[i:i + 2048] for i in range(0, len(dirs), 2048)]
    params = {'viewer': viewer, 'cal': cal, 'j_sun': j_sun, 's_mid': s_mid, 's_ds': s_ds, 'clumps': clumps,
              'f_unres': unresolved_fraction(phi, s_mid, resolved_mag)}
    with mp.Pool(None, _init, (params,)) as pool:
        parts = pool.map(_render_chunk, chunks)
    light = np.concatenate([p[0] for p in parts]).reshape(height, width, 3)
    tau = np.concatenate([p[1] for p in parts]).reshape(height, width)
    return light, tau


def extinction_v(viewer, dirs, dist, cal, steps=64, clumped=True):
    """A_V (mag) from viewer (heliocentric Galactic, kpc) along dirs to dist (kpc).
    Log-spaced from 0.01 pc, so the dust right around the viewer is resolved."""
    start = np.minimum(1e-5 / np.maximum(dist, 1e-9), 0.5)
    edges = dist[:, None] * np.geomspace(start, 1.0, steps + 1).T
    s = np.sqrt(edges[:, 1:] * edges[:, :-1])
    ds = np.diff(edges, axis=1)
    P = viewer[None, None, :] + dirs[:, None, :] * s[:, :, None]
    k = dust_tau(P, cal, np.maximum(s * FOOTPRINT, 1e-6) if clumped else None)
    return 1.0857 * (k * ds).sum(1)


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


def luminance(lin):
    return lin @ np.array([0.2126, 0.7152, 0.0722])


def mean_luminance(lin):
    w = solid_angle_weights(*lin.shape[:2])
    return float((luminance(lin) * w).sum() / w.sum())


# ---------------------------------------------------------------- calibration


def block_mean(img, n):
    h, w = img.shape[0] // n, img.shape[1] // n
    return img[:h * n, :w * n].reshape(h, n, w, n, *img.shape[2:]).mean(axis=(1, 3))


def calibrate(out_dir):
    """Fits, per channel, the factor from rendered light to the real map's units
    (log space, 2.8 deg cells, weighted by solid angle), scanning the disk's dust."""
    w, h = 512, 256
    real = np.maximum(block_mean(load_real_map(4096 // w), 4), 1e-5)
    weight = block_mean(solid_angle_weights(h, w), 4)
    dirs = map_icrf_dirs(w // 4, h // 4)
    for ra, dec, radius, _ in EXTERNAL:  # not part of the model
        weight[dirs @ unit_ra_dec(ra, dec) > math.cos(math.radians(radius * 1.3))] = 0.0
    best = None
    for k in np.geomspace(0.7, 2.8, 13):
        cal = {'dust_scale': float(k)}
        lf = luminosity_function(cal)
        light, _ = render(np.zeros(3), w, h, cal, lf)
        light = block_mean(light, 4)
        logs = np.log(real) - np.log(light)
        cal['map_per_light'] = [float(math.exp((logs[..., c] * weight).sum() / weight.sum())) for c in range(3)]
        model = light * np.array(cal['map_per_light'])
        r = np.log(luminance(model)) - np.log(luminance(real))
        cost = float(np.sqrt((r * r * weight).sum() / weight.sum()))
        print(f'dust scale {k:.3f}: rms log residual {cost:.4f}')
        if best is None or cost < best[0]:
            best = (cost, cal, model)
    cost, cal, model = best
    lf = luminosity_function(cal)
    cal['rms_log_residual'] = cost
    cal['j_sun'] = lf[2]
    cal['a_v_per_kpc_at_sun'] = 1.0857 * float(dust_disk_tau(np.array(R0), np.array(Z_SUN))) * cal['dust_scale']
    # The mean sky brightness from the Sun, rendered as the scenes are.
    light, tau = render(np.zeros(3), 512, 256, cal, lf)
    sun_sky = light * np.array(cal['map_per_light']) + external_galaxies(
        load_real_map(2), np.zeros(3), 512, 256, np.exp(-tau * EXT_RGB[1]))
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


# ---------------------------------------------------------------- point stars


def shell_candidates(viewer, r_min, r_max, j_sun, rng, n):
    """Points between r_min and r_max (kpc) around viewer, log-uniform in radius, with
    weights w such that the emissivity's volume integral is mean(w) * 4 pi ln(r_max/r_min)
    (kpc^3 * L_sun/pc^3) and points drawn with probability ~ w follow the emissivity."""
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    r = r_min * (r_max / r_min) ** rng.uniform(size=n)
    p = viewer + u * r[:, None]
    return p, emissivity(p, j_sun) * r ** 3


def expected_counts(viewer, lf, point_mag, r_min, rng):
    """Expected point stars per magnitude bin (no extinction): the luminosity
    function scaled by the emissivity relative to the local one."""
    phi, _, j_sun = lf
    counts = np.zeros(len(phi))
    for b, lo in enumerate(MAG_EDGES[:-1]):
        D = min(10 ** ((point_mag - lo) / 5 + 1), 15000.0) / 1000.0
        if phi[b] <= 0 or D <= r_min:
            continue
        _, w = shell_candidates(viewer, r_min, D, j_sun, rng, 20000)
        counts[b] = phi[b] * w.mean() * 4 * math.pi * math.log(D / r_min) / j_sun
    return counts


def choose_point_mag(viewer, lf, r_min, rng):
    """NAKED_EYE_MAG, or a brighter limit where that would give more than MAX_POINT_STARS."""
    if expected_counts(viewer, lf, NAKED_EYE_MAG, r_min, rng).sum() <= MAX_POINT_STARS:
        return NAKED_EYE_MAG
    lo, hi = -15.0, NAKED_EYE_MAG
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        if expected_counts(viewer, lf, mid, r_min, rng).sum() > MAX_POINT_STARS:
            hi = mid
        else:
            lo = mid
    return round(lo, 2)


def sample_stars(viewer, cal, lf, point_mag, r_min, rng):
    """Point stars brighter than point_mag seen from viewer (heliocentric Galactic, kpc),
    none closer than r_min: ICRF unit vectors, V and B-V."""
    phi, templates, j_sun = lf
    counts = expected_counts(viewer, lf, point_mag, r_min, rng)
    positions, mags, colors = [], [], []
    for b, lo in enumerate(MAG_EDGES[:-1]):
        n = rng.poisson(counts[b])
        if n == 0 or len(templates[b]) == 0:
            continue
        D = min(10 ** ((point_mag - lo) / 5 + 1), 15000.0) / 1000.0
        p, w = shell_candidates(viewer, r_min, D, j_sun, rng, max(50 * n, 50000))
        pick = rng.choice(len(p), size=n, p=w / w.sum())
        t = templates[b][rng.integers(len(templates[b]), size=n)]
        positions.append(p[pick])
        mags.append(t[:, 0])
        colors.append(t[:, 1])
    p = np.concatenate(positions)
    M = np.concatenate(mags)
    bv0 = np.concatenate(colors)
    rel = p - viewer
    dist = np.linalg.norm(rel, axis=1)
    dirs = rel / dist[:, None]
    av = np.concatenate([extinction_v(viewer, dirs[i:i + 2048], dist[i:i + 2048], cal)
                         for i in range(0, len(dirs), 2048)])
    V = M + 5 * np.log10(dist * 1000.0 / 10.0) + av
    bv = bv0 + av / R_V
    keep = V <= point_mag
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
    """The map (real map units), the point stars and the point/diffuse limit."""
    viewer_icrf, viewer = viewer_from(*viewer_ra_dec_pc)
    lf = luminosity_function(cal)
    near_centre = np.linalg.norm(viewer - NUCLEUS) < 0.5
    # At the centre: stars from 0.2 pc out (closer ones would show parallax as the
    # camera moves around Sgr A*), rays from 0.02 pc.
    r_min = 0.0002 if near_centre else 0.0005
    point_mag = choose_point_mag(viewer, lf, r_min, rng)
    resolved = RESOLVED_MAG if point_mag >= NAKED_EYE_MAG else point_mag
    light, tau = render(viewer, width, height, cal, lf, resolved,
                        s_min=2e-5 if near_centre else 0.001, steps=256 if near_centre else 192)
    model = light * np.array(cal['map_per_light'])
    model += external_galaxies(load_real_map(1), viewer_icrf, width, height, np.exp(-tau * EXT_RGB[1]))
    dirs, V, bv = sample_stars(viewer, cal, lf, point_mag, r_min, rng)
    return model, dirs, V, bv, point_mag


def brightness_for(model, cal):
    """Displayed brightness (Earth's default, with the overall level compressed like
    the stars' fluxes, relative to the model sky at the Sun), the magnitude offset
    with which the scene's stars are shown (> 0 when the eye adapts), and the ratio."""
    ratio = mean_luminance(model) / cal['mean_luminance_sun']
    level = ratio ** MAG_GAMMA
    exposure = min(1.0, SKY_LEVEL_CAP / level)
    mag_offset = max(0.0, -2.5 * math.log10(exposure) / MAG_GAMMA)  # the same factor on the stars' compressed flux
    return DEFAULT_BRIGHTNESS * ratio ** (MAG_GAMMA - 1.0) * exposure, mag_offset, ratio


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
    model, dirs, V, bv, point_mag = sky_for(viewer, width, height, cal, rng)
    white = float(np.percentile(model.max(axis=2), 99.99))
    brightness, mag_offset, ratio = brightness_for(model, cal)
    write_map(os.path.join(ASSETS, sky['milky_way']), model, white)
    header = [
        f'Model sky for {scene.get("name", path)}: stars brighter than V = {point_mag} seen from the viewer.',
        'Generated by tools/sky/make_galaxy_sky.py (see assets/sky/SOURCES.md); statistical, not real stars.',
        f'viewer_ra_dec_distance_pc = {viewer[0]}, {viewer[1]}, {viewer[2]}',
        f'milky_way_brightness = {brightness * white:.4f}',
        f'display_mag_offset = {mag_offset:.3f}',
    ]
    write_stars(os.path.join(ASSETS, sky['stars']), dirs, V, bv, header)
    print(f'{scene.get("name", path)}: {len(V)} stars to V = {point_mag} (brightest V = {V.min():.2f}); '
          f'mean sky {ratio:.3g} x the Sun\'s; milky_way_brightness = {brightness * white:.4f}; '
          f'stars shown {mag_offset:.2f} mag fainter')


def sun(out_dir):
    cal = load_calibration()
    rng = np.random.default_rng(NOISE_SEED)
    model, dirs, V, bv, _ = sky_for((0.0, 0.0, 0.0), 1024, 512, cal, rng)
    b = np.degrees(np.arcsin((dirs @ G2I)[:, 2]))
    real = load_hyg()
    real_b = np.degrees(np.arcsin((unit_ra_dec(real[:, 0], real[:, 1]) @ G2I)[:, 2]))
    real_v = real[:, 3] <= NAKED_EYE_MAG
    print(f'model from the Sun: {len(V)} stars to V = {NAKED_EYE_MAG}, {np.mean(np.abs(b) < 10):.3f} within 10 deg '
          f'of the plane, {np.sum(V < 2):.0f} brighter than V = 2')
    print(f'real (hyg.csv with distances): {real_v.sum()} stars, '
          f'{np.mean(np.abs(real_b[real_v]) < 10):.3f} within 10 deg, {np.sum(real[:, 3] < 2):.0f} brighter than V = 2')
    if out_dir:
        np.save(os.path.join(out_dir, 'model_from_sun.npy'), model)


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
