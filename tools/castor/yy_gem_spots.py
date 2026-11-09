"""Starspots of YY Geminorum from the TESS light curves, for castor.toml.

Usage:
    python -E -P tools/castor/yy_gem_spots.py <TESS SPOC *_lc.fits> ...

Needs numpy; a few minutes. One set of spots per sector, in the
spot model of Torres & Ribas 2002 (circles at latitude +45 deg, 10% cooler than
the photosphere), which the scene draws: only the meridian phase and the size
of each spot are fitted.

The aperture is mostly the light of Castor A and B (71" away, and saturated):
the eclipses are only 1-2% deep instead of 39%, by an amount that changes
between sectors. PDCSAP does not know this (CROWDSAP is about 0.98). The folded
curve is therefore scaled to YY Gem's own light by the depth of its eclipses:
the model's depth over the measured one.

Only the phases between the eclipses are fitted (0.09 from either), with the
ellipsoidal variation taken out (Morris 1985, about 0.5%), and a free level: the
light curve tells how spots differ between phases, not how much is always covered.
It cannot tell which star a spot is on either: spots are given to the two
stars in turn, the largest on Cb, which is the assignment printed.
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'algol'))
from algol_ephemeris import read_light_curve  # noqa: E402

T0, P = 2459690.57515, 0.81428247  # BJD_TDB, d: yy_gem_ephemeris.py
BTJD = 2457000.0
WAVELENGTH_NM = 786.5              # TESS band center (Ricker et al. 2015: 600-1000 nm)
HC_OVER_K = 1.438776877e7          # nm K
T_STAR = 3820.0                    # Torres & Ribas 2002
SPOT_RATIO = 0.90                  # Table 4, fixed
SPOT_LAT_DEG = 45.0                # Table 4, fixed
INCLINATION_DEG = 86.29
LIMB = 0.6                         # as the scene draws the stars
FRACTIONAL_R = 0.1589              # r / a, Torres & Ribas 2002
RINGS, POINTS = 24, 48
NBINS = 80
ECLIPSE_HALF = 0.09
ELLIPSOIDAL = 0.15 * (15.0 + LIMB) * 1.32 / (3.0 - LIMB) * FRACTIONAL_R ** 3 * math.sin(math.radians(INCLINATION_DEG)) ** 2
RADIUS_MIN, RADIUS_MAX = 4.0, 40.0


def planck(t):
    return 1.0 / np.expm1(HC_OVER_K / (WAVELENGTH_NM * t))


CONTRAST = 1.0 - planck(SPOT_RATIO * T_STAR) / planck(T_STAR)


def eclipse_depth():
    """Light lost at minimum by two equal limb-darkened disks (fraction of the pair's)."""
    separation = math.cos(math.radians(INCLINATION_DEG)) / FRACTIONAL_R  # in stellar radii
    n = 800
    covered, total = 0.0, 0.0
    for i in range(n):
        x = (i + 0.5) / n
        mu = math.sqrt(1.0 - x * x)
        weight = x * (1.0 - LIMB * (1.0 - mu))
        total += weight
        if separation + x <= 1.0:
            share = 1.0
        elif separation >= x + 1.0:
            share = 0.0
        else:
            share = math.acos(np.clip((x * x + separation * separation - 1.0) / (2.0 * x * separation), -1.0, 1.0)) / math.pi
        covered += weight * share
    return 0.5 * covered / total


def spot_wave(phases, meridian, radius_deg):
    """Light lost to one spot, as a fraction of the pair's, at the given phases."""
    rho = math.radians(radius_deg)
    lat = math.radians(SPOT_LAT_DEG)
    center = np.array([math.cos(lat), 0.0, math.sin(lat)])     # longitude 0: facing the observer at `meridian`
    u = np.array([0.0, 1.0, 0.0])
    v = np.cross(center, u)
    d_theta = rho / RINGS
    theta = (np.arange(RINGS) + 0.5) * d_theta
    area = (np.cos(np.arange(RINGS) * d_theta) - np.cos((np.arange(RINGS) + 1) * d_theta)) * 2.0 * np.pi / POINTS
    phi = 2.0 * np.pi * (np.arange(POINTS) + 0.5) / POINTS
    t_grid, p_grid = np.meshgrid(theta, phi, indexing='ij')
    normals = (center[None, None, :] * np.cos(t_grid)[..., None]
               + (u[None, None, :] * np.cos(p_grid)[..., None] + v[None, None, :] * np.sin(p_grid)[..., None]) * np.sin(t_grid)[..., None])
    areas = np.repeat(area, POINTS).reshape(RINGS, POINTS)
    beta = math.radians(90.0 - INCLINATION_DEG)                # latitude of the observer on the star
    longitudes = 2.0 * np.pi * (phases - meridian)             # the observer's longitude from the spot's
    observer = np.stack([np.cos(beta) * np.cos(longitudes), np.cos(beta) * np.sin(longitudes), np.full_like(longitudes, np.sin(beta))], axis=1)
    mu = normals.reshape(-1, 3) @ observer.T
    mu = np.clip(mu, 0.0, None)
    seen = (areas.ravel()[None, :] * (mu.T * (1.0 - LIMB * (1.0 - mu.T)))).sum(axis=1)
    return 0.5 * CONTRAST * seen / (np.pi * (1.0 - LIMB / 3.0))


def model(params, phases):
    level, spots = params[0], params[1:]
    light = np.full_like(phases, level)
    for k in range(0, len(spots), 2):
        light = light - spot_wave(phases, spots[k], spots[k + 1])
    return light


def fit(phases, y, sigma, start):
    p = np.array(start, dtype=float)
    damping = 1e-2

    def clamp(q):
        q = q.copy()
        for k in range(2, len(q), 2):
            q[k] = np.clip(q[k], RADIUS_MIN, RADIUS_MAX)
        return q

    def cost(q):
        return np.sum(((y - model(q, phases)) / sigma) ** 2)

    best = cost(p)
    for _ in range(40):
        r = (y - model(p, phases)) / sigma
        jac = np.zeros((len(y), len(p)))
        for j in range(len(p)):
            dq = np.zeros_like(p)
            dq[j] = 1e-3 if j == 0 or j % 2 == 1 else 0.05
            jac[:, j] = (model(clamp(p + dq), phases) - model(clamp(p - dq), phases)) / (2.0 * dq[j] * sigma)
        a = jac.T @ jac
        b = jac.T @ r
        try:
            step = np.linalg.solve(a + damping * np.diag(np.diag(a) + 1e-9), b)
        except np.linalg.LinAlgError:
            break
        trial = clamp(p + step)
        c = cost(trial)
        if c < best:
            p, best = trial, c
            damping = max(damping * 0.4, 1e-6)
        else:
            damping *= 4.0
            if damping > 1e6:
                break
    return best, p


def best_fit(phases, y, sigma, count, rng):
    best = (1e30, None)
    starts = 40
    for _ in range(starts):
        start = [1.0]
        for k in range(count):
            start += [rng.uniform(0.0, 1.0), rng.uniform(10.0, 25.0)]
        c, p = fit(phases, y, sigma, start)
        if c < best[0]:
            best = (c, p)
    return best


def sector_data(path):
    lc = read_light_curve(path, ('TIME', 'PDCSAP_FLUX', 'QUALITY'))
    t, f, q = lc['TIME'] + BTJD, lc['PDCSAP_FLUX'], lc['QUALITY']
    good = np.isfinite(t) & np.isfinite(f) & (q == 0)
    t, f = t[good], f[good]
    phase = ((t - T0) / P) % 1.0
    # A sector's slow trend (Castor's saturated light): a running median of the
    # light between the eclipses, one orbit wide.
    between = (np.abs(phase - 0.5) > ECLIPSE_HALF) & (phase > ECLIPSE_HALF) & (phase < 1.0 - ECLIPSE_HALF)
    knots = np.arange(t[0], t[-1] + P, P)
    centers, levels = [], []
    for a, b in zip(knots[:-1], knots[1:]):
        mask = between & (t >= a) & (t < b)
        if mask.sum() > 50:
            centers.append(0.5 * (a + b))
            levels.append(np.median(f[mask]))
    y = f / np.interp(t, centers, levels)
    primary = (phase < 0.004) | (phase > 0.996)
    secondary = np.abs(phase - 0.5) < 0.004
    depth = 0.5 * (1.0 - np.median(y[primary]) + 1.0 - np.median(y[secondary]))
    dilution = depth / eclipse_depth()
    light = 1.0 + (y - 1.0) / dilution
    light += ELLIPSOIDAL * np.cos(4.0 * np.pi * phase)   # the pair is faintest at conjunction
    edges = np.linspace(0.0, 1.0, NBINS + 1)
    middle, values, errors = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        c = 0.5 * (a + b)
        if abs(c - 0.5) < ECLIPSE_HALF or c < ECLIPSE_HALF or c > 1.0 - ECLIPSE_HALF:
            continue
        v = light[(phase >= a) & (phase < b)]
        spread = 1.4826 * np.median(np.abs(v - np.median(v)))
        v = v[np.abs(v - np.median(v)) < 3.5 * spread]               # flares
        middle.append(c)
        values.append(np.mean(v))
        errors.append(np.std(v) / math.sqrt(max(len(v) / 8.0, 1.0))) # the points are correlated
    jd = 0.5 * (t[0] + t[-1])
    return jd, dilution, np.array(middle), np.array(values), np.array(errors)


def julian_year(jd):
    return 2000.0 + (jd - 2451545.0) / 365.25


def main(paths):
    print('model eclipse depth %.3f, spot contrast in the TESS band %.3f, ellipsoidal %.4f'
          % (eclipse_depth(), CONTRAST, ELLIPSOIDAL))
    rng = np.random.default_rng(1)
    results = []
    for path in paths:
        jd, dilution, phases, y, sigma = sector_data(path)
        name = os.path.basename(path).split('-')[1]
        print('\n%s  year %.2f  dilution %.4f  wave p2p %.4f  error per bin %.4f'
              % (name, julian_year(jd), dilution, np.ptp(y), np.median(sigma)))
        fits = {}
        for count in (1, 2):
            c, p = best_fit(phases, y, sigma, count, rng)
            fits[count] = (c, p)
            bic = c + (1 + 2 * count) * math.log(len(y))
            spots = ', '.join('%.3f %.1f' % (p[1 + 2 * k] % 1.0, p[2 + 2 * k]) for k in range(count))
            print('  %d spots: chi2 %.1f (%d bins), rms %.4f, BIC %.1f : %s'
                  % (count, c, len(y), np.std(y - model(p, phases)), bic, spots))
        results.append((name, jd, phases, y, sigma, fits))

    print('\nFor castor.toml (the count BIC prefers):')
    for i, (name, jd, phases, y, sigma, fits) in enumerate(results):
        bics = {k: fits[k][0] + (1 + 2 * k) * math.log(len(y)) for k in fits}
        count = min(bics, key=bics.get)
        c, p = fits[count]
        spots = sorted([(p[2 + 2 * k], p[1 + 2 * k] % 1.0) for k in range(count)], reverse=True)
        print('# %s (year %.2f): %d spots, rms %.4f' % (name, julian_year(jd), count, np.std(y - model(p, phases))))
        for k, (radius, phase) in enumerate(spots):
            star = 'Cb' if k % 2 == 0 else 'Ca'
            lon = (180.0 - 360.0 * phase) % 360.0 if star == 'Cb' else (-360.0 * phase) % 360.0
            print('  %s: phase %.3f radius %.1f  longitude_deg = %.2f' % (star, phase, radius, lon))
        grid = (phases - 0.0)
        step = max(1, len(grid) // 10)
        sample = ', '.join('{%.4f, %.4f}' % (grid[j], y[j] - p[0]) for j in range(0, len(grid), step))
        print('  points (phase, light - level): ' + sample)


if __name__ == '__main__':
    main(sys.argv[1:])
