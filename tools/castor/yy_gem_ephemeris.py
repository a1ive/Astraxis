"""Linear ephemeris of YY Geminorum (Castor C) from TESS light curves.

Usage:
    python -E -P tools/castor/yy_gem_ephemeris.py <TESS SPOC *_lc.fits> ...

Needs numpy. Reads the SPOC 2-minute light curves of TIC 239187696 (addresses in
assets/ephem/SOURCES.md) with the reader of tools/algol, finds each minimum as
the time of mirror symmetry of the SAP flux, after removing a straight line
through the flux on either side of the eclipse (the spots), and fits a line to
the primary minima. Which minimum is primary comes from the ephemeris of Torres
& Ribas 2002 (the two are almost equally deep). Flares in the eclipse spoil the
symmetry: minima far from the first fit are dropped, and it is fitted again.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'algol'))
from algol_ephemeris import read_light_curve  # noqa: E402

TR02 = (2449345.112327, 0.814282212)  # HJD, d: Torres & Ribas 2002, ApJ 567, 1140, Eq. 1
BTJD = 2457000.0
HALF = 0.05          # the eclipse lasts about 0.08 d; d
SIDE = 0.03          # flux beside the eclipse for the line; d


def minimum_time(t, f, guess):
    near = (np.abs(t - guess) < HALF + SIDE)
    t, f = t[near], f[near]
    if t.size < 0.9 * 2.0 * (HALF + SIDE) * 720.0:
        return None
    side = np.abs(t - guess) > HALF
    line = np.polyfit(t[side] - guess, f[side], 1)
    f = f / np.polyval(line, t - guess)
    grid = np.linspace(guess - 0.01, guess + 0.01, 401)
    cost = []
    tt = np.linspace(0.002, HALF - 0.006, 40)
    for c in grid:
        cost.append(np.sum((np.interp(c + tt, t, f) - np.interp(c - tt, t, f)) ** 2))
    k = int(np.argmin(cost))
    if k == 0 or k == grid.size - 1:
        return None
    a, b, _ = np.polyfit(grid[k - 1:k + 2] - grid[k], np.array(cost[k - 1:k + 2]), 2)
    return grid[k] - b / (2.0 * a)


def main(paths):
    primary, secondary = [], []
    for path in paths:
        lc = read_light_curve(path)
        t, f, q = lc['TIME'], lc['SAP_FLUX'], lc['QUALITY']
        good = (q == 0) & np.isfinite(t) & np.isfinite(f)
        t, f = t[good] + BTJD, f[good]
        first = int(np.ceil((t[0] - TR02[0]) / (0.5 * TR02[1])))
        last = int(np.floor((t[-1] - TR02[0]) / (0.5 * TR02[1])))
        for half_cycle in range(first, last + 1):
            when = minimum_time(t, f, TR02[0] + 0.5 * half_cycle * TR02[1])
            if when is not None:
                (secondary if half_cycle % 2 else primary).append((half_cycle / 2.0, when))
    p = np.array(primary)
    s = np.array(secondary)
    fit = np.polyfit(p[:, 0], p[:, 1], 1)
    for _ in range(4):
        keep = np.abs(p[:, 1] - np.polyval(fit, p[:, 0])) < 2.5 * np.std(p[:, 1] - np.polyval(fit, p[:, 0]))
        p = p[keep]
        fit = np.polyfit(p[:, 0], p[:, 1], 1)
    period = fit[0]
    mean = int(round(p[:, 0].mean()))
    t0 = np.polyval(fit, mean)
    residual = (p[:, 1] - np.polyval(fit, p[:, 0])) * 1440.0
    print('%d primary minima (of %d), rms %.2f min' % (len(p), len(primary), residual.std()))
    print('primary minimum BJD_TDB %.5f, P %.8f d' % (t0, period))
    print('TR02 at that time: O-C %.2f min (HJD against BJD_TDB is about 1.1 min)'
          % ((t0 - TR02[0] - mean * TR02[1]) * 1440.0))
    keep = np.abs(s[:, 1] - np.polyval(fit, s[:, 0])) * 1440.0 < 5.0
    print('secondary minima (%d of %d): mean phase offset from 0.5: %.2f min'
          % (keep.sum(), len(s), ((s[keep, 1] - np.polyval(fit, s[keep, 0])) * 1440.0).mean()))
    for t_low, t_high, name in [(2458800, 2458900, 's20'), (2459500, 2459600, 's46'),
                                (2459900, 2460000, 's60'), (2460250, 2460350, 's72')]:
        one = (p[:, 1] > t_low) & (p[:, 1] < t_high)
        if one.any():
            print('  %s: %d minima, mean residual %+.2f min' % (name, one.sum(), residual[one].mean()))
    # The next primary minima, for the scene's events.
    jd = 2461352.5   # 2026-11-08 00:00 TDB
    cycle = int(np.ceil((jd - t0) / period))
    for k in range(cycle, cycle + 4):
        print('next primary BJD %.5f' % (t0 + (k) * period))


if __name__ == '__main__':
    main(sys.argv[1:])
