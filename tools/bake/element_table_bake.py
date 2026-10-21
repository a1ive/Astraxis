"""Bake tables of slowly varying orbital elements (.elm) for moons close to an
oblate planet, where Kepler-relative ephemerides (horizons_bake.py) need a knot
every few hours: Janus and Epimetheus, Saturn's co-orbital pair.

Usage (numpy is needed; run from the repository root):
    python -E -P tools/bake/element_table_bake.py [name ...]

For each target below, hourly JPL Horizons positions (ICRF, relative to the
planet's center) are downloaded once into tools/bake/cache/ (not in git). In the
reference plane given by its ICRF pole (x: the plane's ascending node on the ICRF
equator, as frame_from_pole in src/core/math.hpp), the orbit is modelled as a
precessing ellipse:

    M = lambda - varpi,  omega = varpi - node,  E - e sin E = M,
    r = Rz(node) Rx(i) Rz(omega) (a (cos E - e), a sqrt(1 - e^2) sin E, 0)

with a free semi-major axis a and a free mean motion n = dlambda/dt (Saturn's
oblateness makes the orbit faster than Kepler's third law says), and the apsidal
and nodal precession rates fixed for each moon (fitted to the whole span). Over a
window of `WINDOW_DAYS` around every `STEP_DAYS`, a, e, varpi, i, node, lambda and
n are fitted by least squares (Gauss-Newton, all windows at once); the series of
fits is then thinned greedily, keeping a knot only where interpolating between
its neighbours would miss the hourly positions by more than `TOLERANCE_KM`. The
scene interpolates a, e, varpi, i and node linearly and lambda as a cubic Hermite
curve through lambda and n (src/ephem/element_table.cpp, mirrored below).

The model leaves out the short-period terms of the oblateness at order e J2 and
the other moons' short-period pull: the fit residual (printed) sets the floor of
the tolerance.

Output format (little endian), see src/ephem/element_table.hpp:
    char[8]  magic "AXELM1\\0\\0"
    int32    target NAIF id, int32 center NAIF id, uint32 knot count, uint32 0
    float64  reference pole RA, Dec (deg), apsidal and nodal rate (rad/s)
    knots:   float64 t (TDB seconds since J2000), a (km), e, varpi, i, node,
             lambda (rad; varpi, node and lambda unwrapped), n (rad/s)
"""

import math
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import horizons_bake as hb  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
OUT = os.path.join(HERE, '..', '..', 'assets', 'ephem')

# [ELEM] (JPL SSD Planetary Satellite Mean Elements, SAT415): Laplace pole of Janus
# and Epimetheus.
POLE_RA_DEC = (40.6, 83.5)
START, STOP = '1972-01-01', '2050-01-01'  # Horizons (sat415) covers 1949-12-27 .. 2050-01-09
STEP_MINUTES = 60
WINDOW_DAYS = 2.0
STEP_DAYS = 0.25
TOLERANCE_KM = 10.0

TARGETS = {
    'janus': {'command': '610', 'naif_id': 610},
    'epimetheus': {'command': '611', 'naif_id': 611},
}
CENTER, CENTER_ID = '500@699', 699
GM = 37940584.8418  # Saturn system [GM], only for the initial guess


# --- Data ---------------------------------------------------------------------

def samples(command):
    """Hourly (t, positions[N, 3], velocities[N, 3]) over START .. STOP, cached."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, f'{command}_{START}_{STOP}_{STEP_MINUTES}m.npz')
    if not os.path.exists(path):
        start_jd, stop_jd = hb.jd_from_text(START), hb.jd_from_text(STOP)
        print(f'    fetching {command} from Horizons ...')
        rows = hb.fetch(command, CENTER, start_jd, stop_jd, STEP_MINUTES)
        t = np.array([r[0] for r in rows])
        p = np.array([r[1] for r in rows])
        v = np.array([r[2] for r in rows])
        np.savez(path, t=t, p=p, v=v)
    d = np.load(path)
    return d['t'], d['p'], d['v']


def reference_frame():
    ra, dec = (math.radians(x) for x in POLE_RA_DEC)
    z = np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)])
    x = np.cross([0.0, 0.0, 1.0], z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return np.array([x, y, z])  # rows: ICRF -> plane


# --- Model --------------------------------------------------------------------

def solve_kepler(m, e):
    big_e = m.copy()
    for _ in range(8):
        big_e -= (big_e - e * np.sin(big_e) - m) / (1.0 - e * np.cos(big_e))
    return big_e


def ellipse(a, e, varpi, i, node, lam):
    """Positions in the reference plane frame (arrays broadcast)."""
    m = lam - varpi
    big_e = solve_kepler(np.remainder(m + math.pi, 2 * math.pi) - math.pi, e)
    xp = a * (np.cos(big_e) - e)
    yp = a * np.sqrt(1.0 - e * e) * np.sin(big_e)
    w = varpi - node
    cw, sw, cn, sn, ci, si = np.cos(w), np.sin(w), np.cos(node), np.sin(node), np.cos(i), np.sin(i)
    x1, y1 = cw * xp - sw * yp, sw * xp + cw * yp  # Rz(omega)
    y2, z2 = ci * y1, si * y1                     # Rx(i)
    return np.stack([cn * x1 - sn * y2, sn * x1 + cn * y2, z2], axis=-1)


def window_model(params, dt, rates):
    """params[..., 7] = a, k = e cos varpi, h = e sin varpi, q = i cos node, p = i sin node,
    lambda, n at the window center; dt[..., S] seconds from it."""
    a, k, h, q, p, lam, n = (params[..., j:j + 1] for j in range(7))
    e = np.hypot(k, h)
    varpi = np.arctan2(h, k) + rates[0] * dt
    i = np.hypot(q, p)
    node = np.arctan2(p, q) + rates[1] * dt
    return ellipse(a, e, varpi, i, node, lam + n * dt)


def fit_chunk(t, pos, vel, rates, centers, half):
    """Gauss-Newton fits for the windows around `centers` (sample indices)."""
    idx = centers[:, None] + np.arange(-half, half + 1)[None, :]
    dt = t[idx] - t[centers][:, None]
    target = pos[idx]

    # Initial guess: osculating elements at the window centers.
    r, v = pos[centers], vel[centers]
    rn = np.linalg.norm(r, axis=1)
    a0 = 1.0 / (2.0 / rn - np.einsum('ij,ij->i', v, v) / GM)
    hvec = np.cross(r, v)
    ecc = np.cross(v, hvec) / GM - r / rn[:, None]
    lam0 = np.arctan2(r[:, 1], r[:, 0])
    varpi0 = np.arctan2(ecc[:, 1], ecc[:, 0])
    e0 = np.hypot(ecc[:, 0], ecc[:, 1])
    hn = hvec / np.linalg.norm(hvec, axis=1)[:, None]
    i0 = np.arccos(np.clip(hn[:, 2], -1, 1))
    node0 = np.arctan2(hn[:, 0], -hn[:, 1])
    n0 = np.full_like(a0, 2 * math.pi / (0.6946 * 86400))
    params = np.stack([a0, e0 * np.cos(varpi0), e0 * np.sin(varpi0), i0 * np.cos(node0), i0 * np.sin(node0),
                       lam0, n0], axis=1)
    scale = np.array([1.0, 1e-4, 1e-4, 1e-5, 1e-5, 1e-5, 1e-10])  # finite-difference steps
    for _ in range(15):
        model = window_model(params, dt, rates)
        res = (target - model).reshape(len(centers), -1)
        jac = np.empty(res.shape + (7,))
        for j in range(7):
            d = params.copy()
            d[:, j] += scale[j]
            jac[..., j] = (window_model(d, dt, rates) - model).reshape(len(centers), -1) / scale[j]
        jtj = np.einsum('wsi,wsj->wij', jac, jac)
        jtr = np.einsum('wsi,ws->wi', jac, res)
        delta = np.linalg.solve(jtj, jtr[..., None])[..., 0]
        params += delta
        if np.max(np.abs(delta[:, 0])) < 1e-6 and np.max(np.abs(delta[:, 5])) < 1e-11:
            break
    model = window_model(params, dt, rates)
    rms = np.sqrt(np.mean(np.sum((target - model) ** 2, axis=2), axis=1))
    return params, rms


def fit_windows(t, pos, vel, rates):
    """Least-squares elements over a window around every STEP_DAYS."""
    step = int(round(STEP_DAYS * 1440 / STEP_MINUTES))
    half = int(round(WINDOW_DAYS / 2 * 1440 / STEP_MINUTES))
    centers = np.arange(half, len(t) - half, step)
    parts = [fit_chunk(t, pos, vel, rates, centers[c:c + 4000], half) for c in range(0, len(centers), 4000)]
    params = np.concatenate([p for p, _ in parts])
    rms = np.concatenate([r for _, r in parts])
    a, k, h, q, p, lam, n = params.T
    out = {
        't': t[centers], 'a': a, 'e': np.hypot(k, h), 'varpi': np.unwrap(np.arctan2(h, k)),
        'i': np.hypot(q, p), 'node': np.unwrap(np.arctan2(p, q)), 'n': n, 'rms': rms,
    }
    # Unwrap lambda by integrating n between the centers.
    lam = np.remainder(lam, 2 * math.pi)
    unwrapped = [lam[0]]
    for j in range(1, len(lam)):
        guess = unwrapped[-1] + 0.5 * (n[j - 1] + n[j]) * (out['t'][j] - out['t'][j - 1])
        unwrapped.append(guess + math.remainder(lam[j] - guess, 2 * math.pi))
    out['lambda'] = np.array(unwrapped)
    return out


# --- Table --------------------------------------------------------------------

def interpolate(knots, t):
    """Elements at times t between knots (dict of arrays), as element_table.cpp."""
    kt = knots['t']
    j = np.clip(np.searchsorted(kt, t, side='right') - 1, 0, len(kt) - 2)
    t0, t1 = kt[j], kt[j + 1]
    u = (t - t0) / (t1 - t0)
    lin = {key: knots[key][j] + u * (knots[key][j + 1] - knots[key][j]) for key in ('a', 'e', 'varpi', 'i', 'node')}
    h = t1 - t0
    u2, u3 = u * u, u * u * u
    lam = ((2 * u3 - 3 * u2 + 1) * knots['lambda'][j] + (u3 - 2 * u2 + u) * h * knots['n'][j] +
           (-2 * u3 + 3 * u2) * knots['lambda'][j + 1] + (u3 - u2) * h * knots['n'][j + 1])
    return ellipse(lin['a'], lin['e'], lin['varpi'], lin['i'], lin['node'], lam)


def subset(series, idx):
    return {key: val[idx] for key, val in series.items()}


def thin(series, t, plane_pos, tol):
    """Greedy: from each kept knot, the farthest next one whose interpolation stays within tol."""
    n = len(series['t'])
    s_idx = np.searchsorted(t, series['t'])  # sample index of every candidate knot

    def ok(a, b):
        lo, hi = s_idx[a], s_idx[b]
        tt = t[lo:hi + 1]
        model = interpolate(subset(series, np.array([a, b])), tt)
        return np.max(np.linalg.norm(model - plane_pos[lo:hi + 1], axis=1)) <= tol

    kept = [0]
    a = 0
    while a < n - 1:
        step = 1
        while a + step * 2 <= n - 1 and ok(a, a + step * 2):
            step *= 2
        lo, hi = a + step, min(n - 1, a + step * 2)
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if ok(a, mid):
                lo = mid
            else:
                hi = mid - 1
        kept.append(lo)
        a = lo
    return np.array(kept)


def write_elm(path, naif_id, knots, rates):
    with open(path, 'wb') as f:
        f.write(b'AXELM1\0\0')
        f.write(struct.pack('<iiII', naif_id, CENTER_ID, len(knots['t']), 0))
        f.write(struct.pack('<4d', POLE_RA_DEC[0], POLE_RA_DEC[1], rates[0], rates[1]))
        for j in range(len(knots['t'])):
            f.write(struct.pack('<8d', *(knots[key][j] for key in ('t', 'a', 'e', 'varpi', 'i', 'node', 'lambda', 'n'))))


def bake(name, spec):
    print(f'{name}: {spec["command"]} @ {CENTER}, {START} .. {STOP}, every {STEP_MINUTES} min')
    t, pos, vel = samples(spec['command'])
    rot = reference_frame()
    plane_pos = pos @ rot.T
    plane_vel = vel @ rot.T
    print(f'    {len(t)} samples')

    # Pass 1 with no precession (short windows); pass 2 with the rates fitted to pass 1.
    rates = (0.0, 0.0)
    for pass_no in (1, 2):
        series = fit_windows(t, plane_pos, plane_vel, rates)
        rates = (np.polyfit(series['t'], series['varpi'], 1)[0], np.polyfit(series['t'], series['node'], 1)[0])
        print(f'    pass {pass_no}: {len(series['t'])} windows, rms {np.median(series["rms"]):.2f} km median, '
              f'{np.max(series["rms"]):.2f} max; apsides {math.degrees(rates[0]) * 86400:+.4f} deg/day, '
              f'node {math.degrees(rates[1]) * 86400:+.4f} deg/day')
    kept = thin(series, t, plane_pos, TOLERANCE_KM)
    knots = subset(series, kept)
    inside = (t >= knots['t'][0]) & (t <= knots['t'][-1])
    err = np.linalg.norm(interpolate(knots, t[inside]) - plane_pos[inside], axis=1)
    path = os.path.join(OUT, f'{name}_sat.elm')
    write_elm(path, spec['naif_id'], knots, rates)
    print(f'    {len(series["t"])} fits -> {len(kept)} knots, max error {np.max(err):.2f} km, '
          f'rms {np.sqrt(np.mean(err ** 2)):.2f} km; wrote {os.path.relpath(path)} '
          f'({os.path.getsize(path) // 1024} KB)')
    return series, knots


def main():
    names = sys.argv[1:] or list(TARGETS)
    for name in names:
        bake(name, TARGETS[name])


if __name__ == '__main__':
    main()
