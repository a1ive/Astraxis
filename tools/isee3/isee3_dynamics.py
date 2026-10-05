"""Force models, integrator and least-squares arc fitting for the ISEE-3/ICE reconstruction.

numpy only. Units: km, s, TDB seconds since J2000; ICRF axes.
"""
import datetime as dt
import os

import numpy as np

J2000 = dt.datetime(2000, 1, 1, 12, 0, 0)
DAY = 86400.0
AU = 149597870.7  # km (IAU 2012)

# GM (km^3/s^2): JPL SSD astrodynamic parameters (DE440), https://ssd.jpl.nasa.gov/astro_par.html
GM = {
    'sun': 1.32712440041279419e11,
    'earth': 398600.435507,
    'moon': 4902.800118,
    'mercury': 22031.868551,
    'venus': 324858.592,
    'mars': 42828.375816,      # system
    'jupiter': 126712764.1,    # system
    'saturn': 37940584.8418,   # system
    'uranus': 5794556.4,       # system
    'neptune': 6836527.10058,  # system
}
EARTH_RE = 6378.1366           # km, NAIF pck00011
EARTH_J2 = 1.08262668e-3       # EGM2008 C20 (unnormalized); only matters on the 1983 perigee passes

# TAI - UTC (s), IERS Bulletin C (USNO tai-utc.dat): (year, month, day, value from that date)
LEAP_SECONDS = [(1978, 1, 1, 17), (1979, 1, 1, 18), (1980, 1, 1, 19), (1981, 7, 1, 20), (1982, 7, 1, 21),
                (1983, 7, 1, 22), (1985, 7, 1, 23), (1988, 1, 1, 24), (1990, 1, 1, 25), (1991, 1, 1, 26),
                (1992, 7, 1, 27), (1993, 7, 1, 28), (1994, 7, 1, 29), (1996, 1, 1, 30), (1997, 7, 1, 31),
                (1999, 1, 1, 32), (2006, 1, 1, 33), (2009, 1, 1, 34), (2012, 7, 1, 35), (2015, 7, 1, 36)]


def utc_to_tdb(d):
    """UTC datetime -> TDB seconds since J2000 (TDB = TT = TAI + 32.184 s, to 2 ms)."""
    leap = 0
    for y, m, dd, n in LEAP_SECONDS:
        if d >= dt.datetime(y, m, dd):
            leap = n
    return (d - J2000).total_seconds() + leap + 32.184


def tdb_to_utc(t):
    d = J2000 + dt.timedelta(seconds=float(t) - 32.184)
    for _ in range(2):   # the leap second count depends on the date itself
        d = d - dt.timedelta(seconds=utc_to_tdb(d) - float(t))
    return d


_ARCSEC = np.pi / 180.0 / 3600.0
_OBLIQUITY_J2000 = np.radians(23.4392911)   # IAU 1980 mean obliquity at J2000


def nutation_vector(t):
    """Small-angle rotation vector (rad, ICRF/J2000 equatorial axes) of the nutation at TDB
    seconds t: the turn from the mean to the true equator and equinox of date. Main terms of
    the IAU 1980 series (Meeus, Astronomical Algorithms, 2nd ed., ch. 22), good to ~0.5".
    SSCWeb's "GEI J2000" ISEE-3 positions turn out to carry it (see isee3_reconstruct.py)."""
    T = np.asarray(t) / (36525.0 * DAY)
    L = np.radians(280.4665 + 36000.7698 * T)        # mean longitude of the Sun
    Lp = np.radians(218.3165 + 481267.8813 * T)      # mean longitude of the Moon
    Om = np.radians(125.04452 - 1934.136261 * T)     # longitude of the Moon's ascending node
    dpsi = (-17.20 * np.sin(Om) - 1.32 * np.sin(2 * L) - 0.23 * np.sin(2 * Lp) + 0.21 * np.sin(2 * Om)) * _ARCSEC
    deps = (9.20 * np.cos(Om) + 0.57 * np.cos(2 * L) + 0.10 * np.cos(2 * Lp) - 0.09 * np.cos(2 * Om)) * _ARCSEC
    # dpsi turns about the ecliptic pole, deps about the equinox (x)
    return np.stack([deps, -dpsi * np.sin(_OBLIQUITY_J2000), dpsi * np.cos(_OBLIQUITY_J2000)], -1)


def remove_nutation(t, v):
    """v (N x 3, true equator of date as SSCWeb delivers it) -> mean equator (J2000)."""
    return v - np.cross(nutation_vector(t), v)


def date(t):
    """TDB seconds -> 'YYYY-MM-DD HH:MM' (TDB)."""
    return (J2000 + dt.timedelta(seconds=float(t))).strftime('%Y-%m-%d %H:%M')


class Track:
    """Cubic Hermite interpolation of Horizons rows (t, x, y, z, vx, vy, vz)."""

    def __init__(self, arr):
        self.t = arr[:, 0]
        self.p = arr[:, 1:4]
        self.v = arr[:, 4:7]

    def _eval(self, t, want_v):
        i = np.clip(np.searchsorted(self.t, t) - 1, 0, len(self.t) - 2)
        t0 = self.t[i]
        h = self.t[i + 1] - t0
        s = (t - t0) / h
        if np.ndim(t):
            h, s = h[:, None], s[:, None]
        p0, p1, v0, v1 = self.p[i], self.p[i + 1], self.v[i] * h, self.v[i + 1] * h
        s2, s3 = s * s, s * s * s
        if not want_v:
            return (2 * s3 - 3 * s2 + 1) * p0 + (s3 - 2 * s2 + s) * v0 + (-2 * s3 + 3 * s2) * p1 + (s3 - s2) * v1
        return ((6 * s2 - 6 * s) * p0 + (3 * s2 - 4 * s + 1) * v0 + (-6 * s2 + 6 * s) * p1 + (3 * s2 - 2 * s) * v1) / h

    def __call__(self, t):
        return self._eval(t, False)

    def velocity(self, t):
        return self._eval(t, True)


class Bodies:
    """Horizons tables of the perturbing bodies (cached .npy files)."""

    def __init__(self, cache_dir):
        def load(name):
            return Track(np.load(os.path.join(cache_dir, f'hzn_{name}.npy')))

        self.moon_geo = load('moon_geo')
        self.sun_geo = load('sun_geo')
        self.helio = {n: load(n + '_hel') for n in
                      ('mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune')}


def _third_body(a, G, r, bodies, want_grad):
    for gm, rk in bodies:
        d = rk - r
        dn = np.linalg.norm(d)
        a = a + gm * (d / dn ** 3 - rk / np.linalg.norm(rk) ** 3)
        if want_grad:
            G = G + gm * (3 * np.outer(d, d) / dn ** 5 - np.eye(3) / dn ** 3)
    return a, G


class GeoModel:
    """Earth-centered: Earth (with J2), Moon, Sun, Venus, Mars, Jupiter; solar radiation pressure."""

    def __init__(self, b):
        self.b = b
        self.center_gm = GM['earth']

    def accel(self, t, r, srp, want_grad):
        rn = np.linalg.norm(r)
        a = -self.center_gm * r / rn ** 3
        G = self.center_gm * (3 * np.outer(r, r) / rn ** 5 - np.eye(3) / rn ** 3) if want_grad else None
        z2 = (r[2] / rn) ** 2      # J2 about the ICRF pole (precession of the equator ignored)
        f = 1.5 * EARTH_J2 * self.center_gm * EARTH_RE ** 2 / rn ** 5
        a = a + f * r * np.array([5 * z2 - 1, 5 * z2 - 1, 5 * z2 - 3])
        sun = self.b.sun_geo(t)
        bodies = [(GM['moon'], self.b.moon_geo(t)), (GM['sun'], sun)]
        bodies += [(GM[n], self.b.helio[n](t) + sun) for n in ('venus', 'mars', 'jupiter')]
        a, G = _third_body(a, G, r, bodies, want_grad)
        u = r - sun
        un = np.linalg.norm(u)
        dsrp = u / un * (AU / un) ** 2
        return a + srp * dsrp, G, dsrp


class HelioModel:
    """Sun-centered: the planets, the Earth and the Moon as third bodies; solar radiation pressure."""

    def __init__(self, b):
        self.b = b
        self.center_gm = GM['sun']

    def accel(self, t, r, srp, want_grad):
        rn = np.linalg.norm(r)
        a = -self.center_gm * r / rn ** 3
        G = self.center_gm * (3 * np.outer(r, r) / rn ** 5 - np.eye(3) / rn ** 3) if want_grad else None
        e = self.b.helio['earth'](t)
        bodies = [(GM['earth'], e), (GM['moon'], e + self.b.moon_geo(t))]
        bodies += [(GM[n], self.b.helio[n](t)) for n in ('mercury', 'venus', 'mars', 'jupiter', 'saturn',
                                                          'uranus', 'neptune')]
        a, G = _third_body(a, G, r, bodies, want_grad)
        dsrp = r / rn * (AU / rn) ** 2
        return a + srp * dsrp, G, dsrp


# --- Integration (Dormand-Prince 5(4), adaptive) ---------------------------------

_C = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1]
_A = [[],
      [1 / 5],
      [3 / 40, 9 / 40],
      [44 / 45, -56 / 15, 32 / 9],
      [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
      [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
      [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]]
_B5 = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0]
_B4 = [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40]


def _rhs(model, t, y, srp, stm):
    a, G, dsrp = model.accel(t, y[0:3], srp, stm)
    out = np.empty_like(y)
    out[0:3] = y[3:6]
    out[3:6] = a
    if stm:
        Phi = y[6:48].reshape(6, 7)       # d(r, v) / d(r0, v0, srp)
        dPhi = np.empty((6, 7))
        dPhi[0:3] = Phi[3:6]
        dPhi[3:6] = G @ Phi[0:3]
        dPhi[3:6, 6] += dsrp
        out[6:48] = dPhi.ravel()
    return out


def integrate(model, t0, y0, t1, srp=0.0, stm=False, rtol=1e-11, h0=60.0, hmax=4 * DAY):
    """Adaptive integration from t0 to t1 (either direction). Returns step times and states
    (6 components, or 48 with the 6x7 state transition matrix appended)."""
    y = np.array(y0, dtype=float)
    if stm and len(y) == 6:
        y = np.concatenate([y, np.hstack([np.eye(6), np.zeros((6, 1))]).ravel()])
    sgn = 1.0 if t1 >= t0 else -1.0
    t = t0
    h = h0 * sgn
    ts, ys = [t], [y.copy()]
    k1 = _rhs(model, t, y, srp, stm)
    while sgn * (t1 - t) > 0:
        if sgn * (t + h - t1) > 0:
            h = t1 - t
        k = [k1]
        for s in range(1, 7):
            k.append(_rhs(model, t + _C[s] * h, y + h * sum(_A[s][j] * k[j] for j in range(s)), srp, stm))
        y5 = y + h * sum(_B5[j] * k[j] for j in range(7))
        y4 = y + h * sum(_B4[j] * k[j] for j in range(7))
        err = np.max(np.abs((y5 - y4)[0:3])) / (rtol * max(np.linalg.norm(y5[0:3]), 1.0))
        if err <= 1.0:
            t, y, k1 = t + h, y5, k[6]
            ts.append(t)
            ys.append(y.copy())
        h *= min(5.0, max(0.2, 0.9 * err ** -0.2 if err > 0 else 5.0))
        if abs(h) > hmax:
            h = hmax * sgn
    return np.array(ts), np.array(ys)


def dense(model, ts, ys, srp, tq):
    """Cubic Hermite through the integrator steps, at the times tq."""
    order = np.argsort(ts)
    ts, ys = ts[order], ys[order]
    stm = ys.shape[1] > 6
    ds = np.array([_rhs(model, t, y, srp, stm) for t, y in zip(ts, ys)])
    i = np.clip(np.searchsorted(ts, tq) - 1, 0, len(ts) - 2)
    h = (ts[i + 1] - ts[i])[:, None]
    s = (tq - ts[i])[:, None] / h
    s2, s3 = s * s, s * s * s
    return ((2 * s3 - 3 * s2 + 1) * ys[i] + (s3 - 2 * s2 + s) * ds[i] * h + (-2 * s3 + 3 * s2) * ys[i + 1]
            + (s3 - s2) * ds[i + 1] * h)


def propagate(model, te, xe, srp, tq, stm=False):
    """States (and STMs) at the times tq of the arc through state xe at te."""
    tq = np.asarray(tq, dtype=float)
    out = np.empty((len(tq), 48 if stm else 6))
    for sel in (tq >= te, tq < te):
        if not sel.any():
            continue
        t_end = tq[sel].max() if tq[sel].max() > te else tq[sel].min()
        ts, ys = integrate(model, te, xe, t_end, srp=srp, stm=stm)
        out[sel] = ys[0] if len(ts) < 2 else dense(model, ts, ys, srp, tq[sel])
    return out


def state_at(model, te, xe, t, srp=0.0):
    return integrate(model, te, xe, t, srp=srp)[1][-1][:6]


# --- Fitting -----------------------------------------------------------------------

def fit_arc(model, tq, obs, te, xe, srp=0.0, iters=6):
    """Least-squares state at te fitting the positions obs (N x 3) at tq, rejecting outliers
    (the source data have isolated spikes of hundreds to thousands of km). Returns the state
    and the residual norms."""
    xe = np.array(xe, float)
    w = np.ones(len(tq), bool)
    for it in range(iters):
        st = propagate(model, te, xe, srp, tq, True)
        res = obs - st[:, 0:3]
        J = st[:, 6:48].reshape(-1, 6, 7)[:, 0:3, 0:6][w].reshape(-1, 6)
        sc = np.linalg.norm(J, axis=0)
        sc[sc == 0] = 1.0
        dx = np.linalg.lstsq(J / sc, res[w].reshape(-1), rcond=None)[0] / sc
        xe = xe + dx
        rn = np.linalg.norm(res, axis=1)
        if it >= 2:
            w = rn < max(5 * np.median(rn), 50.0)
        if np.linalg.norm(dx[:3]) < 1e-2 and it >= 3:
            break
    st = propagate(model, te, xe, srp, tq)
    return xe, np.linalg.norm(obs - st[:, 0:3], axis=1)


def guess(tq, obs, te):
    """Epoch near te and a state from a local quadratic through the samples."""
    i = int(np.argmin(abs(tq - te)))
    j = slice(max(i - 3, 0), i + 4)
    A = np.vstack([np.ones(len(tq[j])), tq[j] - tq[i], (tq[j] - tq[i]) ** 2]).T
    c = np.linalg.lstsq(A, obs[j], rcond=None)[0]
    return tq[i], np.concatenate([c[0], c[1]])
