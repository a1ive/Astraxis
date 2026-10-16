"""Reconstruct Tianwen-1's 2021 path around Mars from the state vectors in its telemetry
(tools/tianwen1/cache/tw1_2021.txt, read by tianwen1_bake.py).

Usage:
    python tianwen1_reconstruct.py <path to the jupyter_notebooks/Tianwen folder>

D. Estevez decoded the orbiter's own position and velocity from its telemetry
(github.com/daniestevez/jupyter_notebooks, Tianwen/orbit/*.script hold the vectors
he used, Tianwen/tm_state_vectors/ some of the raw ones). There are vectors around
each maneuver of 2021 and a few in between, but none between 2021-02-21 and 05-14
and none between 05-14 and 06-03. This tool propagates each vector forward and
backward with its own force model and joins the arcs:

  - Mars orbit insertion (2021-02-10): Estevez's GMAT model of the burn (his
    plane_change.script): from the approach vector of 02-08, a finite burn against the
    velocity at 3 kN (specific impulse 321.6 s, 4,425 kg), 888 s from 11:52:55 UTC. It
    leaves the orbit some 30 km off the capture orbit of the 02-14 vector (its periapsis
    30 km lower), so the two are blended in the hour after the burn.
  - Burns between two vectors (plane change, phasing, deorbit, collision avoidance,
    remote sensing orbit): impulsive, where the forward arc of the vector before and
    the backward arc of the vector after come closest.
  - 2021-02-23, into the 2-sol parking orbit (no vectors from 02-21 to 05-14): Estevez's
    prediction (phasing_orbit.script), -40.9 m/s along the velocity at the periapsis
    after the 02-21 vector; AMSAT-DL's paper: 41 m/s.
  - 2021-05-17, into the 1/3-sol relay orbit (none from 05-14 to 06-03): against the
    velocity at the periapsis of 05-16/17, to the energy of the 06-03 vector's arc.
  - Within an orbit, where the arcs of two vectors meet with small trim maneuvers of
    unknown time in between (at least four in the parking orbit, "some" in the relay
    orbit, AMSAT-DL), the two arcs are blended over one revolution centered on the
    apoapsis where they agree best. The same joins the remote sensing orbit to the
    first MOMAG sample (2021-11-16), where tianwen1_bake.py goes on with MOMAG.

Force model: Mars' gravity field GMM-3 (Genova et al. 2016, PDS
MRO-M-RSS-5-SDP-V1 gmm3_120_sha.tab) to degree and order DEGREE, in the IAU 2000
Mars frame (NAIF pck00008.tpc: the pole GMAT and the telemetry use; GMM-3 is in
Konopliv et al. 2006's frame, which differs by far less than matters here); the Sun,
the Earth-Moon barycenter and Jupiter as third bodies (JPL Horizons, relative to
Mars' center). Estevez's GMAT runs used GMM-2B to degree 50-80 and all planets; the
validation below compares the two.

Axes: the vectors are in Mars-centric body inertial axes (AMSAT-DL: z along Mars'
spin axis at J2000, x along the intersection with the FK5 equator) with the IAU 2000
pole; the output is ICRF, TDB (UTC + 69.184 s), km and km/s, one row a minute.

Downloads (cache/, not in the repository): the gravity field and pck00008 (see
assets/ephem/SOURCES.md), Horizons tables of the third bodies.
"""

import datetime
import hashlib
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'tools', 'bake'))
import horizons_bake as hb  # noqa: E402

CACHE = os.path.join(HERE, 'cache')
GRAVITY_FILE = os.path.join(CACHE, 'gravity', 'gmm3_120_sha.tab')
THIRD_BODY_FILE = os.path.join(CACHE, 'thirdbody_2021.json')
OUT = os.path.join(CACHE, 'tw1_2021.txt')

DEGREE = 10
TDB_MINUS_UTC = 69.184             # 37 leap seconds (since 2017) + 32.184 s
POLE_RA_DEC_IAU2000 = (317.68143, 52.88650)
# IAU 2000 Mars rotation, NAIF pck00008.tpc: RA, Dec (deg, per Julian century), W (deg, per day)
MARS_POLE = ((317.68143, -0.1061), (52.88650, -0.0609))
MARS_PM = (176.630, 350.89198226)
# GM (km^3/s^2): JPL SSD astrodynamic parameters (DE440), https://ssd.jpl.nasa.gov/astro_par.html
GM_SUN = 1.32712440041279419e11
GM_EMB = 398600.435507 + 4902.800118
GM_JUPITER = 126712764.1
GM_MARS_SYSTEM = 42828.375816
STEP = 60.0
G0 = 9.80665

# Mars orbit insertion as Estevez modeled it in GMAT (Tianwen/orbit/plane_change.script):
# start at UTCModJulian 29255.9950810185 (JD - 2430000), 888 s, retrograde along the
# velocity, C1 = 3000 N, K1 = 321.6 s, DryMass 3390 + FuelMass 1035 kg.
MOI_START_GMAT_MJD = 29255.9950810185185185
MOI_DURATION = 888.0
MOI_THRUST_N = 3000.0
MOI_ISP = 321.6
MOI_MASS = 3390.0 + 1035.0
# 2021-02-23 apoapsis lowering as Estevez predicted it (Tianwen/orbit/phasing_orbit.script):
# -0.0409 km/s along the velocity at the periapsis after the 02-21 22:32:58 vector.
PARKING_DV = -0.0409

# Telemetry state vectors: (name, script, UTC epoch as written there).
VECTORS = [
    ('approach', 'plane_change.script', '08 Feb 2021 22:14:25.049'),
    ('capture', 'plane_change_vectors.script', '14 Feb 2021 23:45:31.728'),
    ('polar', 'plane_change_vectors.script', '15 Feb 2021 10:24:29.218'),
    ('pre_phasing', 'phasing_burn_vectors.script', '20 Feb 2021 10:42:16.567'),
    ('phasing', 'phasing_burn_vectors.script', '21 Feb 2021 22:32:58.245'),
    ('pre_deorbit', 'deorbit.script', '14 May 2021 15:52:57.847'),
    ('deorbit', 'deorbit.script', '14 May 2021 18:05:05.349'),
    ('pre_avoidance', 'avoidance.script', '14 May 2021 20:03:29.441'),
    ('avoidance', 'avoidance.script', '14 May 2021 21:36:34.715'),
    ('relay_jun', 'comms_orbit.script', '03 Jun 2021 09:41:25.417'),
    ('relay_sep', 'comms_orbit.script', '07 Sep 2021 10:07:35.999'),   # commented out there
    ('pre_remote', 'remote_sensing_orbit_change.script', '08 Nov 2021 10:21:08.225'),
    ('remote', 'remote_sensing_orbit_change.script', '09 Nov 2021 08:10:41.020'),
]

# Estevez's GMAT reports (propagations of the same vectors), for validation.
REPORTS = [
    ('plane_change/Tianwen1_pre_planechange.txt', 'capture'),
    ('plane_change/Tianwen1_post_planechange.txt', 'polar'),
    ('phasing/Tianwen1_pre_phasing.txt', 'pre_phasing'),
    ('phasing/Tianwen1_post_phasing.txt', 'phasing'),
    ('deorbit/Tianwen1_pre_deorbit.txt', 'pre_deorbit'),
    ('deorbit/Tianwen1_post_deorbit.txt', 'deorbit'),
    ('avoidance/Tianwen1_pre_avoidance.txt', 'pre_avoidance'),
    ('avoidance/Tianwen1_post_avoidance.txt', 'avoidance'),
    ('remote_sensing/Tianwen1_pre_remote_sensing.txt', 'pre_remote'),
    ('remote_sensing/Tianwen1_post_remote_sensing.txt', 'remote'),
]

J2000 = datetime.datetime(2000, 1, 1, 12)


def tdb_from_utc_text(s):
    d = datetime.datetime.strptime(s, '%d %b %Y %H:%M:%S.%f')
    return (d - J2000).total_seconds() + TDB_MINUS_UTC


def tdb_from_gmat_mjd(mjd):
    return (mjd + 2430000.0 - hb.J2000_JD) * 86400.0 + TDB_MINUS_UTC


def date(t):
    return (J2000 + datetime.timedelta(seconds=float(t))).strftime('%Y-%m-%d %H:%M:%S')


# --- Frames ----------------------------------------------------------------------

def body_inertial_matrix():
    """Columns: the body inertial axes in ICRF (IAU 2000 pole at J2000)."""
    ra, dec = (math.radians(a) for a in POLE_RA_DEC_IAU2000)
    z = np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)])
    x = np.array([-math.sin(ra), math.cos(ra), 0.0])
    return np.column_stack([x, np.cross(z, x), z])


BI = body_inertial_matrix()


def mars_fixed_matrix(t):
    """ICRF -> IAU 2000 Mars body-fixed (rows: the body axes in ICRF)."""
    T = t / (36525.0 * 86400.0)
    ra = math.radians(MARS_POLE[0][0] + MARS_POLE[0][1] * T)
    dec = math.radians(MARS_POLE[1][0] + MARS_POLE[1][1] * T)
    w = math.radians((MARS_PM[0] + MARS_PM[1] * t / 86400.0) % 360.0)
    z = np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)])
    node = np.array([-math.sin(ra), math.cos(ra), 0.0])   # RA + 90 deg
    q = np.cross(z, node)
    x = math.cos(w) * node + math.sin(w) * q
    y = np.cross(z, x)
    return np.vstack([x, y, z])


def elements(r, v, gm=GM_MARS_SYSTEM):
    """Periapsis and apoapsis radii (km), period (h), inclination (deg) to the IAU 2000
    equator at J2000."""
    rn, v2 = np.linalg.norm(r), v @ v
    a = 1.0 / (2.0 / rn - v2 / gm)
    h = np.cross(r, v)
    e = np.linalg.norm(np.cross(v, h) / gm - r / rn)
    inc = math.degrees(math.acos((BI.T @ h)[2] / np.linalg.norm(h)))
    period = 2 * math.pi * math.sqrt(a ** 3 / gm) / 3600.0 if a > 0 else 0.0
    return a * (1 - e), a * (1 + e) if a > 0 else 0.0, period, inc


# --- Force model -------------------------------------------------------------------

class Gravity:
    """Spherical harmonics (Montenbruck & Gill, Satellite Orbits, 3.2.4-5: the V/W
    recursions), unnormalized coefficients to degree n_max."""

    def __init__(self, path, n_max):
        with open(path) as f:
            head = [float(c) for c in f.readline().split(',')[:3]]
            self.R, self.gm = head[0], head[1]
            self.n = n_max
            self.C = [[0.0] * (n_max + 1) for _ in range(n_max + 1)]
            self.S = [[0.0] * (n_max + 1) for _ in range(n_max + 1)]
            self.C[0][0] = 1.0
            for line in f:
                c = line.split(',')
                n, m = int(c[0]), int(c[1])
                if n > n_max:
                    break
                k = math.sqrt((2 - (m == 0)) * (2 * n + 1) * math.factorial(n - m) / math.factorial(n + m))
                self.C[n][m] = float(c[2]) * k
                self.S[n][m] = float(c[3]) * k

    def _vw(self, x, y, z):
        N = self.n + 1
        R = self.R
        r2 = x * x + y * y + z * z
        rho = R * R / r2
        x0, y0, z0 = R * x / r2, R * y / r2, R * z / r2
        V = [[0.0] * (N + 1) for _ in range(N + 1)]
        W = [[0.0] * (N + 1) for _ in range(N + 1)]
        V[0][0] = R / math.sqrt(r2)
        V[1][0] = z0 * V[0][0]
        for n in range(2, N + 1):
            V[n][0] = ((2 * n - 1) * z0 * V[n - 1][0] - (n - 1) * rho * V[n - 2][0]) / n
        for m in range(1, N + 1):
            V[m][m] = (2 * m - 1) * (x0 * V[m - 1][m - 1] - y0 * W[m - 1][m - 1])
            W[m][m] = (2 * m - 1) * (x0 * W[m - 1][m - 1] + y0 * V[m - 1][m - 1])
            if m < N:
                V[m + 1][m] = (2 * m + 1) * z0 * V[m][m]
                W[m + 1][m] = (2 * m + 1) * z0 * W[m][m]
            for n in range(m + 2, N + 1):
                V[n][m] = ((2 * n - 1) * z0 * V[n - 1][m] - (n + m - 1) * rho * V[n - 2][m]) / (n - m)
                W[n][m] = ((2 * n - 1) * z0 * W[n - 1][m] - (n + m - 1) * rho * W[n - 2][m]) / (n - m)
        return V, W

    def potential(self, p):
        V, W = self._vw(*p)
        return self.gm / self.R * sum(self.C[n][m] * V[n][m] + self.S[n][m] * W[n][m]
                                      for n in range(self.n + 1) for m in range(n + 1))

    def accel(self, p):
        V, W = self._vw(*p)
        C, S = self.C, self.S
        ax = ay = az = 0.0
        for n in range(self.n + 1):
            c = C[n][0]
            ax -= c * V[n + 1][1]
            ay -= c * W[n + 1][1]
            az -= (n + 1) * c * V[n + 1][0]
            for m in range(1, n + 1):
                c, s = C[n][m], S[n][m]
                if c == 0.0 and s == 0.0:
                    continue
                f = (n - m + 1) * (n - m + 2)
                ax += 0.5 * (-c * V[n + 1][m + 1] - s * W[n + 1][m + 1]
                             + f * (c * V[n + 1][m - 1] + s * W[n + 1][m - 1]))
                ay += 0.5 * (-c * W[n + 1][m + 1] + s * V[n + 1][m + 1]
                             + f * (-c * W[n + 1][m - 1] + s * V[n + 1][m - 1]))
                az += (n - m + 1) * (-c * V[n + 1][m] - s * W[n + 1][m])
        k = self.gm / (self.R * self.R)
        return np.array([ax * k, ay * k, az * k])


class ThirdBodies:
    """The Sun, the Earth-Moon barycenter and Jupiter relative to Mars' center (Horizons,
    6-hourly, cubic Hermite)."""

    BODIES = (('10', GM_SUN), ('3', GM_EMB), ('5', GM_JUPITER))

    def __init__(self, t0, t1):
        tables = {}
        if os.path.exists(THIRD_BODY_FILE):
            tables = json.load(open(THIRD_BODY_FILE))
        jd0, jd1 = t0 / 86400.0 + hb.J2000_JD - 1.0, t1 / 86400.0 + hb.J2000_JD + 1.0
        key = f'{jd0:.1f}-{jd1:.1f}'
        if tables.get('span') != key:
            tables = {'span': key}
            for target, _ in self.BODIES:
                tables[target] = [[t, *p, *v] for t, p, v in hb.fetch(target, '500@499', jd0, jd1, 360)]
            json.dump(tables, open(THIRD_BODY_FILE, 'w'))
        self.tracks = [(np.array(tables[target]), gm) for target, gm in self.BODIES]

    def accel(self, t, r):
        a = np.zeros(3)
        for rows, gm in self.tracks:
            i = min(max(int((t - rows[0, 0]) // (rows[1, 0] - rows[0, 0])), 0), len(rows) - 2)
            t0, h = rows[i, 0], rows[i + 1, 0] - rows[i, 0]
            s = (t - t0) / h
            s2, s3 = s * s, s * s * s
            p = ((2 * s3 - 3 * s2 + 1) * rows[i, 1:4] + (s3 - 2 * s2 + s) * h * rows[i, 4:7]
                 + (-2 * s3 + 3 * s2) * rows[i + 1, 1:4] + (s3 - s2) * h * rows[i + 1, 4:7])
            d = p - r
            a += gm * (d / np.linalg.norm(d) ** 3 - p / np.linalg.norm(p) ** 3)
        return a


class Model:
    def __init__(self, gravity, third):
        self.g = gravity
        self.third = third
        self.burn = None   # (t0, t1, thrust N, mdot kg/s, m0 kg): retrograde finite burn

    def rhs(self, t, y):
        r, v = y[0:3], y[3:6]
        M = mars_fixed_matrix(t)
        a = M.T @ self.g.accel(M @ r) + self.third.accel(t, r)
        if self.burn and self.burn[0] <= t <= self.burn[1]:
            t0, _, thrust, mdot, m0 = self.burn
            a = a - thrust / (m0 - mdot * (t - t0)) / 1000.0 * v / np.linalg.norm(v)
        return np.concatenate([v, a])


# --- Integration (Dormand-Prince 5(4), adaptive; as tools/isee3/isee3_dynamics.py) --------

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
RTOL = 1e-11


class Arc:
    """Integrator steps of one trajectory, cubic Hermite in between."""

    def __init__(self, ts, ys, fs):
        order = np.argsort(ts)
        self.t, self.y, self.f = np.array(ts)[order], np.array(ys)[order], np.array(fs)[order]

    def __call__(self, tq):
        tq = np.atleast_1d(np.asarray(tq, float))
        if tq.min() < self.t[0] - 1e-6 or tq.max() > self.t[-1] + 1e-6:
            raise ValueError(f'outside the arc {date(self.t[0])} .. {date(self.t[-1])}')
        i = np.clip(np.searchsorted(self.t, tq) - 1, 0, len(self.t) - 2)
        h = (self.t[i + 1] - self.t[i])[:, None]
        s = (tq - self.t[i])[:, None] / h
        s2, s3 = s * s, s * s * s
        y0, y1, f0, f1 = self.y[i], self.y[i + 1], self.f[i] * h, self.f[i + 1] * h
        p = (2 * s3 - 3 * s2 + 1) * y0 + (s3 - 2 * s2 + s) * f0 + (-2 * s3 + 3 * s2) * y1 + (s3 - s2) * f1
        dp = ((6 * s2 - 6 * s) * y0 + (3 * s2 - 4 * s + 1) * f0 + (-6 * s2 + 6 * s) * y1
              + (3 * s2 - 2 * s) * f1) / h
        return np.hstack([p[:, 0:3], dp[:, 0:3]])

    def state(self, t):
        return self(t)[0]


def integrate(model, t0, y0, t1, hmax=600.0):
    y = np.array(y0, float)
    sgn = 1.0 if t1 >= t0 else -1.0
    t, h = t0, 10.0 * sgn
    k1 = model.rhs(t, y)
    ts, ys, fs = [t], [y.copy()], [k1]
    while sgn * (t1 - t) > 1e-9:
        if sgn * (t + h - t1) > 0:
            h = t1 - t
        k = [k1]
        for s in range(1, 7):
            k.append(model.rhs(t + _C[s] * h, y + h * sum(_A[s][j] * k[j] for j in range(s))))
        y5 = y + h * sum(_B5[j] * k[j] for j in range(7))
        y4 = y + h * sum(_B4[j] * k[j] for j in range(7))
        err = np.max(np.abs((y5 - y4)[0:3])) / (RTOL * max(np.linalg.norm(y5[0:3]), 1.0))
        if err <= 1.0:
            t, y, k1 = t + h, y5, k[6]
            ts.append(t)
            ys.append(y.copy())
            fs.append(k1)
        h *= min(5.0, max(0.2, 0.9 * err ** -0.2 if err > 0 else 5.0))
        if abs(h) > hmax:
            h = hmax * sgn
    return ts, ys, fs


def propagate(model, te, ye, t_lo, t_hi, stops=()):
    """The arc through state ye at te over [t_lo, t_hi], restarting at the times in
    stops (burn edges)."""
    ts, ys, fs = [], [], []
    for a, b in ((te, t_hi), (te, t_lo)):
        cuts = sorted((s for s in stops if min(a, b) < s < max(a, b)), reverse=b < a)
        t, y = a, ye
        for c in [*cuts, b]:
            if c == t:
                continue
            pts, pys, pfs = integrate(model, t, y, c)
            ts += pts
            ys += pys
            fs += pfs
            t, y = c, pys[-1]
    # Drop duplicate times (restarts); keep the later one in each direction's order.
    order = np.argsort(ts, kind='stable')
    keep = [i for j, i in enumerate(order) if j == len(order) - 1 or ts[order[j + 1]] - ts[i] > 1e-6]
    return Arc([ts[i] for i in keep], [ys[i] for i in keep], [fs[i] for i in keep])


# --- Sources ------------------------------------------------------------------------

def read_vector(orbit_dir, script, epoch):
    """The state vector following the line with this epoch in an Estevez GMAT script
    (body inertial axes) -> (TDB, ICRF state)."""
    lines = open(os.path.join(orbit_dir, script)).read().splitlines()
    start = next(i for i, l in enumerate(lines) if epoch in l and 'Epoch' in l)
    vals = {}
    for line in lines[start + 1:start + 12]:
        m = re.search(r'\.(X|Y|Z|VX|VY|VZ)\s*=\s*([-+0-9.eE]+)', line)
        if m and m.group(1) not in vals:
            vals[m.group(1)] = float(m.group(2))
    p = BI @ np.array([vals['X'], vals['Y'], vals['Z']])
    v = BI @ np.array([vals['VX'], vals['VY'], vals['VZ']])
    return tdb_from_utc_text(epoch), np.concatenate([p, v])


def read_report(path):
    rows = np.loadtxt(path)
    rows = rows[np.argsort(rows[:, 0])]
    t = tdb_from_gmat_mjd(rows[:, 0])
    return t, rows[:, 1:4] @ BI.T


def read_tm(path):
    """Raw telemetry vectors (Tianwen/tm_state_vectors): [timestamp] date time x y z vx vy vz."""
    out = []
    for line in open(path):
        f = line.split()
        d = datetime.datetime.fromisoformat(f[1] + 'T' + f[2])
        out.append(((d - J2000).total_seconds() + TDB_MINUS_UTC,
                    np.concatenate([BI @ np.array([float(c) for c in f[3:6]]),
                                    BI @ np.array([float(c) for c in f[6:9]])])))
    return out


# --- Joins --------------------------------------------------------------------------

def closest(a, b, t_lo, t_hi):
    """Time in [t_lo, t_hi] where arcs a and b come closest (1 s), and the distance."""
    tq = np.arange(t_lo, t_hi, STEP)
    d = np.linalg.norm(a(tq)[:, 0:3] - b(tq)[:, 0:3], axis=1)
    k = int(np.argmin(d))
    tq = np.arange(tq[max(k - 1, 0)], tq[min(k + 1, len(tq) - 1)] + 1.0, 1.0)
    d = np.linalg.norm(a(tq)[:, 0:3] - b(tq)[:, 0:3], axis=1)
    k = int(np.argmin(d))
    return tq[k], d[k]


def apoapses(arc, t_lo, t_hi):
    tq = np.arange(t_lo, t_hi, STEP)
    r = np.linalg.norm(arc(tq)[:, 0:3], axis=1)
    return [tq[i] for i in range(1, len(tq) - 1) if r[i - 1] < r[i] >= r[i + 1]]


def periapsis_after(arc, t, span):
    tq = np.arange(t, t + span, 1.0)
    return tq[int(np.argmin(np.linalg.norm(arc(tq)[:, 0:3], axis=1)))]


def best_apoapsis(a, b, t_lo, t_hi):
    """Among the apoapses in [t_lo, t_hi], the one where arcs a and b agree best."""
    best = None
    for t in apoapses(a, t_lo, t_hi):
        d = np.linalg.norm(a.state(t)[0:3] - b.state(t)[0:3])
        if best is None or d < best[1]:
            best = (t, d)
    return best


class Path:
    """Arcs in time order with the joins between them."""

    def __init__(self):
        self.parts = []   # (arc, blend start, blend end): this arc from blend end on

    def add(self, arc, t_from, t_to=None):
        self.parts.append((arc, t_from, t_from if t_to is None else t_to))

    def states(self, tq):
        out = np.empty((len(tq), 6))
        for i, (arc, b0, b1) in enumerate(self.parts):
            end = self.parts[i + 1][1] if i + 1 < len(self.parts) else math.inf
            sel = (tq >= b0) & (tq < end)
            if not sel.any():
                continue
            out[sel] = arc(tq[sel])
            blend = sel & (tq < b1)
            if b1 > b0 and i > 0 and blend.any():
                mine = arc(tq[blend])
                prev = self.parts[i - 1][0](tq[blend])
                u = ((tq[blend] - b0) / (b1 - b0))[:, None]
                w = u * u * (3 - 2 * u)
                dw = 6 * u * (1 - u) / (b1 - b0)
                mix = prev + w * (mine - prev)
                mix[:, 3:6] += dw * (mine[:, 0:3] - prev[:, 0:3])
                out[blend] = mix
        return out


# --- Main --------------------------------------------------------------------------------

def check_gravity(g):
    p = np.array([2100.0, -2900.0, 1700.0])
    a = g.accel(p)
    num = np.array([(g.potential(p + e) - g.potential(p - e)) / 2e-3 for e in np.eye(3) * 1e-3])
    err = np.linalg.norm(a - num) / np.linalg.norm(a)
    print(f'gravity: degree {g.n}, R {g.R} km, GM {g.gm}; gradient check {err:.1e}')
    assert err < 1e-7


def describe(name, y):
    rp, ra, period, inc = elements(y[0:3], y[3:6])
    if ra:
        return f'{name}: {rp - 3389.5:.0f} x {ra - 3389.5:.0f} km (mean radius), {period:.2f} h, i {inc:.1f} deg'
    return f'{name}: hyperbolic, periapsis {rp - 3389.5:.0f} km, i {inc:.1f} deg'


def energy(y):
    return 0.5 * (y[3:6] @ y[3:6]) - GM_MARS_SYSTEM / np.linalg.norm(y[0:3])


def mso_to_icrf(third, t, p, v):
    """MOMAG's MSO (x to the Sun, z along Mars' orbital angular momentum) -> ICRF, with the
    frame's rotation (as tianwen1_bake.py, from the Sun relative to Mars)."""
    rows = third.tracks[0][0]
    i = min(max(int((t - rows[0, 0]) // (rows[1, 0] - rows[0, 0])), 0), len(rows) - 2)
    u = (t - rows[i, 0]) / (rows[i + 1, 0] - rows[i, 0])
    s = rows[i, 1:4] + u * (rows[i + 1, 1:4] - rows[i, 1:4])
    vs = rows[i, 4:7] + u * (rows[i + 1, 4:7] - rows[i, 4:7])
    x = s / np.linalg.norm(s)
    h = np.cross(s, vs)
    z = h / np.linalg.norm(h)
    m = np.column_stack([x, np.cross(z, x), z])
    r = m @ np.asarray(p)
    return np.concatenate([r, m @ np.asarray(v) + np.cross(h / (s @ s), r)])


def first_momag(third):
    """The first cached MOMAG sample (momag_fetch.py: UTC, MSO position and velocity)."""
    files = sorted(f for f in os.listdir(CACHE) if re.fullmatch(r'\d{8}\.txt', f))
    for line in open(os.path.join(CACHE, files[0]), encoding='ascii'):
        if not line.startswith('#'):
            f = line.split()
            t = (datetime.datetime.fromisoformat(f[0]) - J2000).total_seconds() + TDB_MINUS_UTC
            return t, mso_to_icrf(third, t, [float(c) for c in f[1:4]], [float(c) for c in f[4:7]])


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    tw_dir = sys.argv[1]
    orbit_dir = os.path.join(tw_dir, 'orbit')
    vec = {name: read_vector(orbit_dir, script, epoch) for name, script, epoch in VECTORS}
    t_start = vec['approach'][0]
    t_momag = (datetime.datetime(2021, 11, 16, 12) - J2000).total_seconds()   # within the first MOMAG day

    g = Gravity(GRAVITY_FILE, DEGREE)
    check_gravity(g)
    third = ThirdBodies(t_start, t_momag + 86400.0)
    model = Model(g, third)
    vec['momag'] = first_momag(third)
    t_end = vec['momag'][0]
    day = 86400.0

    def arc(name, lo, hi, stops=(), state=None):
        te, ye = state if state else vec[name]
        # Arcs are cached (cache/arcs/), keyed by everything that goes into them.
        key = repr((te, list(ye), lo, hi, list(stops), model.burn, DEGREE, RTOL)).encode()
        path = os.path.join(CACHE, 'arcs', f'{name}_{hashlib.sha1(key).hexdigest()[:12]}.npz')
        if os.path.exists(path):
            z = np.load(path)
            return Arc(z['t'], z['y'], z['f'])
        print(f'  propagating {name} ({date(te)}) over {date(lo)} .. {date(hi)}', flush=True)
        a = propagate(model, te, ye, lo, hi, stops)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        np.savez(path, t=a.t, y=a.y, f=a.f)
        return a

    print('Arcs:')
    capture = arc('capture', tdb_from_gmat_mjd(MOI_START_GMAT_MJD) - 3600.0, vec['capture'][0] + 12 * 3600.0)
    polar = arc('polar', vec['capture'][0], vec['pre_phasing'][0])
    pre_phasing = arc('pre_phasing', vec['polar'][0] + 3600.0, vec['pre_phasing'][0] + 4 * 3600.0)
    # 280 x 84,600 km: 3.45 days around; the 02-21 vector is near the apoapsis.
    phasing = arc('phasing', vec['pre_phasing'][0], vec['phasing'][0] + 3.5 * day)
    pre_deorbit = arc('pre_deorbit', vec['phasing'][0] + day, vec['pre_deorbit'][0] + 3 * 3600.0)
    deorbit = arc('deorbit', vec['pre_deorbit'][0], vec['pre_avoidance'][0])
    pre_avoidance = arc('pre_avoidance', vec['deorbit'][0], vec['avoidance'][0])
    avoidance = arc('avoidance', vec['pre_avoidance'][0], vec['avoidance'][0] + 3.5 * day)
    relay_jun = arc('relay_jun', vec['avoidance'][0] + 1.5 * day, vec['relay_sep'][0])
    relay_sep = arc('relay_sep', vec['relay_jun'][0], vec['pre_remote'][0])
    pre_remote = arc('pre_remote', vec['relay_sep'][0], vec['pre_remote'][0] + 4 * 3600.0)
    remote = arc('remote', vec['pre_remote'][0], t_end)
    momag = arc('momag', vec['remote'][0], t_end + 3600.0)

    # Mars orbit insertion: Estevez's burn. It leaves the orbit some 30 km off the capture
    # orbit of the 02-14 vector (periapsis 30 km lower); the two are blended in the hour after.
    moi0 = tdb_from_gmat_mjd(MOI_START_GMAT_MJD)
    moi1 = moi0 + MOI_DURATION
    mdot = MOI_THRUST_N / (MOI_ISP * G0)
    model.burn = (moi0, moi1, MOI_THRUST_N, mdot, MOI_MASS)
    approach = arc('approach', t_start, moi1 + 3 * 3600.0, stops=(moi0, moi1))
    model.burn = None
    print('Mars orbit insertion:')
    t, d = closest(approach, capture, moi1, moi1 + 3 * 3600.0)
    dv = MOI_ISP * G0 * math.log(MOI_MASS / (MOI_MASS - mdot * MOI_DURATION))
    print(f'  Estevez\'s burn, {date(moi0)} .. {date(moi1)} ({dv:.0f} m/s): the capture arc within {d:.1f} km '
          f'at {date(t)}, {np.linalg.norm(approach.state(moi1 + 3600.0)[0:3] - capture.state(moi1 + 3600.0)[0:3]):.1f} km '
          f'an hour after the burn')
    print('   ', describe('after the burn', approach.state(moi1 + 600.0)))
    print('   ', describe('the capture orbit (02-14 vector)', capture.state(moi1 + 600.0)))

    print('Validation against Estevez\'s GMAT reports (same vectors, GMM-2B to degree 50-80):')
    for rel, name in REPORTS:
        t, p = read_report(os.path.join(orbit_dir, rel))
        a = {'capture': capture, 'polar': polar, 'pre_phasing': pre_phasing, 'phasing': phasing,
             'pre_deorbit': pre_deorbit, 'deorbit': deorbit, 'pre_avoidance': pre_avoidance,
             'avoidance': avoidance, 'pre_remote': pre_remote, 'remote': remote}[name]
        sel = (t >= a.t[0]) & (t <= a.t[-1])
        d = np.linalg.norm(a(t[sel])[:, 0:3] - p[sel], axis=1)
        print(f'  {rel}: {sel.sum()} rows over {(t[sel][-1] - t[sel][0]) / 3600:.1f} h, max {d.max():.3f} km')
    tm = read_tm(os.path.join(tw_dir, 'tm_state_vectors', '20210208'))
    d = max(np.linalg.norm(approach.state(t)[0:3] - y[0:3]) for t, y in tm)
    print(f'  raw telemetry vectors of 2021-02-08 .. 09 against the approach arc: max {d:.3f} km')

    print('Joins:')
    path = Path()
    path.add(approach, t_start)
    path.add(capture, moi1, moi1 + 3600.0)
    burns = {'moi': [moi0, moi1]}

    def cross(a, b, lo, hi, label):
        t, d = closest(a, b, lo, hi)
        dv = np.linalg.norm(b.state(t)[3:6] - a.state(t)[3:6]) * 1000.0
        print(f'  {label}: {date(t)}, arcs {d:.2f} km apart, {dv:.1f} m/s')
        print('   ', describe('before', a.state(t)))
        print('   ', describe('after', b.state(t)))
        path.add(b, t)
        burns[label] = t
        return t

    def blend(a, b, lo, hi, label):
        """Over the revolution around the apoapsis in [lo, hi] where a and b agree best."""
        for k in range(0, 200, 10):
            t = lo + k * day
            if t < hi:
                d = np.linalg.norm(a.state(t)[0:3] - b.state(t)[0:3])
                print(f'    {date(t)[:10]}: {d:.0f} km apart')
        t_meet, d = best_apoapsis(a, b, lo, hi)
        half = 0.5 * elements(*np.split(a.state(t_meet), 2))[2] * 3600.0
        print(f'  {label}: blended over the revolution around the apoapsis of {date(t_meet)} ({d:.0f} km apart there)')
        path.add(b, t_meet - half, t_meet + half)
        return t_meet

    cross(capture, polar, vec['capture'][0], vec['polar'][0], 'plane change')
    t, d = closest(polar, pre_phasing, vec['polar'][0] + 3600.0, vec['pre_phasing'][0])
    print(f'  02-15 .. 02-20 (no maneuver): the two vectors\' arcs come within {d:.2f} km at {date(t)}')
    path.add(pre_phasing, t)
    cross(pre_phasing, phasing, vec['pre_phasing'][0], vec['pre_phasing'][0] + 4 * 3600.0, 'phasing')

    # Into the parking orbit, at the periapsis after the 02-21 vector: Estevez's prediction.
    # This 2-sol orbit is in resonance with Mars' rotation (its ground track repeats), and the
    # tesseral harmonics make its phase over 80 days depend erratically on the burn: no single
    # burn reaches the 05-14 vector, the spacecraft needed trims. Blended in between instead.
    t_apo = apoapses(phasing, vec['phasing'][0] - 6 * 3600.0, vec['phasing'][0] + 2 * day)[0]
    t_park = periapsis_after(phasing, t_apo, 2.5 * day)
    y = phasing.state(t_park)
    y[3:6] += PARKING_DV * y[3:6] / np.linalg.norm(y[3:6])
    parking = arc('parking', t_park, vec['pre_deorbit'][0], state=(t_park, y))
    print(f'  parking orbit: {date(t_park)}, {abs(PARKING_DV) * 1000:.1f} m/s (Estevez\'s prediction)')
    print('   ', describe('after', parking.state(t_park + 60.0)))
    path.add(parking, t_park)
    burns['parking orbit'] = t_park
    blend(parking, pre_deorbit, t_park + 2 * day, vec['pre_deorbit'][0] - 2 * day, 'parking orbit trims')
    cross(pre_deorbit, deorbit, vec['pre_deorbit'][0], vec['deorbit'][0], 'deorbit')
    t, d = closest(deorbit, pre_avoidance, vec['deorbit'][0], vec['pre_avoidance'][0])
    print(f'  18:05 .. 20:03 (no maneuver): the two vectors\' arcs come within {d:.2f} km at {date(t)}')
    path.add(pre_avoidance, t)
    cross(pre_avoidance, avoidance, vec['pre_avoidance'][0], vec['avoidance'][0], 'collision avoidance')

    # Into the relay orbit, at the periapsis of 05-16/17 (no vectors from 05-14 to 06-03):
    # against the velocity, to the energy of the 06-03 arc there; the two blended in between
    # (a resonant orbit again, with trims).
    t_relay = periapsis_after(avoidance, vec['avoidance'][0] + 1.5 * day, 1.5 * day)
    y = avoidance.state(t_relay)
    v_new = math.sqrt(2.0 * (energy(relay_jun.state(t_relay)) + GM_MARS_SYSTEM / np.linalg.norm(y[0:3])))
    dv = v_new - np.linalg.norm(y[3:6])
    y[3:6] *= v_new / np.linalg.norm(y[3:6])
    relay_may = arc('relay_may', t_relay, vec['relay_jun'][0], state=(t_relay, y))
    print(f'  relay orbit: {date(t_relay)}, {abs(dv) * 1000:.1f} m/s')
    print('   ', describe('before', avoidance.state(t_relay)))
    print('   ', describe('after', relay_may.state(t_relay + 60.0)))
    path.add(relay_may, t_relay)
    burns['relay orbit'] = t_relay
    blend(relay_may, relay_jun, t_relay + day, vec['relay_jun'][0] - day, 'relay orbit trims (05-17 / 06-03)')
    blend(relay_jun, relay_sep, relay_sep.t[0] + day, relay_jun.t[-1] - day, 'relay orbit trims (06-03 / 09-07)')
    blend(relay_sep, pre_remote, pre_remote.t[0] + day, vec['pre_remote'][0] - day, 'relay orbit trims (09-07 / 11-08)')
    cross(pre_remote, remote, vec['pre_remote'][0], vec['pre_remote'][0] + 4 * 3600.0, 'remote sensing orbit')
    blend(remote, momag, vec['remote'][0] + day, t_end - 0.5 * day, 'remote sensing orbit / MOMAG')

    tq = np.arange(t_start, t_end, STEP)
    rows = path.states(tq)
    with open(OUT, 'w') as f:
        f.write('# Tianwen-1 2021, tianwen1_reconstruct.py: TDB s since J2000, Mars-centered ICRF km, km/s\n')
        for t, y in zip(tq, rows):
            f.write(f'{t:.3f} ' + ' '.join(f'{c:.6f}' for c in y[0:3]) + ' '
                    + ' '.join(f'{c:.9f}' for c in y[3:6]) + '\n')
    print(f'wrote {OUT}: {len(tq)} rows, {date(tq[0])} .. {date(tq[-1])}')
    for label, t in burns.items():
        print(f'  burn {label}: ' + (' .. '.join(date(x) for x in t) if isinstance(t, list) else date(t)))


if __name__ == '__main__':
    main()
