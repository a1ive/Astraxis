"""Fit the rotation of comet 67P to ESA's comet attitude kernels.

Usage:
    python c67p_spin.py <directory with the CATT kernels>

writes assets/ephem/c67p_spin.toml, the `spin_table` of 67P in rosetta.toml.

The kernels come from the ESA SPICE Service
(https://spiftp.esac.esa.int/data/SPICE/ROSETTA/kernels/ck/ and former_versions/):
the three that the operational meta-kernel ROS_OPS.TM loads, in its order (a
later kernel wins where they overlap). Each is a navigation solution of the
comet's attitude in the CK-based frame 67P/C-G_CK, relative to J2000 (taken as
ICRF), a quaternion every 10 minutes, propagated to 2017.

The table holds knots of the pole (RA, Dec) and the prime meridian angle W
(from the node of the comet's equator on the ICRF equator, as in the IAU
model) with its rate. Between knots the pole is interpolated linearly and W
by cubic Hermite; before the first knot and after the last, the pole stays
and W turns at the end rate. The knots are fitted by least squares, and more
are added until the whole rotation (the angle between the fitted and the
kernel axes) is within TOL_DEG. The kernels themselves jump by up to half a
degree where one navigation solution hands over to the next, and in 2016 the
pole nods by about 0.1 deg: a smooth fit passes through these.

The module also evaluates the table as the scene does (Table), for
philae_bake.py.

Only the Python standard library is used.
"""

import bisect
import math
import os
import sys
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import ck  # noqa: E402

KERNELS = [  # ROS_OPS.TM (ROS_OPS_V350_20220906_001), lowest priority first
    'CATT_DV_145_02_______00216.BC',  # 2014-01-01 .. (wins until 2015-02-28)
    'CATT_DV_223_02_______00302.BC',  # 2015-02-28 .. (wins until 2015-12-30)
    'CATT_DV_257_03_______00344.BC',  # 2015-12-30 .. 2017-01-01
]
COMET_FRAME_ID = -1000012000  # 67P/C-G_CK
J2000_FRAME_ID = 1
# SCLK -1000012 (ROS_V38.TF): parallel time TDT = 131461084.184 s + ticks / 65536.
# TDB - TDT (under 2 ms) is left out.
SCLK_TDT0 = 131461084.184
SCLK_TICKS_PER_S = 65536.0
J2000_JD = 2451545.0
DAY = 86400.0

TOL_DEG = 0.25
START_SPACING_DAYS = 32.0
MIN_SPACING_DAYS = 1.0
OUT = os.path.join(REPO, 'assets', 'ephem', 'c67p_spin.toml')


# --- Evaluation (as body_rotation in src/scene/scene.cpp) -----------------------------

def pole_axes(ra_deg, dec_deg):
    """Node direction, its quadrature and the pole (ICRF), as iau_pole_frame."""
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    z = (math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec))
    x = (-math.sin(ra), math.cos(ra), 0.0)
    y = (z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0])
    return x, y, z


def body_axes(ra_deg, dec_deg, w_deg):
    """Body x, y, z axes in ICRF."""
    x, y, z = pole_axes(ra_deg, dec_deg)
    w = math.radians(w_deg % 360.0)
    c, s = math.cos(w), math.sin(w)
    return (tuple(c * x[i] + s * y[i] for i in range(3)), tuple(-s * x[i] + c * y[i] for i in range(3)), z)


def hermite(s, h, w0, r0, w1, r1):
    """Value and derivative of the cubic Hermite segment at s in [0, 1]."""
    s2, s3 = s * s, s * s * s
    w = (2 * s3 - 3 * s2 + 1) * w0 + (s3 - 2 * s2 + s) * h * r0 + (-2 * s3 + 3 * s2) * w1 + (s3 - s2) * h * r1
    r = ((6 * s2 - 6 * s) * w0 + (6 * s - 6 * s2) * w1) / h + (3 * s2 - 4 * s + 1) * r0 + (3 * s2 - 2 * s) * r1
    return w, r


class Table:
    """Knots (t TDB s, pole RA deg, pole Dec deg, W deg, dW/dt deg/day)."""

    def __init__(self, knots):
        self.knots = knots
        self.times = [k[0] for k in knots]

    @staticmethod
    def load(path):
        with open(path, 'rb') as f:
            data = tomllib.load(f)
        return Table([((k[0] - J2000_JD) * DAY,) + tuple(k[1:]) for k in data['knots']])

    def state(self, t):
        """(RA, Dec, W, dW/dt) at t TDB."""
        k = self.knots
        if t <= k[0][0] or t >= k[-1][0]:
            e = k[0] if t <= k[0][0] else k[-1]
            return e[1], e[2], e[3] + e[4] * (t - e[0]) / DAY, e[4]
        i = bisect.bisect_right(self.times, t) - 1
        a, b = k[i], k[i + 1]
        h = (b[0] - a[0]) / DAY
        s = (t - a[0]) / (b[0] - a[0])
        w, r = hermite(s, h, a[3], a[4], b[3], b[4])
        return a[1] + (b[1] - a[1]) * s, a[2] + (b[2] - a[2]) * s, w, r

    def axes(self, t):
        ra, dec, w, _ = self.state(t)
        return body_axes(ra, dec, w)

    def spin(self, t):
        """Angular velocity (rad/s, ICRF), leaving out the slow motion of the pole."""
        ra, dec, _, rate = self.state(t)
        z = pole_axes(ra, dec)[2]
        return tuple(math.radians(rate) / DAY * c for c in z)


# --- Fit --------------------------------------------------------------------------------

def read_samples(directory):
    """(t TDB s, rows = body axes in ICRF) every 10 minutes, the later kernel winning."""
    files = []
    for name in KERNELS:
        samples = {}
        for seg in ck.Ck(os.path.join(directory, name)).segments:
            if seg.instrument != COMET_FRAME_ID or seg.frame != J2000_FRAME_ID:
                sys.exit('%s: unexpected segment' % name)
            for e, q in zip(seg.epochs, seg.quaternions):
                # Segments share their end epochs; the later one wins.
                samples[round(SCLK_TDT0 + e / SCLK_TICKS_PER_S, 1)] = q
        files.append(sorted(samples.items()))
    out = []
    for i, f in enumerate(files):
        until = files[i + 1][0][0] if i + 1 < len(files) else math.inf
        out += [(t, ck.rotation(q)) for t, q in f if t < until]
    return out


def solve_banded(a, b, band):
    """Solve the symmetric positive definite system a x = b (a: dense rows, zero outside the band)."""
    n = len(b)
    a = [row[:] for row in a]
    b = b[:]
    for i in range(n):
        for j in range(i + 1, min(n, i + band + 1)):
            f = a[j][i] / a[i][i]
            if f:
                for k in range(i, min(n, i + band + 1)):
                    a[j][k] -= f * a[i][k]
                b[j] -= f * b[i]
    x = [0.0] * n
    for i in reversed(range(n)):
        x[i] = (b[i] - sum(a[i][k] * x[k] for k in range(i + 1, min(n, i + band + 1)))) / a[i][i]
    return x


def interval(times, t):
    return min(max(bisect.bisect_right(times, t) - 1, 0), len(times) - 2)


def fit_linear(times, ts, values):
    """Least-squares piecewise-linear (hat) fit of values at ts, knots at times."""
    n = len(times)
    a = [[0.0] * n for _ in range(n)]
    b = [0.0] * n
    for t, v in zip(ts, values):
        i = interval(times, t)
        s = (t - times[i]) / (times[i + 1] - times[i])
        for p, wp in ((i, 1.0 - s), (i + 1, s)):
            b[p] += wp * v
            for q, wq in ((i, 1.0 - s), (i + 1, s)):
                a[p][q] += wp * wq
    return solve_banded(a, b, 1)


def fit_hermite(times, ts, values):
    """Least-squares cubic Hermite fit (values and slopes per day at the knots)."""
    n = 2 * len(times)
    a = [[0.0] * n for _ in range(n)]
    b = [0.0] * n
    for t, v in zip(ts, values):
        i = interval(times, t)
        h = (times[i + 1] - times[i]) / DAY
        s = (t - times[i]) / (times[i + 1] - times[i])
        s2, s3 = s * s, s * s * s
        terms = ((2 * i, 2 * s3 - 3 * s2 + 1), (2 * i + 1, (s3 - 2 * s2 + s) * h),
                 (2 * i + 2, -2 * s3 + 3 * s2), (2 * i + 3, (s3 - s2) * h))
        for p, wp in terms:
            b[p] += wp * v
            for q, wq in terms:
                a[p][q] += wp * wq
    x = solve_banded(a, b, 3)
    return x[0::2], x[1::2]


def fit(samples, times):
    ts = [t for t, _ in samples]
    poles = [r[2] for _, r in samples]
    ra = fit_linear(times, ts, [math.degrees(math.atan2(z[1], z[0])) for z in poles])
    dec = fit_linear(times, ts, [math.degrees(math.asin(z[2])) for z in poles])
    # W of the kernel's x axis about the fitted pole, unwrapped, less a mean rate
    # (for the conditioning of the fit).
    w = []
    for t, r in samples:
        i = interval(times, t)
        s = (t - times[i]) / (times[i + 1] - times[i])
        x, y, _ = pole_axes(ra[i] + (ra[i + 1] - ra[i]) * s, dec[i] + (dec[i + 1] - dec[i]) * s)
        angle = math.degrees(math.atan2(sum(r[0][k] * y[k] for k in range(3)), sum(r[0][k] * x[k] for k in range(3))))
        if w:
            angle = w[-1] + (angle - w[-1] + 180.0) % 360.0 - 180.0
        w.append(angle)
    rate0 = (w[-1] - w[0]) / ((ts[-1] - ts[0]) / DAY)
    base = [w[0] + rate0 * (t - ts[0]) / DAY for t in ts]
    dw, dr = fit_hermite(times, ts, [a - b for a, b in zip(w, base)])
    # Start with W in [0, 360).
    shift = -360.0 * math.floor((dw[0] + w[0] + rate0 * (times[0] - ts[0]) / DAY) / 360.0)
    knots = [(t, ra[k], dec[k], dw[k] + w[0] + rate0 * (t - ts[0]) / DAY + shift, dr[k] + rate0)
             for k, t in enumerate(times)]
    return Table(knots)


def errors(table, samples):
    """Angle (deg) between the fitted and the kernel axes at each sample."""
    out = []
    for t, r in samples:
        f = table.axes(t)
        trace = sum(f[k][i] * r[k][i] for k in range(3) for i in range(3))
        out.append(math.degrees(math.acos(max(-1.0, min(1.0, (trace - 1.0) / 2.0)))))
    return out


def calendar(t):
    jd = t / DAY + J2000_JD + 0.5
    z = int(jd)
    f = jd - z
    alpha = int((z - 1867216.25) / 36524.25)
    a = z + 1 + alpha - alpha // 4
    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    day = b - d - int(30.6001 * e)
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715
    minutes = round(f * 1440.0)
    return '%04d-%02d-%02d %02d:%02d' % (year, month, day, minutes // 60, minutes % 60)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    samples = read_samples(sys.argv[1])
    ts = [t for t, _ in samples]
    t0, t1 = ts[0], ts[-1]
    print('%d samples, %s .. %s TDB' % (len(samples), calendar(t0), calendar(t1)))
    n = max(1, round((t1 - t0) / DAY / START_SPACING_DAYS))
    times = [t0] + [round((t0 + (t1 - t0) * k / n) / 60.0) * 60.0 for k in range(1, n)] + [t1]
    while True:
        table = fit(samples, times)
        err = errors(table, samples)
        split = set()
        for t, e in zip(ts, err):
            i = interval(times, t)
            if e > TOL_DEG and times[i + 1] - times[i] > 2.0 * MIN_SPACING_DAYS * DAY:
                split.add(i)
        print('%d knots: max %.3f deg, rms %.3f deg; splitting %d' % (
            len(times), max(err), math.sqrt(sum(e * e for e in err) / len(err)), len(split)))
        if not split:
            break
        # Knots on whole minutes.
        times = sorted(times + [round((times[i] + times[i + 1]) / 120.0) * 60.0 for i in split])
    worst = max(range(len(err)), key=err.__getitem__)
    print('worst %.3f deg at %s' % (err[worst], calendar(ts[worst])))
    for k in table.knots[::max(1, len(table.knots) // 12)]:
        print('  %s  pole (%.3f, %.3f)  period %.4f h' % (calendar(k[0]), k[1], k[2], 360.0 / k[4] * 24.0))

    rms = math.sqrt(sum(e * e for e in err) / len(err))
    with open(OUT, 'w', encoding='utf-8', newline='\n') as f:
        f.write('# Rotation of comet 67P/Churyumov-Gerasimenko, fitted to ESA\'s comet attitude kernels\n'
                '# (frame 67P/C-G_CK) by tools/philae/c67p_spin.py; see assets/ephem/SOURCES.md.\n'
                '# Kernels: %s.\n'
                '# Off the kernels by %.3f deg rms, %.3f deg at most (where they jump).\n'
                '#\n'
                '# Knots: [JD TDB, pole RA (deg), pole Dec (deg), W (deg, counted on without wrapping),\n'
                '# dW/dt (deg/day)]. Between knots the pole is linear and W a cubic Hermite; outside them\n'
                '# the pole stays and W turns at the end rate.\n\n'
                % (', '.join(k.split('.')[0] for k in KERNELS), rms, max(err)))
        f.write('knots = [\n')
        for t, ra, dec, w, r in table.knots:
            f.write('    [%.7f, %.4f, %.4f, %.4f, %.6f],\n' % (t / DAY + J2000_JD, ra, dec, w, r))
        f.write(']\n')
    print('wrote', OUT)


if __name__ == '__main__':
    main()
