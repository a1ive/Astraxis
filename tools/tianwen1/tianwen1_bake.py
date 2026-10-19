"""Bake the Tianwen-1 orbiter's path around Mars (assets/ephem/tianwen1_mars.eph).

Usage:
    python tianwen1_bake.py

Sources (see assets/ephem/SOURCES.md):
  - 2021-02-08 .. 2021-11-16: tianwen1_reconstruct.py's path (cache/tw1_2021.txt), the
    state vectors in the spacecraft's telemetry (decoded by D. Estevez) propagated and
    joined at the maneuvers.
  - 2021-11-16 .. 2025-12-17: the MOMAG positions (momag_fetch.py, cache/), MSO
    converted to ICRF with Mars' heliocentric state from JPL Horizons.

Data gaps longer than GAP_S are gaps in the table (the scene hides the orbiter there);
shorter ones are bridged by the interpolation. Kepler-relative (AXEPH2) around the Mars
system GM, as the Horizons tables of mars.toml.
"""

import datetime
import glob
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'tools', 'bake'))
import horizons_bake as hb  # noqa: E402

GM_MARS = 42828.375816            # DE440 Mars system GM (JPL SSD astrodynamic parameters)
TDB_MINUS_UTC = 69.184            # 37 leap seconds (since 2017) + 32.184 s
TOL_KM = 5.0
GAP_S = 7200.0                    # longer data gaps are gaps in the table
MAX_SPAN_S = GAP_S - 300.0        # knot spacing inside a segment stays below GAP_S
DECIMATE_EVERY = 2                # use every 2nd cached (1-minute) MOMAG sample
# Every sample is a knot over these spans (TDB s): the Mars orbit insertion burn (Estevez's
# model, tianwen1_reconstruct.py: 888 s from 2021-02-10 11:52:55 UTC) with a minute on
# either side. Decimated to the tolerance, the whole burn was one interval between two
# knots, whose cubic bent the thrust profile and lost 2% of its delta-v.
MOI_START_TDB = (29255.9950810185185185 + 2430000.0 - 2451545.0) * 86400.0 + TDB_MINUS_UTC  # GMAT UTCModJulian
KEEP_ALL = [(MOI_START_TDB - 60.0, MOI_START_TDB + 888.0 + 60.0)]
OUT = os.path.join(ROOT, 'assets', 'ephem', 'tianwen1_mars.eph')
NAIF_ID = -9901491                # Horizons' (unofficial) Tianwen-1 record
RECONSTRUCTED = os.path.join(HERE, 'cache', 'tw1_2021.txt')

J2000 = datetime.datetime(2000, 1, 1, 12)


def tdb_from_utc(dt):
    return (dt - J2000).total_seconds() + TDB_MINUS_UTC


# --- Frames -------------------------------------------------------------------

def combine(axes, v):
    x, y, z = axes
    return tuple(v[0] * x[i] + v[1] * y[i] + v[2] * z[i] for i in range(3))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    n = math.sqrt(sum(c * c for c in a))
    return tuple(c / n for c in a)


class MarsHelio:
    """Mars' heliocentric state (Horizons 499 @ Sun, ICRF), hourly, linearly interpolated."""

    def __init__(self, t0, t1):
        jd0 = (t0 - 3600.0) / 86400.0 + hb.J2000_JD
        jd1 = (t1 + 3600.0) / 86400.0 + hb.J2000_JD
        self.rows = hb.fetch('499', '500@10', jd0, jd1, 60)
        self.t0 = self.rows[0][0]

    def state(self, t):
        k = min(max(int((t - self.t0) // 3600.0), 0), len(self.rows) - 2)
        (ta, ra, va), (tb, rb, vb) = self.rows[k], self.rows[k + 1]
        u = (t - ta) / (tb - ta)
        return (tuple(ra[i] + u * (rb[i] - ra[i]) for i in range(3)),
                tuple(va[i] + u * (vb[i] - va[i]) for i in range(3)))

    def mso(self, t):
        """MSO axes in ICRF and MSO's angular velocity (rad/s, ICRF)."""
        r, v = self.state(t)
        x = norm(tuple(-c for c in r))
        h = cross(r, v)
        z = norm(h)
        y = cross(z, x)
        r2 = sum(c * c for c in r)
        omega = tuple(c / r2 for c in h)  # the Mars-Sun line turns at h / r^2
        return (x, y, z), omega


# --- Sources ------------------------------------------------------------------

def reconstructed_rows():
    """tianwen1_reconstruct.py's 2021 path (already ICRF, TDB), one row a minute."""
    rows = []
    for line in open(RECONSTRUCTED):
        if line.startswith('#'):
            continue
        f = [float(c) for c in line.split()]
        rows.append((f[0], tuple(f[1:4]), tuple(f[4:7])))
    print(f'2021 reconstruction: {len(rows)} rows, {hb.calendar(rows[0][0])} .. {hb.calendar(rows[-1][0])}')
    return rows


def momag_rows():
    files = sorted(glob.glob(os.path.join(HERE, 'cache', '[0-9]' * 8 + '.txt')))   # one a day
    raw = []
    for path in files:
        for line in open(path, encoding='ascii'):
            if line.startswith('#'):
                continue
            f = line.split()
            t = tdb_from_utc(datetime.datetime.fromisoformat(f[0]))
            raw.append((t, tuple(float(c) for c in f[1:4]), tuple(float(c) for c in f[4:7])))
    raw.sort(key=lambda r: r[0])
    print(f'MOMAG: {len(files)} days, {len(raw)} samples, {hb.calendar(raw[0][0])} .. {hb.calendar(raw[-1][0])}')
    mars = MarsHelio(raw[0][0], raw[-1][0])
    out = []
    for t, p, v in raw:
        axes, omega = mars.mso(t)
        r = combine(axes, p)
        w = combine(axes, v)
        vi = tuple(w[i] + c for i, c in enumerate(cross(omega, r)))  # frame rotation
        out.append((t, r, vi))
    return out


# --- Reduction ------------------------------------------------------------------

def segments(rows, gap_s):
    seg = [rows[0]]
    for r in rows[1:]:
        if r[0] - seg[-1][0] > gap_s:
            yield seg
            seg = []
        seg.append(r)
    yield seg


def decimate_capped(rows, tol_km, interp, max_span_s):
    """horizons_bake.decimate with the knot spacing capped (so that it stays below the
    table's gap threshold)."""
    def ok(a, b):
        if rows[b][0] - rows[a][0] > max_span_s:
            return False
        for j in range(a + 1, b):
            if hb.dist(interp(rows[a], rows[b], rows[j][0]), rows[j][1]) > tol_km:
                return False
        return True

    kept = [0]
    a = 0
    n = len(rows)
    while a < n - 1:
        step = 1
        while a + step * 2 < n and ok(a, a + step * 2):
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
    return [rows[i] for i in kept]


def decimate_keeping(rows, tol_km, interp, max_span_s, keep):
    """decimate_capped, but every row within the `keep` spans (t0, t1) is a knot."""
    out = []
    start = 0
    for t0, t1 in keep:
        inside = [i for i, r in enumerate(rows) if t0 <= r[0] <= t1]
        if not inside or inside[0] < start:
            continue
        lo, hi = max(inside[0] - 1, start), min(inside[-1] + 1, len(rows) - 1)  # samples bracketing the span
        out.extend(decimate_capped(rows[start:lo + 1], tol_km, interp, max_span_s)[:-1] if lo > start else [])
        out.extend(rows[lo:hi])
        start = hi
    out.extend(decimate_capped(rows[start:], tol_km, interp, max_span_s) if start < len(rows) - 1 else rows[start:])
    return out


def main():
    interp = hb.make_interp(GM_MARS)
    rows = reconstructed_rows()
    momag = momag_rows()
    rows.extend(momag[::DECIMATE_EVERY] if DECIMATE_EVERY > 1 else momag)
    rows.sort(key=lambda r: r[0])
    rows = hb.drop_velocity_spikes(rows, max_spacing_s=400.0)
    knots, err, covered = [], 0.0, 0.0
    segs = list(segments(rows, GAP_S))
    for seg in segs:
        if len(seg) < 2:
            continue
        part = decimate_keeping(seg, TOL_KM, interp, MAX_SPAN_S, KEEP_ALL)
        err = max(err, hb.max_error(seg, part, interp))
        knots.extend(part)
        covered += seg[-1][0] - seg[0][0]
    spacing = max(b[0] - a[0] for a, b in zip(knots, knots[1:]) if b[0] - a[0] <= GAP_S)
    assert spacing < GAP_S, spacing
    hb.write_eph(OUT, NAIF_ID, 499, knots, GM_MARS, int(GAP_S))
    print(f'{len(segs)} segments, {covered / 86400:.0f} days covered, {len(rows)} samples -> {len(knots)} knots, '
          f'max error {err:.3f} km; wrote {OUT} ({os.path.getsize(OUT) // 1024} KB)')


if __name__ == '__main__':
    main()
