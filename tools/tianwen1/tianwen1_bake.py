"""Bake the Tianwen-1 orbiter's path around Mars (assets/ephem/tianwen1_mars.eph).

Usage:
    python tianwen1_bake.py [path to the jupyter_notebooks/Tianwen/orbit folder]

Sources (see assets/ephem/SOURCES.md):
  - 2021-02 .. 2021-11: D. Estevez's GMAT reports around the maneuvers (github.com/
    daniestevez/jupyter_notebooks, Tianwen/orbit/*/Tianwen1_pre_*.txt and _post_*.txt):
    trajectories propagated from the state vectors in the spacecraft's telemetry,
    forward from before a burn and backward from after it. They cross at the burn; each
    window takes the pre-burn arc up to there and the post-burn arc after it. Columns:
    UTC modified Julian date (GMAT, JD - 2430000), position and velocity in GMAT's
    "MarsInertial" frame, whose axes are BodyInertial: z along Mars' spin axis at J2000
    (IAU 2000: RA 317.68143, Dec 52.88650 deg, the pole GMAT uses; with it a telemetry
    vector agrees with MOMAG to 6 km, with the IAU 2015 pole to 81 km), x toward the
    ascending node of Mars' equator on the Earth's mean equator of J2000.
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
GMAT_MJD_OFFSET = 2430000.0       # GMAT's UTCModJulian = JD - 2430000
POLE_RA_DEC_IAU2000 = (317.68143, 52.88650)
TOL_KM = 5.0
GAP_S = 7200.0                    # longer data gaps are gaps in the table
MAX_SPAN_S = GAP_S - 300.0        # knot spacing inside a segment stays below GAP_S
DECIMATE_EVERY = 2                # use every 2nd cached (1-minute) MOMAG sample
OUT = os.path.join(ROOT, 'assets', 'ephem', 'tianwen1_mars.eph')
NAIF_ID = -9901491                # Horizons' (unofficial) Tianwen-1 record
WINDOWS = [                        # (folder, pre file, post file)
    ('plane_change', 'Tianwen1_pre_planechange.txt', 'Tianwen1_post_planechange.txt'),
    ('phasing', 'Tianwen1_pre_phasing.txt', 'Tianwen1_post_phasing.txt'),
    ('deorbit', 'Tianwen1_pre_deorbit.txt', 'Tianwen1_post_deorbit.txt'),
    ('avoidance', 'Tianwen1_pre_avoidance.txt', 'Tianwen1_post_avoidance.txt'),
    ('remote_sensing', 'Tianwen1_pre_remote_sensing.txt', 'Tianwen1_post_remote_sensing.txt'),
]

J2000 = datetime.datetime(2000, 1, 1, 12)


def tdb_from_utc(dt):
    return (dt - J2000).total_seconds() + TDB_MINUS_UTC


# --- Frames -------------------------------------------------------------------

def body_inertial_axes(ra_deg, dec_deg):
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    z = (math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec))
    x = (-math.sin(ra), math.cos(ra), 0.0)
    y = (z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0])
    return x, y, z


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

def read_gmat_report(path):
    rows = []
    for line in open(path):
        f = line.split()
        if len(f) < 7:
            continue
        mjd = float(f[0])
        t = (mjd + GMAT_MJD_OFFSET - hb.J2000_JD) * 86400.0 + TDB_MINUS_UTC
        rows.append((t, tuple(float(c) for c in f[1:4]), tuple(float(c) for c in f[4:7])))
    rows.sort(key=lambda r: r[0])
    # Consecutive propagation commands repeat their boundary sample.
    return [r for i, r in enumerate(rows) if i == 0 or r[0] - rows[i - 1][0] > 1e-3]


def window(folder, pre_name, post_name, orbit_dir):
    """Pre-burn arc up to the burn, post-burn arc from it (ICRF). The burn: the time of
    closest approach of the two arcs where both exist."""
    axes = body_inertial_axes(*POLE_RA_DEC_IAU2000)
    pre = read_gmat_report(os.path.join(orbit_dir, folder, pre_name))
    post = read_gmat_report(os.path.join(orbit_dir, folder, post_name))
    t_lo, t_hi = max(pre[0][0], post[0][0]), min(pre[-1][0], post[-1][0])
    post_table = hb.Table(post)
    best = None
    for t, p, _ in pre:
        if t_lo <= t <= t_hi:
            d = hb.dist(p, hb.eval_eph(post_table, t))
            if best is None or d < best[0]:
                best = (d, t)
    if best is None:  # no overlap (the deorbit arcs): join where one ends
        t_burn = pre[-1][0]
        print(f'  {folder}: arcs do not overlap; joined at {hb.calendar(t_burn)}')
    else:
        t_burn = best[1]
        print(f'  {folder}: burn at {hb.calendar(t_burn)} (arcs {best[0]:.2f} km apart there)')
    rows = [r for r in pre if r[0] < t_burn] + [r for r in post if r[0] >= t_burn]
    return [(t, combine(axes, p), combine(axes, v)) for t, p, v in rows]


def momag_rows():
    files = sorted(glob.glob(os.path.join(HERE, 'cache', '*.txt')))
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


def main():
    orbit_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'cache', 'estevez_orbit')
    interp = hb.make_interp(GM_MARS)
    rows = []
    print('Estevez windows:')
    for w in WINDOWS:
        rows.extend(window(*w, orbit_dir))
    momag = momag_rows()
    rows.extend(momag[::DECIMATE_EVERY] if DECIMATE_EVERY > 1 else momag)
    rows.sort(key=lambda r: r[0])
    rows = hb.drop_velocity_spikes(rows, max_spacing_s=400.0)
    knots, err, covered = [], 0.0, 0.0
    segs = list(segments(rows, GAP_S))
    for seg in segs:
        if len(seg) < 2:
            continue
        part = decimate_capped(seg, TOL_KM, interp, MAX_SPAN_S)
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
