"""Reconstructs the trajectory of ISEE-3 / ICE from launch (1978) to its return in 2014.

JPL Horizons only has ISEE-3 for 2014 (its return; baked by tools/bake/ephem.toml as
isee3_2014), so the scene's 1978-2013 trajectory is rebuilt here from public data:

0. SSCWeb's "GEI J2000" ISEE-3 vectors carry the nutation of date (found on its 2008-2014
   track, see 5.); it is removed from all of them first.

1. 1978-08-16 .. 1983-12-31 (L1 halo orbit, geotail passes, five lunar flybys): geocentric
   GEI J2000 positions from NASA GSFC SSCWeb (satellite "isee3", 12-minute resolution,
   https://sscweb.gsfc.nasa.gov/WS/sscr/2). Only the first sample of every hour is used
   (the four after it are extrapolated with a wrong velocity, a sawtooth of up to thousands
   of km). The track is split into ballistic arcs between maneuvers by continuation fitting:
   an arc is extended while the residuals of the newest stretch stay below
   max(50 km, 3 x the scatter of a 4-day fit just before it). Each arc is a least-squares
   orbit through the data (Earth with J2, Moon, Sun, Venus, Mars, Jupiter; outliers rejected).
   The data contain stretches shifted by thousands of km (orbit solutions concatenated, e.g.
   every early January); arcs fitted to such a stretch are dropped when the neighbours meet
   without it. Consecutive arcs are joined where they come closest, and blended over a few
   days where they do not meet (their own data jump).
   The first arc is propagated back to perigee shortly after the launch.
2. 1983-12-31 .. 1985-09-11: ballistic from the last arc (after the 1983-12-22 lunar flyby),
   with two impulses: 1985-06-05, "a major trajectory maneuver", and the last trim three days
   before the encounter, 1985-09-08 (Farquhar 2001, J. Astronaut. Sci. 49, 23, pp. 47-49:
   June 5 ~39 m/s, July 9 ~1 m/s - folded into the first here - and September 8 2.3 m/s).
   Their size is solved so that the path joins the JPL navigation trajectory below.
3. 1985-09-10 .. 14: JPL Navigation Section "save tape" (PDS SBN, ICE-C-PLAWAV-3-RDR-ESP-
   GIACOBIN-ZIN-V1.0, GEOMETRY/TRAJ_ICE.TBL): heliocentric states in the ecliptic and equinox
   of date, rotated to ICRF by matching the file's Earth (heliocentric minus geocentric
   spacecraft) to Horizons' Earth. A ballistic arc through them continues to late 1986.
4. 1986 .. 2013: a ballistic arc through Horizons' 2014-01 .. 06 states (solution s45),
   propagated back and joined to 3. with two impulses (1986 and 1987: "Maneuvers in 1986
   approximately targeted a 2014 August lunar swingby", Dunham, Farquhar et al., IAC-14; the
   2014 firing was the first "since 1987"; the dates are the least-delta-v pair of a scan).
5. Radiation pressure for all heliocentric arcs: fitted to SSCWeb's 2008-2013 track (an
   orbit prediction made before the 2014 recovery, ~2 h along-track from Horizons' s45).
6. Comet 21P/Giacobini-Zinner: Horizons (SAO/1985 elements), moved within a few weeks of the
   encounter onto the navigation file's comet (the two differ by ~270 km).

Outputs (assets/ephem): isee3_geo.eph (geocentric, AXEPH1, launch .. 1984-01-02),
isee3_helio.eph (heliocentric, AXEPH2, 1983-12-31 .. 2014-01-03), gz.eph (comet, AXEPH2). Downloads are cached in tools/isee3/cache (not in git).

Usage: python tools/isee3/isee3_reconstruct.py [--refit]
Needs numpy. Fitting takes ~15 minutes on 8 cores (cached afterwards; --refit redoes it).
"""
import datetime as dt
import json
import multiprocessing
import os
import sys
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'bake'))
import horizons_bake as hb  # noqa: E402
import isee3_dynamics as dyn  # noqa: E402

CACHE = os.path.join(HERE, 'cache')
OUT = os.path.normpath(os.path.join(HERE, '..', '..', 'assets', 'ephem'))
DAY = dyn.DAY

DATA_START = dt.datetime(1978, 8, 16)
DATA_END = dt.datetime(1983, 12, 31, 18)            # later SSCWeb data jump by 10^4 km
GEO_END = dyn.utc_to_tdb(dt.datetime(1984, 1, 2))   # geocentric table end (overlaps the next)
HELIO_START = dyn.utc_to_tdb(dt.datetime(1983, 12, 31))
HELIO_END = dyn.utc_to_tdb(dt.datetime(2014, 1, 3))      # overlaps the Horizons 2014 table
BURNS = [dyn.utc_to_tdb(dt.datetime(1985, 6, 5, 12)), dyn.utc_to_tdb(dt.datetime(1985, 9, 8, 12))]
# After Halley: "Maneuvers in 1986 approximately targeted a 2014 August lunar swingby" (Dunham,
# Farquhar et al., IAC-14.B6.3.4), and the thrusters fired in 2014 "for the first time since
# 1987". Dates: the pair with the least total delta-v in a 30-day scan (big one in 1986).
LATE_BURNS = [dyn.utc_to_tdb(dt.datetime(1986, 10, 8)), dyn.utc_to_tdb(dt.datetime(1987, 5, 1))]
SRP_NOMINAL = 4.0e-11   # km/s^2 at 1 au: Cr 1.3 x 4.56e-6 N/m^2 x 2.8 m^2 (1.74 m x 1.61 m) / 400 kg
TOL_KM = 1.0
ISEE3_ID, COMET_ID = -111, 90000319   # NAIF id; Horizons record of 21P (SAO/1985)


# --- Downloads ------------------------------------------------------------------------

SSC_URL = 'https://sscweb.gsfc.nasa.gov/WS/sscr/2/locations'
SSC_REQUEST = """<?xml version="1.0" encoding="UTF-8"?>
<DataRequest xmlns="http://sscweb.gsfc.nasa.gov/schema">
  <TimeInterval><Start>{start}</Start><End>{end}</End></TimeInterval>
  <Satellites><Id>isee3</Id><ResolutionFactor>1</ResolutionFactor></Satellites>
  <OutputOptions><AllLocationFilters>true</AllLocationFilters>
    <CoordinateOptions><CoordinateSystem>GeiJ2000</CoordinateSystem><Component>X</Component></CoordinateOptions>
    <CoordinateOptions><CoordinateSystem>GeiJ2000</CoordinateSystem><Component>Y</Component></CoordinateOptions>
    <CoordinateOptions><CoordinateSystem>GeiJ2000</CoordinateSystem><Component>Z</Component></CoordinateOptions>
  </OutputOptions>
</DataRequest>"""
PDS_URL = ('https://pdssbn.astro.umd.edu/holdings/ice-c-plawav-3-rdr-esp-giacobin-zin-v1.0/'
           'geometry/traj_ice.tbl')


def cached(name):
    return os.path.join(CACHE, name)


def fetch_ssc(name='ssc_isee3.csv', start=None, end=None):
    """SSCWeb positions over [start, end] (default: the 1978-1983 span), cached as CSV."""
    path = cached(name)
    if os.path.exists(path):
        return path
    start, end = start or DATA_START, end or DATA_END
    rows = {}
    a = start
    for attempt in range(1000):
        if a >= end:
            break
        b = min(a + dt.timedelta(days=60), end)
        body = SSC_REQUEST.format(start=a.strftime('%Y-%m-%dT%H:%M:%S.000Z'),
                                  end=b.strftime('%Y-%m-%dT%H:%M:%S.000Z'))
        req = urllib.request.Request(SSC_URL, data=body.encode(),
                                     headers={'Content-Type': 'application/xml', 'Accept': 'application/json'})
        try:
            d = json.load(urllib.request.urlopen(req, timeout=600))
        except OSError as e:  # transient network errors: retry the chunk
            print(f'  SSCWeb {a:%Y-%m-%d}: {e}, retrying', flush=True)
            continue
        sd = d[1]['Result'][1]['Data'][1][0][1]
        c = sd['Coordinates'][1][0][1]
        for t, x, y, z in zip([x[1] for x in sd['Time'][1]], c['X'][1], c['Y'][1], c['Z'][1]):
            rows[t[:19]] = (x, y, z)
        print(f'  SSCWeb {a:%Y-%m-%d}: {len(rows)} samples', flush=True)
        a = b
    os.makedirs(CACHE, exist_ok=True)
    with open(path, 'w') as f:
        for t in sorted(rows):
            f.write('%s,%.3f,%.3f,%.3f\n' % ((t,) + rows[t]))
    return path


def fetch_bodies():
    jd0, jd1 = hb.jd_from_text('1978-01-01'), hb.jd_from_text('2015-01-01')
    targets = [('moon_geo', '301', '500@399', 120), ('sun_geo', '10', '500@399', 1440)]
    targets += [(n + '_hel', c, '500@10', 1440 if n in ('mercury', 'venus', 'earth') else 5760)
                for n, c in [('mercury', '199'), ('venus', '299'), ('earth', '399'), ('mars', '4'),
                             ('jupiter', '5'), ('saturn', '6'), ('uranus', '7'), ('neptune', '8')]]
    for name, command, center, step in targets:
        path = cached(f'hzn_{name}.npy')
        if not os.path.exists(path):
            rows = hb.fetch(command, center, jd0, jd1, step)
            np.save(path, np.array([[t, *p, *v] for t, p, v in rows]))
            print(f'  Horizons {name}: {len(rows)} rows', flush=True)


def fetch_file(url, name):
    path = cached(name)
    if not os.path.exists(path):
        urllib.request.urlretrieve(url, path)
    return path


def load_obs():
    """Hourly SSCWeb samples (the first of each hour): TDB seconds, GEI J2000 km. SSCWeb's
    "GEI J2000" ISEE-3 vectors carry the nutation of date (up to ~20": ~100 km at L1, ~10^4 km
    at 1 au in 2008-2014), which is removed here."""
    t, p = [], []
    for line in open(fetch_ssc()):
        s, x, y, z = line.strip().split(',')
        if s[14:16] != '08':
            continue
        d = dt.datetime.fromisoformat(s)
        if d > DATA_END:
            continue
        t.append(dyn.utc_to_tdb(d))
        p.append((float(x), float(y), float(z)))
    t = np.array(t)
    return t, dyn.remove_nutation(t, np.array(p))


# --- Segmentation -------------------------------------------------------------------------

_W = {}


def _worker_init():
    b = dyn.Bodies(CACHE)
    _W['model'] = dyn.GeoModel(b)
    _W['t'], _W['p'] = load_obs()
    _W['noise'] = {}


def _fit_window(a, b, te, xe):
    t, p = _W['t'], _W['p']
    s = (t >= a) & (t <= b)
    try:
        xe, rn = dyn.fit_arc(_W['model'], t[s], p[s], te, xe)
    except (FloatingPointError, np.linalg.LinAlgError, ValueError):
        return None, None, None
    return xe, t[s], rn


def _local_noise(b):
    """Median residual of a 4-day ballistic fit ending at b: the scatter of the data there."""
    key = round(b / 3600.0)
    if key not in _W['noise']:
        te, xe = dyn.guess(_W['t'], _W['p'], b - 2 * DAY)
        xe, _, rn = _fit_window(b - 4 * DAY, b, te, xe)
        _W['noise'][key] = np.median(rn) if xe is not None and len(rn) else 1e9
    return _W['noise'][key]


def _segment(span, thresh=50.0):
    t0, t_stop = span
    t, p, model = _W['t'], _W['p'], _W['model']
    t0, t_stop = max(t0, t[0]), min(t_stop, t[-1])
    arcs = []
    a = t0
    while a < t_stop - 2 * DAY:
        b = min(a + 3 * DAY, t_stop)
        te, xe = dyn.guess(t, p, a + 1.5 * DAY)
        xe, _, _ = _fit_window(a, b, te, xe)
        good = (te, xe, b)
        step = 2 * DAY
        noises = []
        while b < t_stop:
            nb = min(b + step, t_stop)
            g_te = 0.5 * (a + nb)
            g_xe = dyn.state_at(model, good[0], good[1], g_te)
            xe2, tq, rn = _fit_window(a, nb, g_te, g_xe)
            ok = xe2 is not None
            if ok:
                new = tq > b
                med_new = np.median(rn[new]) if new.any() else 0.0
                # Scatter just before the new stretch (cannot straddle a maneuver in it).
                noise = _local_noise(b)
                ok = (med_new < max(thresh, 3.0 * noise) and
                      np.median(rn) < max(thresh, 2.0 * np.mean(noises + [noise])))
                if ok:
                    noises.append(noise)
            if ok:
                good = (g_te, xe2, nb)
                b = nb
                step = min(step * 1.5, 8 * DAY)
            elif step > 3 * 3600:
                step *= 0.5
            else:
                break
        te, xe, b = good
        s = (t >= a) & (t <= b)
        _, rn = dyn.fit_arc(model, t[s], p[s], te, xe, iters=2)
        arcs.append(dict(start=a, end=b, epoch=te, state=[float(v) for v in xe], median=float(np.median(rn))))
        print(f'  arc {dyn.date(a)} .. {dyn.date(b)}  median residual {arcs[-1]["median"]:.0f} km', flush=True)
        if b >= t_stop:
            break
        a = b - 3600.0
    return arcs


def segment(refit):
    path = cached('arcs.jsonl')
    if os.path.exists(path) and not refit:
        return [json.loads(l) for l in open(path)]
    years = [dt.datetime(1978, 8, 16), dt.datetime(1979, 6, 1), dt.datetime(1980, 6, 1), dt.datetime(1981, 6, 1),
             dt.datetime(1982, 6, 1), dt.datetime(1983, 1, 1), dt.datetime(1983, 7, 1), DATA_END]
    spans = [(dyn.utc_to_tdb(a), dyn.utc_to_tdb(b)) for a, b in zip(years, years[1:])]
    with multiprocessing.Pool(len(spans), initializer=_worker_init) as pool:
        arcs = [a for part in pool.map(_segment, spans) for a in part]
    arcs.sort(key=lambda a: a['start'])
    with open(path, 'w') as f:
        for a in arcs:
            f.write(json.dumps(a) + '\n')
    return arcs


# --- Joining the arcs ---------------------------------------------------------------------

def closest(model, A, B, after=-1e300):
    """Time (after `after`) where arcs A and B come closest near their common end, the gap
    and the velocity difference."""
    lo, hi = max(A['start'], B['start'] - 3 * DAY), min(A['end'] + 3 * DAY, B['end'])
    if hi <= lo:  # an arc between them was dropped
        lo, hi = A['end'] - DAY, B['start'] + DAY
    lo = max(lo, after)
    tt = np.arange(lo, hi, 600.0)
    sa = dyn.propagate(model, A['epoch'], np.array(A['state']), 0.0, tt)
    sb = dyn.propagate(model, B['epoch'], np.array(B['state']), 0.0, tt)
    d = np.linalg.norm(sa[:, :3] - sb[:, :3], axis=1)
    i = int(np.argmin(d))
    return tt[i], d[i], np.linalg.norm(sa[i, 3:6] - sb[i, 3:6])


def drop_shifted(model, arcs, jump_km=1000.0):
    """Drops short arcs fitted to shifted data: both neighbours jump to them, but not to
    each other."""
    k = 1
    while k < len(arcs) - 1:
        A, M, B = arcs[k - 1], arcs[k], arcs[k + 1]
        if M['end'] - M['start'] < 40 * DAY:
            j1, j2 = closest(model, A, M)[1], closest(model, M, B)[1]
            if min(j1, j2) > jump_km and closest(model, A, B)[1] < 0.25 * min(j1, j2):
                print(f'  dropped {dyn.date(M["start"])} .. {dyn.date(M["end"])} (data shifted by '
                      f'{min(j1, j2):.0f} km)')
                del arcs[k]
                k = max(k - 1, 1)
                continue
        k += 1
    return arcs


def moon_distance(b, t, r):
    return np.linalg.norm(r - b.moon_geo(t), axis=1)


class Composite:
    """The arcs joined at their junctions, blended where they do not meet."""

    def __init__(self, model, bodies, arcs, t_start, t_end):
        self.model, self.arcs = model, arcs
        self.joins = []   # (time, half width, gap km, dv m/s)
        for A, B in zip(arcs, arcs[1:]):
            tj, gap, dv = closest(model, A, B, self.joins[-1][0] + 2 * 3600.0 if self.joins else t_start)
            hw = float(np.clip(gap / 0.04, 6 * 3600.0, 5 * DAY))   # blend at <= ~0.02 km/s
            self.joins.append([tj, hw, gap, dv * 1e3])
        bounds = [t_start] + [j[0] for j in self.joins] + [t_end]
        for k, j in enumerate(self.joins):
            j[1] = min(j[1], 0.45 * (bounds[k + 1] - bounds[k]), 0.45 * (bounds[k + 2] - bounds[k + 1]))
            # keep lunar flybys inside one ballistic arc
            while j[1] > 3600.0:
                tt = np.linspace(j[0] - j[1], j[0] + j[1], 200)
                ra = dyn.propagate(model, arcs[k]['epoch'], np.array(arcs[k]['state']), 0.0, tt)[:, :3]
                if moon_distance(bodies, tt, ra).min() > 7.0e4:
                    break
                j[1] *= 0.5
        self.bounds = bounds

    def states(self, tq):
        tq = np.asarray(tq, float)
        out = np.zeros((len(tq), 6))
        k = np.clip(np.searchsorted(self.bounds, tq, side='right') - 1, 0, len(self.arcs) - 1)
        cache = {}

        def arc_states(i, sel):
            a = self.arcs[i]
            return dyn.propagate(self.model, a['epoch'], np.array(a['state']), 0.0, tq[sel])

        for i in range(len(self.arcs)):
            sel = k == i
            if sel.any():
                out[sel] = arc_states(i, sel)
        for i, (tj, hw, _, _) in enumerate(self.joins):
            sel = (tq > tj - hw) & (tq < tj + hw)
            if not sel.any():
                continue
            sa, sb = arc_states(i, sel), arc_states(i + 1, sel)
            s = (tq[sel] - (tj - hw)) / (2 * hw)
            w = (s * s * (3 - 2 * s))[:, None]
            dw = (6 * s * (1 - s) / (2 * hw))[:, None]
            out[sel, :3] = (1 - w) * sa[:, :3] + w * sb[:, :3]
            out[sel, 3:] = (1 - w) * sa[:, 3:] + w * sb[:, 3:] + dw * (sb[:, :3] - sa[:, :3])
        return out


def sample_times(fn, t0, t1, near):
    """30-minute samples, 1-minute where near(t, states) says so."""
    tt = np.arange(t0, t1, 1800.0)
    tt = np.append(tt, t1)
    st = fn(tt)
    fine = near(tt, st)
    extra = []
    for i in np.where(fine[:-1] | fine[1:])[0]:
        extra.append(np.arange(tt[i] + 60.0, tt[i + 1], 60.0))
    return np.unique(np.concatenate([tt] + extra))


def to_rows(tt, st):
    return [(float(t), tuple(map(float, s[:3])), tuple(map(float, s[3:6]))) for t, s in zip(tt, st)]


def reduce(rows, gm):
    interp = hb.make_interp(gm)
    knots = hb.decimate(rows, TOL_KM, interp)
    return knots, hb.max_error(rows, knots, interp)


# --- JPL navigation file around the comet ---------------------------------------------------

def read_nav(bodies):
    """ICE (heliocentric) and comet-minus-ICE states, ICRF, from TRAJ_ICE.TBL.
    Each record is in the ecliptic and equinox of its own date: the rotation is a fixed one
    times a turn about the ecliptic pole per record, both from the file's Earth."""
    data = open(fetch_file(PDS_URL, 'traj_ice.tbl')).read()
    rec = 2736   # 114 items x 24 characters, no separators (TRAJ_ICE.DOC)
    rows = []
    for k in range(len(data) // rec):
        r = data[k * rec:(k + 1) * rec]
        rows.append([float(r[i * 24:(i + 1) * 24].replace('D', 'E')) if i != 2 else 0.0 for i in range(114)])
    a = np.array(rows)
    t = a[:, 0] - 18262.5 * DAY       # item 1: ET seconds past 1950-01-01 0h
    helio, geo, rel = a[:, 40:46], a[:, 34:40], a[:, 46:52]   # items 41-46, 35-40, 47-52
    earth_nav = helio[:, 0:3] - geo[:, 0:3]
    earth = bodies.helio['earth'](t)

    def rz(x):
        c, s = np.cos(x), np.sin(x)
        return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

    theta = np.zeros(len(t))
    for _ in range(3):
        P = np.array([rz(th) @ e for th, e in zip(theta, earth_nav)])
        U, _, Vt = np.linalg.svd(P.T @ earth)
        R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
        local = earth @ R
        theta = np.arctan2(local[:, 1], local[:, 0]) - np.arctan2(earth_nav[:, 1], earth_nav[:, 0])
        theta -= theta.mean()
    rots = np.array([R @ rz(th) for th in theta])

    def rot(v):
        return np.einsum('nij,nj->ni', rots, v)

    err = np.linalg.norm(rot(earth_nav) - earth, axis=1).max()
    ice = np.hstack([rot(helio[:, 0:3]), rot(helio[:, 3:6])])
    # items 47-52 are "body 1 - S/C" by name, but grow with the spacecraft's motion: S/C - comet
    comet = ice[:, 0:3] - rot(rel[:, 0:3])
    return t, ice, comet, err


# --- After the comets: 1987 .. 2014 -----------------------------------------------------------

def fit_srp(helio, bodies):
    """Solar radiation pressure from SSCWeb's 2008-2013 track (an orbit prediction made before
    the 2014 recovery): one ballistic arc with the pressure free."""
    path = fetch_ssc('ssc_isee3_2008.csv', dt.datetime(2008, 7, 30), dt.datetime(2014, 8, 11))
    t, p = [], []
    for line in open(path):
        s, x, y, z = line.strip().split(',')
        if s[11:16] not in ('00:00', '12:00') or s >= '2014':
            continue
        t.append(dyn.utc_to_tdb(dt.datetime.fromisoformat(s)))
        p.append((float(x), float(y), float(z)))
    t = np.array(t)
    p = dyn.remove_nutation(t, np.array(p)) + bodies.helio['earth'](t)
    i = len(t) // 2
    xe = np.concatenate([p[i], (p[i + 1] - p[i - 1]) / (t[i + 1] - t[i - 1])])
    srp = SRP_NOMINAL
    for _ in range(12):
        st = dyn.propagate(helio, t[i], xe, srp, t, True)
        res = p - st[:, :3]
        J = st[:, 6:48].reshape(-1, 6, 7)[:, 0:3, :].reshape(-1, 7)
        sc = np.linalg.norm(J, axis=0)
        dx = np.linalg.lstsq(J / sc, res.reshape(-1), rcond=None)[0] / sc
        xe, srp = xe + dx[:6], srp + dx[6]
        if np.linalg.norm(dx[:3]) < 1e-3:
            break
    rms = np.sqrt(np.mean(np.sum(res ** 2, axis=1)))
    print(f'  radiation pressure from SSCWeb 2008-2013: {srp:.3e} km/s^2 at 1 au (estimate from size and '
          f'mass {SRP_NOMINAL:.1e}); one ballistic arc fits {len(t)} samples to {rms:.0f} km rms')
    return srp


def arc_2014(helio, srp):
    """A ballistic arc through Horizons' ISEE-3 states of 2014-01 .. 06 (solution s45)."""
    path = cached('hzn_isee3_2014.npy')
    if not os.path.exists(path):
        rows = hb.fetch('-111', '500@10', hb.jd_from_text('2014-01-01'), hb.jd_from_text('2014-06-30'), 1440)
        np.save(path, np.array([[t, *q, *v] for t, q, v in rows]))
    h = np.load(path)
    xe, rn = dyn.fit_arc(helio, h[:, 0], h[:, 1:4], h[0, 0], h[0, 1:7], srp=srp, iters=8)
    return h[0, 0], xe, rn.max()


def shoot(model, t1, x1, t2, r2, srp):
    """Impulse at t1 that takes x1 to the position r2 at t2."""
    dv = np.zeros(3)
    for _ in range(15):
        y = dyn.integrate(model, t1, np.concatenate([x1[:3], x1[3:] + dv]), t2, srp=srp, stm=True)[1][-1]
        miss = r2 - y[:3]
        dv += np.linalg.solve(y[6:48].reshape(6, 7)[0:3, 3:6], miss)
        if np.linalg.norm(miss) < 1e-2:
            break
    return dv


def report_orbit(model, t0, x0, srp):
    """Osculating heliocentric orbit in 1995 (to compare with published figures)."""
    t = dyn.utc_to_tdb(dt.datetime(1995, 1, 1))
    x = dyn.state_at(model, t0, x0, t, srp)
    r, v = x[:3], x[3:]
    mu = dyn.GM['sun']
    a = 1.0 / (2.0 / np.linalg.norm(r) - v @ v / mu)
    h = np.cross(r, v)
    e = np.linalg.norm(np.cross(v, h) / mu - r / np.linalg.norm(r))
    pole = np.array([0.0, -0.397776982902, 0.917482062069])   # ecliptic pole (J2000) in ICRF
    inc = np.degrees(np.arccos(h @ pole / np.linalg.norm(h)))
    period = 2 * np.pi * np.sqrt(a ** 3 / mu) / DAY
    print(f'  orbit in 1995: perihelion {a * (1 - e) / dyn.AU:.3f} au, aphelion {a * (1 + e) / dyn.AU:.3f} au, '
          f'inclination {inc:.2f} deg, period {period:.1f} d')


# --- Main ---------------------------------------------------------------------------------------

def main():
    refit = '--refit' in sys.argv
    os.makedirs(CACHE, exist_ok=True)
    print('Downloading (cached in tools/isee3/cache) ...')
    fetch_ssc()
    fetch_bodies()
    bodies = dyn.Bodies(CACHE)
    geo, helio = dyn.GeoModel(bodies), dyn.HelioModel(bodies)

    print('Fitting ballistic arcs to the SSCWeb track ...')
    arcs = drop_shifted(geo, segment(refit))

    # Launch: the first arc, propagated back, passes perigee an hour after the launch.
    t_launch = dyn.utc_to_tdb(dt.datetime(1978, 8, 12, 15, 12))
    ts, ys = dyn.integrate(geo, arcs[0]['epoch'], np.array(arcs[0]['state']), t_launch + 600.0)
    i = int(np.argmin(np.linalg.norm(ys[:, :3], axis=1)))
    t_start = float(ts[i])
    print(f'  launch perigee {np.linalg.norm(ys[i, :3]) - dyn.EARTH_RE:.0f} km altitude at '
          f'{dyn.tdb_to_utc(t_start):%Y-%m-%d %H:%M} UTC, {np.linalg.norm(ys[i, 3:6]):.2f} km/s')

    while True:   # an arc squeezed to under a day between its junctions adds nothing but a jump
        comp = Composite(geo, bodies, arcs, t_start, GEO_END)
        short = [k for k in range(len(arcs)) if comp.bounds[k + 1] - comp.bounds[k] < DAY]
        if not short:
            break
        print(f'  dropped {dyn.date(arcs[short[0]]["start"])} .. {dyn.date(arcs[short[0]]["end"])} '
              f'(less than a day between its junctions)')
        del arcs[short[0]]
    print(f'  {len(arcs)} arcs; joins (gap > 500 km or dv > 1 m/s):')
    for tj, hw, gap, dv in comp.joins:
        if gap > 500 or dv > 1.0:
            print(f'    {dyn.date(tj)}  gap {gap:7.0f} km  dv {dv:6.2f} m/s  blended over +-{hw / DAY:.1f} d')

    def near_geo(tt, st):
        return (np.linalg.norm(st[:, :3], axis=1) < 8.0e4) | (moon_distance(bodies, tt, st[:, :3]) < 8.0e4)

    tt = sample_times(comp.states, t_start, GEO_END, near_geo)
    knots, err = reduce(to_rows(tt, comp.states(tt)), 0.0)
    hb.write_eph(os.path.join(OUT, 'isee3_geo.eph'), ISEE3_ID, 399, knots)
    print(f'  isee3_geo.eph: {len(tt)} samples -> {len(knots)} knots, max error {err:.2f} km')

    # Heliocentric: escape, two maneuvers, the navigation arc, two more maneuvers, the arc
    # through Horizons' 2014 states.
    print('Heliocentric cruise ...')
    SRP = fit_srp(helio, bodies)
    t_nav, ice, comet_nav, err = read_nav(bodies)
    print(f'  navigation file: Earth matched to {err:.2f} km after the rotation')
    mid = len(t_nav) // 2
    xe, rn = dyn.fit_arc(helio, t_nav, ice[:, :3], t_nav[mid], ice[mid], srp=SRP)
    print(f'  ballistic fit to the navigation states: max residual {rn.max():.3f} km')
    x0 = comp.states([HELIO_START])[0] + np.concatenate(
        [bodies.helio['earth'](HELIO_START), bodies.helio['earth'].velocity(HELIO_START)])
    t_end_conn = BURNS[-1]
    target = dyn.state_at(helio, t_nav[mid], xe, t_end_conn, SRP)
    x_b0 = dyn.state_at(helio, HELIO_START, x0, BURNS[0], SRP)
    dv0 = np.zeros(3)
    for _ in range(10):   # first impulse: reach the navigation arc's position at the second
        y = dyn.integrate(helio, BURNS[0], np.concatenate([x_b0[:3], x_b0[3:] + dv0]), t_end_conn,
                          srp=SRP, stm=True)[1][-1]
        miss = target[:3] - y[:3]
        dv0 += np.linalg.solve(y[6:48].reshape(6, 7)[0:3, 3:6], miss)
        if np.linalg.norm(miss) < 1e-3:
            break
    y = dyn.state_at(helio, BURNS[0], np.concatenate([x_b0[:3], x_b0[3:] + dv0]), t_end_conn, SRP)
    dv1 = target[3:] - y[3:]
    print(f'  maneuvers: {dyn.date(BURNS[0])} {np.linalg.norm(dv0) * 1e3:.1f} m/s, '
          f'{dyn.date(BURNS[1])} {np.linalg.norm(dv1) * 1e3:.1f} m/s')

    t14, x14, err14 = arc_2014(helio, SRP)
    print(f'  arc through Horizons 2014-01 .. 06 (solution s45): max residual {err14:.0f} km')
    x_l0 = dyn.state_at(helio, t_nav[mid], xe, LATE_BURNS[0], SRP)
    late_target = dyn.state_at(helio, t14, x14, LATE_BURNS[1], SRP)
    dv2 = shoot(helio, LATE_BURNS[0], x_l0, LATE_BURNS[1], late_target[:3], SRP)
    y = dyn.state_at(helio, LATE_BURNS[0], np.concatenate([x_l0[:3], x_l0[3:] + dv2]), LATE_BURNS[1], SRP)
    dv3 = late_target[3:] - y[3:]
    print(f'  maneuvers: {dyn.date(LATE_BURNS[0])} {np.linalg.norm(dv2) * 1e3:.1f} m/s, '
          f'{dyn.date(LATE_BURNS[1])} {np.linalg.norm(dv3) * 1e3:.1f} m/s')
    pieces = [(HELIO_START, x0, BURNS[0]),
              (BURNS[0], np.concatenate([x_b0[:3], x_b0[3:] + dv0]), BURNS[1]),
              (BURNS[1], target, LATE_BURNS[0]),
              (LATE_BURNS[0], np.concatenate([x_l0[:3], x_l0[3:] + dv2]), LATE_BURNS[1]),
              (LATE_BURNS[1], late_target, HELIO_END)]
    report_orbit(helio, LATE_BURNS[1], late_target, SRP)
    knots_h = []
    for k, (ta, xa, tb) in enumerate(pieces):
        ts_, ys_ = dyn.integrate(helio, ta, xa, tb, srp=SRP)

        def fn(tq, ts_=ts_, ys_=ys_):
            return dyn.dense(helio, ts_, ys_, SRP, np.asarray(tq, float))[:, :6]

        enc = t_nav[mid]

        def near_h(tq, st):
            e = bodies.helio['earth'](tq)
            return (np.linalg.norm(st[:, :3] - e, axis=1) < 3.0e6) | (abs(tq - enc) < DAY)

        tt = sample_times(fn, ta, tb, near_h)
        kn, err = reduce(to_rows(tt, fn(tt)), dyn.GM['sun'])
        if knots_h:   # impulse: the new piece starts one second later
            t1, p1, v1 = kn[0]
            kn[0] = (t1 + 1.0, tuple(np.array(p1) + np.array(v1)), v1)
        knots_h += kn
        print(f'  piece {k}: {len(tt)} samples -> {len(kn)} knots, max error {err:.2f} km')
    hb.write_eph(os.path.join(OUT, 'isee3_helio.eph'), ISEE3_ID, 10, knots_h, dyn.GM['sun'])

    # Comet: Horizons, moved onto the navigation file's comet near the encounter.
    rows = hb.fetch('DES=21P;CAP<1986', '500@10', hb.jd_from_text('1985-01-01'), hb.jd_from_text('1986-07-01'), 60)
    tc = np.array([r[0] for r in rows])
    pc = np.array([r[1] for r in rows])
    vc = np.array([r[2] for r in rows])
    hz_at_nav = dyn.Track(np.hstack([tc[:, None], pc, vc]))(t_nav)
    offset = (comet_nav - hz_at_nav).mean(axis=0)
    print(f'  comet: navigation file vs Horizons {np.linalg.norm(comet_nav - hz_at_nav, axis=1).min():.0f}'
          f'..{np.linalg.norm(comet_nav - hz_at_nav, axis=1).max():.0f} km; mean offset applied near the encounter')
    x = np.clip((abs(tc - t_nav[mid]) - 4 * DAY) / (11 * DAY), 0.0, 1.0)
    w = (np.cos(np.pi * x) + 1) / 2     # 1 within 4 days of the encounter, 0 beyond 15
    dw = np.gradient(w, tc)
    pc = pc + w[:, None] * offset
    vc = vc + dw[:, None] * offset
    kn, err = reduce(to_rows(tc, np.hstack([pc, vc])), dyn.GM['sun'])
    hb.write_eph(os.path.join(OUT, 'gz.eph'), COMET_ID, 10, kn, dyn.GM['sun'])
    print(f'  gz.eph: {len(tc)} samples -> {len(kn)} knots, max error {err:.2f} km')


if __name__ == '__main__':
    main()
