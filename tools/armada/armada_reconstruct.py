"""Reconstructs four of the Halley Armada spacecraft for 1985-1986: Giotto, Vega 1, Suisei, Sakigake.

JPL Horizons has none of them. Their positions come from the archives of the missions and of
NASA; each source is in its own frame and time scale, and only the Suisei file has a comet
next to it that matches it well enough for the encounter:

  Giotto   ESOC position file for the encounter (PDS SBN, GIO-C-DID-3-RDR-HALLEY-V1.0,
           GEOMETRY/EPGIOTTO.DOC, appendix A1; T. Morley): 1986-03-05 .. 17, 27 a day,
           solar-system-barycentric Giotto, Halley, Earth and Sun, mean equator and equinox of
           B1950, ET (taken as TDB). Rotated with the B1950 -> J2000 matrix, then by a small
           rotation (~0.5": the FK4 equinox) that puts the file's Earth on Horizons' Earth.
  Vega 1   Positions in the MISCHA magnetometer tables (PDS SBN, VEGA1-SW-MISCHA-3-RDR-CRUISE-
           V1.0, DATA/V1CT1SE.TAB; VEGA1-C-MISCHA-3-RDR-HALLEY-V1.0, DATA/V1FBT1SE.TAB,
           V1FBT2SE.TAB): columns XH/YH/ZH_SE 1950, heliocentric, mean ecliptic and equinox of
           B1950, UT; 1984-12-22 .. 1986-06-28 (fill values after it), minutes apart, with gaps
           (one across the Venus flyby). They scatter by ~10^3 km about a smooth orbit.
  Suisei   ISAS ephemeris in the ESP data set (PDS SBN, SUISEI-C-ESP-3-RDR-HALLEY-V1.0,
           GEOMETRY/EPSUISEI.TAB, EPHALLEY.TAB): daily at 0h UT, 1985-09-28 .. 1986-12-31,
           heliocentric km. The frame is not stated: in the mean ecliptic and equinox of date
           the comet table matches Horizons' Halley to 5,000-21,000 km, in the B1950 or J2000
           ecliptic to 10^5-10^6 km.
  Sakigake NASA SPDF HelioWeb (CDAWeb SAKIGAKE_HELIO1HR_POSITION): heliocentric distance,
           solar-ecliptic latitude and longitude, mean of date (HelioWeb documentation), to
           1e-4 au and 0.001 deg. The hourly values are interpolated linearly from daily ones;
           only the 0h samples (the daily values) are used.

Method:
1. Ballistic arcs (the Sun, the planets and the Moon as in tools/isee3; no radiation pressure)
   fitted to chunks of at most 15-30 days of data, split at the data gaps. Consecutive arcs
   are joined where they come closest and blended over hours to days: the joins are at the
   level of the data's own scatter, except where noted below.
2. Vega 1 at Venus (no data 1985-06-02 .. 06-26): two days before the arrival the bus
   released the lander (Wikipedia, Vega 1) and turned off its entry path. The bridge: an
   impulse on 1985-06-09 00:00 and one on 1985-06-25 ("two to four weeks after the Venus
   flyby", ESA SP-1066, the Vega chapter), sized so that the path joins the arcs on both
   sides. Both dates are assumptions; the resulting flyby distance is printed (published:
   39,000 km, Wikipedia).
3. Launch: Vega 1's first arc, propagated back, meets the Earth on the launch day (1984-12-15,
   PDS Vega 1 mission catalog); the table starts where it is 200 km above the surface.
   Suisei's first arc, propagated back, passes the Earth at ~60,000 km on the launch day
   (small early trims are not modeled); the table starts there. Sakigake's first arc misses
   the Earth by 10^6 km (a mid-course correction): its table starts with the data, two weeks
   after the launch.
4. Encounters: within days of the closest approach the spacecraft is placed on the scene's
   Halley (assets/ephem/halley.eph, from Horizons) plus its offset from the comet of its own
   source (ESOC's for Giotto, ISAS's for Suisei; the Vega files have no comet, so for Vega 1
   it is Horizons'), shifted in time and stretched along the miss vector so that the closest
   approach matches the published time and distance where one is given. The direction of the
   miss vector is the source's.

Outputs (assets/ephem, AXEPH2, heliocentric): giotto.eph, vega1.eph, suisei.eph, sakigake.eph.
Downloads are cached in tools/armada/cache (not in git); the planet tables are shared with
tools/isee3/cache.

Usage: python tools/armada/armada_reconstruct.py
Needs numpy.
"""
import datetime as dt
import json
import os
import re
import sys
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'bake'))
sys.path.insert(0, os.path.join(HERE, '..', 'isee3'))
import horizons_bake as hb  # noqa: E402
import isee3_dynamics as dyn  # noqa: E402
import isee3_reconstruct as isee3  # noqa: E402

CACHE = os.path.join(HERE, 'cache')
OUT = os.path.normpath(os.path.join(HERE, '..', '..', 'assets', 'ephem'))
DAY = dyn.DAY
TOL_KM = 1.0
EARTH_RE = dyn.EARTH_RE
SBN = 'https://pdssbn.astro.umd.edu/holdings/'
CDAWEB = ('https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/SAKIGAKE_HELIO1HR_POSITION/data/'
          '{a:%Y%m%dT%H%M%SZ},{b:%Y%m%dT%H%M%SZ}/RAD_AU,SE_LAT,SE_LON?format=text')


def utc(*a):
    return dyn.utc_to_tdb(dt.datetime(*a))


# --- Frames -------------------------------------------------------------------------------

ARCSEC = np.pi / 180.0 / 3600.0


def rx(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def ry(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def precession(T):
    """IAU 1976 precession (Lieske et al. 1977): J2000 mean equator -> mean equator of date,
    T in Julian centuries of TDB from J2000 (v_date = P v_J2000)."""
    zeta = (2306.2181 * T + 0.30188 * T ** 2 + 0.017998 * T ** 3) * ARCSEC
    z = (2306.2181 * T + 1.09468 * T ** 2 + 0.018203 * T ** 3) * ARCSEC
    theta = (2004.3109 * T - 0.42665 * T ** 2 - 0.041833 * T ** 3) * ARCSEC
    return rz(z) @ ry(-theta) @ rz(zeta)


def obliquity(T):
    """IAU 1980 mean obliquity of date (rad)."""
    return (84381.448 - 46.8150 * T - 0.00059 * T ** 2 + 0.001813 * T ** 3) * ARCSEC


T_B1950 = (2433282.4235 - 2451545.0) / 36525.0   # B1950.0 in Julian centuries from J2000
# B1950 mean equator -> J2000 (the precession part of SPICE's B1950 frame, chgirf.f).
B1950_TO_J2000 = precession(T_B1950).T
ECL_B1950_TO_J2000 = B1950_TO_J2000 @ rx(obliquity(T_B1950))


def ecliptic_of_date_to_icrf(t, v):
    """Vectors v (N x 3) in the mean ecliptic and equinox of their dates t (TDB s) -> ICRF."""
    out = np.empty_like(v)
    for i, (ti, vi) in enumerate(zip(t, v)):
        T = ti / (36525.0 * DAY)
        out[i] = precession(T).T @ rx(obliquity(T)) @ vi
    return out


def fit_rotation(src, dst):
    """Rotation R (least squares) with R src ~ dst, for N x 3 arrays."""
    U, _, Vt = np.linalg.svd(src.T @ dst)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    return Vt.T @ D @ U.T


# --- Downloads ----------------------------------------------------------------------------

def fetch(url, name):
    path = os.path.join(CACHE, name)
    if not os.path.exists(path):
        os.makedirs(CACHE, exist_ok=True)
        print(f'  downloading {name}', flush=True)
        urllib.request.urlretrieve(url, path)
    return path


def fetch_sakigake():
    """HelioWeb positions 1985-01-21 .. 1987-01-01 from CDAWeb, a month per request (longer
    requests are refused), cached as one text file."""
    path = os.path.join(CACHE, 'sakigake_helio1hr.txt')
    if not os.path.exists(path):
        headers = {'Accept': 'application/json', 'User-Agent': 'astraxis-tools'}
        lines = []
        a = dt.datetime(1985, 1, 21)
        while a < dt.datetime(1987, 1, 2):
            b = min(dt.datetime(a.year + a.month // 12, a.month % 12 + 1, 1), dt.datetime(1987, 1, 2))
            req = urllib.request.Request(CDAWEB.format(a=a, b=b), headers=headers)
            name = json.load(urllib.request.urlopen(req, timeout=300))['FileDescription'][0]['Name']
            text = urllib.request.urlopen(urllib.request.Request(name, headers=headers), timeout=300).read()
            lines += [l for l in text.decode().splitlines() if re.match(r'\d\d-\d\d-\d{4} ', l)]
            print(f'  CDAWeb {a:%Y-%m}: {len(lines)} rows', flush=True)
            a = b
        os.makedirs(CACHE, exist_ok=True)
        with open(path, 'w') as f:
            f.write('\n'.join(lines) + '\n')
    return path


# --- Sources ------------------------------------------------------------------------------

def load_giotto(bodies):
    """ESOC file: Giotto, Halley (heliocentric ICRF) at TDB t."""
    text = open(fetch(SBN + 'gio-c-did-3-rdr-halley-v1.0/geometry/epgiotto.doc', 'epgiotto.doc')).read()
    lines = text.splitlines()
    start = max(k for k, l in enumerate(lines) if l.strip() == 'The Position Data File')
    rows = [l.rstrip() for l in lines[start + 1:] if l.strip()]
    nums = []
    for l in rows[1:]:   # rows[0]: tape date, record size, record count
        if not re.fullmatch(r'[ \-0-9.D+]+', l):
            break
        nums += [float(l[j:j + 16].replace('D', 'E')) for j in range(0, len(l), 16) if l[j:j + 16].strip()]
    n = int(rows[0].split()[-1])
    rec = np.array(nums[5:5 + 13 * n]).reshape(n, 13)   # ET, Giotto, Halley, Earth, Sun (km)
    t = (rec[:, 0] - 18262.5) * DAY      # MJD counted from 1950-01-01 0h ("14 March 1986 = 13221")
    sun = rec[:, 10:13]
    giotto, comet, earth = [(rec[:, k:k + 3] - sun) @ B1950_TO_J2000.T for k in (1, 4, 7)]
    ref = bodies.helio['earth'](t)
    before = np.linalg.norm(earth - ref, axis=1).max()
    R = fit_rotation(earth, ref)
    after = np.linalg.norm(earth @ R.T - ref, axis=1).max()
    angle = np.degrees(np.arccos((np.trace(R) - 1) / 2)) * 3600
    print(f'  Giotto file: Earth vs Horizons {before:.0f} km -> {after:.1f} km after a {angle:.2f}" rotation')
    return t, giotto @ R.T, comet @ R.T


def load_vega1():
    """MISCHA tables: valid positions (Q_POSITION 0), every 10 minutes at most."""
    files = [(SBN + 'vega1-sw-mischa-3-rdr-cruise-v1.0/data/v1ct1se.tab', 'v1ct1se.tab'),
             (SBN + 'vega1-c-mischa-3-rdr-halley-v1.0/data/v1fbt1se.tab', 'v1fbt1se.tab'),
             (SBN + 'vega1-c-mischa-3-rdr-halley-v1.0/data/v1fbt2se.tab', 'v1fbt2se.tab')]
    rows = {}
    for url, name in files:
        last = -1e300
        for line in open(fetch(url, name)):
            f = line.split()
            # 54 columns: 21-23 XH, YH, ZH_SE 1950 (km), 24 RH_SE; 52 Q_POSITION (0 = good)
            if len(f) != 54 or float(f[51]) != 0.0 or float(f[23]) <= 0.0:
                continue
            t = dyn.utc_to_tdb(dt.datetime.strptime(f[0][:19], '%Y-%m-%dT%H:%M:%S'))
            if t - last < 600.0:
                continue
            last = t
            rows[round(t)] = [float(v) for v in f[20:23]]
    t = np.array(sorted(rows), float)
    return t, np.array([rows[k] for k in sorted(rows)]) @ ECL_B1950_TO_J2000.T


def load_isas(name):
    t, p = [], []
    for line in open(fetch(SBN + 'suisei-c-esp-3-rdr-halley-v1.0/geometry/' + name, name)):
        f = line.split()
        t.append(dyn.utc_to_tdb(dt.datetime(1900 + int(f[0]), int(f[1]), int(f[2]))))
        p.append([float(v) for v in f[3:6]])
    t = np.array(t)
    return t, ecliptic_of_date_to_icrf(t, np.array(p))


def load_sakigake():
    t, p = [], []
    for line in open(fetch_sakigake()):
        if not re.match(r'\d\d-\d\d-\d{4} 00:00', line):
            continue
        d, _, r, lat, lon = line.split()
        t.append(dyn.utc_to_tdb(dt.datetime.strptime(d, '%d-%m-%Y')))
        r, lat, lon = float(r) * dyn.AU, np.radians(float(lat)), np.radians(float(lon))
        p.append([r * np.cos(lat) * np.cos(lon), r * np.cos(lat) * np.sin(lon), r * np.sin(lat)])
    t = np.array(t)
    return t, ecliptic_of_date_to_icrf(t, np.array(p))


# --- Arcs and tracks ----------------------------------------------------------------------

class Arc:
    def __init__(self, model, epoch, state):
        self.model, self.epoch, self.state = model, float(epoch), np.array(state, float)

    def states(self, tq):
        return dyn.propagate(self.model, self.epoch, self.state, 0.0, np.asarray(tq, float))


def fit_chunks(model, t, p, max_days, gap_days=2.0):
    """Arcs fitted to the data in chunks of at most max_days, split at gaps. Returns
    (arc, first, last) triples and prints the residuals."""
    out = []
    for s in np.split(np.arange(len(t)), np.where(np.diff(t) > gap_days * DAY)[0] + 1):
        ts, ps = t[s], p[s]
        n = max(1, int(np.ceil((ts[-1] - ts[0]) / (max_days * DAY))))
        edges = np.linspace(ts[0], ts[-1], n + 1)
        for a, b in zip(edges, edges[1:]):
            k = (ts >= a - 0.5 * DAY) & (ts <= b + 0.5 * DAY)
            if k.sum() < 6:
                continue
            te, xe = dyn.guess(ts[k], ps[k], 0.5 * (a + b))
            xe, rn = dyn.fit_arc(model, ts[k], ps[k], te, xe)
            out.append((Arc(model, te, xe), float(a), float(b)))
            print(f'    {dyn.date(a)} .. {dyn.date(b)}  {k.sum():5d} points, residuals median {np.median(rn):5.0f} '
                  f'max {rn.max():6.0f} km', flush=True)
    return out


def closest(A, B, lo, hi, step=600.0):
    tt = np.arange(lo, hi, step)
    sa, sb = A.states(tt), B.states(tt)
    d = np.linalg.norm(sa[:, :3] - sb[:, :3], axis=1)
    i = int(np.argmin(d))
    return tt[i], d[i], np.linalg.norm(sa[i, 3:] - sb[i, 3:])


class Track:
    """Arcs in sequence. A junction is either an impulse (the arcs meet there) or the point
    where two fitted arcs come closest, blended over a span that keeps the velocity error
    below ~0.02 km/s."""

    def __init__(self, t_start, t_end):
        self.t_start, self.t_end = t_start, t_end
        self.arcs, self.joins = [], []   # joins: [time, half width]

    def add(self, arc, join=None, half_width=0.0):
        if self.arcs:
            self.joins.append([join, half_width])
        self.arcs.append(arc)

    def add_fitted(self, chunks, report=True):
        """Appends fitted arcs, each joined to the previous one where they come closest."""
        for arc, a, b in chunks:
            if self.arcs:
                tj, gap, dv = closest(self.arcs[-1], arc, self.prev_end - DAY, a + DAY)
                hw = float(np.clip(gap / 0.04, 6 * 3600.0, 5 * DAY))
                if report and (gap > 3000.0 or dv > 0.01):
                    print(f'    join {dyn.date(tj)}: gap {gap:6.0f} km, dv {dv * 1e3:5.1f} m/s', flush=True)
                self.add(arc, tj, hw)
            else:
                self.add(arc)
            self.prev_end = b

    def finish(self):
        bounds = [self.t_start] + [j[0] for j in self.joins] + [self.t_end]
        for k, j in enumerate(self.joins):   # blends may not overlap
            j[1] = min(j[1], 0.45 * (bounds[k + 1] - bounds[k]), 0.45 * (bounds[k + 2] - bounds[k + 1]))
        self.bounds = bounds

    def states(self, tq):
        tq = np.asarray(tq, float)
        out = np.zeros((len(tq), 6))
        k = np.clip(np.searchsorted(self.bounds, tq, side='right') - 1, 0, len(self.arcs) - 1)
        for i, arc in enumerate(self.arcs):
            if (k == i).any():
                out[k == i] = arc.states(tq[k == i])
        for i, (tj, hw) in enumerate(self.joins):
            sel = (tq > tj - hw) & (tq < tj + hw)
            if hw <= 0.0 or not sel.any():
                continue
            sa, sb = self.arcs[i].states(tq[sel]), self.arcs[i + 1].states(tq[sel])
            s = (tq[sel] - (tj - hw)) / (2 * hw)
            w = (s * s * (3 - 2 * s))[:, None]
            dw = (6 * s * (1 - s) / (2 * hw))[:, None]
            out[sel, :3] = (1 - w) * sa[:, :3] + w * sb[:, :3]
            out[sel, 3:] = (1 - w) * sa[:, 3:] + w * sb[:, 3:] + dw * (sb[:, :3] - sa[:, :3])
        return out


# --- Encounters -------------------------------------------------------------------------

class Halley:
    """The scene's Halley (assets/ephem/halley.eph): positions, velocities by differencing."""

    def __init__(self):
        self.table = hb.read_eph(os.path.join(OUT, 'halley.eph'))

    def __call__(self, tq):
        return np.array([hb.eval_eph(self.table, t) for t in tq])

    def states(self, tq):
        tq = np.asarray(tq, float)
        p = self(tq)
        v = (self(tq + 1.0) - self(tq - 1.0)) / 2.0
        return np.hstack([p, v])


def sampled_path(t, rel):
    """Cubic Hermite through samples rel (N x 3) at times t, velocities from second-order
    differences: a relative path (spacecraft minus comet) as states."""
    track = dyn.Track(np.hstack([t[:, None], rel, np.gradient(rel, t, axis=0)]))

    def f(tq):
        tq = np.asarray(tq, float)
        return np.hstack([track(tq), track.velocity(tq)])
    return f


def encounter(rel, t_guess, span=3 * 3600.0):
    """Closest approach of the relative path rel(t) (states) within t_guess +- span: time,
    miss vector, speed (coarse-to-fine scans down to 0.01 s)."""
    t = t_guess
    while True:
        tt = np.linspace(t - span, t + span, 401)
        s = rel(tt)
        i = int(np.argmin(np.linalg.norm(s[:, :3], axis=1)))
        t = tt[i]
        if span < 1.0:
            return t, s[i, :3], np.linalg.norm(s[i, 3:])
        span = 4 * span / 400


class Anchored:
    """A track moved near an encounter onto the scene's Halley plus a relative path,
    shifted in time by dt and stretched along the miss vector by `stretch` (km)."""

    def __init__(self, track, halley, rel, t_ca, dt_shift, stretch, inner, outer):
        self.track, self.halley, self.rel = track, halley, rel
        self.t_ca, self.dt, self.stretch = t_ca, dt_shift, np.asarray(stretch, float)
        self.inner, self.outer = inner, outer

    def weight(self, tq):
        x = np.clip((abs(tq - self.t_ca) - self.inner) / (self.outer - self.inner), 0.0, 1.0)
        w = (np.cos(np.pi * x) + 1) / 2
        dw = -np.pi / 2 * np.sin(np.pi * x) * np.sign(tq - self.t_ca) / (self.outer - self.inner)
        dw[(x <= 0.0) | (x >= 1.0)] = 0.0
        return w[:, None], dw[:, None]

    def states(self, tq):
        tq = np.asarray(tq, float)
        base = self.track.states(tq)
        near = abs(tq - self.t_ca) < self.outer
        if near.any():
            t = tq[near]
            target = self.halley.states(t) + self.rel(t - self.dt)
            target[:, :3] += self.stretch
            off = target - base[near]
            w, dw = self.weight(t)
            base[near, :3] += w * off[:, :3]
            base[near, 3:] += w * off[:, 3:] + dw * off[:, :3]
        return base


def anchor(name, track, halley, rel, t_guess, published, inner, outer):
    """Places `track` near the encounter on the scene's Halley plus `rel` (the source's own
    relative path); `published` = (UTC time, km) or None keeps the source's closest approach."""
    t0, b0, v0 = encounter(rel, t_guess)
    if published is None:
        t_ca, dt_shift, stretch = t0, 0.0, np.zeros(3)
    else:
        t_ca = dyn.utc_to_tdb(published[0])
        dt_shift = t_ca - t0
        stretch = (published[1] / np.linalg.norm(b0) - 1.0) * b0
    sun = -halley([t0])[0]
    side = np.degrees(np.arccos(b0 @ sun / np.linalg.norm(b0) / np.linalg.norm(sun)))
    print(f'  {name}: source encounter {dyn.tdb_to_utc(t0):%Y-%m-%d %H:%M:%S} UTC, {np.linalg.norm(b0):.0f} km at '
          f'{v0:.2f} km/s, miss vector {side:.1f} deg from the Sun direction'
          + ('' if published is None else f'; moved by {dt_shift:+.1f} s, {published[1] - np.linalg.norm(b0):+.0f} km'))
    return Anchored(track, halley, rel, t_ca, dt_shift, stretch, inner, outer)


# --- Output -----------------------------------------------------------------------------------

def bake(name, naif_id, fn, t0, t1, fine):
    """30-minute samples, 1 minute where fine(t) says so; knots to TOL_KM (Kepler about the Sun)."""
    tt = np.append(np.arange(t0, t1, 1800.0), t1)
    extra = [np.arange(a + 60.0, b, 60.0) for a, b in zip(tt, tt[1:]) if fine(a) or fine(b)]
    tt = np.unique(np.concatenate([tt] + extra))
    st = fn(tt)
    rows = [(float(t), tuple(map(float, s[:3])), tuple(map(float, s[3:]))) for t, s in zip(tt, st)]
    interp = hb.make_interp(dyn.GM['sun'])
    knots = hb.decimate(rows, TOL_KM, interp)
    err = hb.max_error(rows, knots, interp)
    hb.write_eph(os.path.join(OUT, name), naif_id, 10, knots, dyn.GM['sun'])
    print(f'  {name}: {dyn.tdb_to_utc(t0):%Y-%m-%d %H:%M} .. {dyn.tdb_to_utc(t1):%Y-%m-%d %H:%M} UTC, '
          f'{len(rows)} samples -> {len(knots)} knots, max error {err:.2f} km', flush=True)
    return fn


def report_ca(name, fn, halley, t_guess, half=6 * 3600.0):
    def rel(tq):
        return fn(tq) - halley.states(tq)
    t, b, v = encounter(rel, t_guess, half)
    print(f'  {name} at Halley: {dyn.tdb_to_utc(t):%Y-%m-%d %H:%M:%S} UTC, {np.linalg.norm(b):,.0f} km, '
          f'{v:.2f} km/s', flush=True)


def surface_exit(fn, bodies, t_lo, t_hi, altitude):
    """First time in [t_lo, t_hi] (1 s) when the track is `altitude` km above the Earth's surface."""
    tt = np.arange(t_lo, t_hi, 1.0)
    d = np.linalg.norm(fn(tt)[:, :3] - bodies.helio['earth'](tt), axis=1)
    i = int(np.argmax(d > EARTH_RE + altitude))
    return tt[i]


# --- Spacecraft --------------------------------------------------------------------------------

def giotto(model, bodies, halley):
    print('Giotto (ESOC, 1986-03-05 .. 17):')
    t, g, h = load_giotto(bodies)
    chunks = fit_chunks(model, t, g, 15.0)
    track = Track(t[0], t[-1])
    track.add_fitted(chunks)
    track.finish()
    # Relative path: the file's Giotto minus the file's Halley, 27 samples a day.
    t_guess = utc(1986, 3, 14, 0, 3)
    rel = sampled_path(t, g - h)
    # Closest approach 1986-03-14 00:03:01.84 UTC, 596 km: Curdt et al., A&A 191, L1, from the
    # HMC images, quoted by the file's note (the file itself gives 00:02:59.5 and 629 km).
    fn = anchor('Giotto', track, halley, rel, t_guess, (dt.datetime(1986, 3, 14, 0, 3, 1, 840000), 596.0),
                1.0 * DAY, 3.0 * DAY).states
    bake('giotto.eph', -78, fn, t[0], t[-1], lambda x: abs(x - t_guess) < 2 * 3600.0)
    report_ca('Giotto', fn, halley, t_guess)


def vega1(model, bodies, halley):
    print('Vega 1 (MISCHA positions, 1984-12-22 .. 1986-06-28):')
    t, p = load_vega1()
    chunks = fit_chunks(model, t, p, 15.0)
    venus_gap = (utc(1985, 6, 2, 12), utc(1985, 6, 26))
    pre = [c for c in chunks if c[2] < venus_gap[0]]
    post = [c for c in chunks if c[1] > venus_gap[1]]

    # Launched on 1984-12-15 (PDS Vega 1 mission catalog): the first arc, propagated back to
    # that day, starts where it is 200 km above the surface.
    t_launch = utc(1984, 12, 15)
    first = pre[0][0]
    ts, ys = dyn.integrate(model, first.epoch, first.state, t_launch)
    d = np.linalg.norm(ys[:, :3] - bodies.helio['earth'](ts), axis=1)
    t_min = ts[int(np.argmin(d))]
    print(f'  first arc propagated back: {d.min():.0f} km from the Earth\'s center at '
          f'{dyn.tdb_to_utc(t_min):%Y-%m-%d %H:%M} UTC')
    t_start = surface_exit(first.states, bodies, t_min, t_min + 3 * 3600.0, 200.0)

    # Venus: bridge with two impulses (dates assumed, see the module notes).
    A, B = pre[-1][0], post[0][0]
    t1, t2 = utc(1985, 6, 9), utc(1985, 6, 25)
    x1 = dyn.state_at(model, A.epoch, A.state, t1)
    x2 = dyn.state_at(model, B.epoch, B.state, t2)
    dv1 = isee3.shoot(model, t1, x1, t2, x2[:3], 0.0)
    bridge = Arc(model, t1, np.concatenate([x1[:3], x1[3:] + dv1]))
    dv2 = x2[3:] - bridge.states([t2])[0, 3:]
    tt = np.arange(t1, t2, 10.0)
    dv_ = np.linalg.norm(bridge.states(tt)[:, :3] - bodies.helio['venus'](tt), axis=1)
    t_v = tt[int(np.argmin(dv_))]
    print(f'  Venus: impulses {dyn.date(t1)} {np.linalg.norm(dv1) * 1e3:.0f} m/s, {dyn.date(t2)} '
          f'{np.linalg.norm(dv2) * 1e3:.0f} m/s; flyby {dyn.tdb_to_utc(t_v):%Y-%m-%d %H:%M} UTC at '
          f'{dv_.min():,.0f} km from the center')

    t_end = t[-1]
    track = Track(t_start, t_end)
    track.add_fitted(pre)
    track.add(bridge, t1)
    track.add(B, t2)
    track.prev_end = post[0][2]
    track.add_fitted(post[1:])
    track.finish()

    # Encounter: no comet in the files; the path relative to the scene's Halley is taken from
    # the fitted arcs. Closest approach 1986-03-06 07:20:06 UTC, 8,890 km, 79.2 km/s (PDS SBN
    # Vega 1 mission catalog, "Key data of the Vega mission").
    t_guess = utc(1986, 3, 6, 7, 20)

    def rel(tq):
        return track.states(tq) - halley.states(tq)
    fn = anchor('Vega 1', track, halley, rel, t_guess, (dt.datetime(1986, 3, 6, 7, 20, 6), 8890.0),
                1.0 * DAY, 3.0 * DAY).states

    def fine(x):
        return abs(x - t_guess) < 2 * 3600.0 or abs(x - t_v) < 12 * 3600.0 or x < t_start + 12 * 3600.0
    bake('vega1.eph', -66, fn, t_start, t_end, fine)
    report_ca('Vega 1', fn, halley, t_guess)


def suisei(model, bodies, halley):
    print('Suisei (ISAS ephemeris, 1985-09-28 .. 1986-12-31):')
    t, p = load_isas('epsuisei.tab')
    th, ph = load_isas('ephalley.tab')
    track_halley = {round(x): y for x, y in zip(th, ph)}
    chunks = fit_chunks(model, t, p, 30.0)
    # Launch 1985-08-18 23:33 UTC (JAXA: 08:33 JST on 08-19): the first arc propagated back.
    first = chunks[0][0]
    t_launch = utc(1985, 8, 18, 23, 33)
    ts, ys = dyn.integrate(model, first.epoch, first.state, t_launch)
    d = np.linalg.norm(ys[:, :3] - bodies.helio['earth'](ts), axis=1)
    t_start = float(ts[int(np.argmin(d))])
    print(f'  first arc propagated back: closest to the Earth {dyn.tdb_to_utc(t_start):%Y-%m-%d %H:%M} UTC, '
          f'{d.min():,.0f} km (the table starts there)')
    track = Track(t_start, t[-1])
    track.add_fitted(chunks)
    track.finish()
    # Relative to ISAS's comet (daily samples). JAXA: 151,000 km "on the side facing the Sun"
    # on 1986-03-08; the files' own closest approach is kept.
    t_guess = utc(1986, 3, 8, 13)
    k = np.array([round(x) in track_halley for x in t])
    rel = sampled_path(t[k], p[k] - np.array([track_halley[round(x)] for x in t[k]]))
    fn = anchor('Suisei', track, halley, rel, t_guess, None, 2.0 * DAY, 8.0 * DAY).states
    bake('suisei.eph', -47, fn, t_start, t[-1], lambda x: abs(x - t_guess) < 2 * 3600.0)
    report_ca('Suisei', fn, halley, t_guess)


def sakigake(model, bodies, halley):
    print('Sakigake (HelioWeb, 1985-01-21 .. 1986-12-31):')
    t, p = load_sakigake()
    keep = t < utc(1987, 1, 1)
    t, p = t[keep], p[keep]
    chunks = fit_chunks(model, t, p, 30.0)
    track = Track(t[0], t[-1])
    track.add_fitted(chunks)
    track.finish()
    first = chunks[0][0]
    t_launch = utc(1985, 1, 7, 19, 26)   # JAXA: 04:26 JST on 01-08
    ts, ys = dyn.integrate(model, first.epoch, first.state, t_launch)
    d = np.linalg.norm(ys[:, :3] - bodies.helio['earth'](ts), axis=1)
    print(f'  first arc propagated back to the launch: closest to the Earth {d.min():,.0f} km '
          f'(a mid-course correction; the table starts with the data)')
    fn = bake('sakigake.eph', -46, track.states, t[0], t[-1], lambda x: False)
    report_ca('Sakigake', fn, halley, utc(1986, 3, 11, 4), 3 * DAY)


def main():
    os.makedirs(CACHE, exist_ok=True)
    print('Planet tables (tools/isee3/cache) ...')
    isee3.fetch_bodies()
    bodies = dyn.Bodies(isee3.CACHE)
    model = dyn.HelioModel(bodies)
    halley = Halley()
    which = sys.argv[1:] or ['giotto', 'vega1', 'suisei', 'sakigake']
    for name in which:
        {'giotto': giotto, 'vega1': vega1, 'suisei': suisei, 'sakigake': sakigake}[name](model, bodies, halley)


if __name__ == '__main__':
    main()
