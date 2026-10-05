"""Reconstructs the Halley Armada spacecraft for 1985-1986: Giotto, Vega 1, Vega 2, Suisei,
Sakigake, the Pioneer Venus Orbiter and Pioneer 7 (Vega 2, PVO and Pioneer 7: see their functions).

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

5. Pioneer 7 has no trajectory data after 1968 and made no midcourse maneuvers: the escape
   from the Earth on 1966-08-17 is solved so that the heliocentric orbit after it has the
   published eccentricity and inclination, and the path is integrated to 1987 (1966-1987
   planet tables in tools/armada/cache/p7). Run with the published periods it predicts the
   1986 closest approach to Halley within 4-10 hours and 1% of NASA's; the period is then
   fitted to NASA's time.

Outputs (assets/ephem, AXEPH2, heliocentric): giotto.eph, vega1.eph, vega2.eph, suisei.eph,
sakigake.eph, pioneer7.eph; pvo.eph (Venus-centered). Downloads are cached in tools/armada/cache
(not in git); the 1978-2015 planet tables are shared with tools/isee3/cache.

Usage: python tools/armada/armada_reconstruct.py [giotto vega1 vega2 suisei sakigake pvo pioneer7] [--refit]
Needs numpy. vega2 reads vega1.eph. The PVO fit takes about an hour; its arcs are cached
(tools/armada/cache/pvo_arcs.json, --refit fits them again). pioneer7 takes about a minute.
"""
import datetime as dt
import json
import os
import re
import sys
import urllib.error
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


PVO_URL = 'https://pds-ppi.igpp.ucla.edu/data/pvo-ephem/data-vso-asc/ORB_{:04d}/PVO_EPHEM_VSO_{:04d}.TAB'
VENUS_RADIUS_VSO = 6050.0   # km, the unit of the VSO columns ("Rv = 6050 km")


def load_pvo(bodies, first_orbit, last_orbit, t0, t1):
    """PVO positions relative to Venus (ICRF km) at TDB t, from the VSO files of the orbits
    (one a day; missing ones are gaps in the data), every 2 minutes at most. VSO: x toward the
    Sun, z along the pole of Venus' orbit, y = z x x (against the orbital motion)."""
    t, p = [], []
    for n in range(first_orbit, last_orbit + 1):
        path = os.path.join(CACHE, 'pvo', f'PVO_EPHEM_VSO_{n:04d}.TAB')
        if not os.path.exists(path):
            os.makedirs(os.path.dirname(path), exist_ok=True)
            try:
                urllib.request.urlretrieve(PVO_URL.format(n // 100 * 100, n), path)
            except urllib.error.HTTPError:
                open(path, 'w').close()   # no file for this orbit
        last = -1e300
        for line in open(path):
            f = line.split(',')
            if len(f) < 4:
                continue
            ti = dyn.utc_to_tdb(dt.datetime.fromisoformat(f[0][:19]))
            if ti - last < 120.0 or not t0 <= ti <= t1:
                continue
            x = [float(v) * VENUS_RADIUS_VSO for v in f[1:4]]
            if not 6000.0 < np.linalg.norm(x) < 1.0e5:
                continue   # a few corrupt values (~1e104 km) in the files
            last = ti
            t.append(ti)
            p.append(x)
    t, p = np.array(t), np.array(p)
    o = np.argsort(t)
    t, p = t[o], p[o]
    rv, vv = bodies.helio['venus'](t), bodies.helio['venus'].velocity(t)
    x = -rv / np.linalg.norm(rv, axis=1)[:, None]
    z = np.cross(rv, vv)
    z /= np.linalg.norm(z, axis=1)[:, None]
    y = np.cross(z, x)
    return t, p[:, :1] * x + p[:, 1:2] * y + p[:, 2:3] * z


# --- Arcs and tracks ----------------------------------------------------------------------

class Arc:
    def __init__(self, model, epoch, state):
        self.model, self.epoch, self.state = model, float(epoch), np.array(state, float)

    def states(self, tq):
        return dyn.propagate(self.model, self.epoch, self.state, 0.0, np.asarray(tq, float))


def fit_chunks(model, t, p, max_days, gap_days=2.0, split_above_km=None, cache=None):
    """Arcs fitted to the data in chunks of at most max_days, split at gaps. A chunk whose fit
    fails, or whose median residual exceeds split_above_km, is halved (down to ~7 hours; a
    maneuver or bad data inside it), and left to its neighbours if it still fails. Each fit
    starts from the previous arc (continuation), else from the data around the middle.
    Returns (arc, first, last) triples and prints the residuals; `cache` (a file name in the
    cache directory) keeps them between runs."""
    path = os.path.join(CACHE, cache) if cache else None
    if path and os.path.exists(path):
        rows = json.load(open(path))
        print(f'    {len(rows)} arcs from {cache} (--refit to fit again)')
        return [(Arc(model, r[0], r[1]), r[2], r[3]) for r in rows]
    out = []

    def attempt(ts, ps, te, xe):
        try:
            with np.errstate(invalid='raise', divide='raise', over='raise'):
                return dyn.fit_arc(model, ts, ps, te, xe)
        except (FloatingPointError, np.linalg.LinAlgError, ValueError, RuntimeError):
            return xe, None

    def fit(ts, ps, a, b):
        margin = min(0.5 * DAY, 0.25 * (b - a))
        k = (ts >= a - margin) & (ts <= b + margin)
        if k.sum() < 6:
            return
        te = 0.5 * (a + b)
        guesses = [dyn.guess(ts[k], ps[k], te)]
        if out:
            prev = out[-1][0]
            guesses.insert(0, (te, prev.states([te])[0]))
        for te, xe in guesses:
            xe, rn = attempt(ts[k], ps[k], te, xe)
            ok = rn is not None and (split_above_km is None or np.median(rn) <= split_above_km)
            if ok:
                break
        if not ok and split_above_km is not None and b - a > 0.6 * DAY:
            fit(ts, ps, a, 0.5 * (a + b))
            fit(ts, ps, 0.5 * (a + b), b)
            return
        if not ok:
            print(f'    {dyn.date(a)} .. {dyn.date(b)}  fit failed, left to the neighbours', flush=True)
            return
        out.append((Arc(model, te, xe), float(a), float(b)))
        print(f'    {dyn.date(a)} .. {dyn.date(b)}  {k.sum():5d} points, residuals median {np.median(rn):5.0f} '
              f'max {rn.max():6.0f} km', flush=True)

    for s in np.split(np.arange(len(t)), np.where(np.diff(t) > gap_days * DAY)[0] + 1):
        ts, ps = t[s], p[s]
        n = max(1, int(np.ceil((ts[-1] - ts[0]) / (max_days * DAY))))
        edges = np.linspace(ts[0], ts[-1], n + 1)
        for a, b in zip(edges, edges[1:]):
            fit(ts, ps, a, b)
    if path:
        json.dump([[arc.epoch, arc.state.tolist(), a, b] for arc, a, b in out], open(path, 'w'))
    return out


class SkipModel(dyn.HelioModel):
    """Heliocentric model without the gravity of some bodies: for patched-conic legs that start
    or end at a planet's center."""

    def __init__(self, b, skip):
        super().__init__(b)
        self.skip = set(skip)

    def accel(self, t, r, srp, want_grad):
        rn = np.linalg.norm(r)
        a = -self.center_gm * r / rn ** 3
        G = self.center_gm * (3 * np.outer(r, r) / rn ** 5 - np.eye(3) / rn ** 3) if want_grad else None
        e = self.b.helio['earth'](t)
        bodies = [] if 'earth' in self.skip else [(dyn.GM['earth'], e), (dyn.GM['moon'], e + self.b.moon_geo(t))]
        bodies += [(dyn.GM[n], self.b.helio[n](t)) for n in ('mercury', 'venus', 'mars', 'jupiter', 'saturn',
                                                              'uranus', 'neptune') if n not in self.skip]
        a, G = dyn._third_body(a, G, r, bodies, want_grad)
        dsrp = r / rn * (dyn.AU / rn) ** 2
        return a + srp * dsrp, G, dsrp


class VenusModel:
    """Venus-centered: Venus, the Sun, the Earth and Jupiter (an orbiter of Venus)."""

    def __init__(self, b):
        self.b = b
        self.center_gm = dyn.GM['venus']

    def accel(self, t, r, srp, want_grad):
        rn = np.linalg.norm(r)
        a = -self.center_gm * r / rn ** 3
        G = self.center_gm * (3 * np.outer(r, r) / rn ** 5 - np.eye(3) / rn ** 3) if want_grad else None
        rv = self.b.helio['venus'](t)
        bodies = [(dyn.GM['sun'], -rv), (dyn.GM['earth'], self.b.helio['earth'](t) - rv),
                  (dyn.GM['jupiter'], self.b.helio['jupiter'](t) - rv)]
        a, G = dyn._third_body(a, G, r, bodies, want_grad)
        u = r + rv
        un = np.linalg.norm(u)
        dsrp = u / un * (dyn.AU / un) ** 2
        return a + srp * dsrp, G, dsrp


def shoot(model, t1, x1, t2, r2, dv0=None):
    """Impulse at t1 (Newton from dv0) that takes x1 to the position r2 at t2."""
    dv = np.zeros(3) if dv0 is None else np.array(dv0, float)
    for _ in range(30):
        y = dyn.integrate(model, t1, np.concatenate([x1[:3], x1[3:] + dv]), t2, stm=True)[1][-1]
        miss = r2 - y[:3]
        if np.linalg.norm(miss) < 1e-2:
            return dv
        dv += np.linalg.solve(y[6:48].reshape(6, 7)[0:3, 3:6], miss)
    raise RuntimeError(f'shooting did not converge (miss {np.linalg.norm(miss):.3g} km)')


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
    """A baked table (by default the scene's Halley, assets/ephem/halley.eph): positions, and
    velocities by differencing."""

    def __init__(self, name='halley.eph'):
        self.table = hb.read_eph(os.path.join(OUT, name))

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

def bake(name, naif_id, fn, t0, t1, fine, center=10, gm=dyn.GM['sun'], step=1800.0, fine_step=60.0):
    """Samples every `step` s, every `fine_step` s between two where fine(t, states) says so;
    knots to TOL_KM (Kepler-relative about the center)."""
    tt = np.append(np.arange(t0, t1, step), t1)
    f = fine(tt, fn(tt))
    extra = [np.arange(a + fine_step, b, fine_step) for a, b, fa, fb in zip(tt, tt[1:], f, f[1:]) if fa or fb]
    tt = np.unique(np.concatenate([tt] + extra))
    st = fn(tt)
    rows = [(float(t), tuple(map(float, s[:3])), tuple(map(float, s[3:]))) for t, s in zip(tt, st)]
    interp = hb.make_interp(gm)
    knots = hb.decimate(rows, TOL_KM, interp)
    err = hb.max_error(rows, knots, interp)
    hb.write_eph(os.path.join(OUT, name), naif_id, center, knots, gm)
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
    bake('giotto.eph', -78, fn, t[0], t[-1], lambda x, _: abs(x - t_guess) < 2 * 3600.0)
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

    def fine(x, _):
        return (abs(x - t_guess) < 2 * 3600.0) | (abs(x - t_v) < 12 * 3600.0) | (x < t_start + 12 * 3600.0)
    bake('vega1.eph', -66, fn, t_start, t_end, fine)
    report_ca('Vega 1', fn, halley, t_guess)


def miss_axes(v_rel, to_sun):
    """Unit vectors across the relative velocity: s toward the Sun, w = v x s."""
    v = v_rel / np.linalg.norm(v_rel)
    s = to_sun - (to_sun @ v) * v
    s /= np.linalg.norm(s)
    return s, np.cross(v, s)


def vega_reconstruct(model, bodies, halley, guide, t_launch, t_entry, t_ca, miss_km, miss_angle, t_end):
    """A Vega trajectory from its dates alone: launch, the lander's entry at Venus, the closest
    approach to Halley (time, distance, angle of the miss vector from the Sun's direction
    around the relative velocity). `guide` (a Halley-like table of Vega 1) gives the starting
    guesses. Returns the track and a report."""
    bare = SkipModel(bodies, {'earth', 'venus'})   # patched conics: legs from/to the centers
    no_venus = SkipModel(bodies, {'venus'})
    earth, venus = bodies.helio['earth'], bodies.helio['venus']

    def vs(track, t):
        return np.concatenate([track(np.array([t]))[0], track.velocity(np.array([t]))[0]])

    # Departure 80 minutes after the launch, as Vega 1's first arc shows; the guess for the
    # escape velocity is Vega 1's (2 days out).
    t_dep = t_launch + 80 * 60.0
    g_t = guide.table.start() + 2 * DAY
    v_inf_guess = guide.states([g_t])[0, 3:] - earth.velocity(np.array([g_t]))[0]
    x_dep = vs(earth, t_dep)
    x_dep[3:] += v_inf_guess
    x_dep[3:] += shoot(bare, t_dep, x_dep, t_entry, venus(np.array([t_entry]))[0])
    leg1 = Arc(bare, t_dep, x_dep)
    v_inf_in = leg1.states([t_entry])[0, 3:] - venus.velocity(np.array([t_entry]))[0]

    # Venus -> Halley (patched), aimed at the encounter point; the miss vector depends on the
    # relative velocity, so twice.
    t_v = t_entry
    xv = vs(venus, t_v)
    g_v = guide.states([t_v + 4 * DAY])[0, 3:] - venus.velocity(np.array([t_v + 4 * DAY]))[0]
    h = halley.states([t_ca])[0]
    target = h[:3]
    for _ in range(3):
        x3 = xv.copy()
        x3[3:] += g_v
        x3[3:] += shoot(no_venus, t_v, x3, t_ca, target)
        leg3 = Arc(no_venus, t_v, x3)
        v_rel = leg3.states([t_ca])[0, 3:] - h[3:]
        s, w = miss_axes(v_rel, -h[:3])
        target = h[:3] + miss_km * (np.cos(miss_angle) * s + np.sin(miss_angle) * w)
    v_inf_out = x3[3:] - xv[3:]

    # The flyby that turns v_inf_in into v_inf_out (at the outgoing speed).
    turn = np.arccos(np.clip(v_inf_in @ v_inf_out / np.linalg.norm(v_inf_in) / np.linalg.norm(v_inf_out), -1, 1))
    mu = dyn.GM['venus']
    r_p = mu / (v_inf_out @ v_inf_out) * (1.0 / np.sin(turn / 2) - 1.0)
    r_p_used = max(r_p, 6051.8 + 1000.0)
    peri_dir = v_inf_in / np.linalg.norm(v_inf_in) - v_inf_out / np.linalg.norm(v_inf_out)
    peri_dir /= np.linalg.norm(peri_dir)

    # The bus: leg 1 until it released the lander (2 days before the entry); an impulse onto a
    # flyby that leaves along the patched outbound leg (aimed at its position two weeks later;
    # the first guess passes the periapsis above), and there another impulse onto the
    # encounter point.
    t_d, t_m = t_entry - 2 * DAY, t_v + 14 * DAY
    xd = leg1.states([t_d])[0]
    aim = venus(np.array([t_v]))[0] + r_p_used * peri_dir
    dv1 = shoot(model, t_d, xd, t_v, aim, r_p_used * peri_dir / (t_v - t_d))
    dv1 = shoot(model, t_d, xd, t_m, leg3.states([t_m])[0, :3], dv1)
    flyby = Arc(model, t_d, np.concatenate([xd[:3], xd[3:] + dv1]))
    xm = flyby.states([t_m])[0]
    dv2 = shoot(model, t_m, xm, t_ca, target)
    cruise = Arc(model, t_m, np.concatenate([xm[:3], xm[3:] + dv2]))

    track = Track(t_dep, t_end)
    track.add(leg1)
    track.add(flyby, t_d)
    track.add(cruise, t_m)
    track.finish()
    tt = np.arange(t_v - DAY, t_v + DAY, 10.0)
    dvenus = np.linalg.norm(flyby.states(tt)[:, :3] - venus(tt), axis=1)
    report = dict(v_inf_in=np.linalg.norm(v_inf_in), v_inf_out=np.linalg.norm(v_inf_out),
                  turn=np.degrees(turn), r_p=r_p, flyby_km=dvenus.min(), t_flyby=tt[int(np.argmin(dvenus))],
                  dv1=np.linalg.norm(dv1), dv2=np.linalg.norm(dv2), t_dep=t_dep, t_d=t_d, t_m=t_m,
                  v_rel=np.linalg.norm(cruise.states([t_ca])[0, 3:] - h[3:]))
    return track, report


def print_vega_report(name, r):
    print(f'  {name}: v_inf {r["v_inf_in"]:.2f} km/s at Venus in, {r["v_inf_out"]:.2f} out, turned '
          f'{r["turn"]:.1f} deg -> periapsis {r["r_p"]:,.0f} km from the center; flyby '
          f'{dyn.tdb_to_utc(r["t_flyby"]):%Y-%m-%d %H:%M} UTC at {r["flyby_km"]:,.0f} km; impulses '
          f'{dyn.date(r["t_d"])} {r["dv1"] * 1e3:.0f} m/s, {dyn.date(r["t_m"])} {r["dv2"] * 1e3:.0f} m/s; '
          f'{r["v_rel"]:.2f} km/s at Halley', flush=True)


def vega2(model, bodies, halley):
    print('Vega 2 (reconstructed from its dates; Vega 1 as the guide):')
    guide = Halley('vega1.eph')
    # The angle of Vega 1's miss vector from the Sun's direction, around the relative velocity.
    t1_ca = utc(1986, 3, 6, 7, 20, 6)
    b1 = guide.states([t1_ca])[0] - halley.states([t1_ca])[0]
    s, w = miss_axes(b1[3:], -halley([t1_ca])[0])
    angle = np.arctan2(b1[:3] @ w, b1[:3] @ s)
    print(f'  Vega 1 miss vector: {np.degrees(angle):.1f} deg from the Sun direction')

    # The method on Vega 1 (launch 1984-12-15 09:16:24 UTC, Siddiqi 2018 via Wikipedia; entry
    # 1985-06-11 02:06:10 UTC, Wikipedia), against its track from the data.
    tr1, r1 = vega_reconstruct(model, bodies, halley, guide, utc(1984, 12, 15, 9, 16, 24), utc(1985, 6, 11, 2, 6, 10),
                               t1_ca, 8890.0, angle, utc(1986, 4, 30))
    print_vega_report('method on Vega 1', r1)
    for d in (dt.datetime(1985, 1, 15), dt.datetime(1985, 3, 15), dt.datetime(1985, 5, 15), dt.datetime(1985, 6, 10),
              dt.datetime(1985, 8, 1), dt.datetime(1985, 10, 1), dt.datetime(1985, 12, 1), dt.datetime(1986, 2, 1)):
        t = dyn.utc_to_tdb(d)
        print(f'    {d:%Y-%m-%d}: {np.linalg.norm(tr1.states([t])[0, :3] - guide([t])[0]):,.0f} km from the data')

    # Vega 2: launch 1984-12-21 09:13:52 UTC (Siddiqi 2018 via Wikipedia); the balloon entered
    # at 02:06:04 UT on 1985-06-15 (NASA via Wikipedia); closest approach 1986-03-09 07:20:00
    # UT, 8,030 km, 76.8 km/s (PDS Vega mission catalog); the mission ended 1986-04 (ditto).
    t_ca = utc(1986, 3, 9, 7, 20)
    t_end = utc(1986, 4, 30)
    track, r2 = vega_reconstruct(model, bodies, halley, guide, utc(1984, 12, 21, 9, 13, 52),
                                 utc(1985, 6, 15, 2, 6, 4), t_ca, 8030.0, angle, t_end)
    print_vega_report('Vega 2', r2)
    # The table starts 20,000 km out (the escape from the parking orbit is not modeled).
    t_start = surface_exit(track.states, bodies, r2['t_dep'], r2['t_dep'] + 6 * 3600.0, 20000.0 - EARTH_RE)

    def fine(x, _):
        return (abs(x - t_ca) < 2 * 3600.0) | (abs(x - r2['t_flyby']) < 12 * 3600.0)
    fn = bake('vega2.eph', -67, track.states, t_start, t_end, fine)
    report_ca('Vega 2', fn, halley, t_ca)


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
    bake('suisei.eph', -47, fn, t_start, t[-1], lambda x, _: abs(x - t_guess) < 2 * 3600.0)
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
    fn = bake('sakigake.eph', -46, track.states, t[0], t[-1], lambda x, _: np.zeros(len(x), bool))
    report_ca('Sakigake', fn, halley, utc(1986, 3, 11, 4), 3 * DAY)


def pvo(model, bodies, halley):
    """The Pioneer Venus Orbiter around Venus, 1985-1986, from its VSO positions."""
    print('Pioneer Venus Orbiter (PDS PPI VSO positions, 1985-1986):')
    t0, t1 = utc(1985, 1, 1), utc(1987, 1, 1)
    t, p = load_pvo(bodies, 2215, 2955, t0, t1)
    print(f'  {len(t)} positions, {dyn.date(t[0])} .. {dyn.date(t[-1])}')
    venus_model = VenusModel(bodies)
    chunks = fit_chunks(venus_model, t, p, 2.0, gap_days=30.0, split_above_km=30.0, cache='pvo_arcs.json')
    track = Track(t[0], t[-1])
    track.add_fitted(chunks, report=False)
    track.finish()
    gaps = [(a, b) for a, b in zip(t, t[1:]) if b - a > 2 * DAY]
    print('  gaps over 2 days (bridged by the arcs): '
          + ', '.join(f'{dyn.date(a)} .. {dyn.date(b)}' for a, b in gaps))
    bake('pvo.eph', -12, track.states, t[0], t[-1],
         lambda x, st: np.linalg.norm(st[:, :3], axis=1) < 15000.0,
         center=299, gm=dyn.GM['venus'], step=600.0, fine_step=60.0)


# --- Pioneer 7 (1966 .. 1986, no trajectory data) -----------------------------------------------

P7_CACHE = os.path.join(CACHE, 'p7')
P7_LAUNCH = (1966, 8, 17, 15, 20, 17)   # UT (NASA, science.nasa.gov/mission/pioneer-7)
P7_Q = (1.01, 1.125)                    # perihelion, aphelion (au): NASA; NSSDC prints 1.0100, 1.1250
P7_PERIOD = (402.9, 402.95)             # days: NSSDC 1975 catalog, NASA
P7_INC = 0.09767                        # deg to the ecliptic (NSSDC 1975 catalog)
P7_CA = (dt.datetime(1986, 3, 20, 23, 36), 12.1e6)   # NASA: "23:36 UT March 20, 1986 ... 12.1 million km"
ECL_POLE = np.array([0.0, -0.397776982902, 0.917482062069])   # J2000 ecliptic pole in ICRF


def fetch_bodies_1966():
    """Planet tables for 1966-07 .. 1987-01 (tools/isee3/cache only has 1978 .. 2015), same
    names and steps as tools/isee3."""
    jd0, jd1 = hb.jd_from_text('1966-07-01'), hb.jd_from_text('1987-01-02')
    targets = [('moon_geo', '301', '500@399', 720), ('sun_geo', '10', '500@399', 1440)]
    targets += [(n + '_hel', c, '500@10', 1440 if n in ('mercury', 'venus', 'earth') else 5760)
                for n, c in [('mercury', '199'), ('venus', '299'), ('earth', '399'), ('mars', '4'),
                             ('jupiter', '5'), ('saturn', '6'), ('uranus', '7'), ('neptune', '8')]]
    os.makedirs(P7_CACHE, exist_ok=True)
    for name, command, center, step in targets:
        path = os.path.join(P7_CACHE, f'hzn_{name}.npy')
        if not os.path.exists(path):
            rows = hb.fetch(command, center, jd0, jd1, step)
            np.save(path, np.array([[t, *p, *v] for t, p, v in rows]))
            print(f'  Horizons {name} 1966-1987: {len(rows)} rows', flush=True)
    return dyn.Bodies(P7_CACHE)


def helio_elements(x, t_dep_dir):
    """Osculating heliocentric period (d), eccentricity and inclination to the J2000 ecliptic
    (deg; negative if the spacecraft left the Earth going south, i.e. near the descending node)."""
    mu = dyn.GM['sun']
    r, v = x[:3], x[3:6]
    a = 1.0 / (2.0 / np.linalg.norm(r) - v @ v / mu)
    h = np.cross(r, v)
    e = np.linalg.norm(np.cross(v, h) / mu - r / np.linalg.norm(r))
    inc = np.degrees(np.arccos(h @ ECL_POLE / np.linalg.norm(h)))
    node = np.cross(ECL_POLE, h)
    return np.array([2 * np.pi * np.sqrt(a ** 3 / mu) / DAY, e, inc * np.sign(node @ t_dep_dir)])


def escape_state(bodies, t_inj, vinf, r_p=EARTH_RE + 200.0):
    """Heliocentric state at perigee (t_inj) of the escape hyperbola with excess velocity vinf
    (ICRF km/s): the plane holding the asymptote that is closest to the equator (an eastward
    launch), perigee r_p."""
    mu = dyn.GM['earth']
    vi = np.linalg.norm(vinf)
    s = vinf / vi
    n = np.array([0.0, 0.0, 1.0]) - s[2] * s     # normal: the Earth's pole, across the asymptote
    n /= np.linalg.norm(n)
    w = np.cross(n, s)
    e_h = 1.0 + r_p * vi * vi / mu
    nu = np.arccos(-1.0 / e_h)                   # perigee to the outgoing asymptote
    p_hat = np.cos(nu) * s - np.sin(nu) * w
    q_hat = np.cross(n, p_hat)
    geo = np.concatenate([r_p * p_hat, np.sqrt(vi * vi + 2 * mu / r_p) * q_hat])
    t = np.array([t_inj])
    return geo + np.concatenate([bodies.helio['earth'](t)[0], bodies.helio['earth'].velocity(t)[0]])


def p7_departure(model, bodies, t_inj, target, vinf0, t_osc=150 * DAY):
    """Excess velocity such that the heliocentric osculating period, eccentricity and signed
    inclination 150 days after the injection (the Earth then ~2e7 km away) are `target`.
    Newton with a finite-difference Jacobian."""
    vinf = np.array(vinf0, float)
    dep = bodies.helio['earth'](np.array([t_inj]))[0]   # the node line passes near the departure point
    dep /= np.linalg.norm(dep)

    def f(v):
        x = dyn.state_at(model, t_inj, escape_state(bodies, t_inj, v), t_inj + t_osc)
        return helio_elements(x, dep)
    # Converged at 1e-7 d in the period (the 1986 encounter moves by ~5e5 s per day of it).
    scale = np.array([1e-4, 1e-6, 1e-4])          # days, -, deg
    for _ in range(12):
        r = (f(vinf) - target) / scale
        if np.max(np.abs(r)) < 1e-3:
            return vinf
        J = np.column_stack([((f(vinf + dv) - target) / scale - r) / 1e-5 for dv in np.eye(3) * 1e-5])
        vinf = vinf - np.linalg.solve(J, r)
    raise RuntimeError('Pioneer 7 departure did not converge')


def p7_arc(model, bodies, t_inj, period, ecc, sign, vinf0):
    vinf = p7_departure(model, bodies, t_inj, np.array([period, ecc, sign * P7_INC]), vinf0)
    return vinf, Arc(model, t_inj, escape_state(bodies, t_inj, vinf))


def p7_conic_vinf(bodies, t_inj, period, ecc):
    """First guess: the heliocentric conic of the given period and eccentricity through the
    Earth's position, leaving outward, in the ecliptic (patched conics)."""
    mu = dyn.GM['sun']
    t = np.array([t_inj])
    re, ve = bodies.helio['earth'](t)[0], bodies.helio['earth'].velocity(t)[0]
    a = (mu * (period * DAY / (2 * np.pi)) ** 2) ** (1 / 3)
    r = np.linalg.norm(re)
    p = a * (1 - ecc * ecc)
    nu = np.arccos((p / r - 1) / ecc)
    h = np.sqrt(mu * p)
    u = re / r
    w = np.cross(ECL_POLE, u)
    w /= np.linalg.norm(w)
    return mu / h * ecc * np.sin(nu) * u + h / r * w - ve


def pioneer7(model, bodies, halley):
    """Pioneer 7 from its launch and published orbit, with no maneuvers ("the Pioneers had no
    midcourse maneuvers"), integrated to 1986. The free choices are fitted to NASA's Halley
    encounter: see the notes in isee3.toml."""
    print('Pioneer 7 (launch 1966-08-17 and published orbit, integrated to 1986):')
    b66 = fetch_bodies_1966()
    m66 = dyn.HelioModel(b66)
    t_inj = utc(*P7_LAUNCH) + 30 * 60.0     # third-stage burnout: assumed 30 min after the launch
    t_ca_pub = dyn.utc_to_tdb(P7_CA[0])
    q, Q = P7_Q
    e_pub = (Q - q) / (Q + q)

    def ca(arc):
        t0 = utc(1986, 2, 1)
        x0 = arc.states([t0])[0]
        near = Arc(m66, t0, x0)
        tt = np.arange(utc(1986, 3, 1), utc(1986, 4, 20), 3 * 3600.0)
        d = np.linalg.norm(near.states(tt)[:, :3] - halley(tt), axis=1)
        tg = tt[int(np.argmin(d))]

        def rel(tq):
            return near.states(tq) - halley.states(tq)
        t, b, v = encounter(rel, tg, 2 * DAY)
        return t, np.linalg.norm(b), b, near

    # The published elements as they are: a prediction of the 1986 encounter.
    vinfs = {}
    for sign in (+1, -1):
        for P in P7_PERIOD:
            vinf, arc = p7_arc(m66, b66, t_inj, P, e_pub, sign, p7_conic_vinf(b66, t_inj, P, e_pub))
            t, d, b, _ = ca(arc)
            vinfs[(sign, P)] = vinf
            print(f'  published elements, P {P} d, {"north" if sign > 0 else "south"}bound: v_inf '
                  f'{np.linalg.norm(vinf):.3f} km/s; closest to Halley {dyn.tdb_to_utc(t):%Y-%m-%d %H:%M} UTC, '
                  f'{d / 1e6:.2f} million km', flush=True)

    # The fit: the published shape (e), northbound (nearer the published distance), and the
    # period such that the closest approach falls at NASA's time. The two sources differ by
    # 0.05 d, and an osculating period this soon after the escape still depends on when it is
    # taken by a few hundredths of a day (the Earth's potential, ~2e7 km away). The distance is
    # left as it comes.
    sign = +1
    vinf = vinfs[(sign, P7_PERIOD[1])]
    ps, ts = [], []
    P = P7_PERIOD[1]
    for it in range(8):
        vinf, arc = p7_arc(m66, b66, t_inj, P, e_pub, sign, vinf)
        t, d, b, near = ca(arc)
        ps.append(P)
        ts.append(t)
        print(f'    P {P:.4f} d: {dyn.tdb_to_utc(t):%Y-%m-%d %H:%M:%S} UTC, {d / 1e6:.3f} million km', flush=True)
        if abs(t - t_ca_pub) < 10.0:
            break
        P = P + 0.03 if len(ps) == 1 else ps[-1] - (ts[-1] - t_ca_pub) * (ps[-1] - ps[-2]) / (ts[-1] - ts[-2])
    sun = -halley([t])[0]
    side = np.degrees(np.arccos(b @ sun / np.linalg.norm(b) / np.linalg.norm(sun)))
    print(f'  fitted: period {P:.3f} d (published {P7_PERIOD[0]}, {P7_PERIOD[1]}); v_inf {np.linalg.norm(vinf):.3f} km/s; '
          f'closest to Halley {dyn.tdb_to_utc(t):%Y-%m-%d %H:%M} UTC, {d / 1e6:.3f} million km (NASA '
          f'{P7_CA[1] / 1e6} million km), {side:.0f} deg from the Sun direction')

    # Along the way: the Earth in 1977, the osculating orbit in 1985.
    tt = np.arange(utc(1967, 1, 1), utc(1985, 1, 1), DAY)
    st = arc.states(tt)
    de = np.linalg.norm(st[:, :3] - b66.helio['earth'](tt), axis=1)
    k = int(np.argmin(de))
    x85 = arc.states([utc(1985, 1, 1)])[0]
    el = helio_elements(x85, x85[:3] / np.linalg.norm(x85[:3]))
    print(f'  closest to the Earth after 1966: {dyn.tdb_to_utc(tt[k]):%Y-%m-%d}, {de[k] / 1e6:.2f} million km; '
          f'osculating orbit 1985-01-01: period {el[0]:.2f} d, e {el[1]:.5f}, i {abs(el[2]):.4f} deg')

    t0, t1 = utc(1985, 1, 1), utc(1987, 1, 1)
    out = Arc(model, t0, x85)      # the 1978-2015 planet tables from here on
    fn = bake('pioneer7.eph', -7, out.states, t0, t1, lambda x, _: np.zeros(len(x), bool))
    report_ca('Pioneer 7', fn, halley, t_ca_pub, 2 * DAY)


def p7_semi_major_axis(period):
    return (dyn.GM['sun'] * (period * DAY / (2 * np.pi)) ** 2) ** (1 / 3)


def main():
    os.makedirs(CACHE, exist_ok=True)
    print('Planet tables (tools/isee3/cache) ...')
    isee3.fetch_bodies()
    bodies = dyn.Bodies(isee3.CACHE)
    model = dyn.HelioModel(bodies)
    halley = Halley()
    if '--refit' in sys.argv and os.path.exists(os.path.join(CACHE, 'pvo_arcs.json')):
        os.remove(os.path.join(CACHE, 'pvo_arcs.json'))
    which = [a for a in sys.argv[1:] if a != '--refit'] or ['giotto', 'vega1', 'vega2', 'suisei', 'sakigake', 'pvo',
                                                            'pioneer7']
    for name in which:
        {'giotto': giotto, 'vega1': vega1, 'vega2': vega2, 'suisei': suisei, 'sakigake': sakigake,
         'pvo': pvo, 'pioneer7': pioneer7}[name](model, bodies, halley)


if __name__ == '__main__':
    main()
