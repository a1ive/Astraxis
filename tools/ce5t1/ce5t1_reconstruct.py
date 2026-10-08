"""Reconstructs the extended mission of the Chang'e-5 T1 service module (2014-10-31 .. 2015-01-13).

No public ephemeris exists. [L16] Liu & Li 2016, "CHANG'E-5T1 extended mission: The first lunar
libration point flight via a lunar swing-by", Adv. Space Res. 58, 1257 (doi:10.1016/j.asr.2016.05.015),
give the state after the separation from the reentry capsule (Table 1, with the dynamics used),
the 12 maneuvers flown (Table 5: times and delta-v magnitudes, not directions) and plots of the
flight (Fig. 5). The directions are solved for here, each from what the paper says the maneuver
did; the rest is the force model (Earth with J2, Moon, Sun, planets, radiation pressure with
Table 1's coefficient, area and mass; isee3_dynamics.GeoModel) and impulsive burns.

- Table 1 is the state the extended mission was designed from: the design's first burn (Table 3,
  trajectory a, 73.2 m/s) put retrograde gives apogee an hour from the design apogee maneuver,
  533,600 km away (Table 4: 540,000). The flown first burn was 26.0 m/s, so the flown orbit had
  less energy: the speed at the epoch is lowered (11 m/s) until apogee falls on the apogee
  maneuver. That orbit's next perigee then falls within an hour of the perigee maneuver by
  itself (but inside the Earth: the apogee maneuver, prograde, raises it to 650 km).
- Burn 1 (in-plane angle, besides retrograde) and burn 3 (perigee, near retrograde) aim the
  lunar swing-by; burn 4 (the 1.9 m/s correction) is left along the velocity. Burn 5 brakes
  anti-velocity at perilune ([L16] 3.2.1). To stay at EML2 the swing-by must deliver the craft
  onto the stable manifold of an L2 orbit; the arrivals that do so with small out-of-plane
  motion pass low (perilune 200-450 km). Burn 3 is bisected onto the manifold; burn 1 was chosen
  so that the craft is near Fig. 5b's point 6 at burn 6 (the figure read to ~500 km).
- Burns 6 and 7 (EML2 injection and maintenance, 2.2 and 2.7 m/s) are bisected onto the manifold
  again. The L2 orbit is balanced on a knife edge: the solution only holds for the integration
  path of Flight (one propagation call between burns); differences of 1e-9 grow to 10^4 km in a
  month.
- Burns 8 and 9: brake 1 at perilune (the paper's "perilune arrival time"), and brakes 2 and 3
  at perilunes too (Fig. 5c). The two-body orbit through the brakes of Table 5 alone fixes the
  perilune: 1940 km from the center at 2.210 km/s, so that brake 2 comes 3 revolutions after
  brake 1 and brake 3 7 after brake 2 (brake_perilune); the final orbit is then ~200 km
  circular, Eq. 2's design orbit. The design inclination (45 deg) cannot be met as well (the
  arrival energy is set by the L2 orbit): it comes out 39.5 deg.

Checks against [L16] are printed (report) and tested in tests/mission_tests.cpp.

Usage: python tools/ce5t1/ce5t1_reconstruct.py [--solve [steps]]
Needs numpy. Without --solve, bakes from tools/ce5t1/solution.json (seconds); --solve re-derives
the directions from those in the file (steps 1-4 of solve; about an hour). Downloads are cached
in tools/ce5t1/cache (not in git).
"""
import datetime as dt
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'bake'))
sys.path.insert(0, os.path.join(HERE, '..', 'isee3'))
import horizons_bake as hb  # noqa: E402
import isee3_dynamics as dyn  # noqa: E402

CACHE = os.path.join(HERE, 'cache')
OUT = os.path.normpath(os.path.join(HERE, '..', '..', 'assets', 'ephem'))
DAY = dyn.DAY
HOUR = 3600.0
MOON_R = 1737.4      # km, mean radius (JPL SSD)


def utc(*a):
    return dyn.utc_to_tdb(dt.datetime(*a))


def utc_text(t):
    return f'{dyn.tdb_to_utc(t):%Y-%m-%d %H:%M:%S}'


# [L16] Table 1: Earth-centric J2000.0, epoch 2014-10-31 22:18:08 UTC
T0 = utc(2014, 10, 31, 22, 18, 8)
X0 = np.array([1598.667, 6093.421, 1636.044, -8.465, 0.330, 7.044])
# Table 1: SRP coefficient 1.24, area 21.7 m^2, mass 1965 kg; 4.56e-6 N/m^2 at 1 au -> km/s^2
SRP = 1.24 * 4.56e-6 * 21.7 / 1965.0 * 1e-3

# [L16] Table 5: (UTC, delta-v m/s, event)
BURNS = [
    (utc(2014, 11, 1, 8, 13, 0), 26.0, 'maneuver 10 h after separation'),
    (utc(2014, 11, 9, 6, 42, 0), 17.4, 'apogee maneuver'),
    (utc(2014, 11, 17, 15, 37, 30), 26.4, 'perigee maneuver'),
    (utc(2014, 11, 21, 8, 0, 0), 1.9, 'trajectory correction'),
    (utc(2014, 11, 23, 15, 11, 41), 233.6, 'perilune maneuver'),
    (utc(2014, 11, 28, 9, 0, 20), 2.2, 'EML2 orbit injection'),
    (utc(2014, 12, 10, 16, 10, 20), 2.7, 'EML2 orbit maintenance 1'),
    (utc(2014, 12, 26, 8, 30, 20), 9.6, 'EML2 orbit maintenance 2'),
    (utc(2015, 1, 4, 15, 0, 0), 41.3, 'return from the EML2 orbit'),
    (utc(2015, 1, 10, 19, 4, 33), 219.7, 'lunar brake 1'),
    (utc(2015, 1, 11, 17, 35, 39), 171.9, 'lunar brake 2'),
    (utc(2015, 1, 12, 19, 33, 21), 231.2, 'lunar brake 3'),
]
TB = [b[0] for b in BURNS]
DV = [b[1] * 1e-3 for b in BURNS]

# EML2 in the circular restricted problem with the DE440 mass ratio (as in earth_moon.toml)
MU = 4902.800118 / (398600.435507 + 4902.800118)


def _l2_gamma(mu):
    g = (mu / 3) ** (1 / 3)
    for _ in range(50):
        f = g ** 5 + (3 - mu) * g ** 4 + (3 - 2 * mu) * g ** 3 - mu * g ** 2 - 2 * mu * g - mu
        df = 5 * g ** 4 + 4 * (3 - mu) * g ** 3 + 3 * (3 - 2 * mu) * g ** 2 - 2 * mu * g - 2 * mu
        g -= f / df
    return g


GAMMA2 = _l2_gamma(MU)   # L2 beyond the Moon, as a fraction of the Earth-Moon distance

SPAN = (hb.jd_from_text('2014-10-30'), hb.jd_from_text('2015-02-01'))


def fetch_bodies():
    """The perturbing bodies of isee3_dynamics.GeoModel over the mission (cached)."""
    os.makedirs(CACHE, exist_ok=True)
    targets = [('moon_geo', '301', '500@399', 30), ('sun_geo', '10', '500@399', 720)]
    targets += [(n + '_hel', c, '500@10', 1440) for n, c in
                [('mercury', '199'), ('venus', '299'), ('earth', '399'), ('mars', '4'),
                 ('jupiter', '5'), ('saturn', '6'), ('uranus', '7'), ('neptune', '8')]]
    for name, command, center, step in targets:
        path = os.path.join(CACHE, f'hzn_{name}.npy')
        if not os.path.exists(path):
            rows = hb.fetch(command, center, SPAN[0], SPAN[1], step)
            np.save(path, np.array([[t, *p, *v] for t, p, v in rows]))
            print(f'  Horizons {name}: {len(rows)} rows', flush=True)
    return dyn.Bodies(CACHE)


class Mission:
    def __init__(self):
        self.b = fetch_bodies()
        self.model = dyn.GeoModel(self.b)

    # --- propagation ---------------------------------------------------------------------

    def prop(self, t0, x, t1):
        return dyn.state_at(self.model, t0, x, t1, SRP)

    def track(self, t0, x, t1):
        return dyn.integrate(self.model, t0, x, t1, srp=SRP)

    def moon(self, t):
        return np.concatenate([self.b.moon_geo(t), self.b.moon_geo.velocity(t)])

    def kick(self, t, x, dv, al, be, lunar=False):
        """Impulse dv (km/s) at in-plane angle al from the velocity (towards the radial-out side
        at +90 deg) and out-of-plane angle be (towards the orbit normal); Earth- or Moon-relative."""
        rel = x - self.moon(t) if lunar else x
        T = rel[3:] / np.linalg.norm(rel[3:])
        N = np.cross(rel[:3], rel[3:])
        N /= np.linalg.norm(N)
        B = np.cross(N, T)
        u = np.cos(be) * (np.cos(al) * T + np.sin(al) * B) + np.sin(be) * N
        y = np.array(x, float)
        y[3:] += dv * u
        return y

    # --- events --------------------------------------------------------------------------

    def _refine(self, t0, x0, ta, tb, g):
        """Root of g(t, state) in [ta, tb] (bisection on fresh propagations)."""
        xa = self.prop(t0, x0, ta)
        ga = g(ta, xa)
        for _ in range(40):
            tm = 0.5 * (ta + tb)
            xm = self.prop(ta, xa, tm)
            gm = g(tm, xm)
            if (gm > 0) == (ga > 0):
                ta, xa, ga = tm, xm, gm
            else:
                tb = tm
            if tb - ta < 0.01:
                break
        return ta, xa

    def event(self, t0, x0, g, t_end, nth=0, sign=0):
        """nth root of g(t, state) after t0 (sign: +1 rising, -1 falling, 0 either)."""
        ts, ys = self.track(t0, x0, t_end)
        gs = np.array([g(t, y) for t, y in zip(ts, ys)])
        k = 0
        for i in range(len(ts) - 1):
            if (gs[i] > 0) != (gs[i + 1] > 0) and (sign == 0 or sign * (gs[i + 1] - gs[i]) > 0):
                if k == nth:
                    return self._refine(ts[i], ys[i], ts[i], ts[i + 1], g)
                k += 1
        return None

    def apsis(self, t0, x0, kind, t_end):
        return self.event(t0, x0, lambda t, y: y[:3] @ y[3:], t_end, sign=-1 if kind == 'apo' else 1)

    def perilune(self, t0, x0, t_end):
        def g(t, y):
            d = y - self.moon(t)
            return d[:3] @ d[3:]
        return self.event(t0, x0, g, t_end, sign=1)

    # --- the Earth-Moon synodic frame centered on L2 (Liu & Li Fig. 5b) ------------------

    def synodic(self, t, x):
        m = self.moon(t)
        X = m[:3] / np.linalg.norm(m[:3])
        Z = np.cross(m[:3], m[3:])
        w = Z / (m[:3] @ m[:3])
        Z /= np.linalg.norm(Z)
        Y = np.cross(Z, X)
        R = np.array([X, Y, Z])
        r = x[:3] - (1 + GAMMA2) * m[:3]
        v = x[3:] - (1 + GAMMA2) * m[3:] - np.cross(w, r)
        return np.concatenate([R @ r, R @ v])


def lunar_elements(rel):
    """Perilune radius and inclination (deg, to the ICRF equator) of a Moon-relative state."""
    gm = dyn.GM['moon']
    r, v = rel[:3], rel[3:]
    h = np.cross(r, v)
    e = np.cross(v, h) / gm - r / np.linalg.norm(r)
    p = h @ h / gm
    return p / (1 + np.linalg.norm(e)), np.degrees(np.arccos(h[2] / np.linalg.norm(h)))


def newton(f, p, scale, tol, iters=20, verbose=False):
    """Minimum-norm Newton on the residuals f(p) = 0 (fewer residuals than parameters stay
    near the starting point). scale: typical parameter steps for the finite differences."""
    p = np.array(p, float)
    scale = np.asarray(scale, float)
    r = f(p)
    for it in range(iters):
        if verbose:
            print('    newton', it, np.round(r, 4), flush=True)
        if np.all(np.abs(r) < tol):
            return p, r
        J = np.empty((len(r), len(p)))
        for j in range(len(p)):
            q = p.copy()
            q[j] += 1e-3 * scale[j]
            J[:, j] = (f(q) - r) / (1e-3 * scale[j])
        Js = J * scale
        step = np.linalg.lstsq(Js, -r, rcond=None)[0] * scale
        lam = 1.0
        while True:
            rn = f(p + lam * step)
            if np.linalg.norm(rn) < np.linalg.norm(r):
                p, r = p + lam * step, rn
                break
            lam *= 0.5
            if lam < 1e-2:
                # Newton direction fails: damped (Levenberg-Marquardt) steps
                g = np.linalg.norm(Js.T @ Js)
                for mu in (1e-4, 1e-3, 1e-2, 1e-1, 1.0):
                    dmp = np.linalg.solve(Js.T @ Js + mu * g * np.eye(len(p)), -Js.T @ r) * scale
                    rn = f(p + dmp)
                    if np.linalg.norm(rn) < np.linalg.norm(r):
                        p, r = p + dmp, rn
                        break
                else:
                    return p, r
                break
    return p, r



def osculating_periapsis(rel, gm):
    """Two-body time to the next periapsis (s), its radius and the inclination (deg, to the ICRF
    equator) of an inbound relative state, elliptic or hyperbolic."""
    r, v = rel[:3], rel[3:]
    rn = np.linalg.norm(r)
    h = np.cross(r, v)
    e = np.linalg.norm(np.cross(v, h) / gm - r / rn)
    rp = h @ h / gm / (1 + e)
    inc = np.degrees(np.arccos(h[2] / np.linalg.norm(h)))
    a = 1.0 / (2.0 / rn - v @ v / gm)
    s = np.sign(r @ v)
    if a > 0:
        E = np.arccos(np.clip((1.0 - rn / a) / e, -1.0, 1.0)) * s
        n = np.sqrt(gm / a ** 3)
        dt = -(E - e * np.sin(E)) / n
        if dt < 0:
            dt += 2 * np.pi / n
    else:
        F = np.arccosh(max((1.0 - rn / a) / e, 1.0)) * s
        dt = -(e * np.sinh(F) - F) / np.sqrt(gm / (-a) ** 3)
    return dt, rp, inc


def brake_perilune():
    """Perilune radius and speed at brake 1 for which the three anti-velocity brakes of Table 5,
    at one perilune (Fig. 5c), make brake 2 fall 3 revolutions after brake 1 and brake 3 7
    revolutions after brake 2 (two-body)."""
    gm = dyn.GM['moon']
    want = np.array([(TB[10] - TB[9]) / 3, (TB[11] - TB[10]) / 7])

    def periods(x):
        v, out = x[1], []
        for k in (9, 10):
            v -= DV[k]
            out.append(2 * np.pi * np.sqrt((1 / (2 / x[0] - v * v / gm)) ** 3 / gm))
        return np.array(out)
    x = np.array([1937.0, 2.214])
    for _ in range(30):
        r = periods(x) - want
        J = np.column_stack([(periods(x + d) - periods(x)) / np.linalg.norm(d) for d in ([1e-3, 0], [0, 1e-7])])
        x = x - np.linalg.solve(J, r)
    return x


def departure(mission, t0, x0, span, radius=5.5e4):
    """Which way an L2 orbit falls off (sign of the synodic x when it leaves `radius` km of L2;
    0 if it stays) and when."""
    ts, ys = mission.track(t0, x0, t0 + span)
    for t, y in zip(ts, ys):
        s = mission.synodic(t, y)
        if np.linalg.norm(s[:3]) > radius:
            return np.sign(s[0]), t
    return 0.0, ts[-1]


# --- the flown trajectory --------------------------------------------------------------------

SOLUTION = os.path.join(HERE, 'solution.json')
END = TB[11] + DAY                       # a day in the final orbit; nothing is known after it
SEL_START = TB[9] - 2 * DAY              # Moon-centered table from two days before brake 1
# [L16] Fig. 5b, read off the plot (+-500 km): positions at burns 6-9 relative to EML2
FIG5B = {5: (800, -21800), 6: (3300, 16800), 7: (2700, 3800), 8: (-13600, 22800)}


def secant(f, x0, x1, tol, maxstep):
    f0, f1 = f(x0), f(x1)
    for _ in range(40):
        if abs(f1) < tol:
            return x1
        dx = np.clip(-f1 * (x1 - x0) / (f1 - f0), -maxstep, maxstep)
        x0, x1, f0 = x1, x1 + dx, f1
        f1 = f(x1)
    raise RuntimeError('secant did not converge')


class Flight:
    """The trajectory for a solution: Table 1's state with its speed changed by dv0, then the
    12 impulses of Table 5 with the solution's directions."""

    def __init__(self, mission, sol):
        self.m, self.sol = mission, sol

    def start(self):
        x = X0.copy()
        x[3:] *= 1 + self.sol['dv0_m_s'] * 1e-3 / np.linalg.norm(x[3:])
        return x

    def burn(self, k, x):
        b = self.sol['burns'][k]
        return self.m.kick(TB[k], x, DV[k], b['al'], b['be'], lunar=b['frame'] == 'moon')

    def after(self, k):
        """State just after burn k (k = -1: the start)."""
        t, x = T0, self.start()
        for j in range(k + 1):
            x = self.burn(j, self.m.prop(t, x, TB[j]))
            t = TB[j]
        return x

    def segments(self):
        """(start time, state, end time) of the ballistic arcs between the burns."""
        out, t, x = [], T0, self.start()
        for j in range(len(TB)):
            out.append((t, x, TB[j]))
            x = self.burn(j, self.m.prop(t, x, TB[j]))
            t = TB[j]
        out.append((t, x, END))
        return out


def solve(mission, sol, verbose=True, steps=(1, 2, 3, 4)):
    """Re-derives the burn directions, starting from those in `sol`:
    1. dv0 puts apogee at the apogee maneuver (burn 1 direction fixed), the in-plane angle of
       burn 2 the next perigee at the perigee maneuver.
    2. Burn 3's in-plane angle: bisection onto the stable manifold of the EML2 orbits (the arrival
       neither falls back towards the Moon nor away), with burn 5 anti-velocity at the flyby.
    3. Burns 6 and 7: in-plane angle bisected onto the stable manifold again.
    4. Burns 8 and 9: brakes 1, 2 and 3 at perilunes (Fig. 5c): first aimed with the two-body
       orbit 12 h before brake 1 at the perilune and energy of brake_perilune, then refined."""
    m, f = mission, Flight(mission, sol)
    b = sol['burns']

    def say(*a):
        if verbose:
            print(*a, flush=True)

    def apogee_err(dv0):
        sol['dv0_m_s'] = dv0
        return (m.apsis(TB[0], f.after(0), 'apo', TB[0] + 12 * DAY)[0] - TB[1]) / HOUR

    def perigee_err(a2):
        b[1]['al'] = a2
        return (m.apsis(TB[1], f.burn(1, x2), 'peri', TB[1] + 12 * DAY)[0] - TB[2]) / HOUR
    if 1 in steps:
        sol['dv0_m_s'] = secant(apogee_err, sol['dv0_m_s'], sol['dv0_m_s'] + 0.01, 1e-5, 0.5)
        x2 = m.prop(TB[0], f.after(0), TB[1])
        b[1]['al'] = secant(perigee_err, b[1]['al'], b[1]['al'] + 1e-3, 1e-5, np.radians(3))
        say(f'  dv0 {sol["dv0_m_s"]:.3f} m/s, burn 2 in-plane {np.degrees(b[1]["al"]):.3f} deg')

    def bisect(k, x_before, t_check, span, half):
        """Bisect burn k's in-plane angle within +-half of its value onto the stable manifold."""
        def side(al):
            b[k]['al'] = al
            x = f.burn(k, x_before)
            for j in range(k + 1, 5):           # through the flyby when k = 2
                x = f.burn(j, m.prop(TB[j - 1], x, TB[j]))
            t0 = TB[max(k, 4)]
            return departure(m, t_check, m.prop(t0, x, t_check), span)
        al0 = b[k]['al']
        lo, hi = al0 - half, al0 + half
        slo, shi = side(lo)[0], side(hi)[0]
        if slo == shi:
            raise RuntimeError(f'burn {k + 1}: no stable-manifold crossing within +-{np.degrees(half)} deg')
        best = None
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            s, tl = side(mid)
            if best is None or tl > best[1]:
                best = (mid, tl)
            if s == 0.0:
                break
            if s == slo:
                lo = mid
            else:
                hi = mid
        b[k]['al'] = best[0]
        return best[1]

    if 2 in steps:
        x3 = m.prop(TB[1], f.after(1), TB[2])
        tl = bisect(2, x3, TB[5], 60 * DAY, np.radians(0.02))
        say(f'  burn 3 in-plane {np.degrees(b[2]["al"]):.10f} deg: ballistic in the EML2 region until {utc_text(tl)}')
    if 3 in steps:
        for k in (5, 6):
            xk = m.prop(TB[k - 1], f.after(k - 1), TB[k])
            tl = bisect(k, xk, TB[k], TB[k + 1] - TB[k] + 20 * DAY, np.radians(1.0))
            say(f'  burn {k + 1} in-plane {np.degrees(b[k]["al"]):.6f} deg: stays until {utc_text(tl)}')
    if 4 not in steps:
        return sol

    x8 = m.prop(TB[6], f.after(6), TB[7])

    def from_perilune(t, y):
        """Hours from the nearest perilune of the osculating lunar orbit at t (> 0: past it)."""
        gm = dyn.GM['moon']
        rel = y - m.moon(t)
        r, v = rel[:3], rel[3:]
        rn = np.linalg.norm(r)
        a = 1.0 / (2.0 / rn - v @ v / gm)
        E = np.arctan2(r @ v / np.sqrt(gm * a), 1.0 - rn / a)
        e = np.hypot(r @ v / np.sqrt(gm * a), 1.0 - rn / a)
        return (E - e * np.sin(E)) / np.sqrt(gm / a ** 3) / HOUR

    def lunar_residuals(p):
        # brakes at perilunes (Fig. 5c), from the osculating orbit at each brake: smooth, unlike
        # counting perilunes between them
        b[7]['al'], b[7]['be'], b[8]['al'], b[8]['be'] = p
        y = m.prop(TB[8], f.burn(8, m.prop(TB[7], f.burn(7, x8), TB[8])), TB[9])
        if np.linalg.norm(y[:3] - m.moon(TB[9])[:3]) > 2e4:
            return np.full(4, 1e3)
        res = [from_perilune(TB[9], y)]
        for k in (9, 10):
            y = m.prop(TB[k], f.burn(k, y), TB[k + 1])
            res.append(from_perilune(TB[k + 1], y))
        return np.array(res)

    rp_want, vp_want = brake_perilune()
    e_want = vp_want ** 2 / 2 - dyn.GM['moon'] / rp_want
    tgt = np.array([0.0, rp_want, e_want])

    def osculating_residuals(p):
        # aim with the two-body orbit 12 h before brake 1 (smooth); the targets are corrected
        # for the perturbations below
        b[7]['al'], b[7]['be'], b[8]['al'], b[8]['be'] = p
        tref = TB[9] - 12 * HOUR
        x = m.prop(TB[8], f.burn(8, m.prop(TB[7], f.burn(7, x8), TB[8])), tref)
        rel = x - m.moon(tref)
        dt, rp, inc = osculating_periapsis(rel, dyn.GM['moon'])
        en = rel[3:] @ rel[3:] / 2 - dyn.GM['moon'] / np.linalg.norm(rel[:3])
        return np.array([(tref + dt - TB[9]) / HOUR - tgt[0], (rp - tgt[1]) / 100, (en - tgt[2]) * 100])
    p = np.array([b[7]['al'], b[7]['be'], b[8]['al'], b[8]['be']])
    for it in range(4):
        p, r = newton(osculating_residuals, p, [0.02] * 4, [1e-3] * 3, iters=30, verbose=verbose)
        osculating_residuals(p)
        y = m.prop(TB[8], f.burn(8, m.prop(TB[7], f.burn(7, x8), TB[8])), TB[9])
        rel = y - m.moon(TB[9])
        rp, _ = lunar_elements(rel)
        en = rel[3:] @ rel[3:] / 2 - dyn.GM['moon'] / np.linalg.norm(rel[:3])
        miss = np.array([-from_perilune(TB[9], y), rp - rp_want, en - e_want])
        say(f'  aim {it}: perilune {miss[0] * 60:+.1f} min, {miss[1]:+.1f} km, energy {miss[2]:+.5f} km^2/s^2')
        if abs(miss[0]) < 0.01 and abs(miss[1]) < 1.0 and abs(miss[2]) < 2e-5:
            break
        tgt = tgt - miss
    if lunar_residuals(p)[0] >= 1e3:
        raise RuntimeError('burns 8, 9: no lunar arrival from the starting directions')
    p, r = newton(lunar_residuals, p, [0.02] * 4, [2e-3] * 3, iters=30, verbose=verbose)
    lunar_residuals(p)
    y = m.prop(TB[8], f.after(8), TB[9])
    say(f'  burns 8, 9: brakes {r[0] * 60:+.1f}, {r[1] * 60:+.1f}, {r[2] * 60:+.1f} min after perilunes, '
        f'inclination {lunar_elements(y - m.moon(TB[9]))[1]:.2f} deg')
    return sol


def load_solution():
    with open(SOLUTION) as fh:
        return json.load(fh)


def save_solution(sol):
    with open(SOLUTION, 'w') as fh:
        json.dump(sol, fh, indent=1)
        fh.write('\n')


# --- tables ------------------------------------------------------------------------------------

TOL_KM = 1.0


def bake(mission, flight, name, center, t0, t1):
    """Knots of the flight over [t0, t1], Earth- or Moon-centered, one ballistic arc at a time
    (each arc after an impulse starts a second later, as in isee3_reconstruct.py)."""
    m = mission
    gm = dyn.GM['earth'] if center == 399 else dyn.GM['moon']
    interp = hb.make_interp(gm)
    knots, worst = [], 0.0
    for ta, xa, tb in flight.segments():
        a, b = max(ta, t0), min(tb, t1)
        if b <= a:
            continue
        ts, ys = dyn.integrate(m.model, ta, xa, b, srp=SRP)
        tq = np.append(np.arange(a, b, 600.0), b)
        st = dyn.dense(m.model, ts, ys, SRP, tq)[:, :6]
        near = ((np.linalg.norm(st[:, :3], axis=1) < 5.0e4) |
                (np.linalg.norm(st[:, :3] - m.b.moon_geo(tq), axis=1) < 4.0e4))
        extra = [np.arange(tq[i] + 60.0, tq[i + 1], 60.0) for i in np.where(near[:-1] | near[1:])[0]]
        tq = np.unique(np.concatenate([tq] + extra))
        st = dyn.dense(m.model, ts, ys, SRP, tq)[:, :6]
        if center == 301:
            st = st - np.hstack([m.b.moon_geo(tq), m.b.moon_geo.velocity(tq)])
        rows = [(float(t), tuple(map(float, s[:3])), tuple(map(float, s[3:]))) for t, s in zip(tq, st)]
        kn = hb.decimate(rows, TOL_KM, interp)
        worst = max(worst, hb.max_error(rows, kn, interp))
        if knots:
            t, p, v = kn[0]
            kn[0] = (t + 1.0, tuple(np.array(p) + np.array(v)), v)
        knots += kn
    hb.write_eph(os.path.join(OUT, name), 0, center, knots, gm)   # no NAIF id
    print(f'  {name}: {len(knots)} knots, max error {worst:.2f} km')


# --- report ------------------------------------------------------------------------------------

def report(mission, flight):
    """The events of the reconstruction, to compare with the paper."""
    m = mission
    x = flight.start()
    print(f'  initial speed {flight.sol["dv0_m_s"]:+.3f} m/s from Table 1')
    e = m.apsis(TB[0], flight.after(0), 'apo', TB[1] + DAY)
    print(f'  apogee {utc_text(e[0])} UTC, {np.linalg.norm(e[1][:3]):.0f} km from the Earth')
    e = m.apsis(TB[1], flight.after(1), 'peri', TB[2] + DAY)
    print(f'  perigee {utc_text(e[0])} UTC, {np.linalg.norm(e[1][:3]) - 6378.1366:.0f} km up')
    e = m.perilune(TB[3], flight.after(3), TB[4] + HOUR)
    rel = e[1] - m.moon(e[0])
    print(f'  flyby perilune {utc_text(e[0])} UTC, {np.linalg.norm(rel[:3]) - MOON_R:.0f} km up, '
          f'{np.linalg.norm(rel[3:]):.3f} km/s')
    for k, (fx, fy) in FIG5B.items():
        s = m.synodic(TB[k], m.prop(TB[k - 1], flight.after(k - 1), TB[k]))
        print(f'  burn {k + 1}: ({s[0]:7.0f}, {s[1]:7.0f}, {s[2]:7.0f}) km from EML2; Fig. 5b ({fx}, {fy})')
    far = 0.0
    for ta, xa, tb in flight.segments():
        if tb <= TB[8]:
            ts, ys = m.track(ta, xa, tb)
            far = max(far, np.linalg.norm(ys[:, :3], axis=1).max())
    print(f'  farthest from the Earth {far:.0f} km')
    e = m.perilune(TB[8], flight.after(8), TB[9] + 6 * HOUR)
    rel = e[1] - m.moon(e[0])
    print(f'  lunar arrival perilune {utc_text(e[0])} UTC, {np.linalg.norm(rel[:3]) - MOON_R:.0f} km up, '
          f'{np.linalg.norm(rel[3:]):.4f} km/s, inclination {lunar_elements(rel)[1]:.2f} deg')
    gm = dyn.GM['moon']
    for k in (9, 10, 11):
        rel = flight.after(k) - m.moon(TB[k])
        a = 1 / (2 / np.linalg.norm(rel[:3]) - rel[3:] @ rel[3:] / gm)
        rp, inc = lunar_elements(rel)
        print(f'  after brake {k - 8}: {rp - MOON_R:.0f} x {2 * a - rp - MOON_R:.0f} km, '
              f'{2 * np.pi * np.sqrt(a ** 3 / gm) / HOUR:.2f} h, inclination {inc:.2f} deg')


def main():
    m = Mission()
    sol = load_solution()
    if '--solve' in sys.argv:
        steps = tuple(int(a) for a in sys.argv[sys.argv.index('--solve') + 1:] if a.isdigit()) or (1, 2, 3, 4)
        print(f'Solving the burn directions (steps {steps}) ...')
        save_solution(solve(m, sol, steps=steps))
    flight = Flight(m, sol)
    print('Events:')
    report(m, flight)
    print('Tables:')
    bake(m, flight, 'ce5t1_geo.eph', 399, T0, SEL_START)
    bake(m, flight, 'ce5t1_sel.eph', 301, SEL_START, END)


if __name__ == '__main__':
    main()
