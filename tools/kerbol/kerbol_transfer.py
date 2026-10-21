"""Bake a schematic Kerbin -> Mun Hohmann transfer for the Kerbol scene.

    python -I kerbol_transfer.py [output.eph]

A DESIGNED trajectory, not a recorded flight: a craft in a circular parking orbit
80 km above Kerbin (just above the 70 km atmosphere) fires prograde at UT BURN_UT
onto a Hohmann ellipse whose apoapsis touches the Mun's orbit, coasts to the Mun,
and fires retrograde at the lowest point of its pass to stay in a circular orbit
CAPTURE_ALT_KM above the Mun. The motion is patched conics, exactly as KSP flies
it: around Kerbin alone until the craft enters the Mun's sphere of influence, then
around the Mun alone. All bodies' values (GMs, radii, the Mun's orbit and SOI)
come from the kOS dump kerbol_bodies.txt (KSP 1.12.5) through kerbol_scene.py.

The parking orbit's phase is chosen so that the burn at BURN_UT is the right one:
the burn point is 180 deg from where the Mun will be when the craft reaches
apoapsis. The craft arrives moving along the Mun's orbit and the Mun overtakes it
from behind, so the closest approach is set by the apoapsis radius, not by the
timing: the burn aims the apoapsis slightly inside the Mun's orbit, which makes
the craft pass on the Mun's inner side and the capture orbit prograde
(counterclockwise seen from the north, like every orbit in the game), at the
radius that puts the pass's lowest point CAPTURE_ALT_KM above the Mun.

Output: an AXEPH2 table relative to Kerbin (assets/ephem/kerbol_mun_transfer.eph
by default), from UT 0 to a day after the capture; after each impulse the next
arc starts a second later, as in tools/ce5t1/. Standard library only.
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "bake"))

import horizons_bake as hb  # noqa: E402
import kerbol_scene as ks  # noqa: E402

PARKING_ALT_KM = 80.0
CAPTURE_ALT_KM = 30.0
BURN_UT = 7200.0            # Year 1, Day 1, 2:00
END_AFTER_CAPTURE_S = 86400.0
TOL_KM = 0.05


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(a, s):
    return tuple(x * s for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


class Setup:
    def __init__(self, dump):
        by = ks.parse(dump)
        kerbin, mun = by["Kerbin"], by["Mun"]
        self.gm_kerbin = kerbin["mu_m3_s2"] / 1e9
        self.gm_mun = mun["mu_m3_s2"] / 1e9
        self.r_kerbin = kerbin["radius_m"] / 1e3
        self.r_mun = mun["radius_m"] / 1e3
        self.soi_mun = mun["soi_radius_m"] / 1e3
        # The Mun's orbit is a circle in the reference plane (e = i = node = arg_peri = 0).
        assert mun["eccentricity"] == 0 and mun["inclination_deg"] == 0 and mun["lan_deg"] == 0
        self.a_mun = mun["semi_major_axis_m"] / 1e3
        self.n_mun = 2.0 * math.pi / mun["period_s"]
        self.l0_mun = math.radians(mun["mean_anomaly_at_epoch_deg"] + mun["arg_periapsis_deg"])

    def mun_state(self, t):
        lam = self.l0_mun + self.n_mun * t
        v = self.a_mun * self.n_mun
        return ((self.a_mun * math.cos(lam), self.a_mun * math.sin(lam), 0.0),
                (-v * math.sin(lam), v * math.cos(lam), 0.0))


def circular_state(r, gm, lam):
    v = math.sqrt(gm / r)
    return (r * math.cos(lam), r * math.sin(lam), 0.0), (-v * math.sin(lam), v * math.cos(lam), 0.0)


def propagate(state, gm, dt):
    if dt == 0.0:
        return state
    return hb.kepler_propagate(state[0], state[1], gm, dt)


def transfer(s, apo_offset):
    """Patched-conic flight with the apoapsis apo_offset km from the Mun's orbit:
    the post-burn state, the burn's delta-v and, if it reaches the Mun's SOI, the
    entry time and the Mun-relative state there."""
    r0 = s.r_kerbin + PARKING_ALT_KM
    a_t = 0.5 * (r0 + s.a_mun + apo_offset)
    t_half = math.pi * math.sqrt(a_t ** 3 / s.gm_kerbin)
    arrival = s.l0_mun + s.n_mun * (BURN_UT + t_half)
    lam_burn = arrival - math.pi
    p, v = circular_state(r0, s.gm_kerbin, lam_burn)
    v_burn = math.sqrt(s.gm_kerbin * (2.0 / r0 - 1.0 / a_t))
    dv = v_burn - norm(v)
    post = (p, scale(v, v_burn / norm(v)))
    # Coast until the Mun's SOI (step, then bisect).
    t, step = BURN_UT, 60.0
    def rel(tt):
        c = propagate(post, s.gm_kerbin, tt - BURN_UT)
        m = s.mun_state(tt)
        return sub(c[0], m[0]), sub(c[1], m[1])
    while t < BURN_UT + 2.0 * t_half:
        if norm(rel(t + step)[0]) < s.soi_mun:
            lo, hi = t, t + step
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                lo, hi = (lo, mid) if norm(rel(mid)[0]) < s.soi_mun else (mid, hi)
            return post, dv, hi, rel(hi)
        t += step
    return post, dv, None, None


def periapsis(state, gm):
    """Time (from the state) and radius of the next periapsis, and the angular momentum's z."""
    p, v = state
    h = (p[1] * v[2] - p[2] * v[1], p[2] * v[0] - p[0] * v[2], p[0] * v[1] - p[1] * v[0])
    t, step = 0.0, 10.0
    while dot(*propagate(state, gm, t + step)) < 0.0:
        t += step
        if t > 1e6:
            raise RuntimeError("no periapsis")
    lo, hi = t, t + step
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if dot(*propagate(state, gm, mid)) < 0.0 else (lo, mid)
    return lo, norm(propagate(state, gm, lo)[0]), h[2]


def solve(s):
    """Apoapsis offset (km) whose pass has its periapsis CAPTURE_ALT_KM above the Mun, prograde."""
    target = s.r_mun + CAPTURE_ALT_KM
    def miss(offset):
        _, _, t_soi, st = transfer(s, offset)
        if t_soi is None:
            return None
        _, rp, hz = periapsis(st, s.gm_mun)
        return rp if hz > 0 else -rp
    grid = [-2000.0 + 10.0 * k for k in range(401)]
    values = [miss(o) for o in grid]
    for (o0, f0), (o1, f1) in zip(zip(grid, values), zip(grid[1:], values[1:])):
        if f0 is None or f1 is None or f0 <= 0 or f1 <= 0:
            continue
        if (f0 - target) * (f1 - target) <= 0:
            lo, hi, flo = o0, o1, f0 - target
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                fm = miss(mid)
                if fm is None or fm <= 0:
                    raise RuntimeError("root left the prograde branch")
                if (fm - target) * flo <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fm - target
            found = 0.5 * (lo + hi)
            assert abs(miss(found) - target) < 1e-3, "converged on a jump, not a root"
            return found
    raise RuntimeError("no prograde pass at the capture altitude")


def sample_arc(t0, t1, state_at, step):
    ts = [t0 + k * step for k in range(int((t1 - t0) / step))] + [t1]
    return [(t, *state_at(t)) for t in ts]


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "assets", "ephem",
                                                             "kerbol_mun_transfer.eph")
    s = Setup(os.path.join(HERE, "kerbol_bodies.txt"))
    apo_offset = solve(s)
    post, dv1, t_soi, st_soi = transfer(s, apo_offset)
    dt_p, rp, _ = periapsis(st_soi, s.gm_mun)
    t_cap = t_soi + dt_p
    p_cap, v_cap = propagate(st_soi, s.gm_mun, dt_p)
    v_circ = math.sqrt(s.gm_mun / rp)
    dv2 = norm(v_cap) - v_circ
    captured = (p_cap, scale(v_cap, v_circ / norm(v_cap)))
    t_end = t_cap + END_AFTER_CAPTURE_S

    r0 = s.r_kerbin + PARKING_ALT_KM
    parking0 = propagate((post[0], scale(post[1], math.sqrt(s.gm_kerbin / r0) / norm(post[1]))),
                         s.gm_kerbin, -BURN_UT)

    def in_mun(state0, t0):
        def f(t):
            c = propagate(state0, s.gm_mun, t - t0)
            m = s.mun_state(t)
            return add(c[0], m[0]), add(c[1], m[1])
        return f

    arcs = [
        sample_arc(0.0, BURN_UT, lambda t: propagate(parking0, s.gm_kerbin, t), 30.0),
        sample_arc(BURN_UT, t_soi, lambda t: propagate(post, s.gm_kerbin, t - BURN_UT), 30.0) +
        sample_arc(t_soi, t_cap, in_mun(st_soi, t_soi), 5.0)[1:],
        sample_arc(t_cap, t_end, in_mun(captured, t_cap), 10.0),
    ]
    interp = hb.make_interp(s.gm_kerbin)
    knots, worst = [], 0.0
    for rows in arcs:
        kn = hb.decimate(rows, TOL_KM, interp)
        worst = max(worst, hb.max_error(rows, kn, interp))
        if knots:  # the impulse: the new arc starts a second later
            t, p, v = kn[0]
            kn[0] = (t + 1.0, add(p, v), v)
        knots += kn
    hb.write_eph(out, 0, 0, knots, s.gm_kerbin)

    fmt = lambda t: f"{t:.1f} s ({ks_date(t)})"
    print(f"apoapsis {apo_offset:+.3f} km from the Mun's orbit ({s.a_mun + apo_offset:.3f} km from Kerbin)")
    print(f"burn 1 (prograde) at UT {fmt(BURN_UT)}: {dv1 * 1000:.1f} m/s, "
          f"{math.sqrt(s.gm_kerbin / r0) * 1000:.1f} -> {norm(post[1]) * 1000:.1f} m/s")
    print(f"enters the Mun's SOI at UT {fmt(t_soi)}, {(t_soi - BURN_UT) / 3600:.3f} h after the burn")
    print(f"burn 2 (retrograde) at UT {fmt(t_cap)}: {dv2 * 1000:.1f} m/s, "
          f"{norm(v_cap) * 1000:.1f} -> {v_circ * 1000:.1f} m/s, {rp - s.r_mun:.3f} km above the Mun")
    print(f"total {(dv1 + dv2) * 1000:.1f} m/s; table ends at UT {fmt(t_end)}")
    print(f"{out}: {len(knots)} knots, max error {worst * 1000:.1f} m")


def ks_date(t):
    year, rest = divmod(int(t), 9201600)
    day, rest = divmod(rest, 21600)
    return f"Year {year + 1}, Day {day + 1}, {rest // 3600}:{rest // 60 % 60:02d}:{rest % 60:02d}"


if __name__ == "__main__":
    main()
