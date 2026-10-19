"""Join a burn in a baked table to the arc that follows a seam after it.

Usage:
    python rejoin_burn.py <config.toml> <target name> <burn start> <burn end> <seam>

(times TDB, 'YYYY-MM-DD HH:MM[:SS]'; <seam> within a few minutes of the step).
Some Horizons trajectories are chains of arcs fitted to tracking ("tag-up fits");
the arc that carries an orbit insertion burn can hold the planned burn, and the
next arc, fitted after it, the orbit actually reached (Hope: 2,400 km and 126 m/s
apart at 2021-02-10 00:04 TDB, eight hours after the burn). This takes the later
arc as the truth and rebuilds the stretch from the burn to the seam:

1. 1-minute Horizons samples from the burn's end to just past the seam; the step
   is the minute with the largest velocity change two-body gravity does not
   explain. Over each minute before it, that change is the perturbation (the
   planet's oblateness near it, ~0.1 m/s a minute just after Hope's burn).
2. The first sample after the step is integrated back to the burn's end under
   two-body gravity plus those perturbations (the two arcs are close where they
   matter). The same back-integration of the last sample before the step must
   land on the burn's end: its miss is printed as a check.
3. The burn is integrated again from its start with the planned thrust (the
   acceleration two-body gravity does not explain between the table's knots,
   which through the burn must be every 1-minute sample: ephem.toml keep_all)
   plus an extra acceleration A + B u (u from 0 to 1 over the burn) scaled by the
   thrust's profile, so that it acts only while thrusting; A and B are solved for
   so that the burn ends on the later arc.
4. The table's knots through the burn are replaced by the integrated states, those
   between the burn and the step by the back-integrated path, reduced to the
   target's tolerance. All other knots stay. Run it after (re)baking the target.

Only the Python standard library is used.
"""

import math
import os
import sys
import tomllib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import horizons_bake as hb  # noqa: E402


def add(a, b, s=1.0):
    return tuple(a[i] + s * b[i] for i in range(3))


def solve(m, y):
    """Gaussian elimination with partial pivoting."""
    n = len(y)
    a = [row[:] + [y[i]] for i, row in enumerate(m)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(a[r][c]))
        a[c], a[p] = a[p], a[c]
        for r in range(n):
            if r != c:
                f = a[r][c] / a[c][c]
                a[r] = [a[r][k] - f * a[c][k] for k in range(n + 1)]
    return [a[i][n] / a[i][i] for i in range(n)]


def rk4(p, v, t, h, acc):
    """One RK4 step of x'' = acc(t, x) (h may be negative)."""
    k1v, k1p = acc(t, p), v
    k2v, k2p = acc(t + h / 2, add(p, k1p, h / 2)), add(v, k1v, h / 2)
    k3v, k3p = acc(t + h / 2, add(p, k2p, h / 2)), add(v, k2v, h / 2)
    k4v, k4p = acc(t + h, add(p, k3p, h)), add(v, k3v, h)
    return (tuple(p[k] + h / 6 * (k1p[k] + 2 * k2p[k] + 2 * k3p[k] + k4p[k]) for k in range(3)),
            tuple(v[k] + h / 6 * (k1v[k] + 2 * k2v[k] + 2 * k3v[k] + k4v[k]) for k in range(3)))


def unexplained(rows, gm):
    """Mean acceleration two-body gravity does not explain over each interval."""
    out = []
    for (ta, pa, va), (tb, _, vb) in zip(rows, rows[1:]):
        _, cv = hb.kepler_propagate(pa, va, gm, tb - ta)
        out.append(tuple((vb[k] - cv[k]) / (tb - ta) for k in range(3)))
    return out


def integrate(p, v, times, extra_acc, gm, reverse=False):
    """States at `times` (ascending) from the state at times[0] (or times[-1] if
    reverse), under two-body gravity plus extra_acc(n, t) on interval n."""
    def gravity(pp):
        r3 = math.sqrt(sum(x * x for x in pp)) ** 3
        return tuple(-gm * pp[k] / r3 for k in range(3))

    order = range(len(times) - 1)
    if reverse:
        order = reversed(order)
    out = {(len(times) - 1 if reverse else 0): (p, v)}
    for n in order:
        ta, tb = times[n], times[n + 1]
        steps = max(1, round(tb - ta))
        h = (tb - ta) / steps
        t = tb if reverse else ta
        for _ in range(steps):
            p, v = rk4(p, v, t, -h if reverse else h,
                       lambda tt, pp, n=n: add(gravity(pp), extra_acc(n, tt)))
            t += -h if reverse else h
        out[n if reverse else n + 1] = (p, v)
    return [(times[i],) + out[i] for i in range(len(times))]


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    with open(sys.argv[1], 'rb') as f:
        cfg = tomllib.load(f)
    tgt = next(t for t in cfg['target'] if t['name'] == sys.argv[2])
    root = os.path.dirname(os.path.abspath(sys.argv[1]))
    out = os.path.normpath(os.path.join(root, tgt.get('output', os.path.join(cfg['defaults'].get('output_dir', '.'),
                                                                             tgt['name'] + '.eph'))))
    tol = tgt.get('tolerance_km', cfg['defaults'].get('tolerance_km', 1.0))
    t_a, t_b, t_seam = ((hb.jd_from_text(s) - hb.J2000_JD) * hb.DAY for s in sys.argv[3:6])
    table = hb.read_eph(out)
    gm = table.gm
    if not gm:
        sys.exit('needs a Kepler-relative (AXEPH2) table')
    knots = table.knots
    burn = [i for i, k in enumerate(knots) if t_a - 1.0 <= k[0] <= t_b + 1.0]
    if len(burn) < 3 or any(knots[j][0] - knots[i][0] > 61.0 for i, j in zip(burn, burn[1:])):
        sys.exit('the burn must have a knot every minute (keep_all)')
    times = [knots[i][0] for i in burn]
    t0, t1 = times[0], times[-1]

    # 1. The stretch from the burn's end past the seam, and the step in it.
    rows = hb.fetch(tgt['command'], tgt['center'], t1 / hb.DAY + hb.J2000_JD, (t_seam + 300.0) / hb.DAY + hb.J2000_JD,
                    hb.MIN_STEP_MINUTES)
    pert = unexplained(rows, gm)
    jump = max(range(len(pert)), key=lambda n: sum(x * x for x in pert[n]))
    if math.sqrt(sum(x * x for x in pert[jump])) * (rows[jump + 1][0] - rows[jump][0]) < 1e-3:
        sys.exit('no step of more than 1 m/s between the burn and the seam')
    pert[jump] = pert[jump - 1] if jump > 0 else (0.0, 0.0, 0.0)
    coast = [r[0] for r in rows[:jump + 2]]

    # 2. Back to the burn's end along the perturbations of the arc before the step.
    def perturbation(n, t):
        return pert[n]
    t_post, p_post, v_post = rows[jump + 1]
    path = integrate(p_post, v_post, coast, perturbation, gm, reverse=True)
    check = integrate(rows[jump][1], rows[jump][2], coast[:-1], perturbation, gm, reverse=True)[0]
    print(f'{out}: step of {math.dist(hb.kepler_propagate(rows[jump][1], rows[jump][2], gm, t_post - rows[jump][0])[0], p_post):.1f} km '
          f'between {hb.calendar(rows[jump][0])} and {hb.calendar(t_post)}')
    print(f'    the arc before it, integrated back to {hb.calendar(t1)}: {math.dist(check[1], rows[0][1]) * 1e3:.0f} m and '
          f'{math.dist(check[2], rows[0][2]) * 1e3:.3f} m/s off the burn\'s end (the check)')
    target_p, target_v = path[0][1], path[0][2]
    print(f'    the arc after it there: {math.dist(target_p, rows[0][1]):.1f} km and '
          f'{math.dist(target_v, rows[0][2]) * 1e3:.1f} m/s off')

    # 3. The burn, again, ending on the later arc.
    thrust = unexplained([knots[i] for i in burn], gm)
    peak = max(math.sqrt(sum(x * x for x in a)) for a in thrust)
    shape = [math.sqrt(sum(x * x for x in a)) / peak for a in thrust]

    def burn_path(extra):
        def acc(n, t):
            u = (t - t0) / (t1 - t0)
            return tuple(thrust[n][k] + shape[n] * (extra[k] + extra[3 + k] * u) for k in range(3))
        return integrate(knots[burn[0]][1], knots[burn[0]][2], times, acc, gm)

    def miss(extra):
        _, p, v = burn_path(extra)[-1]
        return [x - y for x, y in zip(list(p) + list(v), list(target_p) + list(target_v))]

    extra = [0.0] * 6
    for _ in range(4):
        m0 = miss(extra)
        jac = []
        for k in range(6):
            e = extra[:]
            e[k] += 1e-6  # km/s^2
            jac.append([(a - b) / 1e-6 for a, b in zip(miss(e), m0)])
        step = solve([[jac[c][r] for c in range(6)] for r in range(6)], [-x for x in m0])
        extra = [x + s for x, s in zip(extra, step)]
    states = burn_path(extra)
    end = miss(extra)
    planned = sum(math.sqrt(sum(x * x for x in a)) * (tb - ta) for a, ta, tb in zip(thrust, times, times[1:]))
    total = sum(math.sqrt(sum((a[k] + w * (extra[k] + extra[3 + k] * ((ta + tb) / 2 - t0) / (t1 - t0))) ** 2
                              for k in range(3))) * (tb - ta)
                for a, w, ta, tb in zip(thrust, shape, times, times[1:]))
    moved = max(math.dist(s[1], knots[i][1]) for s, i in zip(states, burn))
    print(f'    burn delta-v {planned * 1e3:.1f} -> {total * 1e3:.1f} m/s; the path through it moves up to '
          f'{moved:.1f} km; ends on the later arc to {math.sqrt(sum(x * x for x in end[:3])) * 1e3:.3f} m')

    # 4. The new knots.
    after = next(i for i, k in enumerate(knots) if k[0] >= t_post - 1.0)
    if abs(knots[after][0] - t_post) > 1.0:
        sys.exit(f'the table has no knot at {hb.calendar(t_post)}, the first sample after the step')
    stretch = hb.decimate(path[:-1] + [knots[after]], tol, table.interp)
    print(f'    {after - burn[-1] - 1} knots between the burn and the step -> {len(stretch) - 2}, max error '
          f'{hb.max_error(path, stretch, table.interp):.3f} km')
    new = knots[:burn[0]] + states + stretch[1:-1] + knots[after:]
    hb.write_eph(out, int(tgt['naif_id']), int(tgt['center_naif_id']), new, gm, table.max_gap)
    print(f'    wrote {out} ({len(new)} knots)')


if __name__ == '__main__':
    main()
