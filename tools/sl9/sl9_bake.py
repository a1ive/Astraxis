"""Bake comet Shoemaker-Levy 9 about Jupiter: assets/ephem/sl9_*.eph, one table per fragment
and one for the parent comet. Prints the bodies for assets/scenes/jupiter.toml and how each
impact compares with Chodas & Yeomans 1996.

Usage (needs the network; run outside the sandbox):
    python -I tools/sl9/sl9_bake.py

The orbits are JPL Horizons' small-body records D/1993 F2-A .. -W, the final solutions of
[CY96] P. W. Chodas & D. K. Yeomans 1996, in "The Collision of Comet Shoemaker-Levy 9 and
Jupiter" (Cambridge), pp. 1-30: Table 2's orbit references (A38, G52, Q63, ...) are the
records' `soln ref`. Jupiter-centered, one sample every two hours, five minutes near Jupiter,
one minute in the last six hours.

Each table ends where the fragment reaches the 100-mbar level, [CY96]'s definition of impact:
Horizons does not stop there, it carries the orbit on through the planet. The level is taken
50 km above the 1-bar ellipsoid (the scene's radii; an assumption, it moves the time by a
second).

The records say they are valid from 1992-07-15, and the orbits, run backward, do not meet at
one point: the perijove times differ by 45 minutes ([CY96] Fig. 2), which at 50 km/s puts the
fragments up to 94,000 km apart two hours after perijove, when they were still one body.
So each fragment leaves its parent:
  - The parent comet is fragment K's path up to the effective breakup, 2.5 h after the perijove
    of 1992 July 7.860364 TDB ([SEK] Z. Sekanina 1996, same book, pp. 55-80, Sec. 4.2). [CY96]
    Sec. 6 takes K as the parent's center of mass.
  - From its splitting time to the time given below, a fragment's position is
    parent + w (orbit - parent), with w a smooth step from 0 to 1; its velocity follows.
    After that, it is the JPL solution unchanged.
  - Q2 left Q1 on 1993 Apr 12 and P2 left Q on 1992 Dec 29 ([SEK] Sec. 4.3; Q is Q1 here).
    Their blends end at the first observation that separates them (1993 Jul 1, HST at 0.3
    arcsec, and 1993 Mar 27).
  - P1: [SEK] gives no date. It is taken to leave with P2 and join its orbit by the start of
    its data, 1993 Jul 1 ([CY96] Table 2): an assumption.
  - The other off-train fragments (B, F, N, T, U, V; [WEA] Table 2) also split later, "most" of
    them between July 1992 and the start of 1993 ([SEK] Sec. 4.3), but no dates are published:
    they leave K with the on-train ones.
J and M, lost by January 1994, and G2 have no orbit solutions and are left out.

Sizes and brightness: [WEA] H. A. Weaver et al. 1995, Science 267, 1282, Table 2. Relative size
times 1.4 is an upper limit on the diameter in km. The coma magnitudes share out the comet's
total m1 of 14 at discovery (IAU Circular 5725, 1993 Mar 24.355; 4.44 au from the Earth) in
proportion to [WEA]'s relative intensities of May 1994.

Only the Python standard library is used. Downloads are cached in tools/sl9/cache (not in git).
"""

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'tools', 'bake'))
import horizons_bake as hb  # noqa: E402

EPHEM = os.path.join(REPO, 'assets', 'ephem')
CACHE = os.path.join(HERE, 'cache')
CENTER = '500@599'
CENTER_ID = 599
GM = 126712764.1          # Jupiter system (DE440), as the scene's Jupiter
TOLERANCE_KM = 5.0
RADIUS_KM = 71492.0
POLAR_KM = 66854.0
LEVEL_KM = 50.0           # the 100-mbar level above 1 bar: see above
POLE_RA, POLE_DEC = math.radians(268.056595), math.radians(64.495303)  # pck00011
POLE = (math.cos(POLE_DEC) * math.cos(POLE_RA), math.cos(POLE_DEC) * math.sin(POLE_RA), math.sin(POLE_DEC))

PARENT_START = '1991-01-01'
TOTAL_M1 = 14.0           # IAUC 5725
DISCOVERY_DELTA_AU = 4.44 # Earth-comet distance on 1993 Mar 24 (Horizons)


def tdb(text):
    return (hb.jd_from_text(text) - hb.J2000_JD) * hb.DAY


PERIJOVE = tdb('1992-07-07') + 0.860364 * hb.DAY  # [SEK]
BREAKUP = PERIJOVE + 2.5 * 3600.0
UNTIL = tdb('1992-07-15')   # the records' validity starts here

# name, Horizons record, parent, split time, end of blend, [WEA] relative intensity and size,
# [CY96] Table 5: July day, accepted impact time (UT at the Earth), latitude, W longitude, 1 sigma
FRAGMENTS = [
    ('A', 90002227, 'K', BREAKUP, UNTIL, 1.0, 1.0, 16, '20:10:40', -43.35, 184),
    ('B', 90002228, 'K', BREAKUP, UNTIL, 1.5, 1.2, 17, '02:50:00', -43.22, 67),
    ('C', 90002229, 'K', BREAKUP, UNTIL, 2.2, 1.5, 17, '07:10:50', -43.47, 222),
    ('D', 90002230, 'K', BREAKUP, UNTIL, 1.0, 1.0, 17, '11:52:30', -43.53, 33),
    ('E', 90002231, 'K', BREAKUP, UNTIL, 4.2, 2.0, 17, '15:11:40', -43.54, 153),
    ('F', 90002232, 'K', BREAKUP, UNTIL, 2.3, 1.5, 18, '00:35:45', -43.68, 135),
    ('G', 90002233, 'K', BREAKUP, UNTIL, 8.2, 2.9, 18, '07:33:33', -43.66, 26),
    ('H', 90002234, 'K', BREAKUP, UNTIL, 4.9, 2.2, 18, '19:31:59', -43.79, 99),
    ('K', 90002235, None, BREAKUP, BREAKUP, 7.3, 2.7, 19, '10:24:17', -43.86, 278),
    ('L', 90002236, 'K', BREAKUP, UNTIL, 6.1, 2.5, 19, '22:16:49', -43.96, 348),
    ('N', 90002237, 'K', BREAKUP, UNTIL, 1.0, 1.0, 20, '10:29:20', -44.31, 71),
    ('Q1', 90002239, 'K', BREAKUP, UNTIL, 8.6, 2.9, 20, '20:13:53', -44.00, 63),
    ('P2', 90002246, 'Q1', tdb('1992-12-29'), tdb('1993-03-27'), 1.8, 1.4, 20, '15:21:11', -44.69, 249),
    ('P1', 90002238, 'P2', tdb('1992-12-29'), tdb('1993-07-01'), 0.38, 0.61, 20, '16:32:35', -45.02, 293),
    ('Q2', 90002247, 'Q1', tdb('1993-04-12'), tdb('1993-07-01'), 4.8, 2.2, 20, '19:44:00', -44.32, 46),
    ('R', 90002240, 'K', BREAKUP, UNTIL, 3.3, 1.8, 21, '05:34:57', -44.10, 42),
    ('S', 90002241, 'K', BREAKUP, UNTIL, 4.3, 2.1, 21, '15:16:30', -44.22, 33),
    ('T', 90002242, 'K', BREAKUP, UNTIL, 0.21, 0.46, 21, '18:09:56', -45.01, 141),
    ('U', 90002243, 'K', BREAKUP, UNTIL, 0.42, 0.65, 21, '22:00:02', -44.48, 278),
    ('V', 90002244, 'K', BREAKUP, UNTIL, 0.95, 0.97, 22, '04:23:20', -44.47, 149),
    ('W', 90002245, 'K', BREAKUP, UNTIL, 3.0, 1.7, 22, '08:06:16', -44.13, 283),
]
ORDER = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'K', 'L', 'N', 'P1', 'P2', 'Q1', 'Q2', 'R', 'S', 'T', 'U', 'V', 'W']


def sampled(command, a, b, minutes):
    """Rows in [a, b], seconds past J2000, cached by the query."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, f'{command.strip(";")}_{a:.3f}_{b:.3f}_{minutes}.json')
    if os.path.exists(path):
        with open(path, encoding='utf-8') as f:
            return [(r[0], tuple(r[1]), tuple(r[2])) for r in json.load(f)]
    rows = hb.fetch(command, CENTER, a / hb.DAY + hb.J2000_JD, b / hb.DAY + hb.J2000_JD, minutes)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(rows, f)
    return rows


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def level_km(p):
    z = dot(p, POLE) / norm(p)
    return 1.0 / math.sqrt((1.0 - z * z) / (RADIUS_KM + LEVEL_KM) ** 2 + z * z / (POLAR_KM + LEVEL_KM) ** 2)


def impact_utc_seconds(day, hms):
    """[CY96]'s time, UT at the Earth, as seconds past J2000 TDB at Jupiter (rough: 43 min of light)."""
    return tdb(f'1994-07-{day:02d} {hms[:5]}') + float(hms[6:]) + 61.2 - 43.0 * 60.0


def hermite(a, b, t):
    h = b[0] - a[0]
    s = (t - a[0]) / h
    h00, h10, h01, h11 = 2 * s**3 - 3 * s**2 + 1, s**3 - 2 * s**2 + s, -2 * s**3 + 3 * s**2, s**3 - s**2
    return tuple(h00 * a[1][i] + h10 * h * a[2][i] + h01 * b[1][i] + h11 * h * b[2][i] for i in range(3))


def rows_to_impact(command, start, day, hms):
    guess = impact_utc_seconds(day, hms)
    fast_until = min(UNTIL, guess)
    rows = []
    edges = [(start, 5, fast_until), (fast_until, 120, guess - 3 * hb.DAY),
             (guess - 3 * hb.DAY, 10, guess - 6 * 3600.0), (guess - 6 * 3600.0, 1, guess + 900.0)]
    for a, minutes, b in edges:
        if b <= a:
            continue
        part = sampled(command, a, b, minutes)
        rows.extend(r for r in part if not rows or r[0] > rows[-1][0])
    for r0, r1 in zip(rows, rows[1:]):
        if norm(r1[1]) <= level_km(r1[1]) and dot(r0[1], r0[2]) < 0.0:
            lo, hi = r0[0], r1[0]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if norm(hermite(r0, r1, mid)) > level_km(hermite(r0, r1, mid)):
                    lo = mid
                else:
                    hi = mid
            f = (lo - r0[0]) / (r1[0] - r0[0])
            velocity = tuple(r0[2][i] + f * (r1[2][i] - r0[2][i]) for i in range(3))
            return [r for r in rows if r[0] < lo] + [(lo, hermite(r0, r1, lo), velocity)]
    raise RuntimeError(f'{command}: no crossing of the 100-mbar level')


def state_of(table, t):
    """Position and velocity, by a central difference inside the table."""
    margin = 1e-3
    def at(x):
        return hb.eval_eph(table, min(max(x, table.start() + margin), table.end() - margin))
    dt = 30.0
    a, b = max(t - dt, table.start() + margin), min(t + dt, table.end() - margin)
    pa, pb = at(a), at(b)
    return at(t), tuple((pb[i] - pa[i]) / (b - a) for i in range(3))


def blended(rows, parent, split, until):
    out = []
    for t, p, v in rows:
        if t < split:
            continue
        if t >= until or parent is None:
            out.append((t, p, v))
            continue
        x = (t - split) / (until - split)
        w = x * x * (3.0 - 2.0 * x)
        dw = 6.0 * x * (1.0 - x) / (until - split)
        q, u = state_of(parent, t)
        out.append((t, tuple(q[i] + w * (p[i] - q[i]) for i in range(3)),
                    tuple(u[i] + w * (v[i] - u[i]) + dw * (p[i] - q[i]) for i in range(3))))
    return out


def write(name, naif_id, rows):
    interp = hb.make_interp(GM)
    knots = hb.decimate(rows, TOLERANCE_KM, interp)
    hb.write_eph(os.path.join(EPHEM, name + '.eph'), naif_id, CENTER_ID, knots, GM)
    table = hb.Table(knots, GM)
    print(f'{name}: {len(rows)} samples -> {len(knots)} knots, max error '
          f'{hb.max_error(rows, knots, interp):.2f} km, {hb.calendar(rows[0][0])} .. {hb.calendar(rows[-1][0])}')
    return table


def west_longitude(p, t):
    w = math.radians(284.95 + 870.536 * t / hb.DAY)  # pck00011, System III
    n = (-math.sin(POLE_RA), math.cos(POLE_RA), 0.0)
    e = (POLE[1] * n[2] - POLE[2] * n[1], POLE[2] * n[0] - POLE[0] * n[2], POLE[0] * n[1] - POLE[1] * n[0])
    x = tuple(math.cos(w) * n[i] + math.sin(w) * e[i] for i in range(3))
    y = tuple(-math.sin(w) * n[i] + math.cos(w) * e[i] for i in range(3))
    return (-math.degrees(math.atan2(dot(p, y), dot(p, x)))) % 360.0


def main():
    tables = {}
    k = next(f for f in FRAGMENTS if f[0] == 'K')
    # The parent comet: K up to the breakup, every minute around the perijove.
    command = f'{k[1]};'
    around = [(tdb(PARENT_START), 120, PERIJOVE - hb.DAY), (PERIJOVE - hb.DAY, 1, BREAKUP)]
    rows = []
    for a, b_minutes, b in around:
        rows.extend(r for r in sampled(command, a, b, b_minutes) if not rows or r[0] > rows[-1][0])
    rows = [r for r in rows if r[0] <= BREAKUP]
    tail = sampled(command, BREAKUP - 120.0, BREAKUP + 120.0, 1)
    last = min(tail, key=lambda r: abs(r[0] - BREAKUP))
    if last[0] != rows[-1][0]:
        rows = [r for r in rows if r[0] < last[0]] + [last]
    write('sl9_parent', k[1], rows)
    breakup = last[0]
    print(f'breakup {hb.calendar(breakup)}, perijove of the parent: '
          f'{min(norm(r[1]) for r in rows):.0f} km')

    fragments = sorted(FRAGMENTS, key=lambda f: 0 if f[2] is None else 1)
    report = {}
    for name, record, parent, split, until, _, _, day, hms, lat, lon in fragments:
        split = breakup if split == BREAKUP else split
        rows = rows_to_impact(f'{record};', split, day, hms)
        rows = blended(rows, tables.get(parent), split, until)
        tables[name] = write('sl9_' + name.lower(), record, rows)
        end = rows[-1]
        report[name] = (end, day, hms, lat, lon)

    print('\nimpacts (ours / [CY96] Table 5):')
    for name in ORDER:
        end, day, hms, lat, lon = report[name]
        p = end[1]
        ours_lat = math.degrees(math.asin(dot(p, POLE) / norm(p)))
        print(f'{name:3s} {hb.calendar(end[0])} TDB, {end[0] - impact_utc_seconds(day, hms):+6.0f} s, '
              f'lat {ours_lat:7.2f} / {lat}, W lon {west_longitude(p, end[0]):6.1f} / {lon}, '
              f'{norm(end[2]):.2f} km/s')

    total = sum(f[5] for f in FRAGMENTS)
    print('\n# --- the bodies for assets/scenes/jupiter.toml ---')
    for f in sorted(FRAGMENTS, key=lambda f: ORDER.index(f[0])):
        name, intensity, size = f[0], f[5], f[6]
        m1 = TOTAL_M1 - 5.0 * math.log10(DISCOVERY_DELTA_AU) + 2.5 * math.log10(total / intensity)
        print(f'''
[[bodies]]
name = "SL9 {name}"
parent = "Jupiter"
color = "#d9cfc2"
trail = "none"
radii_km = [{0.7 * size:.2f}]
[bodies.comet]
dust_only = true
dust_speed_km_s = 0.066
magnitude_m1_k1 = [{m1:.2f}, 0.0]
[bodies.orbit]
type = "ephemeris"
file = "ephem/sl9_{name.lower()}.eph"''')


if __name__ == '__main__':
    main()
