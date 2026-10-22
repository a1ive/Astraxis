"""Bake DART's last days and LICIACube's flight past the Didymos system, relative to the
system barycenter: assets/ephem/dart_didymos.eph and assets/ephem/liciacube_didymos.eph.

Usage (after s547_bake.py, whose tables place Dimorphos):
    python approach_bake.py

JPL Horizons, DART (-135) and LICIACube (-210) relative to @20065803, the Didymos system
barycenter, one sample a minute, knots reduced to 5 m (cubic Hermite; Didymos' pull
bends these paths by less than a metre).

DART: Horizons joins the navigation reconstruction (to 2022-09-26 23:14:45.3 TDB) to a
48-second segment ending at the impact ("spc"), which starts 6.8 km off the
reconstruction and 1.3 degrees off its direction: a jump the scene would show
right before the impact. The reconstruction alone, carried on in a straight line,
passes 0.14 km from Dimorphos' center (solution 547, as the scene places it) 0.3 s
before the impact time, 1.9 km short. The table is the reconstruction up to that
segment, shifted by that constant offset so that it reaches Dimorphos' center at the
impact time (the center's position from assets/ephem/didymos_bary.eph and
dimorphos.eph), and ends there.

LICIACube: from its release on 2022-09-11 to a day after the flyby, as Horizons gives
it; it passes 57.6 km from Dimorphos' center (Horizons notes: 57.7 km).

Only the Python standard library is used.
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'tools', 'bake'))
import horizons_bake as hb  # noqa: E402

EPHEM = os.path.join(REPO, 'assets', 'ephem')
IMPACT = 717506133.3653613         # 2022-09-26 23:15:33.365 TDB
SEAM = IMPACT - 48.0               # the "spc" segment starts 48.07 s before the impact
BARYCENTER = 20065803
TOLERANCE_KM = 0.005


def jd(t):
    return t / hb.DAY + hb.J2000_JD


def dimorphos_center(t):
    did = hb.read_eph(os.path.join(EPHEM, 'didymos_bary.eph'))
    dim = hb.read_eph(os.path.join(EPHEM, 'dimorphos.eph'))
    a, b = hb.eval_eph(did, t), hb.eval_eph(dim, t)
    return tuple(a[i] + b[i] for i in range(3))


def write(name, naif_id, rows):
    knots = hb.decimate(rows, TOLERANCE_KM, hb.hermite)
    err = hb.max_error(rows, knots, hb.hermite)
    path = os.path.join(EPHEM, name + '.eph')
    hb.write_eph(path, naif_id, BARYCENTER, knots)
    print(f'{name}: {len(rows)} samples -> {len(knots)} knots, max error {err * 1e3:.1f} m, '
          f'{hb.calendar(rows[0][0])} .. {hb.calendar(rows[-1][0])} -> {path}')


def main():
    rows = [r for r in hb.fetch('-135', '@20065803', jd(IMPACT - 3 * hb.DAY), jd(SEAM - 2.0), 1) if r[0] < SEAM]
    t, p, v = rows[-1]
    end = tuple(p[i] + v[i] * (IMPACT - t) for i in range(3))
    target = dimorphos_center(IMPACT)
    shift = tuple(target[i] - end[i] for i in range(3))
    speed = math.sqrt(sum(x * x for x in v))
    along = sum(shift[i] * v[i] for i in range(3)) / speed
    print(f'DART: shifted {math.sqrt(sum(x * x for x in shift)):.3f} km ({along:.3f} km along the path, '
          f'{along / speed:.2f} s), from {hb.calendar(rows[0][0])}')
    rows = [(t, tuple(p[i] + shift[i] for i in range(3)), v) for t, p, v in rows]
    rows.append((IMPACT, target, v))
    write('dart_didymos', -135, rows)

    release = (hb.jd_from_text('2022-09-11 23:16') - hb.J2000_JD) * hb.DAY
    write('liciacube_didymos', -210, hb.fetch('-210', '@20065803', jd(release), jd(IMPACT + hb.DAY), 1))


if __name__ == '__main__':
    main()
