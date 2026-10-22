"""Bake Philae's descent, bounces and resting place into assets/ephem/philae_67p.eph.

Usage:
    python philae_bake.py <directory with the ESA kernels>

The kernels come from the ESA SPICE Service
(https://spiftp.esac.esa.int/data/SPICE/ROSETTA/kernels/spk/):
  LORB_SUN_J2000_SDL_V1_1.BSP     separation, descent and first touchdown
                                  (SONC-FD), Philae relative to the Sun, J2000
  LORB_C_G_FIXED_RBD_7_V2_0.BSP   the bounces to the resting place (Grieger and
                                  the 7Dwarfs group), relative to 67P in its
                                  body-fixed frame 67P/C-G_CK

The table is in 67P's body-fixed frame (the scene's `body_fixed = true`), so that
Philae touches down and rests on the drawn shape whatever the comet's rotation:
  - descent: Philae minus 67P (assets/ephem/c67p.eph, the ESA comet solution
    from Horizons that the scene uses), turned into the body frame with the
    scene's rotation of 67P (its spin_table, fitted by c67p_spin.py to the
    attitude kernels of the same frame as the bounces). The small offset left
    from the 7Dwarfs touchdown point is taken up linearly over the descent
    (from nothing at the release to all of it at touchdown);
  - bounces: the 7Dwarfs states;
  - rest: the last 7Dwarfs position, held (zero velocity) until 2100.
Knots are reduced to 0.5 m (cubic Hermite, as at run time).

Only the Python standard library is used.
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'tools', 'bake'))
sys.path.insert(0, os.path.join(REPO, 'tools'))
sys.path.insert(0, HERE)
import c67p_spin  # noqa: E402
import horizons_bake as hb  # noqa: E402
import spk  # noqa: E402
from scene_templates import load_scene  # noqa: E402

LANDER_ID = -226800
COMET_ID = 1000012
DESCENT_STEP_S = 10.0
TOL_KM = 0.0005
HOLD_UNTIL_TDB = (2488069.5 - 2451545.0) * 86400.0  # 2100-01-01 TDB


def comet_rotation(scene_path):
    """The scene's spin table of 67P."""
    scene = load_scene(scene_path)
    body = next(b for b in scene['bodies'] if b['name'] == '67P')
    return c67p_spin.Table.load(os.path.join(REPO, 'assets', body['spin_table']))


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    kernels = sys.argv[1]
    sdl = spk.Spk(os.path.join(kernels, 'LORB_SUN_J2000_SDL_V1_1.BSP'))
    rbd = spk.Spk(os.path.join(kernels, 'LORB_C_G_FIXED_RBD_7_V2_0.BSP'))
    comet = hb.read_eph(os.path.join(REPO, 'assets', 'ephem', 'c67p.eph'))
    rotation = comet_rotation(os.path.join(REPO, 'assets', 'scenes', 'rosetta.toml'))
    for seg in rbd.segments:
        if seg.target != LANDER_ID or seg.center != COMET_ID:
            sys.exit('unexpected segment in the bounce kernel')

    def comet_state(t):
        p = hb.eval_eph(comet, t)
        a, b = hb.eval_eph(comet, t - 1.0), hb.eval_eph(comet, t + 1.0)
        return p, tuple((b[i] - a[i]) / 2.0 for i in range(3))

    def descent_body(t):
        lp, lv = sdl.state(t)
        cp, cv = comet_state(t)
        r = tuple(lp[i] - cp[i] for i in range(3))
        v = tuple(lv[i] - cv[i] - cross(rotation.spin(t), r)[i] for i in range(3))
        ax = rotation.axes(t)
        return tuple(dot(a, r) for a in ax), tuple(dot(a, v) for a in ax)

    release, _ = sdl.span()
    touchdown, rest = rbd.span()
    end_p, _ = descent_body(touchdown)
    td_p, _ = rbd.state(touchdown)
    offset = tuple(td_p[i] - end_p[i] for i in range(3))
    print('release %.3f, touchdown %.3f, rest %.3f (TDB s)' % (release, touchdown, rest))
    print('descent ends %.1f m from the 7Dwarfs touchdown point' % (1000.0 * math.sqrt(dot(offset, offset))))

    span = touchdown - release
    descent = []
    n = int(math.ceil(span / DESCENT_STEP_S))
    for k in range(n):
        t = release + span * k / n
        p, v = descent_body(t)
        f = (t - release) / span
        descent.append((t, tuple(p[i] + offset[i] * f for i in range(3)), tuple(v[i] + offset[i] / span for i in range(3))))
    bounce = [(t, tuple(s[0:3]), tuple(s[3:6])) for seg in rbd.segments for t, s in zip(seg.epochs, seg.states)]
    rest_p = bounce[-1][1]
    bounce[-1] = (bounce[-1][0], rest_p, (0.0, 0.0, 0.0))

    knots = hb.decimate(descent + bounce[:1], TOL_KM, hb.hermite)[:-1] + hb.decimate(bounce, TOL_KM, hb.hermite)
    knots.append((HOLD_UNTIL_TDB, rest_p, (0.0, 0.0, 0.0)))
    print('descent %d + bounce %d samples -> %d knots; error %.2f m / %.2f m' % (
        len(descent), len(bounce), len(knots),
        1000.0 * hb.max_error(descent, knots, hb.hermite), 1000.0 * hb.max_error(bounce, knots, hb.hermite)))
    print('rest at %s km (body frame), %.4f km from the center' % (
        tuple(round(x, 5) for x in rest_p), math.sqrt(dot(rest_p, rest_p))))
    out = os.path.join(REPO, 'assets', 'ephem', 'philae_67p.eph')
    hb.write_eph(out, LANDER_ID, COMET_ID, knots)
    print('wrote', out)


if __name__ == '__main__':
    main()
