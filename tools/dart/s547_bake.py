"""Bake Didymos and Dimorphos from JPL's orbit solution 547 (Naidu et al. 2024) into
assets/ephem/didymos_bary.eph (the primary around the system barycenter) and
assets/ephem/dimorphos.eph (the satellite around the primary).

Usage:
    python s547_bake.py <dimorphos_s547.bsp>

The SPK is on figshare (Naidu & Chesley 2023, "Dimorphos orbit solution JPL 547",
doi:10.6084/m9.figshare.24256069, CC BY 4.0; see assets/ephem/SOURCES.md). It has
type 5 segments (discrete states every 1728 s, in ECLIPJ2000) of 120065803 Dimorphos
relative to 920065803 Didymos and of Didymos relative to 20065803, the system
barycenter, each split at the DART impact. Horizons serves the same solution, but
(2026-10) answers nothing for these bodies before 2026-07.

The SPK's own states are the samples (cubic Hermite between them is good to about a
centimetre on the 11.9 h orbit). At the impact the velocity jumps by 2.6 mm/s: the
states a second either side (from the segments before and after) are kept as knots,
so the scene reads the impulse from the table (unexplained_delta_v) and the cubic
between them is off the true path by under a millimetre. Knots are reduced to 0.5 m.

Only the Python standard library is used.
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'tools', 'bake'))
sys.path.insert(0, os.path.join(REPO, 'tools', 'philae'))
import horizons_bake as hb  # noqa: E402
import spk  # noqa: E402

OUT_DIR = os.path.join(REPO, 'assets', 'ephem')
DIMORPHOS, DIDYMOS, BARYCENTER = 120065803, 920065803, 20065803
ECLIPJ2000 = 17
OBLIQUITY = math.radians(84381.448 / 3600.0)  # NAIF's ECLIPJ2000 (IAU 1976)
# The reference GMs only shape the Kepler-relative interpolation (AXEPH2); the knots are
# exact samples. Dimorphos: n^2 a^3 of the circular orbit before the impact (the GM in
# the trailer of the SPK's Dimorphos segments), 1.7% above the solution's GM of the
# system (3.544e-8; Didymos' J2 speeds the orbit up), so the two-body arcs keep pace.
GM_ORBIT = 3.604713774874992e-08
# The primary's circle around the barycenter (~10 m) as a two-body orbit: n^2 r^3, the
# GM in the trailer of the SPK's Didymos segments.
GM_PRIMARY_REFLEX = 2.7299184854083606e-14
IMPACT = 717506133.3653613        # 2022-09-26 23:15:33.365 TDB, where the segments meet
SPAN = (hb.jd_from_text('2021-11-24'), hb.jd_from_text('2023-04-01'))
TOLERANCE_KM = 0.0005             # 0.5 m
MARGIN_S = 1.0                    # the knots either side of the impact


def to_icrf(v):
    c, s = math.cos(OBLIQUITY), math.sin(OBLIQUITY)
    return (v[0], c * v[1] - s * v[2], s * v[1] + c * v[2])


def segments(kernel, target, center):
    out = [s for s in kernel.segments if s.target == target and s.center == center]
    if not out or any(s.frame != ECLIPJ2000 or s.kind != 5 for s in out):
        sys.exit(f'{target} relative to {center}: expected type 5 segments in ECLIPJ2000')
    return sorted(out, key=lambda s: s.start())


def state(seg, t):
    """Cubic Hermite between the segment's states (as spk.Spk.state)."""
    one = spk.Spk.__new__(spk.Spk)
    one.segments = [seg]
    p, v = one.state(t)
    return to_icrf(p), to_icrf(v)


def rows_for(kernel, target, center):
    before, after = segments(kernel, target, center)
    if abs(before.end() - IMPACT) > 1e-3 or abs(after.start() - IMPACT) > 1e-3:
        sys.exit('the segments do not meet at the impact')
    t0 = (SPAN[0] - hb.J2000_JD) * hb.DAY
    t1 = (SPAN[1] - hb.J2000_JD) * hb.DAY
    rows = []
    for seg, lo, hi in ((before, t0, IMPACT - MARGIN_S), (after, IMPACT + MARGIN_S, t1)):
        times = [t for t in seg.epochs if lo < t < hi]
        for t in [lo] + times + [hi]:
            rows.append((t,) + state(seg, t))
    return rows


def bake(kernel, name, target, center, gm):
    rows = rows_for(kernel, target, center)
    interp = hb.make_interp(gm)
    keep = [(IMPACT - MARGIN_S, IMPACT + MARGIN_S)]
    knots = hb.decimate(rows, TOLERANCE_KM, interp, keep)
    err = hb.max_error(rows, knots, interp)
    path = os.path.join(OUT_DIR, name + '.eph')
    hb.write_eph(path, target, center, knots, gm)
    print(f'{name}: {len(rows)} samples -> {len(knots)} knots, max error {err * 1e3:.3f} m -> {path}')
    return rows


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    kernel = spk.Spk(sys.argv[1])
    bake(kernel, 'didymos_bary', DIDYMOS, BARYCENTER, GM_PRIMARY_REFLEX)
    rows = bake(kernel, 'dimorphos', DIMORPHOS, DIDYMOS, GM_ORBIT)

    # The impulse: the velocity either side, carried across the minutes by two-body
    # motion (the J2 of Didymos changes it by much less over two minutes).
    i = next(k for k, r in enumerate(rows) if r[0] > IMPACT)
    (ta, pa, va), (tb, pb, vb) = rows[i - 1], rows[i]
    _, va_at = hb.kepler_propagate(pa, va, GM_ORBIT, IMPACT - ta)
    _, vb_at = hb.kepler_propagate(pb, vb, GM_ORBIT, IMPACT - tb)
    dv = tuple(vb_at[k] - va_at[k] for k in range(3))
    r = math.sqrt(sum(x * x for x in pb))
    print(f'impulse at {hb.calendar(IMPACT)}: {math.sqrt(sum(x * x for x in dv)) * 1e6:.3f} mm/s, '
          f'{sum(dv[k] * va_at[k] for k in range(3)) / math.sqrt(sum(x * x for x in va_at)) * 1e6:.3f} mm/s '
          f'along the motion; Dimorphos {r * 1e3:.1f} m from Didymos')


if __name__ == '__main__':
    main()
