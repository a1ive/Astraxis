"""Check the Saturn small-moon tables against independently timed Horizons vectors.

Run from the repository root (numpy required):
    python -E -P tools/bake/verify_saturn_small_moons.py [name ...]

Samples start at 00:30 TDB and recur every 17 hours across the baked span. They
are cached alongside the hourly fit samples. Daphnis' source seams get separate
minute samples. This checks interpolation against the source ephemeris; it does
not measure the ephemeris' observational uncertainty. Prints C++ fixture rows
for tests/planet_scene_tests.cpp at three dates per moon.
"""

import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import element_table_bake as bake  # noqa: E402

NAMES = ('prometheus', 'pandora', 'pan', 'daphnis', 'telesto', 'calypso', 'helene', 'polydeuces')


def verify(name):
    spec = bake.TARGETS[name]
    path = os.path.join(bake.OUT, f'{name}_sat.elm')
    with open(path, 'rb') as f:
        raw = f.read()
    assert raw[:8] == b'AXELM1\0\0'
    target, center, count, _ = struct.unpack_from('<iiII', raw, 8)
    assert (target, center) == (spec['naif_id'], 699)
    assert struct.unpack_from('<2d', raw, 24) == bake.POLE_RA_DEC
    data = np.frombuffer(raw, '<f8', offset=56).reshape(count, 8)
    knots = dict(zip(('t', 'a', 'e', 'varpi', 'i', 'node', 'lambda', 'n'), data.T))
    start, stop = spec.get('start', bake.START), spec.get('stop', bake.STOP)
    t, p, _ = bake.samples(spec['command'], start + ' 00:30', stop, 17 * 60)
    inside = (t >= knots['t'][0]) & (t <= knots['t'][-1])
    t, p = t[inside], p[inside]
    rot = bake.reference_frame()
    error = np.linalg.norm(bake.interpolate(knots, t) @ rot - p, axis=1)
    limit = spec.get('tolerance_km', bake.TOLERANCE_KM) * 1.1
    print(f'{name}: {len(t)} withheld vectors, max {np.max(error):.4f} km, '
          f'rms {np.sqrt(np.mean(error**2)):.4f} km', flush=True)
    if not np.all(np.isfinite(error)) or np.max(error) > limit:
        raise RuntimeError(f'{name}: withheld error exceeds {limit} km')
    for date in (('1991-08-25', '2006-01-01', '2017-01-16') if name == 'daphnis' else
                 ('1981-08-25', '2006-01-01', '2047-10-17')):
        j = np.argmin(abs(t - (bake.hb.jd_from_text(date) - bake.hb.J2000_JD) * 86400))
        xyz = ', '.join(f'{x:.6f}' for x in p[j])
        print(f'    {{"{name.title()}", {t[j]:.6f}, {{{xyz}}}, {limit}}},')
    for date in spec.get('breaks', ()):
        jd = bake.hb.jd_from_text(date)
        rows = bake.hb.horizons_vectors(spec['command'], bake.CENTER, jd - 2/1440, jd + 2/1440, 1)
        tt = np.array([r[0] for r in rows])
        pp = np.array([r[1] for r in rows])
        err = np.linalg.norm(bake.interpolate(knots, tt) @ rot - pp, axis=1)
        print(f'    {date} seam: max {np.max(err):.4f} km', flush=True)
        if np.max(err) > limit:
            raise RuntimeError(f'{name}: seam error exceeds {limit} km')


if __name__ == '__main__':
    for name in sys.argv[1:] or NAMES:
        verify(name)
