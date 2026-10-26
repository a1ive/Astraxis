"""Bake Proteus from a local NEP098 SPK excerpt (numpy and jplephem required).

See assets/ephem/SOURCES.md for download and reproduction commands.
The reference basis is Neptune's IAU pole (pck00011.tpc). Inclination and
node are fitted in that fixed basis; JPL's NEP098 mean-elements table uses
an equatorial frame and supplies no Laplace pole for Proteus.
"""

import sys

import numpy as np
from jplephem.spk import SPK

import element_table_bake as et


def state(kernel, t):
    """ICRF km and km/s, Proteus minus Neptune, at TDB seconds since J2000."""
    jd = et.hb.J2000_JD + t / et.hb.DAY
    result = np.zeros((2, 3, len(t)))
    for target, sign in ((808, 1), (899, -1)):
        covered = np.zeros(len(t), dtype=bool)
        for segment in kernel.segments:
            if segment.target != target or segment.center != 8:
                continue
            if segment.frame != 1:
                raise ValueError('expected ICRF segments')
            # jplephem excerpts can label each segment with the whole requested
            # span. Intersect it with the actual Chebyshev record coverage.
            first, step, coefficients = segment.load_array()
            last = first + step * coefficients.shape[1]
            mask = (jd >= max(segment.start_jd, first)) & (jd <= min(segment.end_jd, last)) & ~covered
            if mask.any():
                p, v = segment.compute_and_differentiate(jd[mask])
                result[0][:, mask] += sign * p
                result[1][:, mask] += sign * v / et.hb.DAY
                covered |= mask
        if not covered.all():
            raise ValueError(f'SPK does not cover target {target} at all requested times')
    return result[0].T, result[1].T


def main():
    if len(sys.argv) != 2:
        sys.exit('usage: python -E tools/bake/proteus_bake.py <nep098 excerpt.bsp>')
    et.CENTER, et.CENTER_ID = '500@899', 899
    et.POLE_RA_DEC = (299.36, 43.46)
    et.GM = 6835099.97  # Neptune alone, JPL SSD physical parameters (NEP097).
    start, stop = [(et.hb.jd_from_text(s) - et.hb.J2000_JD) * et.hb.DAY
                   for s in ('1972-01-01', '2100-01-01')]
    with SPK.open(sys.argv[1]) as kernel:
        t = np.arange(start, stop + 1, 3600.0)
        p, v = state(kernel, t)
        rot = et.reference_frame()
        print(f'{len(t)} hourly NEP098 samples', flush=True)
        series, rates = et.fit_segment(t, p @ rot.T, v @ rot.T)
        series = et.extend_endpoint(series, start, rates, True)
        series = et.extend_endpoint(series, stop, rates, False)
        knots = et.subset(series, et.thin(series, t, p @ rot.T, 10.0))
        # Fresh times at 4-hour intervals, offset by 37 seconds from the fit grid.
        check_t = np.arange(start + 37, stop, 14400.0)
        check_p, _ = state(kernel, check_t)
        error = np.linalg.norm(et.interpolate(knots, check_t) @ rot - check_p, axis=1)
        print(f'{len(check_t)} independent samples: max {error.max():.4f} km, '
              f'rms {np.sqrt(np.mean(error**2)):.4f} km', flush=True)
        if not np.isfinite(error).all() or error.max() > 10.0:
            raise RuntimeError('Proteus independent validation exceeds 10 km')
        et.write_elm('assets/ephem/proteus_nep.elm', 808, knots, rates)
        print(f'wrote {len(knots["t"])} knots', flush=True)
        for jd in (2447763.5, 2461041.5, 2488067.5):
            p, _ = state(kernel, np.array([(jd - et.hb.J2000_JD) * et.hb.DAY]))
            print(f'JD {jd}: {p[0].tolist()}')


if __name__ == '__main__':
    main()
