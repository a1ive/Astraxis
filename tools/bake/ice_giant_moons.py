"""Bake Neptune and Uranus moon orbits from JPL SPK excerpts.

See assets/ephem/SOURCES.md for the excerpts and reproduction commands.
"""

import argparse
from functools import partial
from pathlib import Path

import numpy as np
from jplephem.spk import SPK

import element_table_bake as et
from proteus_bake import barycentric_state

hb = et.hb
SYSTEMS = {
    'nep': (8, 899, 6835099.97, (299.36, 43.46), 'nep098', 'nep105',
            {'naiad': 803, 'thalassa': 804, 'despina': 805, 'galatea': 806,
             'larissa': 807, 'hippocamp': 814},
            {'nereid': (802, 25.0)}),
    'ura': (7, 799, 5793951.3, (77.311, 15.175), 'ura184', 'ura117',
            {'cordelia': 706, 'ophelia': 707, 'portia': 712, 'belinda': 714,
             'puck': 715, 'mab': 726},
            {'caliban': (716, 15.0), 'sycorax': (717, 25.0)}),
}


def relative_state(t, *, kernel, planet_kernel, barycenter, planet, target):
    return (barycentric_state(kernel, t, target, barycenter)
            - barycentric_state(planet_kernel, t, planet, barycenter))


def validate(name, predicted, actual, tolerance):
    error = np.linalg.norm(predicted - actual, axis=1)
    print(f'{name}: {len(error)} check samples, max {error.max():.3f} km, '
          f'RMS {np.sqrt(np.mean(error**2)):.3f} km')
    if not np.isfinite(error).all() or error.max() > tolerance:
        raise ValueError(f'{name} exceeds {tolerance} km')


def bake_inner(path, target, state, start, stop, tolerance):
    t = np.arange(start, stop + 1, 3600.0)
    p, v = state(t)
    rot = et.reference_frame()
    plane_pos = p @ rot.T
    series, rates = et.fit_segment(t, plane_pos, v @ rot.T)
    series = et.extend_endpoint(series, start, rates, True)
    series = et.extend_endpoint(series, stop, rates, False)
    knots = et.subset(series, et.thin(series, t, plane_pos, tolerance * 0.8))
    validate(path.stem, et.interpolate(knots, t) @ rot, p, tolerance)
    check = np.arange(start + 1837.0, stop, 10800.0)
    validate(path.stem, et.interpolate(knots, check) @ rot, state(check)[0], tolerance)
    et.write_elm(path, target, knots, rates)
    print(f'wrote {path} ({len(knots["t"])} knots)')


def bake_outer(path, target, state, start, stop, tolerance):
    step = 8 * 3600.0
    t = np.r_[np.arange(start, stop, step), stop]
    p, v = state(t)
    rows = list(zip(t, p, v))
    interp = hb.make_interp(et.GM)
    knots = hb.decimate(rows, tolerance * 0.8, interp)
    check = np.arange(start + 907.0, stop, 3 * 3600.0)
    table = hb.Table(knots, et.GM)
    predicted = np.array([hb.eval_eph(table, t) for t in check])
    validate(path.stem, predicted, state(check)[0], tolerance)
    hb.write_eph(path, target, et.CENTER_ID, knots, et.GM)
    print(f'wrote {path} ({len(knots)} knots)')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kernel_dir', type=Path)
    parser.add_argument('system', choices=SYSTEMS)
    parser.add_argument('moons', nargs='*')
    args = parser.parse_args()

    barycenter, et.CENTER_ID, et.GM, et.POLE_RA_DEC, inner_spk, outer_spk, inner, outer = SYSTEMS[args.system]
    moons = args.moons or list(inner) + list(outer)
    if unknown := set(moons) - (inner.keys() | outer.keys()):
        parser.error(f'unknown moons: {", ".join(sorted(unknown))}')
    start, stop = [(hb.jd_from_text(date) - hb.J2000_JD) * hb.DAY
                   for date in ('1972-01-01', '2100-01-01')]
    et.STEP_DAYS = 2.0
    with (SPK.open(str(args.kernel_dir / f'{inner_spk}.bsp')) as planet_kernel,
          SPK.open(str(args.kernel_dir / f'{outer_spk}.bsp')) as outer_kernel):
        for name in moons:
            path = Path(et.OUT) / f'{name}_{args.system}'
            print(name, flush=True)
            kernel = planet_kernel if name in inner else outer_kernel
            target = inner[name] if name in inner else outer[name][0]

            state = partial(relative_state, kernel=kernel, planet_kernel=planet_kernel,
                            barycenter=barycenter, planet=et.CENTER_ID, target=target)
            if name in inner:
                tolerance = 15.0 if name in ('naiad', 'hippocamp') else 10.0
                bake_inner(path.with_suffix('.elm'), target, state, start, stop, tolerance)
            else:
                bake_outer(path.with_suffix('.eph'), target, state, start, stop, outer[name][1])


if __name__ == '__main__':
    main()
