"""Bake ESA's Hera approach arc and the contemporaneous Didymos system.

Download the three SPKs below from https://spiftp.esac.esa.int/data/SPICE/HERA/kernels/spk/
into one directory. With SpiceyPy installed, run:
    python -E -P tools/dart/hera_bake.py <kernel-directory>

SPICE evaluates the type 19 Hera and type 13 asteroid kernels in J2000.
The historical DART tables remain separate from these approach tables.
"""

import argparse
import math
from pathlib import Path
import sys

import spiceypy as spice

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'bake'))
import horizons_bake as hb

HERA, PRIMARY, MOON, SUN = -91, -658030, -658031, 10
BARYCENTER = 2065803
MOON_FRACTION = 3.0268e-10 / 3.5521178e-8  # Same GMs as the scene.
HERA_SPK = 'hera_fsr_000008_261028_261107_v01.bsp'
ASTEROID_SPKS = ('didymos_fsp_000005_260701_270701_v01.bsp',
                 'dimorphos_fsp_000005_260701_270701_v01.bsp')


def state(target, center, t):
    values, _ = spice.spkgeo(target, t, 'J2000', center)
    return tuple(float(value) for value in values)


def sample(t):
    primary = state(PRIMARY, SUN, t)
    moon = state(MOON, PRIMARY, t)
    # ESA publishes the two bodies heliocentrically. Resolve the same hierarchy
    # as the DART scene using their mutual separation and the scene's mass ratio.
    reflex = tuple(-MOON_FRACTION * value for value in moon)
    barycenter = tuple(a - b for a, b in zip(primary, reflex))
    return barycenter, reflex, moon, state(HERA, PRIMARY, t)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kernels', type=Path)
    args = parser.parse_args()
    hera_kernel = str(args.kernels / HERA_SPK)
    for kernel in (*ASTEROID_SPKS, HERA_SPK):
        spice.furnsh(str(args.kernels / kernel))
    start, end = spice.spkcov(hera_kernel, HERA)
    regular = range(math.floor(start) + 1, math.ceil(end), 30)
    # The October 28 braking maneuver drops the relative speed from 2.83 to
    # 0.31 m/s in a few seconds. Sample the source densely across the impulse.
    burn = round((hb.jd_from_text('2026-10-28 02:01') - hb.J2000_JD) * hb.DAY)
    fine = (burn + .5 * i for i in range(241))
    times = sorted({start, end, *regular, *fine})
    samples = [sample(t) for t in times]

    tables = (
        ('didymos_hera', BARYCENTER, SUN, .010),
        ('didymos_bary_hera', PRIMARY, BARYCENTER, .001),
        ('dimorphos_hera', MOON, PRIMARY, .001),
        ('hera_didymos', HERA, PRIMARY, .002),
    )
    outputs = []
    for index, (name, target, center, tolerance) in enumerate(tables):
        rows = [(t, states[index][:3], states[index][3:]) for t, states in zip(times, samples, strict=True)]
        # Hold out every other sample to check the compressed ephemeris independently.
        knots = hb.decimate(rows[:-1:2] + [rows[-1]], tolerance * .75, hb.hermite)
        error = hb.max_error(rows, knots, hb.hermite)
        if error > tolerance:
            raise ValueError(f'{name}: {error * 1000:.3f} m exceeds {tolerance * 1000:.1f} m')
        outputs.append((name, target, center, knots, error))

    for name, target, center, knots, error in outputs:
        path = ROOT / 'assets' / 'ephem' / f'{name}.eph'
        hb.write_eph(str(path), target, center, knots)
        print(f'{name}: {len(knots)} knots, checked error {error * 1000:.3f} m')
    print(f'Coverage (TDB): {spice.etcal(start)} .. {spice.etcal(end)}')


if __name__ == '__main__':
    main()
