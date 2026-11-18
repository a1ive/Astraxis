"""Bakes MSWIM2D's modelled solar wind at a planet into an hourly record like the Earth's
(assets/geomag/solar_wind.bin), for the planet's magnetopause and bow shock.

    python -I tools/geomag/mswim2d_bake.py <dir with the body's YYYY.csv> <dir with the l1 YYYY.csv>
        assets/geomag/solar_wind.bin <out.bin>

Input: MSWIM2D (Keebler et al. 2022, ApJS 260, 43), a 2-D MHD model of the solar wind in the
ecliptic out to 75 au, driven by OMNI (1985 to July 2004) and then by MIDL at L1 and STEREO,
from https://csem.engin.umich.edu/MSWIM2D/:
    precomputed_trajectories/chunks/<Body>/<year>.csv: the model sampled at the body every
        hour (UT, on the hour): datetime, r_AU, phi_deg (HGI), rho (amu/cm^3), ux, uy, uz
        (km/s), ...
    MSWIM2D_Data_New/Satellite_Data/chunks/l1/<year>.csv: the observations at L1 that drive
        it; only the last time is read, the end of the data.
and the Earth's record (OMNI2, solar_wind.bin), to tell where the model ran on OMNI's gaps.

Each hour of the output is the model at its middle, between the samples on either side
(none if they are more than 6 hours apart). Samples of all zeros, and single samples whose
speed is over 30% off both neighbours (a few in late 1985), are faults and dropped. Left out
(as no measurement):
  - the model's start: hours whose wind left 1 au within a solar rotation (27.2753 days) of
    the run's first time, while the model's starting state still filled the space between
    (some 45 days at Jupiter, 55 at Saturn, 85 at Uranus);
  - hours whose wind left 1 au, by corotation from the Earth's longitude, at a time with no
    OMNI pressure within a day: MSWIM2D bridges OMNI's gaps with straight lines (only
    before July 2004, where OMNI drives it; about a fifth of 1985-94 at Jupiter);
  - everything from the first hour whose wind left 1 au after the end of the L1 data: there
    the model is a forecast, the last solar rotation repeated.
The wind at hour t left 1 au at s = t - (r - 1 au) / v (radial at the local speed) at the
planet's longitude, and the Sun turns that longitude to the Earth's (27.2753 days a turn,
synodic) after ((lon_Earth - lon_planet) mod 360) / 360 of a turn; the Earth's record is
read at s plus that.

Output (little-endian), version 2 of the Earth's format:
    char[4]  "AXSW"
    uint16   version (2)
    uint16   first year (the record starts on its January 1 at 00:00 UT)
    uint32   count of hours
    float32  pressure step, nPa: the power of ten that holds the largest pressure in 65534
             steps (1e-4 at Jupiter, 1e-5 at Saturn and Uranus)
    per hour:
      uint16 dynamic pressure rho v^2 in pressure steps (at least 1); 65535 = none
      int8   Bz: always -128 = none (the model's field is not in the planet's axes)
      uint8  speed, 150 + 8 n km/s (n up to 254); 255 = none

Standard library only.
"""

import bisect
import csv
import datetime
import hashlib
import math
import os
import re
import struct
import sys

AU_KM = 1.495978707e8
AMU_KG = 1.66053907e-27
NPA_PER_AMU_CM3_KM2_S2 = AMU_KG * 1e6 * 1e6 * 1e9  # rho v^2: amu/cm^3 (km/s)^2 -> nPa
NO_PRESSURE, NO_BZ, NO_SPEED = 65535, -128, 255
MAX_SAMPLE_GAP = datetime.timedelta(hours=6)
SYNODIC_DAYS = 27.2753
SPIN_UP = datetime.timedelta(days=SYNODIC_DAYS)
OMNI_DRIVEN_UNTIL = datetime.datetime(2004, 7, 1)
OMNI_REACH_HOURS = 24
HGI_NODE_DEG = 75.76  # the solar equator's ascending node on the J2000 ecliptic


def earth_hgi_longitude(t):
    """The Earth's heliocentric longitude in HGI, degrees (the Sun's low-precision longitude,
    Astronomical Almanac; about a degree)."""
    d = (t - datetime.datetime(2000, 1, 1, 12)).total_seconds() / 86400.0
    mean = 280.460 + 0.9856474 * d
    anomaly = math.radians(357.528 + 0.9856003 * d)
    sun = mean + 1.915 * math.sin(anomaly) + 0.020 * math.sin(2.0 * anomaly)
    # From the equinox of date to J2000, then from the equinox to the node.
    return (sun + 180.0 - 1.397 * d / 36525.0 - HGI_NODE_DEG) % 360.0


def read_years(directory):
    years = sorted(int(m.group(1)) for name in os.listdir(directory)
                   if (m := re.fullmatch(r'(\d{4})\.csv', name)))
    if not years or years != list(range(years[0], years[-1] + 1)):
        sys.exit(f'expected consecutive YYYY.csv files in {directory}, found {years}')
    return years


def read_model(directory):
    """(time, r_AU, phi_deg, pressure nPa, speed km/s) for each valid sample, in order."""
    samples = []
    for year in read_years(directory):
        with open(os.path.join(directory, f'{year}.csv'), newline='', encoding='ascii') as f:
            rows = csv.reader(f)
            header = next(rows)
            if header[:7] != ['datetime', 'r_AU', 'phi_deg', 'rho', 'ux', 'uy', 'uz']:
                sys.exit(f'{year}.csv: unexpected columns {header[:7]}')
            for row in rows:
                t = datetime.datetime.fromisoformat(row[0])
                r, phi, rho = float(row[1]), float(row[2]), float(row[3])
                speed = math.sqrt(sum(float(x) ** 2 for x in row[4:7]))
                if samples and t <= samples[-1][0]:
                    sys.exit(f'{year}.csv: {t} out of order')
                if rho > 0.0 and speed > 0.0:
                    samples.append((t, r, phi, rho * speed * speed * NPA_PER_AMU_CM3_KM2_S2, speed))
    # Single samples far off both neighbours that agree with each other are faults in the
    # files, not the wind (a shock stays): in late 1985, at the same hours at every planet,
    # some 1000 km/s between two of 400.
    def spike(k):
        a, m, c = samples[k - 1][4], samples[k][4], samples[k + 1][4]
        return abs(m - a) > 0.3 * a and abs(m - c) > 0.3 * c and abs(c - a) < 0.1 * a
    faults = {k for k in range(1, len(samples) - 1) if spike(k)}
    for k in sorted(faults):
        print(f'left out a fault: {samples[k][0]}, {samples[k][4]:.0f} km/s between '
              f'{samples[k - 1][4]:.0f} and {samples[k + 1][4]:.0f}')
    return [s for k, s in enumerate(samples) if k not in faults]


def data_end(l1_dir):
    last = read_years(l1_dir)[-1]
    with open(os.path.join(l1_dir, f'{last}.csv'), newline='', encoding='ascii') as f:
        rows = list(csv.reader(f))
    return datetime.datetime.fromisoformat(rows[-1][0]) + datetime.timedelta(hours=1)


def read_omni(path):
    """The start of the Earth's record and, per hour, whether it has a pressure."""
    with open(path, 'rb') as f:
        data = f.read()
    if data[:4] != b'AXSW':
        sys.exit(f'{path}: not a solar wind record')
    version, first_year, count = struct.unpack_from('<HHI', data, 4)
    if version != 1:
        sys.exit(f'{path}: expected the Earth\'s record (version 1)')
    measured = [struct.unpack_from('<H', data, 12 + 4 * i)[0] != NO_PRESSURE for i in range(count)]
    return datetime.datetime(first_year, 1, 1), measured


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    model_dir, l1_dir, omni_path, out_path = sys.argv[1:]
    samples = read_model(model_dir)
    end_of_data = data_end(l1_dir)
    omni_start, omni = read_omni(omni_path)
    # Running count of OMNI's measured hours, for the share within a span.
    omni_sum = [0]
    for m in omni:
        omni_sum.append(omni_sum[-1] + m)

    def omni_has(t):
        i = math.floor((t - omni_start).total_seconds() / 3600.0)
        lo, hi = max(0, i - OMNI_REACH_HOURS), min(len(omni), i + OMNI_REACH_HOURS)
        return hi > lo and omni_sum[hi] > omni_sum[lo]

    times = [s[0] for s in samples]
    run_start = times[0]
    first_year = run_start.year
    start = datetime.datetime(first_year, 1, 1)
    values = []  # (pressure nPa, speed km/s) or None, hour by hour
    counts = {'spin_up': 0, 'sampling': 0, 'omni_gap': 0}
    gaps_by_year = {}
    t = start
    while True:
        middle = t + datetime.timedelta(minutes=30)
        t += datetime.timedelta(hours=1)
        j = bisect.bisect_right(times, middle)
        if j >= len(samples):
            break
        if j == 0:
            counts['spin_up'] += 1
            values.append(None)
            continue
        a, b = samples[j - 1], samples[j]
        if b[0] - a[0] > MAX_SAMPLE_GAP:
            counts['sampling'] += 1
            values.append(None)
            continue
        f = (middle - a[0]) / (b[0] - a[0])
        r = a[1] + f * (b[1] - a[1])
        phi = a[2] + f * ((b[2] - a[2] + 180.0) % 360.0 - 180.0)
        pressure = a[3] + f * (b[3] - a[3])
        speed = a[4] + f * (b[4] - a[4])
        left = middle - datetime.timedelta(seconds=(r - 1.0) * AU_KM / speed)
        if left < run_start + SPIN_UP:
            counts['spin_up'] += 1
            values.append(None)
            continue
        turn = ((earth_hgi_longitude(left) - phi) % 360.0) / 360.0
        at_earth = left + datetime.timedelta(days=turn * SYNODIC_DAYS)
        if at_earth > end_of_data:
            break  # a forecast from here on
        if at_earth < OMNI_DRIVEN_UNTIL and not omni_has(at_earth):
            counts['omni_gap'] += 1
            gaps_by_year[middle.year] = gaps_by_year.get(middle.year, 0) + 1
            values.append(None)
            continue
        values.append((pressure, speed))
    while values and values[-1] is None:
        values.pop()

    largest = max(v[0] for v in values if v)
    step = 10.0 ** math.ceil(math.log10(largest / (NO_PRESSURE - 1)))
    hours = []
    for v in values:
        if v is None:
            hours.append((NO_PRESSURE, NO_BZ, NO_SPEED))
        else:
            p = max(1, min(NO_PRESSURE - 1, round(v[0] / step)))
            s = max(0, min(254, round((v[1] - 150.0) / 8.0)))
            hours.append((p, NO_BZ, s))

    payload = b''.join(struct.pack('<HbB', *h) for h in hours)
    with open(out_path, 'wb') as f:
        f.write(b'AXSW')
        f.write(struct.pack('<HHIf', 2, first_year, len(hours), step))
        f.write(payload)
    end = start + datetime.timedelta(hours=len(hours))
    print(f'{out_path}: {len(hours)} hours, {start:%Y-%m-%d} .. {end:%Y-%m-%d %H:%M} UT '
          f'(L1 data end {end_of_data:%Y-%m-%d %H:%M}); pressure step {step:g} nPa')
    print(f'left out: {counts["spin_up"]} hours of spin-up, {counts["sampling"]} between sparse samples, '
          f'{counts["omni_gap"]} on OMNI\'s gaps')
    for year in sorted(gaps_by_year):
        print(f'  {year}: {gaps_by_year[year]} hours on OMNI\'s gaps')
    measured = sorted(h[0] * step for h in hours if h[0] != NO_PRESSURE)
    speeds = sorted(150 + 8 * h[2] for h in hours if h[2] != NO_SPEED)
    print(f'pressure: median {measured[len(measured) // 2]:.5f} nPa, 5% {measured[len(measured) // 20]:.5f}, '
          f'95% {measured[19 * len(measured) // 20]:.5f}, max {measured[-1]:.5f}; speed median {speeds[len(speeds) // 2]} km/s')
    with open(out_path, 'rb') as f:
        print(f'SHA-256 {hashlib.sha256(f.read()).hexdigest()}')


if __name__ == '__main__':
    main()
