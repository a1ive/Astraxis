"""Bakes the hourly solar wind at the Earth into assets/geomag/solar_wind.bin (the Earth's
magnetopause and bow shock).

    python -I tools/geomag/solar_wind_bake.py <dir with omni2_YYYY.dat> assets/geomag/solar_wind.bin

Input: the yearly files of the hourly OMNI2 set (omni2_1963.dat .. ), from
https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/ (NASA GSFC SPDF; King & Papitashvili
2005, JGR 110, A02104; format in omni2.text there). One line per hour, words separated by
blanks: 1 year, 2 day of the year, 3 hour (the average over that hour UT, shifted to the
Earth's bow shock), 17 Bz in GSM (nT; 999.9 = none), 25 flow speed (km/s; 9999. = none),
29 flow pressure (nPa; 99.99 = none).

Output (little-endian):
    char[4]  "AXSW"
    uint16   version (1)
    uint16   first year (the record starts on its January 1 at 00:00 UT)
    uint32   count of hours
    per hour:
      uint16 flow pressure in 0.01 nPa; 65535 = none
      int8   Bz (GSM) in 0.5 nT, held within -63.5 .. 63.5; -128 = none
      uint8  speed, 150 + 8 n km/s (n up to 254); 255 = none

Standard library only.
"""

import datetime
import os
import re
import struct
import sys

NO_PRESSURE, NO_BZ, NO_SPEED = 65535, -128, 255


def encode(words):
    bz, speed, pressure = float(words[16]), float(words[24]), float(words[28])
    p = NO_PRESSURE if pressure >= 99.985 else min(round(pressure * 100), NO_PRESSURE - 1)
    b = NO_BZ if abs(bz) >= 999.85 else max(-127, min(127, round(bz * 2)))
    s = NO_SPEED if speed >= 9998.5 else max(0, min(254, round((speed - 150) / 8)))
    return p, b, s


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src_dir, out_path = sys.argv[1], sys.argv[2]
    years = sorted(int(m.group(1)) for name in os.listdir(src_dir)
                   if (m := re.fullmatch(r'omni2_(\d{4})\.dat', name)))
    if not years or years != list(range(years[0], years[-1] + 1)):
        sys.exit(f'expected consecutive omni2_YYYY.dat files in {src_dir}, found {years}')

    hours = []
    for year in years:
        expected = datetime.datetime(year, 1, 1)
        with open(os.path.join(src_dir, f'omni2_{year}.dat'), encoding='ascii') as f:
            for line in f:
                words = line.split()
                if len(words) < 29:
                    continue
                at = datetime.datetime(int(words[0]), 1, 1) + datetime.timedelta(days=int(words[1]) - 1, hours=int(words[2]))
                if at != expected:
                    sys.exit(f'omni2_{year}.dat: expected {expected}, found {at}')
                hours.append(encode(words))
                expected += datetime.timedelta(hours=1)
        if expected != datetime.datetime(year + 1, 1, 1):
            sys.exit(f'omni2_{year}.dat ends at {expected}')
    # The current year is filled with blanks to its end.
    while hours and hours[-1] == (NO_PRESSURE, NO_BZ, NO_SPEED):
        hours.pop()

    with open(out_path, 'wb') as f:
        f.write(b'AXSW')
        f.write(struct.pack('<HHI', 1, years[0], len(hours)))
        f.write(b''.join(struct.pack('<HbB', *h) for h in hours))
    start = datetime.datetime(years[0], 1, 1)
    end = start + datetime.timedelta(hours=len(hours))
    print(f'{out_path}: {len(hours)} hours, {start:%Y-%m-%d} .. {end:%Y-%m-%d %H:%M} UT')

    def median(values):
        values = sorted(values)
        return values[len(values) // 2]
    pressures = [h[0] / 100 for h in hours if h[0] != NO_PRESSURE]
    speeds = [150 + 8 * h[2] for h in hours if h[2] != NO_SPEED]
    print(f'  medians: {median(pressures):.2f} nPa, {median(speeds)} km/s; largest pressure {max(pressures):.2f} nPa, '
          f'Bz held at the limit in {sum(abs(h[1]) == 127 for h in hours)} hours')
    print('  hours with a pressure / with Bz, by decade:')
    for decade in range(years[0] // 10 * 10, years[-1] + 1, 10):
        first = max((datetime.datetime(decade, 1, 1) - start) // datetime.timedelta(hours=1), 0)
        last = min((datetime.datetime(decade + 10, 1, 1) - start) // datetime.timedelta(hours=1), len(hours))
        part = hours[first:last]
        if part:
            print(f'    {decade}s: {sum(h[0] != NO_PRESSURE for h in part) / len(part):5.1%} / '
                  f'{sum(h[1] != NO_BZ for h in part) / len(part):5.1%}')


if __name__ == '__main__':
    main()
