"""Bakes the planetary Kp index into assets/geomag/kp.bin (the aurora's activity).

    python kp_bake.py <dir with kpYYYY.wdc> <out kp.bin>

Input: GFZ's yearly files in WDC format (kp1932.wdc .. ), from
https://ftp.gwdg.de/pub/geophys/kp-ap/wdc/yearly/ (a mirror of
ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/wdc/yearly/; format in wdc_fmt.txt
there). Columns 13-28 hold the day's eight 3-hourly Kp values (0-3, 3-6, ...
21-24 UT) as two-digit codes in tenths: 0, 3, 7, 10, 13, 17, ... for
0o, 0+, 1-, 1o, 1+, 2-, ... (thirds of a unit). Only Kp is kept: the files'
sunspot numbers have another licence (SILSO, CC BY-NC).

Output (little-endian):
    char[4]  "AXKP"
    uint16   version (1)
    uint16   first year (the record starts on its January 1 at 00:00 UT)
    uint32   count of 3-hour intervals
    uint8    Kp in thirds of a unit (0 .. 27) for each interval; 255 = missing

Standard library only.
"""

import datetime
import os
import re
import struct
import sys

MISSING = 255


def read_year(path):
    """{date: [8 Kp codes]} from one yearly WDC file."""
    days = {}
    with open(path, encoding='ascii') as f:
        for line in f:
            if len(line.rstrip('\n')) < 28 or not line[:6].strip():
                continue
            yy, mm, dd = int(line[0:2]), int(line[2:4]), int(line[4:6])
            codes = [int(line[12 + 2 * k:14 + 2 * k]) for k in range(8)]
            days[(yy, mm, dd)] = codes
    return days


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src_dir, out_path = sys.argv[1], sys.argv[2]
    years = sorted(int(m.group(1)) for name in os.listdir(src_dir)
                   if (m := re.fullmatch(r'kp(\d{4})\.wdc', name)))
    if not years or years != list(range(years[0], years[-1] + 1)):
        sys.exit(f'expected consecutive kpYYYY.wdc files in {src_dir}, found {years}')

    values = []
    last_valid = None
    for year in years:
        days = read_year(os.path.join(src_dir, f'kp{year}.wdc'))
        date = datetime.date(year, 1, 1)
        while date.year == year:
            codes = days.get((year % 100, date.month, date.day))
            if codes is None:
                values.extend([MISSING] * 8)
            else:
                for code in codes:
                    thirds = round(code * 3 / 10)
                    if not 0 <= thirds <= 27 or abs(thirds * 10 / 3 - code) > 0.5:
                        sys.exit(f'{year}-{date.month:02}-{date.day:02}: bad Kp code {code}')
                    values.append(thirds)
                last_valid = len(values)
            date += datetime.timedelta(days=1)
    values = values[:last_valid]  # the current year ends at the last day in the file
    missing = sum(v == MISSING for v in values)

    with open(out_path, 'wb') as f:
        f.write(b'AXKP')
        f.write(struct.pack('<HHI', 1, years[0], len(values)))
        f.write(bytes(values))
    end = datetime.datetime(years[0], 1, 1) + datetime.timedelta(hours=3 * len(values))
    print(f'{out_path}: {len(values)} intervals, {years[0]}-01-01 .. {end:%Y-%m-%d %H:%M} UT, {missing} missing')


if __name__ == '__main__':
    main()
