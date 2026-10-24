"""Bracket Spitzer's monthly Horizons navigation seams in its baked table.

Run after horizons_bake.py tools/bake/ephem.toml spitzer_geo:
    python tools/bake/spitzer_seams.py

Horizons -79 has position discontinuities at some month boundaries during the
tracking era. A coarse spline spreads them over hours. Fetch states one second
on either side of every boundary, using TLIST (ICRF, geometric, TDB), and insert
them into spitzer_geo.eph. The two-second bridge represents a source-data seam,
not a spacecraft maneuver. See assets/ephem/SOURCES.md.
"""

import math
from pathlib import Path
import struct
import time
import urllib.parse
import urllib.request

import horizons_bake as h


def vectors_at(dates):
    params = {
        'format': 'text', 'COMMAND': "'-79'", 'CENTER': "'500@399'",
        'OBJ_DATA': "'NO'", 'MAKE_EPHEM': "'YES'", 'EPHEM_TYPE': "'VECTORS'",
        'TLIST': ' '.join(f"'{jd:.12f}'" for jd in dates),
        'TLIST_TYPE': "'JD'", 'TIME_TYPE': "'TDB'", 'VEC_TABLE': "'2'",
        'REF_SYSTEM': "'ICRF'", 'REF_PLANE': "'FRAME'", 'OUT_UNITS': "'KM-S'",
        'VEC_CORR': "'NONE'", 'CSV_FORMAT': "'YES'",
    }
    url = h.API + '?' + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                text = response.read().decode('utf-8')
            break
        except OSError:
            if attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))
    if '$$SOE' not in text:
        raise RuntimeError(text[-2000:])
    rows = []
    for line in text.split('$$SOE', 1)[1].split('$$EOE', 1)[0].strip().splitlines():
        f = line.split(',')
        rows.append(((float(f[0]) - h.J2000_JD) * h.DAY,
                     tuple(float(x) for x in f[2:5]), tuple(float(x) for x in f[5:8])))
    if len(rows) != len(dates):
        raise RuntimeError('Horizons returned an unexpected number of seam samples')
    for row, jd in zip(rows, dates):
        if abs(row[0] - (jd - h.J2000_JD) * h.DAY) > 0.001 or not all(
                math.isfinite(x) for part in row[1:] for x in part):
            raise RuntimeError('Invalid seam state or timestamp')
    return rows


def main():
    path = Path(__file__).resolve().parents[2] / 'assets/ephem/spitzer_geo.eph'
    header = path.read_bytes()[:24]
    if header[:8] != b'AXEPH1\0\0' or struct.unpack_from('<ii', header, 8) != (-79, 399):
        raise ValueError('Expected the AXEPH1 Spitzer geocentric table')
    table = h.read_eph(path)
    dates = []
    for year in range(2003, 2021):
        for month in range(1, 13):
            if (year, month) < (2003, 9) or (year, month) > (2020, 2):
                continue
            jd = h.jd_from_text(f'{year:04}-{month:02}-01')
            dates.extend([jd - 1.0 / h.DAY, jd + 1.0 / h.DAY])
    rows = []
    for begin in range(0, len(dates), 16):
        rows.extend(vectors_at(dates[begin:begin + 16]))
        print(f'Fetched {len(rows)}/{len(dates)} seam states', flush=True)
    if rows[0][0] <= table.start() or rows[-1][0] >= table.end():
        raise ValueError('Seam samples are outside the baked span')
    # Remove any old knot within a bridge, including a previous run's anchors.
    windows = [(rows[i][0] - 0.001, rows[i + 1][0] + 0.001) for i in range(0, len(rows), 2)]
    kept = [row for row in table.knots if not h.in_windows(row[0], windows)]
    jumps = []
    for a, b in zip(rows[::2], rows[1::2]):
        dt = b[0] - a[0]
        smooth = tuple(a[1][j] + 0.5 * (a[2][j] + b[2][j]) * dt for j in range(3))
        jump = h.dist(smooth, b[1])
        if jump > 1.0:
            jumps.append(jump)
    h.write_eph(path, -79, 399, sorted(kept + rows), max_gap=table.max_gap)
    print(f'{len(windows)} month boundaries, {len(jumps)} jumps > 1 km; '
          f'largest {max(jumps, default=0.0):.3f} km; {len(kept) + len(rows)} knots')


if __name__ == '__main__':
    main()
