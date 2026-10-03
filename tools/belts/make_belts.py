"""Fetch orbital elements of the main asteroid belt and the trans-Neptunian objects
from the JPL Small-Body Database and write them as point clouds for assets/belts/.

Usage (from the repository root):
    python tools/belts/make_belts.py [output_dir]

Queries the SBDB Query API (https://ssd-api.jpl.nasa.gov/sbdb_query.api):
  main_belt.bin    classes MBA, IMB, OMB (main belt, inner, outer) with H < 15
  kuiper_belt.bin  class TNO (all trans-Neptunian objects)
Jupiter Trojans (TJN) are not included.

The elements are heliocentric, osculating, J2000 ecliptic, each at its own epoch;
the mean anomaly is moved to a common epoch with the two-body mean motion
(GM of the Sun, DE440), so that the renderer can propagate every object from one
epoch. Objects are sorted by absolute magnitude H (brightest first).

Output format (little endian), see src/scene/scene_loader.cpp:
    char[8]  magic "AXBELT1\\0"
    uint32   object count
    uint32   reserved (0)
    float64  epoch (JD TDB)
    objects: float32 a (au), e, i, node, arg_peri, mean_anomaly (rad, at the epoch), H

Only the Python standard library is used.
"""

import json
import math
import os
import struct
import sys
import time
import urllib.parse
import urllib.request

API = 'https://ssd-api.jpl.nasa.gov/sbdb_query.api'
FIELDS = 'spkid,a,e,i,om,w,ma,epoch,H'
PAGE = 20000
EPOCH_JD = 2461200.5  # current SBDB standard epoch (2026-06-09)
GM_SUN = 1.32712440041279419e11  # km^3/s^2, DE440 (JPL SSD astrodynamic parameters)
AU_KM = 149597870.7
DAY = 86400.0

BELTS = [
    ('main_belt.bin', 'MBA,IMB,OMB', {'AND': ['H|LT|15']}),
    ('kuiper_belt.bin', 'TNO', None),
]


def query(sb_class, cdata, offset):
    params = {'fields': FIELDS, 'sb-class': sb_class, 'limit': str(PAGE), 'limit-from': str(offset)}
    if cdata:
        params['sb-cdata'] = json.dumps(cdata)
    url = API + '?' + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=300) as resp:
                return json.load(resp)
        except OSError as exc:
            print(f'  retry {attempt + 1}: {exc}')
            time.sleep(5 * (attempt + 1))
    sys.exit(f'query failed: {url}')


def fetch(sb_class, cdata):
    rows, offset = [], 0
    while True:
        d = query(sb_class, cdata, offset)
        fields = d['fields']
        page = [dict(zip(fields, r)) for r in d.get('data', [])]
        rows.extend(page)
        print(f'  {sb_class}: {len(rows)} / {d["count"]}')
        if len(page) < PAGE:
            return rows
        offset += PAGE


def mean_motion_rad_per_day(a_au):
    a_km = a_au * AU_KM
    return math.sqrt(GM_SUN / a_km ** 3) * DAY


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', '..', 'assets', 'belts')
    os.makedirs(out_dir, exist_ok=True)
    for name, sb_class, cdata in BELTS:
        rows = fetch(sb_class, cdata)
        objects, skipped = [], 0
        for r in rows:
            try:
                a, e = float(r['a']), float(r['e'])
                i, om, w, ma = (math.radians(float(r[k])) for k in ('i', 'om', 'w', 'ma'))
                epoch, h = float(r['epoch']), float(r['H'])
            except (TypeError, ValueError):
                skipped += 1  # missing elements or H
                continue
            if not (a > 0.0 and 0.0 <= e < 1.0):
                skipped += 1
                continue
            m = math.remainder(ma + mean_motion_rad_per_day(a) * (EPOCH_JD - epoch), 2.0 * math.pi) % (2.0 * math.pi)
            objects.append((h, a, e, i, om, w, m))
        objects.sort()
        path = os.path.normpath(os.path.join(out_dir, name))
        with open(path, 'wb') as f:
            f.write(b'AXBELT1\0')
            f.write(struct.pack('<IId', len(objects), 0, EPOCH_JD))
            for h, a, e, i, om, w, m in objects:
                f.write(struct.pack('<7f', a, e, i, om, w, m, h))
        print(f'{name}: {len(objects)} objects ({skipped} skipped) -> {path} ({os.path.getsize(path) // 1024} KB)')


if __name__ == '__main__':
    main()
