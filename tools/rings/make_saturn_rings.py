"""Builds Saturn's ring optical-depth profile (assets/rings/saturn_rss.ring)
from a Cassini Radio Science ring occultation (PDS volume CORSS_8001).

Usage:
    python make_saturn_rings.py <output.ring> [samples]

Source (see assets/rings/SOURCES.md): Rev 7 egress, X band (3.6 cm), DSN 43,
2005-05-03, 10 km resolution, sampled every 2.5 km:
    RSS_2005_123_X43_E_TAU_10KM.TAB
Columns used: 1 ring radius, 2 and 3 radius corrections (improved pole,
timing offset), 6 normalized signal power, 7 normal optical depth, 9 normal
optical depth threshold. Where the signal is lost (power <= 0) or the optical
depth reaches the threshold (~5.06, the B ring core), the threshold is used:
the ring is opaque there either way. Negative noise in empty gaps becomes 0.
The profile is averaged into `samples` equal bins (default 8192).

Output format (little endian), read by src/scene/scene_loader.cpp:
    char[8]  "AXRING1\\0"
    float64  inner radius, outer radius (km)
    uint32   sample count, uint32 reserved (0)
    float32  normal optical depth per bin, from inner to outer

Only the Python standard library is used.
"""

import os
import struct
import sys
import urllib.request

URL = ('https://pds-rings.seti.org/holdings/volumes/CORSS_8xxx/CORSS_8001/data/Rev007/Rev007E/'
       'Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_TAU_10KM.TAB')


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    out = sys.argv[1]
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 8192

    print(f'fetching {URL}')
    request = urllib.request.Request(URL, headers={'User-Agent': 'astraxis-tools'})  # the default UA is refused
    with urllib.request.urlopen(request, timeout=300) as resp:
        text = resp.read().decode('ascii')

    rows = []
    for line in text.splitlines():
        f = [x.strip() for x in line.split(',')]
        if len(f) < 9:
            continue
        r = float(f[0]) + float(f[1]) + float(f[2])
        power, tau, threshold = float(f[5]), float(f[6]), float(f[8])
        if power <= 0.0 or tau >= threshold:
            tau = threshold
        rows.append((r, max(0.0, tau)))
    rows.sort()
    inner, outer = rows[0][0], rows[-1][0]
    print(f'{len(rows)} samples, {inner:.1f} .. {outer:.1f} km')

    sums = [0.0] * count
    counts = [0] * count
    for r, tau in rows:
        k = min(count - 1, int((r - inner) / (outer - inner) * count))
        sums[k] += tau
        counts[k] += 1
    profile = []
    for k in range(count):
        if counts[k]:
            profile.append(sums[k] / counts[k])
        else:  # no sample in this bin: interpolate from the previous one
            profile.append(profile[-1] if profile else 0.0)

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'wb') as f:
        f.write(b'AXRING1\0')
        f.write(struct.pack('<ddII', inner, outer, count, 0))
        f.write(struct.pack(f'<{count}f', *profile))
    print(f'wrote {out} ({os.path.getsize(out) // 1024} KB), {(outer - inner) / count:.2f} km per bin')

    # Landmarks, for the record (bin averages).
    def tau_at(r):
        return profile[min(count - 1, max(0, int((r - inner) / (outer - inner) * count)))]
    for name, r in [('C ring (80,000 km)', 80000.0), ('B ring core (110,000 km)', 110000.0),
                    ('Cassini Division (119,000 km)', 119000.0), ('A ring (130,000 km)', 130000.0),
                    ('Encke gap (133,580 km)', 133580.0)]:
        print(f'  tau at {name}: {tau_at(r):.3f}')


if __name__ == '__main__':
    main()
