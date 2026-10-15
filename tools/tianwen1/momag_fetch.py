"""Download the Tianwen-1 orbiter positions from the MOMAG magnetometer data files.

Usage:
    python momag_fetch.py [year ...]

The Tianwen-1/MOMAG team (USTC, http://space.ustc.edu.cn/dreams/tw1_momag/) publishes
calibrated 1 Hz magnetic field data, one file per day, each line with the orbiter's
position in MSO coordinates (in Mars radii of 3389.5 km with three decimals in 2021,
later in km). Only the positions are kept: one sample a minute, position and velocity
from a least-squares cubic over the 1 Hz samples within +-2 minutes (which smooths the
3.4 km steps of the 2021 files), written to cache/<YYYYMMDD>.txt as

    UTC (ISO)   x y z (km, MSO)   vx vy vz (km/s, MSO)   quality flag

The raw files (~12.5 MB each, ~15 GB in all) are not kept. The server sends whole
files only (no range requests, no compression), so the download is slow; files
already in the cache are skipped, so the script can be stopped and rerun.

MSO: origin at Mars' center, x toward the Sun, z along Mars' orbital angular
momentum, y = z x x. tianwen1_bake.py converts to Mars-centered ICRF.

Data use: the MOMAG team asks for an acknowledgement and citations of Wang et al.
2023 (Earth Planet. Phys. 7, 1), Zou et al. 2023 (Sci. China Tech. Sci. 66, 2396) and
Wang et al. 2024 (JGR Space Physics 129, e2023JA031757); see assets/ephem/SOURCES.md.
Only the Python standard library is used.
"""

import datetime
import os
import re
import sys
import time
import urllib.request

BASE = 'https://space.ustc.edu.cn/dreams/tw1_momag/'
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
USER_AGENT = 'Astraxis-tianwen1-ephemeris/1.0 (non-commercial; positions only)'
YEARS = range(2021, 2027)
PAUSE_S = 2.0


def http_get(url, attempts=5):
    req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
    for k in range(attempts):
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                return resp.read()
        except Exception as e:  # network hiccup: back off and retry
            if k == attempts - 1:
                raise
            print(f'    retry ({e})')
            time.sleep(10 * (k + 1))


def list_files(year):
    page = http_get(f'{BASE}?magdata=cal&sr=1&year={year}').decode('utf-8', errors='replace')
    names = sorted(set(re.findall(r'TW1_MOMAG_MSO_01Hz_(\d{8})_2C_v(\d+)\.dat', page)))
    latest = {}
    for day, ver in names:  # keep the highest version of each day
        if day not in latest or int(ver) > int(latest[day]):
            latest[day] = ver
    return [(day, f'TW1_MOMAG_MSO_01Hz_{day}_2C_v{ver}.dat') for day, ver in sorted(latest.items())]


def parse(raw):
    """([(datetime, (x, y, z) km, flag)] in time order, unit note) from a MOMAG file.
    The files of 2021 give positions in Mars radii (header: "Position_MSO (x,y,z) (R_M)",
    "R_M: Mars mean radius = 3389.5 km") with three decimals (3.4 km steps); later
    ones in km with six decimals."""
    text = raw.decode('latin-1')
    scale, unit = 1.0, 'km'
    columns = re.search(r'##Columns = .*Position_MSO \(x,y,z\) \(([^)]+)\)', text)
    if columns and columns.group(1).strip() == 'R_M':
        radius = re.search(r'##R_M: .*?=\s*([0-9.]+)\s*km', text)
        if not radius:
            raise RuntimeError('positions in R_M without the radius in the header')
        scale, unit = float(radius.group(1)), f'R_M = {radius.group(1)} km'
    elif not columns or columns.group(1).strip() != 'km':
        raise RuntimeError('unknown position unit: ' + (columns.group(1) if columns else 'no Columns line'))
    rows = []
    for line in text.splitlines():
        if not line or line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 12:
            continue
        try:
            t = datetime.datetime.strptime(f[0].rstrip('Z')[:26], '%Y-%m-%dT%H:%M:%S.%f')
            p = (float(f[5]) * scale, float(f[6]) * scale, float(f[7]) * scale)
        except ValueError:
            continue
        if p == (0.0, 0.0, 0.0) or any(v != v for v in p):
            continue
        rows.append((t, p, f[11]))
    rows.sort(key=lambda r: r[0])
    return rows, unit


HALF_WINDOW_S = 120.0  # cubic fit over +-2 min: smooths the 3.4 km steps of the R_M files


def cubic_fit(ts, ps):
    """Least-squares cubic through (ts, ps) (ts centered on the fit epoch): value and slope at 0."""
    # Normal equations for [1, t, t^2, t^3], scaled by the half window for conditioning.
    s = HALF_WINDOW_S
    m = [[0.0] * 4 for _ in range(4)]
    rhs = [[0.0] * 4 for _ in range(3)]
    for t, p in zip(ts, ps):
        u = t / s
        b = (1.0, u, u * u, u * u * u)
        for i in range(4):
            for j in range(4):
                m[i][j] += b[i] * b[j]
            for k in range(3):
                rhs[k][i] += b[i] * p[k]
    out = []
    for k in range(3):  # Gaussian elimination (4x4, well conditioned)
        a = [row[:] + [rhs[k][i]] for i, row in enumerate(m)]
        for c in range(4):
            piv = max(range(c, 4), key=lambda r: abs(a[r][c]))
            a[c], a[piv] = a[piv], a[c]
            for r in range(4):
                if r != c:
                    f = a[r][c] / a[c][c]
                    for cc in range(c, 5):
                        a[r][cc] -= f * a[c][cc]
        coef = [a[i][4] / a[i][i] for i in range(4)]
        out.append((coef[0], coef[1] / s))
    return tuple(o[0] for o in out), tuple(o[1] for o in out)


def reduce(rows):
    """One sample per minute: the first in each minute whose +-2 min window is nearly
    complete (at least 200 of 241 samples, spanning both sides), with position and
    velocity from a least-squares cubic over the window."""
    out = []
    last_minute = None
    times = [r[0] for r in rows]
    lo = 0
    for i, (t, _, flag) in enumerate(rows):
        minute = t.replace(second=0, microsecond=0)
        if minute == last_minute:
            continue
        while (t - times[lo]).total_seconds() > HALF_WINDOW_S:
            lo += 1
        hi = i
        while hi + 1 < len(rows) and (times[hi + 1] - t).total_seconds() <= HALF_WINDOW_S:
            hi += 1
        if hi - lo + 1 < 200 or (t - times[lo]).total_seconds() < 0.9 * HALF_WINDOW_S \
                or (times[hi] - t).total_seconds() < 0.9 * HALF_WINDOW_S:
            continue
        ts = [(times[j] - t).total_seconds() for j in range(lo, hi + 1)]
        p, v = cubic_fit(ts, [rows[j][1] for j in range(lo, hi + 1)])
        out.append((t, p, v, flag))
        last_minute = minute
    return out


def main():
    years = [int(y) for y in sys.argv[1:]] or list(YEARS)
    os.makedirs(CACHE, exist_ok=True)
    files = []
    for y in years:
        listed = list_files(y)
        print(f'{y}: {len(listed)} daily files')
        files.extend(listed)
    done = 0
    for day, name in files:
        out = os.path.join(CACHE, day + '.txt')
        if os.path.exists(out):
            continue
        t_start = time.time()
        for attempt in range(3):  # the server now and then sends a truncated file
            raw = http_get(BASE + 'fetch.php?datafile=' + name)
            try:
                rows, unit = parse(raw)
                break
            except RuntimeError as e:
                print(f'{day}: {e} ({len(raw)} bytes), attempt {attempt + 1}')
                time.sleep(10)
        else:
            print(f'{day}: skipped')
            continue
        samples = reduce(rows)
        tmp = out + '.part'
        with open(tmp, 'w', encoding='ascii') as f:
            f.write(f'# {name}, {len(raw)} bytes, positions in {unit}; one sample a minute, cubic fit over +-2 min\n')
            for t, p, v, flag in samples:
                f.write(f'{t.isoformat()} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f} {v[0]:.7f} {v[1]:.7f} {v[2]:.7f} {flag}\n')
        os.replace(tmp, out)
        done += 1
        print(f'{day}: {len(samples)} samples ({len(raw) / 1e6:.1f} MB in {time.time() - t_start:.1f} s)', flush=True)
        time.sleep(PAUSE_S)
    print(f'{done} new days')


if __name__ == '__main__':
    main()
