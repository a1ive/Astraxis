"""Bake JPL Horizons state vectors into compact ephemeris files (.eph).

Usage:
    python horizons_bake.py <config.toml> [target_name ...]

Each [[target]] in the config is fetched from the Horizons API
(https://ssd.jpl.nasa.gov/api/horizons.api) as geometric ICRF state vectors in
TDB, then:
  1. Adaptive refinement: the cubic Hermite interpolation error at every knot is
     estimated by a leave-one-out test (interpolate across the knot from its
     neighbours, which have twice the spacing; the error scales as h^4, so the
     error at the real spacing is about 1/16 of that). Windows whose estimate
     exceeds the tolerance are re-fetched with a step 8x smaller, repeatedly.
     This automatically resolves planetary flybys down to minutes.
  2. Knot reduction: knots are dropped greedily as long as every original
     sample is still reproduced within the tolerance.

Output format (little endian), see src/ephem/ephemeris.hpp:
    char[8]  magic "AXEPH1\\0\\0"
    int32    target NAIF id, int32 center NAIF id
    uint32   knot count, uint32 reserved (0)
    knots:   float64 t (TDB seconds since J2000), x, y, z (km), vx, vy, vz (km/s)

The script also prints close approaches between each spacecraft and the
targets listed in its `encounters` (for scene event lists).
Only the Python standard library is used.
"""

import math
import os
import struct
import sys
import time
import tomllib
import urllib.parse
import urllib.request

API = 'https://ssd.jpl.nasa.gov/api/horizons.api'
J2000_JD = 2451545.0
DAY = 86400.0
MAX_ROWS_PER_QUERY = 40000
MIN_STEP_MINUTES = 1


# --- Horizons ---------------------------------------------------------------

def horizons_vectors(command, center, start_jd, stop_jd, step_minutes):
    """Returns a list of (t_tdb_seconds, (x, y, z), (vx, vy, vz))."""
    params = {
        'format': 'text',
        'COMMAND': f"'{command}'",
        'OBJ_DATA': "'NO'",
        'MAKE_EPHEM': "'YES'",
        'EPHEM_TYPE': "'VECTORS'",
        'CENTER': f"'{center}'",
        'START_TIME': f"'JD{start_jd:.9f}'",
        'STOP_TIME': f"'JD{stop_jd:.9f}'",
        'STEP_SIZE': f"'{int(step_minutes)} m'",
        'VEC_TABLE': "'2'",
        'REF_SYSTEM': "'ICRF'",
        'REF_PLANE': "'FRAME'",
        'OUT_UNITS': "'KM-S'",
        'VEC_CORR': "'NONE'",
        'CSV_FORMAT': "'YES'",
        'TIME_TYPE': "'TDB'",
    }
    url = API + '?' + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                text = resp.read().decode('utf-8', errors='replace')
            break
        except Exception as e:  # network hiccup: back off and retry
            if attempt == 4:
                raise
            print(f'    retry ({e})')
            time.sleep(5 * (attempt + 1))

    if '$$SOE' not in text:
        raise RuntimeError(f'Horizons returned no data for {command}:\n{text[-2000:]}')
    body = text.split('$$SOE', 1)[1].split('$$EOE', 1)[0]
    rows = []
    for line in body.strip().splitlines():
        f = [s.strip() for s in line.split(',')]
        jd = float(f[0])
        p = (float(f[2]), float(f[3]), float(f[4]))
        v = (float(f[5]), float(f[6]), float(f[7]))
        rows.append(((jd - J2000_JD) * DAY, p, v))
    return rows


def fetch(command, center, start_jd, stop_jd, step_minutes):
    """Fetches [start, stop] in chunks small enough for the API."""
    rows = []
    span_minutes = (stop_jd - start_jd) * 1440.0
    chunks = max(1, math.ceil(span_minutes / step_minutes / MAX_ROWS_PER_QUERY))
    for c in range(chunks):
        a = start_jd + (stop_jd - start_jd) * c / chunks
        b = start_jd + (stop_jd - start_jd) * (c + 1) / chunks
        part = horizons_vectors(command, center, a, b, step_minutes)
        if rows and part and part[0][0] <= rows[-1][0]:
            part = [r for r in part if r[0] > rows[-1][0]]
        rows.extend(part)
    # Horizons only emits multiples of the step; make sure the stop time is a knot.
    if rows and rows[-1][0] < (stop_jd - J2000_JD) * DAY - 1.0:
        rows.extend(horizons_vectors(command, center, stop_jd, stop_jd + 1e-6, 1))
    return rows


# --- Hermite ----------------------------------------------------------------

def hermite(k0, k1, t):
    t0, p0, v0 = k0
    t1, p1, v1 = k1
    h = t1 - t0
    s = (t - t0) / h
    s2 = s * s
    s3 = s2 * s
    h00 = 2 * s3 - 3 * s2 + 1
    h10 = s3 - 2 * s2 + s
    h01 = -2 * s3 + 3 * s2
    h11 = s3 - s2
    return tuple(h00 * p0[i] + h10 * h * v0[i] + h01 * p1[i] + h11 * h * v1[i] for i in range(3))


def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def leave_one_out_errors(rows):
    """Estimated interpolation error (km) of the interval around each knot."""
    err = [0.0] * len(rows)
    for i in range(1, len(rows) - 1):
        e2h = dist(hermite(rows[i - 1], rows[i + 1], rows[i][0]), rows[i][1])
        err[i] = e2h / 16.0
    if len(rows) > 2:
        err[0] = err[1]
        err[-1] = err[-2]
    return err


def refine(command, center, rows, tol_km, step_minutes, depth):
    for level in range(depth):
        err = leave_one_out_errors(rows)
        bad = [i for i, e in enumerate(err) if e > tol_km]
        if not bad or step_minutes <= MIN_STEP_MINUTES:
            break
        new_step = max(MIN_STEP_MINUTES, step_minutes // 8)
        # Merge bad knots into windows, padded by one knot on each side.
        windows = []
        for i in bad:
            a, b = max(0, i - 1), min(len(rows) - 1, i + 1)
            if windows and a <= windows[-1][1]:
                windows[-1][1] = max(windows[-1][1], b)
            else:
                windows.append([a, b])
        print(f'    level {level + 1}: {len(bad)} knots over tolerance -> {len(windows)} windows at {new_step} min')
        out = []
        cursor = 0
        for a, b in windows:
            out.extend(rows[cursor:a])
            ja = rows[a][0] / DAY + J2000_JD
            jb = rows[b][0] / DAY + J2000_JD
            fine = fetch(command, center, ja, jb, new_step)
            out.extend(fine)
            cursor = b + 1
            while cursor < len(rows) and rows[cursor][0] <= fine[-1][0]:
                cursor += 1
        out.extend(rows[cursor:])
        rows = out
        step_minutes = new_step
    return rows


def decimate(rows, tol_km):
    """Greedy knot reduction: keep a knot only when skipping it would break the tolerance."""
    def ok(a, b):
        for j in range(a + 1, b):
            if dist(hermite(rows[a], rows[b], rows[j][0]), rows[j][1]) > tol_km:
                return False
        return True

    kept = [0]
    a = 0
    n = len(rows)
    while a < n - 1:
        # Exponential search for the reach, then binary search.
        step = 1
        while a + step * 2 < n and ok(a, a + step * 2):
            step *= 2
        lo, hi = a + step, min(n - 1, a + step * 2)
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if ok(a, mid):
                lo = mid
            else:
                hi = mid - 1
        kept.append(lo)
        a = lo
    return [rows[i] for i in kept]


def max_error(rows, kept):
    """Max deviation of the reduced ephemeris from all fetched samples (km)."""
    worst = 0.0
    k = 0
    for t, p, _ in rows:
        while k < len(kept) - 2 and kept[k + 1][0] < t:
            k += 1
        worst = max(worst, dist(hermite(kept[k], kept[k + 1], t), p))
    return worst


# --- Output -------------------------------------------------------------------

def write_eph(path, target_id, center_id, knots):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'wb') as f:
        f.write(b'AXEPH1\0\0')
        f.write(struct.pack('<iiII', target_id, center_id, len(knots), 0))
        for t, p, v in knots:
            f.write(struct.pack('<7d', t, *p, *v))


def read_eph(path):
    with open(path, 'rb') as f:
        data = f.read()
    _, _, n, _ = struct.unpack_from('<iiII', data, 8)
    knots = []
    for i in range(n):
        vals = struct.unpack_from('<7d', data, 24 + i * 56)
        knots.append((vals[0], vals[1:4], vals[4:7]))
    return knots


def eval_eph(knots, t):
    lo, hi = 0, len(knots) - 1
    if t <= knots[0][0] or t >= knots[-1][0]:
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if knots[mid][0] <= t:
            lo = mid
        else:
            hi = mid
    return hermite(knots[lo], knots[hi], t)


def calendar(t):
    jd = t / DAY + J2000_JD
    return time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime((jd - 2440587.5) * DAY)) + ' TDB'


def report_encounters(name, knots, others, radii):
    """Closest approaches to other baked bodies (same center), within 0.1 AU."""
    for other_name, other in others.items():
        t0 = max(knots[0][0], other[0][0])
        t1 = min(knots[-1][0], other[-1][0])
        t = t0 + DAY
        prev = None
        candidates = []
        while t < t1 - DAY:
            a, b = eval_eph(knots, t), eval_eph(other, t)
            d = dist(a, b)
            if prev is not None and d < 1.5e7:
                candidates.append((d, t))
            prev = d
            t += 0.25 * DAY
        if not candidates:
            continue
        # Refine the global minimum of each separate close-approach episode.
        candidates.sort(key=lambda c: c[1])
        episodes = []
        for d, t in candidates:
            if episodes and t - episodes[-1][-1][1] < 5 * DAY:
                episodes[-1].append((d, t))
            else:
                episodes.append([(d, t)])
        for ep in episodes:
            d, tc = min(ep)
            lo, hi = tc - 0.5 * DAY, tc + 0.5 * DAY
            for _ in range(60):  # golden-section-like shrink
                m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
                f1 = dist(eval_eph(knots, m1), eval_eph(other, m1))
                f2 = dist(eval_eph(knots, m2), eval_eph(other, m2))
                if f1 < f2:
                    hi = m2
                else:
                    lo = m1
            tc = (lo + hi) / 2
            dc = dist(eval_eph(knots, tc), eval_eph(other, tc))
            r = radii.get(other_name)
            alt = f', altitude {dc - r:,.0f} km' if r else ''
            print(f'  encounter {name} - {other_name}: {calendar(tc)}  distance {dc:,.0f} km{alt}')


# --- Main ---------------------------------------------------------------------

def jd_from_text(s):
    """'YYYY-MM-DD[ HH:MM]' (TDB) -> JD."""
    date, _, clock = s.partition(' ')
    y, m, d = (int(x) for x in date.split('-'))
    hh, mm = (int(x) for x in clock.split(':')) if clock else (0, 0)
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + b - 1524.5 + (hh + mm / 60) / 24


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    config_path = sys.argv[1]
    only = set(sys.argv[2:])
    with open(config_path, 'rb') as f:
        cfg = tomllib.load(f)
    root = os.path.dirname(os.path.abspath(config_path))
    defaults = cfg.get('defaults', {})
    default_tol = defaults.get('tolerance_km', 1.0)
    output_dir = defaults.get('output_dir', '.')

    baked = {}
    for tgt in cfg['target']:
        name = tgt['name']
        out = os.path.normpath(os.path.join(root, tgt.get('output', os.path.join(output_dir, name + '.eph'))))
        if only and name not in only:
            if os.path.exists(out):
                baked[name] = read_eph(out)
            continue
        tol = tgt.get('tolerance_km', default_tol)
        step = int(tgt['base_step_minutes'])
        start, stop = jd_from_text(tgt['start']), jd_from_text(tgt['stop'])
        print(f'{name}: {tgt["command"]} @ {tgt["center"]}, {tgt["start"]} .. {tgt["stop"]}, step {step} min, tol {tol} km')
        rows = fetch(tgt['command'], tgt['center'], start, stop, step)
        print(f'    fetched {len(rows)} samples')
        rows = refine(tgt['command'], tgt['center'], rows, tol, step, tgt.get('max_refine', 4))
        knots = decimate(rows, tol)
        err = max_error(rows, knots)
        print(f'    {len(rows)} samples -> {len(knots)} knots, max error {err:.3f} km')
        write_eph(out, int(tgt['naif_id']), int(tgt['center_naif_id']), knots)
        print(f'    wrote {out} ({os.path.getsize(out) // 1024} KB)')
        baked[name] = knots

    radii = {t['name']: t.get('radius_km') for t in cfg['target']}
    for tgt in cfg['target']:
        enc = tgt.get('encounters')
        if enc and tgt['name'] in baked:
            others = {o: baked[o] for o in enc if o in baked}
            report_encounters(tgt['name'], baked[tgt['name']], others, radii)


if __name__ == '__main__':
    main()
