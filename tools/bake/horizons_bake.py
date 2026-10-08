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

With `kepler_gm_km3_s2` (the center's GM), the spline interpolates only the
deviation from the two-body arc starting at the earlier knot (AXEPH2). For
orbits around a planet this needs several times fewer knots; the GM only
shapes the interpolation, the knots themselves are exact samples.

Output format (little endian), see src/ephem/ephemeris.hpp:
    char[8]  magic "AXEPH1\\0\\0" (plain Hermite) or "AXEPH2\\0\\0" (Kepler-relative)
    int32    target NAIF id, int32 center NAIF id
    uint32   knot count, uint32 max knot spacing in seconds (0: none; longer
             intervals are gaps without data, see `windows_near`)
    AXEPH2:  float64 reference GM (km^3/s^2)
    knots:   float64 t (TDB seconds since J2000), x, y, z (km), vx, vy, vz (km/s)

The script also prints close approaches between each spacecraft and the
targets listed in its `encounters` (local distance minima below `encounter_km`,
sampled every `encounter_step_minutes`), and with `periapsides = true` the
periapsides around the center (for scene event lists).

`refine_near = [targets]` re-fetches at 1 minute wherever the target comes
within `refine_near_km` of those (already baked) targets, e.g. moon flybys that
are too brief for the base step to notice.

`verify_step_minutes = n` checks the reduced table against fresh samples every n
minutes between the base samples, re-fetches 1-minute samples around any over
the tolerance (maneuvers, seams between trajectory files) and reduces again
(at most twice). Samples with a velocity far off their neighbours' (bogus
states during burns) are dropped before every reduction.

`windows_near = "<Horizons command>"` keeps data only where the target comes
within `window_km` of that object (e.g. a moon near a spacecraft), padded by
`window_pad_days` on each side; windows closer than 2 pads merge. The rest is
a gap (the scene falls back to mean elements there). Each window is refined
and reduced on its own.
Only the Python standard library is used.
"""

import bisect
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
    # (Ask for the last two minutes: a zero-length span is rejected right at the end of
    # coverage, and a one-minute span can round to just under the one-minute step.)
    if rows and rows[-1][0] < (stop_jd - J2000_JD) * DAY - 1.0:
        tail = horizons_vectors(command, center, stop_jd - 2.0 / 1440.0, stop_jd, 1)
        rows.extend(r for r in tail if r[0] > rows[-1][0])
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


def stumpff_c(z):
    if z > 1e-6:
        return (1.0 - math.cos(math.sqrt(z))) / z
    if z < -1e-6:
        return (math.cosh(math.sqrt(-z)) - 1.0) / -z
    return 0.5 - z / 24.0 + z * z / 720.0


def stumpff_s(z):
    if z > 1e-6:
        r = math.sqrt(z)
        return (r - math.sin(r)) / (r * r * r)
    if z < -1e-6:
        r = math.sqrt(-z)
        return (math.sinh(r) - r) / (r * r * r)
    return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0


def kepler_propagate(p, v, gm, dt):
    """Two-body state after dt (universal variables, Curtis Algorithms 3.3/3.4).
    Mirrors propagate_kepler in src/ephem/kepler.cpp."""
    r0 = math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
    if dt == 0.0:
        return p, v
    sqrt_mu = math.sqrt(gm)
    vr0 = (p[0] * v[0] + p[1] * v[1] + p[2] * v[2]) / r0
    alpha = 2.0 / r0 - (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / gm
    x = sqrt_mu * abs(alpha) * dt if abs(alpha) >= 1e-12 else sqrt_mu * dt / r0
    # f(x) increases monotonically (df/dx = r >= q, the periapsis distance): steps
    # leaving the bracket [0, sqrt(mu) dt / q] fall back to bisection (plain Newton
    # can cycle on eccentric orbits).
    hx = (p[1] * v[2] - p[2] * v[1], p[2] * v[0] - p[0] * v[2], p[0] * v[1] - p[1] * v[0])
    h2 = hx[0] * hx[0] + hx[1] * hx[1] + hx[2] * hx[2]
    q = h2 / (gm * (1.0 + math.sqrt(max(0.0, 1.0 - alpha * h2 / gm))))
    lo, hi = 0.0, 0.0
    if q > 0.0:
        lo, hi = sorted((0.0, sqrt_mu * dt / q))
    for _ in range(100):
        x2 = x * x
        z = alpha * x2
        c, sz = stumpff_c(z), stumpff_s(z)
        f = r0 * vr0 / sqrt_mu * x2 * c + (1.0 - alpha * r0) * x2 * x * sz + r0 * x - sqrt_mu * dt
        df = r0 * vr0 / sqrt_mu * x * (1.0 - z * sz) + (1.0 - alpha * r0) * x2 * c + r0
        nxt = x - f / df
        if q > 0.0:
            if f < 0.0:
                lo = x
            else:
                hi = x
            if not lo <= nxt <= hi:
                nxt = 0.5 * (lo + hi)
        dx = nxt - x
        x = nxt
        if abs(dx) <= 1e-12 * max(1.0, abs(x)):
            break
    x2 = x * x
    z = alpha * x2
    c, sz = stumpff_c(z), stumpff_s(z)
    f = 1.0 - x2 / r0 * c
    g = dt - x2 * x / sqrt_mu * sz
    r = tuple(f * p[i] + g * v[i] for i in range(3))
    rn = math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
    fdot = sqrt_mu / (rn * r0) * (z * x * sz - x)
    gdot = 1.0 - x2 / rn * c
    return r, tuple(fdot * p[i] + gdot * v[i] for i in range(3))


def make_interp(gm):
    """Position interpolator between two knots: plain Hermite, or Kepler-relative."""
    if not gm:
        return hermite

    def kepler_hermite(k0, k1, t):
        t0, p0, v0 = k0
        t1, p1, v1 = k1
        arc_t, _ = kepler_propagate(p0, v0, gm, t - t0)
        arc_1p, arc_1v = kepler_propagate(p0, v0, gm, t1 - t0)
        zero = (0.0, 0.0, 0.0)
        d = hermite((t0, zero, zero),
                    (t1, tuple(p1[i] - arc_1p[i] for i in range(3)), tuple(v1[i] - arc_1v[i] for i in range(3))), t)
        return tuple(arc_t[i] + d[i] for i in range(3))
    return kepler_hermite


def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def leave_one_out_errors(rows, interp):
    """Estimated interpolation error (km) of the interval around each knot."""
    err = [0.0] * len(rows)
    for i in range(1, len(rows) - 1):
        try:
            e2h = dist(interp(rows[i - 1], rows[i + 1], rows[i][0]), rows[i][1])
        except (OverflowError, ValueError):  # a bogus velocity (see drop_velocity_spikes)
            e2h = math.inf
        err[i] = e2h / 16.0
    if len(rows) > 2:
        err[0] = err[1]
        err[-1] = err[-2]
    return err


def refine(command, center, rows, tol_km, step_minutes, depth, interp):
    for level in range(depth):
        err = leave_one_out_errors(rows, interp)
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


def drop_velocity_spikes(rows, threshold_kms=1.0, max_spacing_s=150.0):
    """Drops samples whose velocity is far from the median of their neighbours
    (three on each side, all at most `max_spacing_s` apart). Horizons can return
    bogus states for a few minutes during short burns (THEMIS-B 2008-10-18 03:22
    TDB: up to 120 km/s instead of 1 km/s, positions ~1,700 km off a smooth path);
    the table then bridges the burn with the knots around it. Refinement narrows
    such stretches down to minute samples first. Position jumps at trajectory-file
    seams (velocities continuous) keep their samples."""
    keep = []
    for i, (t, _, v) in enumerate(rows):
        window = rows[max(0, i - 3):i + 4]
        if len(window) < 5 or any(b[0] - a[0] > max_spacing_s for a, b in zip(window, window[1:])):
            keep.append(rows[i])
            continue
        median = [sorted(w[2][k] for w in window)[len(window) // 2] for k in range(3)]
        if dist(v, median) <= threshold_kms:
            keep.append(rows[i])
        else:
            print(f'    dropped {calendar(t)}: velocity {math.sqrt(sum(x * x for x in v)):.3f} km/s, '
                  f'neighbours {math.sqrt(sum(x * x for x in median)):.3f} km/s')
    return keep


def verify(command, center, rows, knots, gm, tol_km, step_minutes, pad_minutes=15):
    """Checks the reduced table against fresh samples off the base grid (every
    `step_minutes`, shifted by 37 s). The leave-one-out estimate assumes smooth
    motion, so a jump between two base samples (a maneuver, or a seam between
    trajectory files: ARTEMIS-P1 2012-07-05 20:01 TDB, 44 km within a minute)
    escapes it. Returns `rows` with 1-minute samples `pad_minutes` around every
    sample over the tolerance, or None when all are within it."""
    table = Table(knots, gm)
    t0, t1 = rows[0][0], rows[-1][0]
    check = fetch(command, center, t0 / DAY + J2000_JD + 37.0 / DAY, t1 / DAY + J2000_JD - 1.0 / 1440, step_minutes)
    pad = pad_minutes * 60.0
    windows = []
    for t, p, _ in check:
        q = eval_eph(table, t)
        if q is None or dist(p, q) <= tol_km:
            continue
        a, b = max(t0, t - pad), min(t1, t + pad)
        if windows and a <= windows[-1][1]:
            windows[-1][1] = max(windows[-1][1], b)
        else:
            windows.append([a, b])
    print(f'    verified at {len(check)} samples every {step_minutes} min: {len(windows)} windows over tolerance')
    if not windows:
        return None
    fine = []
    for a, b in windows:
        fine.extend(fetch(command, center, a / DAY + J2000_JD, b / DAY + J2000_JD, MIN_STEP_MINUTES))
    starts = [a for a, _ in windows]
    kept = []
    for r in rows:
        i = bisect.bisect_right(starts, r[0]) - 1
        if i < 0 or r[0] > windows[i][1]:
            kept.append(r)
    return sorted(kept + fine, key=lambda r: r[0])


def decimate(rows, tol_km, interp):
    """Greedy knot reduction: keep a knot only when skipping it would break the tolerance."""
    def ok(a, b):
        for j in range(a + 1, b):
            if dist(interp(rows[a], rows[b], rows[j][0]), rows[j][1]) > tol_km:
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


def max_error(rows, kept, interp):
    """Max deviation of the reduced ephemeris from all fetched samples (km)."""
    worst = 0.0
    k = 0
    for t, p, _ in rows:
        while k < len(kept) - 2 and kept[k + 1][0] < t:
            k += 1
        worst = max(worst, dist(interp(kept[k], kept[k + 1], t), p))
    return worst


# --- Output -------------------------------------------------------------------

def write_eph(path, target_id, center_id, knots, gm=0.0, max_gap=0):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'wb') as f:
        f.write(b'AXEPH2\0\0' if gm else b'AXEPH1\0\0')
        f.write(struct.pack('<iiII', target_id, center_id, len(knots), int(max_gap)))
        if gm:
            f.write(struct.pack('<d', gm))
        for t, p, v in knots:
            f.write(struct.pack('<7d', t, *p, *v))


class Table:
    """Baked knots with their interpolator."""

    def __init__(self, knots, gm=0.0, max_gap=0):
        self.knots = knots
        self.gm = gm
        self.max_gap = max_gap  # seconds; longer intervals between knots have no data
        self.interp = make_interp(gm)

    def start(self):
        return self.knots[0][0]

    def end(self):
        return self.knots[-1][0]


def read_eph(path):
    with open(path, 'rb') as f:
        data = f.read()
    _, _, n, max_gap = struct.unpack_from('<iiII', data, 8)
    offset = 24
    gm = 0.0
    if data[:6] == b'AXEPH2':
        gm = struct.unpack_from('<d', data, offset)[0]
        offset += 8
    knots = []
    for i in range(n):
        vals = struct.unpack_from('<7d', data, offset + i * 56)
        knots.append((vals[0], vals[1:4], vals[4:7]))
    return Table(knots, gm, max_gap)


def eval_eph(table, t):
    knots = table.knots
    lo, hi = 0, len(knots) - 1
    if t <= knots[0][0] or t >= knots[-1][0]:
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if knots[mid][0] <= t:
            lo = mid
        else:
            hi = mid
    if table.max_gap and knots[hi][0] - knots[lo][0] > table.max_gap:
        return None  # in a gap
    return table.interp(knots[lo], knots[hi], t)


def calendar(t):
    jd = t / DAY + J2000_JD
    return time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime((jd - 2440587.5) * DAY)) + ' TDB'


def distance_samples(table, other, t0, t1, step):
    """(t, distance) every `step` seconds over [t0, t1] (inf where either table has no data)."""
    out = []
    t = t0
    while t <= t1:
        a, b = eval_eph(table, t), eval_eph(other, t)
        out.append((t, dist(a, b) if a is not None and b is not None else math.inf))
        t += step
    return out


def close_windows(table, other, t0, t1, step, threshold):
    """Intervals where the two bodies are closer than `threshold`, padded by one step."""
    windows = []
    inside = None
    samples = distance_samples(table, other, t0, t1, step)
    for t, d in samples:
        if d < threshold and inside is None:
            inside = t
        elif d >= threshold and inside is not None:
            windows.append((max(t0, inside - step), t))
            inside = None
    if inside is not None:
        windows.append((max(t0, inside - step), t1))
    return windows


def closest_approaches(table, other, t0, t1, step, threshold):
    """Local distance minima below `threshold`, refined by ternary search."""
    samples = distance_samples(table, other, t0, t1, step)
    out = []
    for i in range(1, len(samples) - 1):
        (ta, da), (_, db), (tc, dc) = samples[i - 1], samples[i], samples[i + 1]
        if not (db <= da and db < dc and db < threshold):
            continue
        lo, hi = ta, tc
        for _ in range(80):
            m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
            if dist(eval_eph(table, m1), eval_eph(other, m1)) < dist(eval_eph(table, m2), eval_eph(other, m2)):
                hi = m2
            else:
                lo = m1
        tm = 0.5 * (lo + hi)
        out.append((tm, dist(eval_eph(table, tm), eval_eph(other, tm))))
    return out


def refine_near(command, center, rows, gm, others, threshold, scan_step):
    """Re-fetches at 1 minute wherever the target comes within `threshold` of
    one of `others` (e.g. moon flybys): a flyby can bend the path within a few
    minutes, which hours-apart samples never see."""
    coarse = Table(rows, gm)
    t0, t1 = coarse.start(), coarse.end()
    windows = []
    for other in others:
        windows.extend(close_windows(coarse, other, max(t0, other.start()), min(t1, other.end()), scan_step,
                                     threshold))
    windows.sort()
    merged = []
    for a, b in windows:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    if not merged:
        return rows
    fine = []
    for a, b in merged:
        fine.extend(fetch(command, center, a / DAY + J2000_JD, b / DAY + J2000_JD, MIN_STEP_MINUTES))
    print(f'    near other bodies: {len(merged)} windows, {len(fine)} samples at {MIN_STEP_MINUTES} min')
    kept = [r for r in rows if not any(a <= r[0] <= b for a, b in merged)]
    return sorted(kept + fine, key=lambda r: r[0])


def near_windows(rows, command, center, start_jd, stop_jd, step_minutes, threshold, pad):
    """Windows (t0, t1) where `command` comes within `threshold` km of the target
    (both sampled every step), padded by `pad` seconds; merged when closer than 2 pads."""
    other = {round(r[0]): r[1] for r in fetch(command, center, start_jd, stop_jd, step_minutes)}
    windows = []
    for t, p, _ in rows:
        q = other.get(round(t))
        if q is None or dist(p, q) >= threshold:
            continue
        a, b = t - pad, t + pad
        if windows and a <= windows[-1][1] + 2.0 * pad:
            windows[-1][1] = max(windows[-1][1], b)
        else:
            windows.append([a, b])
    return windows


def report_encounters(name, table, others, radii, step, threshold):
    """Closest approaches to other baked bodies (same center)."""
    for other_name, other in others.items():
        t0 = max(table.start(), other.start())
        t1 = min(table.end(), other.end())
        for tc, dc in closest_approaches(table, other, t0, t1, step, threshold):
            r = radii.get(other_name)
            alt = f', altitude {dc - r:,.0f} km' if r else ''
            print(f'  encounter {name} - {other_name}: {calendar(tc)}  distance {dc:,.0f} km{alt}')


def report_periapsides(name, table):
    """Distance minima from the center: r.v changes sign from - to + between knots."""
    knots = table.knots

    def radial(t):
        lo, hi = 0, len(knots) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if knots[mid][0] <= t:
                lo = mid
            else:
                hi = mid
        # Velocity of the segment, by central difference of positions.
        p = table.interp(knots[lo], knots[hi], t)
        dt = 1.0
        a = table.interp(knots[lo], knots[hi], t - dt)
        b = table.interp(knots[lo], knots[hi], t + dt)
        return sum(p[i] * (b[i] - a[i]) for i in range(3))

    def sub_steps(k):
        """Knot interval k split to ~1/16 revolution (Kepler-relative knots can be
        whole orbits apart, hiding periapsides between them)."""
        (t0, p0, v0), (t1, _, _) = knots[k], knots[k + 1]
        if not table.gm:
            return [t0, t1]
        r = math.sqrt(sum(x * x for x in p0))
        energy = 0.5 * sum(x * x for x in v0) - table.gm / r
        pieces = 1
        if energy < 0.0:
            sma = -table.gm / (2.0 * energy)
            n = math.sqrt(table.gm / sma ** 3)
            pieces = max(1, min(256, math.ceil(n * (t1 - t0) / (2.0 * math.pi / 16.0))))
        return [t0 + (t1 - t0) * j / pieces for j in range(pieces + 1)]

    for k in range(len(knots) - 1):
        steps = sub_steps(k)
        for t0, t1 in zip(steps, steps[1:]):
            if not radial(t0) < 0.0 <= radial(t1):
                continue
            lo, hi = t0, t1
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if radial(mid) < 0.0:
                    lo = mid
                else:
                    hi = mid
            tc = 0.5 * (lo + hi)
            r = math.sqrt(sum(x * x for x in eval_eph(table, tc) or knots[k][1]))
            print(f'  periapsis {name}: {calendar(tc)}  distance {r:,.0f} km')


# --- Main ---------------------------------------------------------------------

def jd_from_text(s):
    """'YYYY-MM-DD[ HH:MM[:SS]]' (TDB) -> JD."""
    date, _, clock = s.partition(' ')
    y, m, d = (int(x) for x in date.split('-'))
    hh, mm, ss = (float(x) for x in (clock.split(':') + ['0'])[:3]) if clock else (0, 0, 0)
    mm += ss / 60
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
        gm = tgt.get('kepler_gm_km3_s2', 0.0)
        interp = make_interp(gm)
        step = int(tgt['base_step_minutes'])
        start, stop = jd_from_text(tgt['start']), jd_from_text(tgt['stop'])
        print(f'{name}: {tgt["command"]} @ {tgt["center"]}, {tgt["start"]} .. {tgt["stop"]}, step {step} min, tol {tol} km')
        rows = fetch(tgt['command'], tgt['center'], start, stop, step)
        print(f'    fetched {len(rows)} samples')
        max_gap = 0
        if tgt.get('windows_near'):
            pad = tgt.get('window_pad_days', 1.0) * DAY
            windows = near_windows(rows, tgt['windows_near'], tgt['center'], start, stop, step, tgt['window_km'], pad)
            knots, all_rows, err = [], [], 0.0
            for a, b in windows:
                part = [r for r in rows if a <= r[0] <= b]
                if len(part) < 3:
                    continue
                part = refine(tgt['command'], tgt['center'], part, tol, step, tgt.get('max_refine', 4), interp)
                part = drop_velocity_spikes(part)
                seg = decimate(part, tol, interp)
                err = max(err, max_error(part, seg, interp))
                knots.extend(seg)
                all_rows.extend(part)
            # Knots within a window are far closer than the windows are apart (>= 2 pads).
            spacing = max(b[0] - a[0] for a, b in zip(knots, knots[1:]) if b[0] - a[0] < pad)
            max_gap = int(pad)
            assert spacing < max_gap, f'{name}: knot spacing {spacing} s exceeds the gap threshold'
            covered = sum(b - a for a, b in windows) / (stop - start) / DAY
            print(f'    {len(windows)} windows near {tgt["windows_near"]} ({covered:.1%} of the span), '
                  f'{len(all_rows)} samples -> {len(knots)} knots, max error {err:.3f} km')
            write_eph(out, int(tgt['naif_id']), int(tgt['center_naif_id']), knots, gm, max_gap)
            print(f'    wrote {out} ({os.path.getsize(out) // 1024} KB)')
            baked[name] = Table(knots, gm, max_gap)
            continue
        rows = refine(tgt['command'], tgt['center'], rows, tol, step, tgt.get('max_refine', 4), interp)
        if tgt.get('refine_near'):
            missing = [o for o in tgt['refine_near'] if o not in baked]
            if missing:
                sys.exit(f'{name}: refine_near needs {missing} baked first (earlier in the config)')
            rows = refine_near(tgt['command'], tgt['center'], rows, gm, [baked[o] for o in tgt['refine_near']],
                               tgt['refine_near_km'], tgt.get('refine_near_step_minutes', 20) * 60.0)
        rows = drop_velocity_spikes(rows)
        knots = decimate(rows, tol, interp)
        if tgt.get('verify_step_minutes'):
            for _ in range(2):
                checked = verify(tgt['command'], tgt['center'], rows, knots, gm, tol, tgt['verify_step_minutes'])
                if checked is None:
                    break
                rows = drop_velocity_spikes(checked)
                knots = decimate(rows, tol, interp)
        err = max_error(rows, knots, interp)
        kind = 'Kepler-relative' if gm else 'Hermite'
        print(f'    {len(rows)} samples -> {len(knots)} knots ({kind}), max error {err:.3f} km')
        write_eph(out, int(tgt['naif_id']), int(tgt['center_naif_id']), knots, gm)
        print(f'    wrote {out} ({os.path.getsize(out) // 1024} KB)')
        baked[name] = Table(knots, gm)

    radii = {t['name']: t.get('radius_km') for t in cfg['target']}
    for tgt in cfg['target']:
        if only and tgt['name'] not in only:
            continue
        enc = tgt.get('encounters')
        if enc and tgt['name'] in baked:
            others = {o: baked[o] for o in enc if o in baked}
            report_encounters(tgt['name'], baked[tgt['name']], others, radii,
                              tgt.get('encounter_step_minutes', 360) * 60.0, tgt.get('encounter_km', 1.5e7))
        if tgt.get('periapsides') and tgt['name'] in baked:
            report_periapsides(tgt['name'], baked[tgt['name']])


if __name__ == '__main__':
    main()
