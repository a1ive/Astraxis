"""Bake geometric Earth/Sun and Moon/Earth vectors around solar and lunar eclipses.

Run from any directory: python -E -P tools/corona/eclipse_bake.py
Uses horizons_bake.py's ICRF/TDB vector query and AXEPH1 writer. Cached JSON
responses live in tools/bake/cache/corona/. Off-grid samples independently
check the 30-minute Hermite tables. Prints designed observer offsets for a
space view 1000 km above Earth (solar eclipse) or 200 km above the lunar near
side (lunar eclipse), at each event's reference time.
These fixed ICRF offsets are not ground observatories or spacecraft orbits.
"""

import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'bake'))
import horizons_bake as h


def cached(command, center, start, stop, step, name):
    cache = ROOT / 'tools' / 'bake' / 'cache' / 'corona' / (name + '.json')
    if cache.exists():
        return json.loads(cache.read_text(encoding='utf-8'))
    rows = h.fetch(command, center, start, stop, step)
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(rows), encoding='utf-8')
    return rows


def main():
    windows = [('2017-08-20', '2017-08-23', '2017-08-21 18:26:40'),
               ('2022-11-07', '2022-11-10', '2022-11-08 10:59:00'),
               ('2024-04-07', '2024-04-10', '2024-04-08 18:17:20')]
    all_rows = {'earth': [], 'moon': []}
    for start, stop, utc in windows:
        a, b = h.jd_from_text(start), h.jd_from_text(stop)
        tables = {}
        for name, command, center, target_id, center_id in [
                ('earth', '399', '500@10', 399, 10), ('moon', '301', '500@399', 301, 399)]:
            rows = cached(command, center, a, b, 30, name + start)
            table = h.Table(rows)
            mid = cached(command, center, a + 15 / 1440, b - 15 / 1440, 30, name + start + '-verify')
            error = max(h.dist(h.eval_eph(table, r[0]), r[1]) for r in mid)
            if error > 0.01:
                raise RuntimeError(f'{name}: {error} km exceeds 10 m interpolation budget')
            print(f'{name} {start} .. {stop}: {len(rows)} knots, off-grid max {error:.9f} km', flush=True)
            all_rows[name].extend(rows)
            tables[name] = table
        # UTC to TDB: TT-UTC=69.184 s on these dates; the project's <2 ms
        # TDB-TT approximation also applies here (core/time.hpp).
        t = (h.jd_from_text(utc) - h.J2000_JD) * h.DAY + 69.184
        earth = h.eval_eph(tables['earth'], t)
        moon = h.eval_eph(tables['moon'], t)
        if start == '2022-11-07':
            # A designed point 200 km above the lunar near side at this instant.
            # It follows the Moon at a fixed ICRF offset, not a landed observer.
            norm = math.sqrt(sum(v*v for v in moon))
            observer = [-v / norm * (1737.4 + 200.0) for v in moon]
            print(f'{utc} UTC: lunar observer_offset_km = [{", ".join(f"{v:.9f}" for v in observer)}]', flush=True)
            print(f'  geocentric lunar distance {norm:.6f} km', flush=True)
            continue
        direction = [moon[i] + earth[i] for i in range(3)]
        norm = math.sqrt(sum(v*v for v in direction))
        direction = [v / norm for v in direction]
        dot = sum(moon[i] * direction[i] for i in range(3))
        side = [moon[i] - dot * direction[i] for i in range(3)]
        radius = 6378.1366 + 1000.0
        discriminant = radius*radius - sum(v*v for v in side)
        if discriminant <= 0:
            raise RuntimeError('shadow axis misses the observer altitude')
        distance = -dot - math.sqrt(discriminant)
        observer = [moon[i] + distance * direction[i] for i in range(3)]
        sun_distance = math.sqrt(sum((earth[i]+observer[i])**2 for i in range(3)))
        ratio = (1737.4/distance) / (695700.0/sun_distance)
        print(f'{utc} UTC: observer_offset_km = [{", ".join(f"{v:.9f}" for v in observer)}]', flush=True)
        print(f'  Moon distance {distance:.6f} km, apparent Moon/Sun radius ratio {ratio:.9f}', flush=True)
    for name, target, center in [('earth', 399, 10), ('moon', 301, 399)]:
        h.write_eph(str(ROOT / 'assets' / 'ephem' / f'eclipse_{name}.eph'), target, center,
                    all_rows[name], max_gap=3600)


if __name__ == '__main__':
    main()
