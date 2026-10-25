"""Check WMAP/Chandra tables against independent Horizons samples.

Run after horizons_bake.py: python tools/bake/verify_wmap_chandra.py [target ...]
Queries a grid offset by 17 seconds, every 173 minutes, plus minute samples
after WMAP launch, at its lunar encounter and in Chandra's initial orbit. Downloads
are cached in tools/bake/cache; remove the corresponding JSON to refresh.
The 5.1 km check allows a small margin around the bake's 5 km knot tolerance.
"""

import json
from pathlib import Path
import sys

import horizons_bake as h


ROOT = Path(__file__).resolve().parents[2]


def main():
    failed = False
    selected = set(sys.argv[1:]) or {'wmap_geo', 'chandra_geo'}
    unknown = selected - {'wmap_geo', 'chandra_geo'}
    if unknown:
        sys.exit(f'Unknown targets: {sorted(unknown)}')
    for name, command, detail_start, detail_stop in (
        ('wmap_geo', '-165', '2001-07-30', '2001-07-31'),
        ('chandra_geo', '-151', '1999-08-08', '1999-08-14'),
    ):
        if name not in selected:
            continue
        table = h.read_eph(ROOT / 'assets/ephem' / (name + '.eph'))
        ranges = [
            (table.start() / h.DAY + h.J2000_JD + 17 / h.DAY,
             table.end() / h.DAY + h.J2000_JD - 2 / 1440, 173),
            (h.jd_from_text(detail_start) + 17 / h.DAY,
             h.jd_from_text(detail_stop), 1),
        ]
        if name == 'wmap_geo':
            # Resolve the fast motion just after launch as well as the lunar flyby.
            ranges.append((table.start() / h.DAY + h.J2000_JD + 29 / h.DAY,
                           h.jd_from_text('2001-07-02'), 1))
        samples = []
        for start, stop, step in ranges:
            cache = ROOT / 'tools/bake/cache' / f'{name}_{start:.9f}_{stop:.9f}_{step}.json'
            cache.parent.mkdir(exist_ok=True)
            if cache.exists():
                rows = json.loads(cache.read_text(encoding='utf-8'))
            else:
                rows = h.fetch(command, '500@399', start, stop, step)
                cache.write_text(json.dumps(rows), encoding='utf-8')
            samples.extend(rows)
        errors = sorted((h.dist(h.eval_eph(table, t), p), t) for t, p, _ in samples)
        over = [(error, t) for error, t in errors if error > 5.1]
        print(f'{name}: {len(samples)} samples, max {errors[-1][0]:.6f} km, '
              f'{len(over)} over 5.1 km', flush=True)
        for error, t in errors[-10:]:
            print(f'  {h.calendar(t)}: {error:.6f} km', flush=True)
        failed |= bool(over)
    return int(failed)


if __name__ == '__main__':
    sys.exit(main())
