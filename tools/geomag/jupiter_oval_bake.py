"""Bake Head (2024) Juno-UVS reference ovals (CC BY 4.0).

Download UVSRefOval_N.csv and UVSRefOval_S.csv from
https://zenodo.org/records/12067692 into the input directory.
Usage: python tools/geomag/jupiter_oval_bake.py INPUT_DIR OUTPUT_FILE
Only the Python standard library is needed.
"""

import bisect
import csv
import hashlib
import math
from pathlib import Path
import sys

COUNT = 256
RADIUS_KM = 71492.0
HASHES = {
    "UVSRefOval_N.csv": "96023e6b935814acbb76c6dfa694a12e9b64731e338b2e18a22c7590e364e033",
    "UVSRefOval_S.csv": "5ae737ae83515b1e6d4563209424680af82259b501f9fed670ea519149bb95a4",
}


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def direction(lon, lat):
    lon, lat = math.radians(lon), math.radians(lat)
    return math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)


def polar(lon_west, lat, center_west, center_lat):
    center = direction(-center_west, center_lat)
    lon = math.radians(-center_west)
    latitude = math.radians(center_lat)
    east = (-math.sin(lon), math.cos(lon), 0.0)
    north = (-math.sin(latitude) * math.cos(lon), -math.sin(latitude) * math.sin(lon), math.cos(latitude))
    p = direction(-lon_west, lat)
    return (math.atan2(dot(p, north), dot(p, east)) % math.tau,
            math.acos(max(-1.0, min(1.0, dot(p, center)))))


def sample(points, angle):
    i = bisect.bisect_right(points, (angle, math.inf))
    a, b = points[i - 1], points[i % len(points)]
    x0 = a[0] if i else a[0] - math.tau
    x1 = b[0] if i < len(points) else b[0] + math.tau
    return a[1] + (b[1] - a[1]) * (angle - x0) / (x1 - x0)


def bake(path, center_west, center_lat):
    if hashlib.sha256(path.read_bytes()).hexdigest() != HASHES[path.name]:
        raise ValueError(f"Unexpected source checksum: {path.name}")
    with path.open(encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
    points = sorted(set(polar(float(r["longitude"]), float(r["latitude"]), center_west, center_lat) for r in rows))
    gaps = [(points[(i + 1) % len(points)][0] - p[0]) % math.tau for i, p in enumerate(points)]
    if max(gaps) > 0.02:
        raise ValueError("Incomplete angular coverage")
    result = [(math.tau * i / COUNT, sample(points, math.tau * i / COUNT)) for i in range(COUNT)]
    # Check every original point, including those between the baked samples.
    errors = [abs(sample(result, a) - r) * RADIUS_KM for a, r in points]
    worst = max(errors)
    if worst > 100.0:
        raise ValueError(f"Resampling error {worst:.3f} km exceeds 100 km")
    print(f"{path.name}: {len(rows)} points, max radial error {worst:.3f} km")
    return [r for _, r in result]


def main():
    source, output = map(Path, sys.argv[1:])
    # Head et al. (2024), section 4.1: System III WEST longitude,
    # planetocentric latitude. The tangent frame uses east, north, center.
    north = bake(source / "UVSRefOval_N.csv", 185.0, 74.0)
    south = bake(source / "UVSRefOval_S.csv", 32.0, -82.0)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as target:
        target.write(f"AXJOVIAN1 {COUNT}\n")
        for n, s in zip(north, south):
            target.write(f"{n:.10f} {s:.10f}\n")


if __name__ == "__main__":
    main()
