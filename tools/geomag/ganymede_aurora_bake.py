"""Bake Marzok et al. 2022 Fig. 5 (CC BY 4.0) into a zenith brightness map.

Usage: python tools/geomag/ganymede_aurora_bake.py Fig5_BrightnessMap.txt output.txt
Source and display assumptions: assets/geomag/SOURCES.md.
"""
import hashlib
import math
from pathlib import Path
import sys


def bake(source: Path, output: Path) -> None:
    raw = source.read_bytes()
    expected = "f0dd34cdb5fd0f43353631eb26c488866c78e4081cf829879bd419028646a5d5"
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError("Unexpected Fig5_BrightnessMap.txt revision")
    lines = raw.decode("ascii").splitlines()
    longitude = list(map(float, lines[6].split()))
    latitude = list(map(float, lines[8].split()))
    rows = [list(map(float, line.split())) for line in lines[10:]]
    assert longitude == [2.25 + 4.5 * i for i in range(80)]
    assert latitude == [-87.75 + 4.5 * i for i in range(40)]
    assert len(rows) == 40 and all(len(row) == 80 for row in rows)
    assert all(math.isfinite(v) and -100 < v < 200 for row in rows for v in row)
    # Keep west longitude and south-to-north rows. Negative background residuals
    # cannot emit light; zero-valued missing cells remain dark.
    text = "AXGANYMEDE1 80 40\n"
    text += "".join(" ".join(f"{max(v, 0):.2f}" for v in row) + "\n" for row in rows)
    output.write_text(text, encoding="ascii", newline="\n")


if __name__ == "__main__":
    bake(Path(sys.argv[1]), Path(sys.argv[2]))
