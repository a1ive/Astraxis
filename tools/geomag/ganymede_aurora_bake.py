"""Bake Marzok et al. 2022 Figs. 5 and 6 (CC BY 4.0) into Ganymede's brightness maps.

Usage: python -E -P tools/geomag/ganymede_aurora_bake.py <directory with Fig5_*, Fig6_*> output.txt

Four layers of 80 x 40 zenith Rayleigh: the all-exposure map, then Ganymede above, inside
and below Jupiter's current sheet. The exposures of a group cover less of the moon than
all of them do (between 180 and 225 W nearly only "above" looked), and each has a third
of the photons. So a group is not drawn from its own cells: it is the main map times a
smooth ratio of the group to the main map over the cells it did observe, which goes to 1
where it has no data. The ratios are then scaled so that, over the time Ganymede spends in
each group, the layers average to the main map in every cell. Source and display
assumptions: assets/geomag/SOURCES.md.
"""
import hashlib
import math
from pathlib import Path
import sys

import numpy as np

SOURCES = {
    "Fig5_BrightnessMap.txt": "f0dd34cdb5fd0f43353631eb26c488866c78e4081cf829879bd419028646a5d5",
    "Fig6_AboveCurrentSheet.txt": "905b24c0c211834681943802f829146d1ee38a95bf35435f708f615ff0c937d0",
    "Fig6_InsideCurrentSheet.txt": "22cc0e9b5cdcefa70161e6ebb09f8b89d599ca0d62ea27d3f4c238e35463e0c4",
    "Fig6_BelowCurrentSheet.txt": "1bc357c1d8a3044949facccb77647eb2f1b9a3c9cd637e6962ac0b6c516961b5",
}
EXPOSURES = (14.0, 19.0, 13.0)  # above, inside, below: Marzok section 3.2
SMOOTHING_DEG = 8.0  # keeps the 31 degree wide ovals' shift, loses single noisy cells
PRIOR_RAYLEIGH = 5.0  # a kernel's mean brightness that weighs as much as "no change"
RANGE_RAYLEIGH = 260  # the layers' largest value, rounded up

# Marzok's table 1 is a dipole tilted 9.5 degrees (test_ganymede_aurora), the groups
# split at +-6 degrees of magnetic latitude, section 3.2.
_INSIDE = 2.0 / math.pi * math.asin(math.sin(math.radians(6.0)) / math.sin(math.radians(9.5)))
TIME_SHARES = np.array([(1.0 - _INSIDE) / 2.0, _INSIDE, (1.0 - _INSIDE) / 2.0])


def read_map(path: Path) -> np.ndarray:
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SOURCES[path.name]:
        raise ValueError(f"Unexpected {path.name} revision")
    lines = raw.decode("ascii").splitlines()
    assert list(map(float, lines[6].split())) == [2.25 + 4.5 * i for i in range(80)]
    assert list(map(float, lines[8].split())) == [-87.75 + 4.5 * i for i in range(40)]
    rows = np.array([list(map(float, line.split())) for line in lines[10:]])
    assert rows.shape == (40, 80) and np.all(np.isfinite(rows)) and np.all(np.abs(rows) < 300.0)
    return rows


def oval_band_mean(cells: np.ndarray, west_low: float, west_high: float) -> float:
    """Marzok's oval bands: Eq. 8 of the paper plus three cells on each side, Table 2."""
    total = 0.0
    count = 0
    for i in range(80):
        west = 2.25 + 4.5 * i
        if not west_low <= west <= west_high:
            continue
        for mean, amplitude, phase in ((32.3, -16.9, 1.5), (-29.4, 11.6, 8.8)):
            latitude = mean + amplitude * math.sin(math.radians(west + phase))
            center = round((latitude + 87.75) / 4.5)
            for j in range(center - 3, center + 4):
                if cells[j, i] != 0.0:
                    total += cells[j, i]
                    count += 1
    return total / count


def inside_over_outside(groups, west_low: float, west_high: float) -> float:
    above, inside, below = (oval_band_mean(g, west_low, west_high) for g in groups)
    return inside / ((EXPOSURES[0] * above + EXPOSURES[2] * below) / (EXPOSURES[0] + EXPOSURES[2]))


def bake(directory: Path, output: Path) -> None:
    main = read_map(directory / "Fig5_BrightnessMap.txt")
    groups = [read_map(directory / f"Fig6_{name}CurrentSheet.txt") for name in ("Above", "Inside", "Below")]

    latitude = np.radians(-87.75 + 4.5 * np.arange(40))
    west = np.radians(2.25 + 4.5 * np.arange(80))
    lat_grid, west_grid = np.meshgrid(latitude, west, indexing="ij")
    points = np.stack([np.cos(lat_grid) * np.cos(west_grid), np.cos(lat_grid) * np.sin(west_grid), np.sin(lat_grid)], -1).reshape(-1, 3)
    area = np.cos(lat_grid).ravel()
    angle = np.arccos(np.clip(points @ points.T, -1.0, 1.0))
    kernel = np.exp(-0.5 * (angle / math.radians(SMOOTHING_DEG)) ** 2)
    prior = PRIOR_RAYLEIGH * (kernel * area).sum(1)

    ratios = []
    for group in groups:
        # The negative residuals stay in here: clipping them would bias a faint region up.
        seen = ((group != 0.0) & (main != 0.0)).ravel()
        ratios.append(((kernel * (area * group.ravel() * seen)).sum(1) + prior) /
                      ((kernel * (area * main.ravel() * seen)).sum(1) + prior))
    ratios = np.array(ratios)
    ratios /= (TIME_SHARES[:, None] * ratios).sum(0)
    # As the map has always been: negative residuals emit nothing, missing cells stay dark.
    base = np.maximum(main, 0.0)
    layers = [base] + [base * r.reshape(40, 80) for r in ratios]
    largest = max(layer.max() for layer in layers)
    assert largest < RANGE_RAYLEIGH, largest

    # Marzok section 3.3: inside the sheet the downstream aurora brightens by 1.3 +- 0.31
    # and the upstream one dims by 0.78 +- 0.19.
    filled = layers[1:]
    for label, band, value, error in (("downstream", (40.0, 150.0), 1.3, 0.31), ("upstream", (240.0, 340.0), 0.78, 0.19)):
        raw = inside_over_outside(groups, *band)
        baked = inside_over_outside(filled, *band)
        print(f"{label}: inside/outside raw {raw:.3f}, baked {baked:.3f}, Marzok {value} +- {error}")
        assert abs(raw - value) < error and abs(baked - value) < error
    print(f"time shares {TIME_SHARES}, ratios {ratios.min():.3f} .. {ratios.max():.3f}, largest {largest:.1f} R")
    mean = sum(share * layer for share, layer in zip(TIME_SHARES, filled))
    assert np.abs(mean - base).max() < 1e-9

    text = f"AXAURORAMAP4 80 40 {RANGE_RAYLEIGH}\n"
    for layer in layers:
        text += "".join(" ".join(f"{v:.2f}" for v in row) + "\n" for row in layer)
    output.write_text(text, encoding="ascii", newline="\n")


if __name__ == "__main__":
    bake(Path(sys.argv[1]), Path(sys.argv[2]))
