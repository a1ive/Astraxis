"""Bake the radial field maps of [bodies.magnetic_map]; see assets/geomag/SOURCES.md.

    python -E -P tools/geomag/field_map_bake.py <Morschhauser2014.txt> assets/geomag

Needs numpy. The output maps the body's ellipsoid like its texture: rows from north
to south by parametric latitude (the body mesh's), columns from 180 W eastward. A
header line "AXFIELDMAP1 width height scale_nT" is followed by width * height signed
bytes, 127 * sign(Br) * sqrt(|Br| / scale): the square root keeps the weak fields of
a map whose peak is a hundred times its typical value.
"""

import sys
from pathlib import Path

import numpy as np


def load_gauss(path):
    """g, h by [n][m] from "g n m value" lines, or "n m value" (m < 0: h) lines."""
    rows = []
    for line in Path(path).read_text(encoding="ascii").splitlines():
        parts = line.split()
        if not parts or parts[0].startswith("#") or parts[0] == "degree":
            continue
        if parts[0] in ("g", "h"):
            rows.append((parts[0], int(parts[1]), int(parts[2]), float(parts[3])))
        else:
            n, m = int(parts[0]), int(parts[1])
            rows.append(("h" if m < 0 else "g", n, abs(m), float(parts[2])))
    degree = max(r[1] for r in rows)
    g = np.zeros((degree + 1, degree + 1))
    h = np.zeros((degree + 1, degree + 1))
    for kind, n, m, value in rows:
        (g if kind == "g" else h)[n, m] = value
    return g, h


def radial_field(g, h, cos_theta, sin_theta, radius, width):
    """Br (nT) on rows of constant colatitude and radius (in model radii), width longitudes
    from -180 eastward, each at the middle of its cell."""
    degree = g.shape[0] - 1
    rows = len(cos_theta)
    even = np.zeros((rows, degree + 1))  # cos(m phi) terms
    odd = np.zeros((rows, degree + 1))   # sin(m phi) terms
    before = {}
    last = {0: np.ones(rows)}
    diagonal = np.ones(rows)
    for n in range(1, degree + 1):
        current = {}
        for m in range(n):
            above = np.sqrt((n - 1) ** 2 - m * m) * before[m] if m in before else 0.0
            current[m] = ((2 * n - 1) * cos_theta * last[m] - above) / np.sqrt(n * n - m * m)
        diagonal = diagonal * sin_theta * (1.0 if n == 1 else np.sqrt((2 * n - 1) / (2 * n)))
        current[n] = diagonal
        factor = (n + 1) * radius ** -(n + 2)
        for m in range(n + 1):
            even[:, m] += factor * g[n, m] * current[m]
            odd[:, m] += factor * h[n, m] * current[m]
        before, last = last, current
    longitude = (np.arange(width) + 0.5) / width * 2.0 * np.pi - np.pi
    out = np.zeros((rows, width))
    for m in range(degree + 1):
        out += np.outer(even[:, m], np.cos(m * longitude)) + np.outer(odd[:, m], np.sin(m * longitude))
    return out


def ellipsoid_rows(height, equatorial_km, polar_km, altitude_km, reference_km):
    """Colatitude and radius (in model radii) of the mesh's rows, altitude_km above them."""
    beta = np.pi / 2 - (np.arange(height) + 0.5) / height * np.pi
    x = equatorial_km * np.cos(beta)
    z = polar_km * np.sin(beta)
    r = np.hypot(x, z)
    return z / r, x / r, (r + altitude_km) / reference_km


def write_map(path, br):
    scale = float(np.abs(br).max())
    encoded = np.rint(127.0 * np.sign(br) * np.sqrt(np.abs(br) / scale)).astype(np.int8)
    height, width = br.shape
    path.write_bytes(f"AXFIELDMAP1 {width} {height} {scale:.3f}\n".encode("ascii") + encoded.tobytes())
    print(f"{path}: {width} x {height}, Br {br.min():.0f} .. {br.max():.0f} nT, rms {np.sqrt(np.mean(br ** 2)):.0f}")


def report_extreme(br, band_deg, what):
    height, width = br.shape
    latitude = 90.0 - (np.arange(height) + 0.5) * 180.0 / height
    longitude = -180.0 + (np.arange(width) + 0.5) * 360.0 / width
    rows = np.abs(latitude) <= band_deg
    sub = br[rows]
    i, j = np.unravel_index(np.argmin(sub), sub.shape)
    print(f"  {what}: minimum within {band_deg} deg of the equator {sub.min():.0f} nT, "
          f"parametric latitude {latitude[rows][i]:+.1f}, east longitude {longitude[j] % 360:.1f}")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    mars_model, out = sys.argv[1], Path(sys.argv[2])

    # Mars 150 km up: where the model is measured, and where it does not ring.
    g, h = load_gauss(mars_model)
    cos_t, sin_t, radius = ellipsoid_rows(360, 3396.19, 3376.20, 150.0, 3393.5)
    mars = radial_field(g, h, cos_t, sin_t, radius, 720)
    write_map(out / "mars_br.fieldmap", mars)
    i, j = np.unravel_index(np.argmax(np.abs(mars)), mars.shape)
    print(f"  peak at parametric latitude {89.75 - 0.5 * i:+.1f}, east longitude {(-179.75 + 0.5 * j) % 360:.1f}")

    # Jupiter at the 1-bar surface, and at the top of the dynamo: same directions from the
    # center, so that the two maps line up on the mesh.
    g, h = load_gauss(out / "field_jupiter_jrm33.txt")
    cos_t, sin_t, radius = ellipsoid_rows(180, 71492.0, 66854.0, 0.0, 71492.0)
    surface = radial_field(g, h, cos_t, sin_t, radius, 360)
    write_map(out / "jupiter_br.fieldmap", surface)
    report_extreme(surface, 15.0, "cloud tops")
    deep = radial_field(g, h, cos_t, sin_t, np.full_like(radius, 0.81), 360)
    write_map(out / "jupiter_br_deep.fieldmap", deep)
    report_extreme(deep, 15.0, "0.81 RJ")


if __name__ == "__main__":
    main()
