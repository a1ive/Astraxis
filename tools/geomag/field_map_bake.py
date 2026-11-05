"""Bake the field maps of [bodies.magnetic_map]; see assets/geomag/SOURCES.md.

    python -E -P tools/geomag/field_map_bake.py <Morschhauser2014.txt> <T2015_449.sh.gz> assets/geomag

Needs numpy. The output maps the body's ellipsoid like its texture: rows from north
to south by parametric latitude (the body mesh's), columns from 180 W eastward. Radial
maps: a header line "AXFIELDMAP1 width height scale_nT", then width * height signed
bytes, 127 * sign(Br) * sqrt(|Br| / scale): the square root keeps the weak fields of
a map whose peak is a hundred times its typical value. Strength maps: "AXFIELDMAG1
width height min_nT max_nT", then bytes 127 * (2 t - 1), t = (|B| - min) / (max - min).
"""

import gzip
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


def load_pyshtools(path, nmax):
    """g, h in nT from a gzipped pyshtools file ("radius, degree" then "n m c s" in tesla)."""
    g = np.zeros((nmax + 1, nmax + 1))
    h = np.zeros((nmax + 1, nmax + 1))
    with gzip.open(path, "rt", encoding="ascii") as lines:
        next(lines)
        for line in lines:
            parts = line.split()
            n, m = int(parts[0]), int(parts[1])
            if 1 <= n <= nmax:
                g[n, m] = float(parts[2]) * 1e9
                h[n, m] = float(parts[3]) * 1e9
    return g, h


def total_field(g, h, cos_theta, sin_theta, radius, width):
    """|B| (nT), as radial_field: from the three components, for a low degree."""
    degree = g.shape[0] - 1
    longitude = (np.arange(width) + 0.5) / width * 2.0 * np.pi - np.pi
    b_r = np.zeros((len(cos_theta), width))
    b_t = np.zeros_like(b_r)
    b_p = np.zeros_like(b_r)
    before = {}
    last = {0: np.ones(len(cos_theta))}
    diagonal = np.ones(len(cos_theta))
    for n in range(1, degree + 1):
        current = {}
        for m in range(n):
            above = np.sqrt((n - 1) ** 2 - m * m) * before[m] if m in before else 0.0
            current[m] = ((2 * n - 1) * cos_theta * last[m] - above) / np.sqrt(n * n - m * m)
        diagonal = diagonal * sin_theta * (1.0 if n == 1 else np.sqrt((2 * n - 1) / (2 * n)))
        current[n] = diagonal
        factor = radius ** -(n + 2)
        for m in range(n + 1):
            below = last[m] if m < n else 0.0
            derivative = (n * cos_theta * current[m] - np.sqrt(n * n - m * m) * below) / sin_theta
            cos_m, sin_m = np.cos(m * longitude), np.sin(m * longitude)
            along = np.outer(np.ones(len(cos_theta)), g[n, m] * cos_m + h[n, m] * sin_m)
            across = np.outer(np.ones(len(cos_theta)), g[n, m] * sin_m - h[n, m] * cos_m)
            b_r += ((n + 1) * factor * current[m])[:, None] * along
            b_t -= (factor * derivative)[:, None] * along
            b_p += (m * factor * current[m] / sin_theta)[:, None] * across
        before, last = last, current
    return np.sqrt(b_r ** 2 + b_t ** 2 + b_p ** 2)


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


def write_magnitude_map(path, b):
    low, high = float(b.min()), float(b.max())
    encoded = np.rint(127.0 * (2.0 * (b - low) / (high - low) - 1.0)).astype(np.int8)
    height, width = b.shape
    path.write_bytes(f"AXFIELDMAG1 {width} {height} {low:.3f} {high:.3f}\n".encode("ascii") + encoded.tobytes())
    i, j = np.unravel_index(np.argmin(b), b.shape)
    print(f"{path}: {width} x {height}, |B| {low:.0f} .. {high:.0f} nT, weakest at parametric latitude "
          f"{90.0 - (i + 0.5) * 180.0 / height:+.1f}, east longitude {(-180.0 + (j + 0.5) * 360.0 / width) % 360:.1f}")


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
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    mars_model, moon_model, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])

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

    # The Moon 30 km up, near the altitude of the low orbits that measured it. At the surface
    # the expansion of the gridded model shows north-south stripes around the anomalies; going
    # up removes them, as it does the small scales the data do not hold. Degrees above 240
    # are then below 1% of their surface value.
    g, h = load_pyshtools(moon_model, 240)
    cos_t, sin_t, radius = ellipsoid_rows(512, 1737.4, 1737.4, 30.0, 1737.4)
    moon = radial_field(g, h, cos_t, sin_t, radius, 1024)
    write_map(out / "moon_br.fieldmap", moon)
    i, j = np.unravel_index(np.argmax(np.abs(moon)), moon.shape)
    print(f"  peak at latitude {90.0 - (i + 0.5) * 180.0 / 512:+.1f}, east longitude {(-180.0 + (j + 0.5) * 360.0 / 1024) % 360:.1f}")

    # The Earth: the strength, not Br, which only shows the dipole's two halves.
    g, h = load_gauss(out / "field_earth_igrf14.txt")
    cos_t, sin_t, radius = ellipsoid_rows(180, 6378.1366, 6356.7519, 0.0, 6371.2)
    write_magnitude_map(out / "earth_b.fieldmap", total_field(g, h, cos_t, sin_t, radius, 360))


if __name__ == "__main__":
    main()
