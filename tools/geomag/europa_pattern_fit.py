"""Fit the pattern of Europa's aurora (scene/satellite_aurora.hpp) to Roth et al. 2016.

Usage: python -E -P tools/geomag/europa_pattern_fit.py

Forward model: an optically thin shell (the scene's height profile) seen from outside at
the geometries of the paper's images, through a small PSF, measured over the regions of
its Figures 7 and 8. Prints the fitted parameters and the ratios against the paper's.
Needs numpy. Source and assumptions: assets/geomag/SOURCES.md.

Body axes: x toward Jupiter, y toward the trailing hemisphere, z north. Sunlight comes from
the observer (HST's phase angle is under 11 degrees).
"""
import math

import numpy as np

SCALE_KM = 150.0  # Roth 2016 section 4.2.4: the favored height for OI 1356 A
TOP_KM = 600.0
RADIUS_KM = 1560.8
PSF_SIGMA_RE = 0.08  # about the core of STIS's PSF over a one arcsecond disk
GRID = 48
EXTENT = 1.5  # image half width, RE

# The dipole of Roth's Figure 7 is VIPAL's 9.8 degrees; Europa's largest tilt of the field
# is two of these on the observer's sky.
TILT_DEG = 9.5


def field_north_end(psi_deg, azimuth_deg):
    """The north end of the dipole line through Europa, at magnetic latitude psi. In the
    body axes the moon at System III longitude sees the dipole axis tilted by the rest."""
    # The dipole axis's component toward the moon is sin(psi).
    tilt = math.radians(TILT_DEG)
    psi = math.radians(psi_deg)
    toward_moon = math.sin(psi)
    trailing = math.sqrt(max(math.sin(tilt) ** 2 - toward_moon ** 2, 0.0)) * azimuth_deg
    # m has -toward_moon along x (x points at Jupiter); b = m - 3 sin(psi) r.
    b = np.array([-toward_moon + 3.0 * toward_moon, trailing, math.cos(tilt)])
    return b / np.linalg.norm(b)


def pattern(n, b, sheet, dusk, p):
    """europa_pattern in scene/satellite_aurora.hpp, before its normalization."""
    across, sheet_gain, upstream, tilt = p
    return (1.0 + across * (n @ dusk)) * np.exp(sheet * n[..., 2] * (sheet_gain + upstream * n[..., 1]) +
                                                tilt * ((n @ b) ** 2 - n[..., 2] ** 2))


def image(view, b, sheet, p):
    """The brightness over the sky plane, north up."""
    w = np.array(view, float)
    up = np.array([0.0, 0.0, 1.0])
    right = np.cross(up, w)
    axis = (np.arange(GRID) + 0.5) / GRID * 2.0 * EXTENT - EXTENT
    u, v = np.meshgrid(axis, axis, indexing="xy")
    radius_top = 1.0 + TOP_KM / RADIUS_KM
    r2 = u * u + v * v
    depth_steps = 32
    # Toward the observer from the far side of the shell, stopping at the surface.
    far = np.where(r2 < 1.0, np.sqrt(np.maximum(1.0 - r2, 0.0)), -np.sqrt(np.maximum(radius_top ** 2 - r2, 0.0)))
    near = np.sqrt(np.maximum(radius_top ** 2 - r2, 0.0))
    fractions = (np.arange(depth_steps) + 0.5) / depth_steps
    depth = far[None] + (near - far)[None] * fractions[:, None, None]
    points = u[None, ..., None] * right + v[None, ..., None] * up + depth[..., None] * w
    distance = np.linalg.norm(points, axis=-1)
    n = points / distance[..., None]
    height = (distance - 1.0) * RADIUS_KM
    emission = np.exp(-np.maximum(height, 0.0) / SCALE_KM) * pattern(n, b, sheet, right, p)
    emission[height > TOP_KM] = 0.0
    picture = emission.sum(0) * (near - far) * RADIUS_KM / depth_steps / SCALE_KM
    picture[r2 >= radius_top ** 2] = 0.0
    # The PSF, by FFT.
    step = 2.0 * EXTENT / GRID
    frequency = np.fft.fftfreq(GRID, step)
    fx, fy = np.meshgrid(frequency, frequency, indexing="xy")
    blur = np.exp(-2.0 * (math.pi * PSF_SIGMA_RE) ** 2 * (fx * fx + fy * fy))
    return np.fft.ifft2(np.fft.fft2(picture) * blur).real, u, v


def measures(view, b, sheet, p):
    picture, u, v = image(view, b, sheet, p)
    inside = u * u + v * v < 1.25 ** 2
    mean = lambda mask: picture[mask & inside].mean()
    north, south = v > 0.5, v < -0.5
    # Figure 7a's quadrants 1 are upper left and lower right.
    quadrant_1 = ((u < 0) & (v > 0)) | ((u > 0) & (v < 0))
    return {
        "north_south": mean(north) / mean(south),
        "right_left": mean(u > 0) / mean(u < 0),
        "poles_equator": mean(north | south) / mean(~(north | south)),
        "quadrants": mean(quadrant_1) / mean(~quadrant_1),
        "disk": picture[inside].sum() * (2.0 * EXTENT / GRID) ** 2 / math.pi,
    }


TRAILING, LEADING, SUBJOVIAN = (0, 1, 0), (0, -1, 0), (1, 0, 0)
# Above the sheet at psi = +9.5, field tilted; sheet = -1 makes the south brighter.
ABOVE = (field_north_end(TILT_DEG, 0.0), -1.0)
BELOW = (field_north_end(-TILT_DEG, 0.0), 1.0)
CROSSING = (field_north_end(0.0, 1.0), 0.0)


def residuals(p):
    out = []
    # Roth's numbers, read from the figures. The bins' spread is large: ln weights say how
    # far to trust each, and the eclipse panel has eight points of short exposures.
    def add(value, target, weight):
        out.append(weight * math.log(value / target))
    # Figure 7c: trailing, 0-90 bins about 3.6, 180-270 about 0.26.
    add(measures(TRAILING, *BELOW, p)["north_south"], 3.6, 1.0)
    add(measures(TRAILING, *ABOVE, p)["north_south"], 0.26, 1.0)
    # Figure 7d: leading, about 1.8 and 0.7.
    add(measures(LEADING, *BELOW, p)["north_south"], 1.8, 1.0)
    add(measures(LEADING, *ABOVE, p)["north_south"], 0.7, 1.0)
    # Figure 7e: sub-Jovian, inverted, about 0.33 and 3.
    add(measures(SUBJOVIAN, *BELOW, p)["north_south"], 0.33, 0.25)
    add(measures(SUBJOVIAN, *ABOVE, p)["north_south"], 3.0, 0.25)
    # Figure 7a: |beta| = 20 degrees, about 1.6 (trailing and leading images together).
    for view in (TRAILING, LEADING):
        q = measures(view, *BELOW, p)["quadrants"]
        add(max(q, 1.0 / q), 1.6, 1.0)
    # Figure 8a: 1.59 (0.33), all visits out of eclipse.
    for view in (TRAILING, LEADING):
        add(measures(view, *CROSSING, p)["right_left"], 1.59, 1.5)
    return np.array(out)


def fit(start):
    p = np.array(start, float)
    step = np.full(len(p), 0.2)
    cost = (residuals(p) ** 2).sum()
    for _ in range(25):
        improved = False
        for k in range(len(p)):
            for sign in (1.0, -1.0):
                trial = p.copy()
                trial[k] = max(trial[k] + sign * step[k], 0.0)
                c = (residuals(trial) ** 2).sum()
                if c < cost:
                    p, cost, improved = trial, c, True
                    break
        if not improved:
            step *= 0.5
            if step.max() < 0.01:
                break
    return p, cost


def report(p):
    print(f"cost {(residuals(p) ** 2).sum():.3f}")
    print("across %.3f sheet %.3f upstream %.3f tilt %.3f" % tuple(p))
    for label, view in (("trailing", TRAILING), ("leading", LEADING), ("sub-Jovian", SUBJOVIAN)):
        for name, state in (("below", BELOW), ("above", ABOVE), ("crossing", CROSSING)):
            m = measures(view, *state, p)
            print(f"  {label:10} {name:8} N/S {m['north_south']:.2f} R/L {m['right_left']:.2f} "
                  f"poles/eq {m['poles_equator']:.2f} quadrants {m['quadrants']:.2f} disk {m['disk']:.3f}")


def main():
    p, _ = fit([0.45, 0.6, 0.3, 1.0])
    report(p)


if __name__ == "__main__":
    main()
