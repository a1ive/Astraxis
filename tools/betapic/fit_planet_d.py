"""Fit a representative Jacobi orbit for Beta Pic d (Sutlieff + Gibbs 2026).

Run with NumPy and SciPy. Fix b and c to the scene's Lacour 2021 orbits;
subtract the star's reflex motion before projecting the new orbit onto the sky.
The semi-major axis is fixed to 32 au, in the overlap between Sutlieff Table 5
and the region favored by Gibbs's stability tests. This is a conditional
least-squares fit, not a new estimate of the orbital posterior.
"""

from datetime import date
from pathlib import Path
import tomllib

import numpy as np
from scipy.optimize import least_squares


AU_KM = 149597870.7
DAY_S = 86400.0
GM_D = 2.4 * 126686531.9
EPOCH = date(2000, 1, 1)
PARALLAX = 51.44  # mas [Lacour et al. 2021]
A_D = 32.0


def mjd(day):
    return 51544 + (date.fromisoformat(day) - EPOCH).days


def state(a, e, inc, node, arg, mean, gm):
    """Kepler state in north, east and away axes, in au and km/s."""
    anomaly = np.remainder(mean + np.pi, 2 * np.pi) - np.pi
    ecc = anomaly.copy()
    for _ in range(8):
        ecc -= (ecc - e * np.sin(ecc) - anomaly) / (1 - e * np.cos(ecc))
    cos_e, sin_e = np.cos(ecc), np.sin(ecc)
    cos_n, sin_n = np.cos(node), np.sin(node)
    p = np.array([cos_n, sin_n, 0.0])
    q = np.array([-sin_n * np.cos(inc), cos_n * np.cos(inc), np.sin(inc)])
    peri = np.cos(arg) * p + np.sin(arg) * q
    side = -np.sin(arg) * p + np.cos(arg) * q
    b = np.sqrt(1 - e * e)
    position = a * ((cos_e - e)[:, None] * peri + (b * sin_e)[:, None] * side)
    speed = np.sqrt(gm / (a * AU_KM)) / (1 - e * cos_e)
    velocity = speed[:, None] * (-sin_e[:, None] * peri + (b * cos_e)[:, None] * side)
    return position, velocity


def reflex(orbit, days):
    a = orbit["a_au"]
    gm = orbit["gm_km3_s2"]
    n = np.sqrt(gm / (a * AU_KM)**3) * DAY_S
    mean = n * (days + 2400000.5 - orbit["t_peri_jd_tdb"])
    angles = np.radians([orbit[key] for key in ("i_deg", "node_deg", "arg_peri_deg")])
    return np.array(state(a, orbit["e"], *angles, mean, gm)) * orbit["scale"]


def main():
    scene = Path(__file__).resolve().parents[2] / "assets/scenes/beta_pictoris.toml"
    bodies = {body["name"]: body for body in tomllib.loads(scene.read_text(encoding="utf-8"))["bodies"]}
    inner = bodies["Beta Pic A+c"]["orbit"]
    star = bodies["Beta Pic"]["orbit"]
    gm = inner["gm_km3_s2"] + GM_D
    n = np.sqrt(gm / (A_D * AU_KM)**3) * DAY_S
    ref_mjd = mjd("2025-12-03")

    # Sutlieff Table 1: separation, position angle, and their independent errors.
    imaging = [
        ("2014-12-08", 375, 208.0, 20, 2.0),
        ("2019-03-11", 719, 209.6, 12, 1.0),
        ("2020-02-08", 784, 209.4, 15, 1.0),
        ("2023-03-18", 1012, 210.4, 20, 1.0),
        ("2025-03-22", 1104, 211.1, 20, 1.0),
        ("2025-12-03", 1126.9, 210.35, 3.9, 0.22),
    ]
    # Gibbs Table 1: astrometry and absolute RV; subtract the stellar 20 km/s
    # and add its 0.7 km/s error in quadrature. Use symmetric RV errors.
    spectra = [
        ("2025-11-29", -1000.0, -528.0, 22.2, 10.7, -5.3, np.hypot(2.0, 0.7)),
        ("2026-03-29", -1001.7, -525.5, 22.7, 30.8, -7.4, np.hypot(2.7, 0.7)),
    ]
    dates = [r[0] for r in imaging + spectra]
    days = np.array([mjd(day) for day in dates], dtype=float)
    offsets, speeds = reflex(inner, days) + reflex(star, days)
    sep, pa, sep_err, pa_err = np.array([row[1:] for row in imaging]).T
    north, east, north_err, east_err, rv, rv_err = np.array([row[1:] for row in spectra]).T
    split = len(imaging)

    def prediction(elements):
        e, inc, node, arg, phase = elements
        position, velocity = state(A_D, e, inc, node, arg, phase + n * (days - ref_mjd), gm)
        return (position - offsets) * PARALLAX, velocity - speeds

    def residuals(elements):
        position, velocity = prediction(elements)
        x, y = position[:split, :2].T
        angle = np.degrees(np.arctan2(y, x)) % 360
        return np.concatenate(((np.hypot(x, y) - sep) / sep_err, (angle - pa) / pa_err,
                               (position[split:, 0] - north) / north_err,
                               (position[split:, 1] - east) / east_err,
                               (velocity[split:, 2] - rv) / rv_err))

    fits = []
    # Imaging alone has a 180-degree ambiguity in the ascending node.
    for node in (31, 211):
        bounds = ([0, np.radians(85), np.radians(node - 10), -np.inf, -np.inf],
                  [0.4, np.radians(93), np.radians(node + 10), np.inf, np.inf])
        for arg in np.linspace(0, 2 * np.pi, 8, endpoint=False):
            for phase in (-1.0, 0.0, 1.0):
                fits.append(least_squares(residuals, [0.15, *np.radians([89, node]), arg, phase],
                                         bounds=bounds))
    fit = min(fits, key=lambda result: result.cost)
    e, inc, node, arg, phase = fit.x
    print(f"a_au={A_D:g}, e={e:.10f}, i_deg={np.degrees(inc):.8f}, node_deg={np.degrees(node):.8f}")
    print(f"arg_peri_deg={np.degrees(arg) % 360:.8f}, t_peri_jd_tdb={ref_mjd + 2400000.5 - (phase % (2*np.pi))/n:.6f}")
    print(f"gm_km3_s2={gm:.4f}, period={2*np.pi/n/365.25:.3f} yr, chi-squared={2*fit.cost:.2f}")
    position, velocity = prediction(fit.x)
    for day, (x, y, _), (_, _, v) in zip(dates, position, velocity):
        print(f"{day}: north={x:8.2f} mas east={y:7.2f} mas RV={v:+.2f} km/s")


if __name__ == "__main__":
    main()
