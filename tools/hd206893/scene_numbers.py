#!/usr/bin/env python3
"""The derived numbers of assets/scenes/hd206893.toml.

The disk: Marino et al. 2020 (MNRAS 498, 1319, arXiv:2010.12582), Table 2: surface
density ~ r^0.8 from 28 to 114 au, r^-3.9 beyond it (their fixed outer edge at 250
au), times a Gaussian gap at 74 au, 27 au FWHM, empty at its center. Their inner edge
is a sharp cut at 28 au, poorly known (14 to 40 au), so here it falls off as a Gaussian
of 3 au, and likewise their fixed outer edge at 250 au, by 12 au. Prints radial_au and column, peak 1, and how that looks at 0.88 mm against
Sect. 2.1.1.

The orbits: Kral et al. 2025 (arXiv:2511.20091), Table 2, Jacobi coordinates: c around
the star, B around the star and c. tau in periods after MJD 59000: this is not stated,
but it alone fits their Table 1 (c to 0.3 mas, B to 1.7 mas). Prints the GMs, the
barycenter's scales and the periastrons. Also the offset of the observer that stands
toward Earth.

Standard library only: python -I tools/hd206893/scene_numbers.py
"""

import math

GM_SUN = 1.32712440041279419e11  # JPL DE440
GM_JUP = 126686531.9  # JUP365
AU_KM = 149597870.7
YEAR_S = 365.25 * 86400.0
RA, DEC = 326.34126, -12.78335  # assets/stars/hyg.csv, HIP 107412

STAR, B, C, MOON = 1.31 * GM_SUN, 19.5 * GM_JUP, 11.1 * GM_JUP, 0.4 * GM_JUP


def period_days(a_au, gm):
    return 2.0 * math.pi * math.sqrt((a_au * AU_KM) ** 3 / gm) / 86400.0


def column(r):
    r_min, r_c, r_out = 28.0, 114.0, 250.0
    edge = 1.0
    if r > r_out:
        edge = math.exp(-0.5 * ((r - r_out) / 12.0) ** 2)
    if r < r_min:
        edge = math.exp(-0.5 * ((r_min - r) / 3.0) ** 2)
        r = r_min
    power = edge * (r / r_c) ** (0.8 if r < r_c else -3.9)
    sigma = 27.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    return max(power * (1.0 - math.exp(-0.5 * ((r - 74.0) / sigma) ** 2)), 1e-6)


def main():
    radii = [5.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0]
    r = 30.0
    while r < 250.0:
        radii.append(r)
        r += 2.0 if 62.0 <= r < 86.0 else 4.0 if r < 130.0 else 12.0
    radii += [250.0, 256.0, 262.0, 268.0, 274.0, 280.0, 290.0, 300.0]
    print("radial_au = [" + ", ".join(f"{x:g}" for x in radii) + "]")
    print("column = [" + ", ".join(f"{column(x):.4g}" for x in radii) + "]")
    # Brightness at 0.88 mm: column * T, T ~ r^-1/2.
    profile = [(x, column(x) * x**-0.5) for x in [20 + i for i in range(200)]]
    peak = max(v for _, v in profile)
    print("0.88 mm peak at", max(profile, key=lambda p: p[1])[0], "au;", end=" ")
    inner = max((p for p in profile if p[0] < 60), key=lambda p: p[1])
    print("inner peak at", inner[0], "au,", round(inner[1] / peak, 2), "of it")
    # The dust's blackbody temperature: 278 K L^(1/4) / sqrt(r), L = 2.9 L_sun [MA20].
    print("T at 100 au", 278.3 * 2.9**0.25 / 10.0)

    pair = STAR + C
    total = pair + B
    p_b = period_days(10.75, total)
    p_c = period_days(3.74, pair)
    print(f"star {STAR:.8e}, B {B:.8e}, c {C:.8e}, moon {MOON:.8e}")
    print(f"M_star + M_c {pair:.8e}, all {total:.8e}, B + moon {B + MOON:.8e}")
    print(f"scales: A+c {-B / total:.9f}, B {pair / total:.9f}, star {-C / pair:.9f}, c {STAR / pair:.9f}")
    print(f"P_B {p_b / 365.25:.3f} yr, P_c {p_c / 365.25:.3f} yr, moon {period_days(0.22, B + MOON) / 365.25:.3f} yr")
    print(f"t_peri B JD {2459000.5 + 0.133 * p_b:.3f}, c JD {2459000.5 + 0.694 * p_c:.3f}")

    # Toward Earth, 4000 au.
    ra, dec = math.radians(RA), math.radians(DEC)
    away = (math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec))
    d = 4000.0 * AU_KM
    print("observer_offset_km = [" + ", ".join(f"{-c * d:.1f}" for c in away) + "]")


if __name__ == "__main__":
    main()
