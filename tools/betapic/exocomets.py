#!/usr/bin/env python3
"""Schematic exocomets for assets/scenes/beta_pictoris.toml, from TESS's transits.

Lecavelier des Etangs et al. 2022 (Sci. Rep. 12, 5855, arXiv:2204.13618) list 30
transits (Extended Data Table 2: the time t0, the duration dt, the absorption depth
AD) and turn them into a comet each:
  periastron    dt = 13 h sqrt(q / au)
  nucleus       R  = 1.5 km sqrt(AD / 1e-4)
What a transit does not say is assumed here, so the orbits are schematic: e = 0.99
(Heller 2024, A&A 691, arXiv:2407.20829, highly eccentric, a little less than the
hyperbolas he tries); the orbit in the disk's plane, within a degree or so, tilted
just enough to cross the star at the impact parameter below; and where on its orbit
the comet crosses the line of sight (Lecavelier: the depth does not depend on it
over a wide range).

The transits chosen are the run of Dec 2018 - Jan 2019 in TESS's Sector 6, with the
deepest one. The total magnitude is scaled from the comet Hale-Bopp as they scale the
dust production, M ~ R^2 L, then dimmed and flattened for the display (a comet this
near an A star would otherwise fill the view).

Prints the [[bodies]] blocks. Standard library only: python -I tools/betapic/exocomets.py
"""

import math

GM_STAR = 2.32246770e11  # km^3 s^-2, 1.75 M_sun [LA21]
AU = 149597870.7  # km
R_STAR_AU = 1069987.0 / AU  # 1.538 R_sun [ZW19]
DISK_NODE = 29.3  # deg, the main disk [MA19]
E = 0.99

# Hale-Bopp, for the magnitudes [LE22]: 30 km; JPL SBDB's M1 of about -1. Display: dimmed
# by DISPLAY_DIM, K1 smaller than a real comet's 10.
HALE_BOPP_R_KM, HALE_BOPP_M1 = 30.0, -1.0
DISPLAY_DIM, K1 = 8.0, 5.0
L_STAR = 8.7  # [LE22]

# [LE22] Extended Data Table 2: ID, t0 (BJD - 2457000), dt (days), AD (1e-4). Then the
# choices of this tool: the true anomaly at transit (deg) and the impact parameter.
TRANSITS = [
    (5, 1458.80, 0.34, 5.65, -40.0, 0.2),
    (7, 1468.49, 0.27, 1.97, 25.0, -0.5),
    (8, 1474.58, 0.12, 1.77, -15.0, 0.3),
    (9, 1476.34, 0.28, 3.76, 35.0, 0.6),
    (10, 1484.31, 0.20, 1.63, -30.0, -0.2),
    (11, 1485.72, 0.62, 19.63, 10.0, 0.1),
]


def main():
    for ident, t0, dt, depth, nu_deg, impact in TRANSITS:
        q = (dt * 24.0 / 13.0) ** 2
        radius = 1.5 * math.sqrt(depth)
        a = q / (1.0 - E)
        nu = math.radians(nu_deg)
        r_transit = a * (1.0 - E * E) / (1.0 + E * math.cos(nu))
        # On the near side the on-sky offset is r cos(i); i < 90 like the planets.
        i_deg = 90.0 - math.degrees(math.asin(impact * R_STAR_AU / r_transit))
        # On the near side argument = omega + nu = 270 deg.
        omega = (270.0 - nu_deg) % 360.0
        ecc = 2.0 * math.atan2(math.sqrt(1.0 - E) * math.sin(nu / 2.0), math.sqrt(1.0 + E) * math.cos(nu / 2.0))
        mean = ecc - E * math.sin(ecc)
        n = math.sqrt(GM_STAR / (a * AU) ** 3)
        t_peri = 2457000.0 + t0 - mean / n / 86400.0
        m1 = HALE_BOPP_M1 + DISPLAY_DIM - 2.5 * math.log10((radius / HALE_BOPP_R_KM) ** 2 * L_STAR)
        print(f"""[[bodies]]
name = "Exocomet {ident}"
parent = "Beta Pic"
radii_km = [{radius:.2f}]                 # [LE22] from the depth {depth}e-4
color = "#d8d0c0"
orbit_color = "#b8c8e8"
trail = "orbit"
[bodies.comet]
magnitude_m1_k1 = [{m1:.2f}, {K1}]        # scaled from Hale-Bopp, see the tool
ion_tail_diameter_km = 2.0e5               # illustrative
ion_tail_length_km = 2.0e7                 # illustrative
solar_wind_km_s = 400.0                    # illustrative
dust_beta_range = [2.0e-3, 0.1]            # illustrative
dust_tail_days = 15.0                      # illustrative
dust_speed_km_s = 0.3                      # illustrative
[bodies.orbit]                             # q = {q:.3f} au from dt = {dt} d [LE22]; e and the angles are assumed
type = "visual_orbit"
gm_km3_s2 = {GM_STAR:.8e}
a_au = {a:.4f}
e = {E}
i_deg = {i_deg:.4f}
node_deg = {DISK_NODE}
arg_peri_deg = {omega:.1f}
t_peri_jd_tdb = {t_peri:.4f}               # crossing the star at BJD 2457000 + {t0}
sky_ra_dec_deg = [86.82118, -51.06651]
""")


if __name__ == "__main__":
    main()
