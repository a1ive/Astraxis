"""Two flares of YY Gem in March 1971, for castor.toml (standard library only).

Usage:
    python -I tools/castor/yy_gem_flares.py <directory of assets/stars>

Readings by eye of the figures of Moffett & Bopp 1971 (ApJ 168, L117; the ADS
scans), as thousands of detected photons in 2-second integrations:
  Fig. 1, 1971 March 24, U band: the peak is 1.17 mag above the star (their text);
  Fig. 3, 1971 March 28, B, V and U bands.

Writes for each the light of the flaring star relative to its quiescent light in
V, which is what the scene plots. YY Gem's stars are equal, so a flare adding s
of the pair's light adds 2 s of the star's. March 28 has V itself (between the
peaks and after them; the first peak is scaled from B). March 24 has only U: V
takes the ratio of the excesses in V and U measured at the main peak of March 28.

Bopp 1974 (ApJ 193, 389, Table 3 and Sect. V) puts these flares on the primary, in
one region at most 60 deg wide around longitude 80 deg (phase 0.222 at the meridian).
The patch: a 9000 K blackbody on a 3820 K star, the area giving the V excess,
seen at the angle the site has at the peak.
"""

import math
import sys

HC_OVER_K = 1.438776877e7  # nm K
LIMB = 0.6
T0, P = 2459690.57515, 0.81428247
SITE_MERIDIAN_PHASE = 80.0 / 360.0
SITE_LAT_DEG = 45.0            # assumed, the latitude of the spots of Torres & Ribas 2002
INCLINATION_DEG = 86.29

# Fig. 1, minutes from 04:03 UT, thousands of photons.
MARCH_24 = [
    (0.0, 3.15), (0.2, 3.55), (1.0, 3.60), (1.36, 3.69), (2.06, 4.67), (2.43, 5.19), (2.7, 4.96),
    (2.96, 5.94), (3.29, 6.67), (3.6, 6.40), (3.82, 7.20), (4.25, 6.98), (4.62, 6.57), (5.0, 6.17),
    (5.53, 5.63), (6.06, 5.77), (7.13, 5.30), (8.47, 5.21), (9.37, 4.84), (10.6, 4.87), (11.94, 4.61),
    (12.7, 4.70), (13.8, 4.32), (14.9, 4.59), (15.95, 4.40), (16.75, 3.98), (18.1, 3.75), (19.6, 3.52),
    (21.7, 3.31), (22.2, 3.24), (23.4, 3.19), (25.8, 3.05), (27.4, 2.97), (29.4, 3.17), (30.8, 3.10),
    (32.6, 2.88), (34.1, 2.83), (38.6, 2.74), (42.7, 2.76), (50.0, 2.72),
]
SKY_U = 0.5                    # the "sky" mark of Fig. 1
PEAK_U_MAG = 1.17              # text


def base_march_24(t):          # the star before and after the flare: 2.45 at 3:51, 2.72 at the end
    return 2.45 + (t + 12.0) * 0.27 / 62.0


# Fig. 3, minutes from 04:23 UT, thousands of photons.
MARCH_28_V = [(1.0, 54.06), (2.4, 53.5), (6.7, 51.7), (10.1, 52.9), (13.6, 54.5), (16.5, 53.5), (24.2, 51.7),
              (31.0, 51.3), (35.4, 51.0), (41.8, 50.7), (47.7, 50.4)]
BASE_V = 50.4                  # the end of the night's V curve
SKY_V = 3.3                    # as the "B sky" mark
BASE_B, B_FIRST_PEAK, B_MAIN_PEAK = 50.5, 61.0, 61.5
V_MAIN_PEAK = 54.5
U_MAIN_PEAK, BASE_U, SKY_U_28 = 14.7, 8.0, 1.0


def planck(t, wavelength_nm):
    return 1.0 / math.expm1(HC_OVER_K / (wavelength_nm * t))


def site_mu(phase):
    lat = math.radians(SITE_LAT_DEG)
    beta = math.radians(90.0 - INCLINATION_DEG)
    delta = 2.0 * math.pi * (phase - SITE_MERIDIAN_PHASE)
    return math.sin(lat) * math.sin(beta) + math.cos(lat) * math.cos(beta) * math.cos(delta)


def phase_of(jd):
    return ((jd - T0) / P) % 1.0


def write(path, comment, rows):
    with open(path, 'w', encoding='ascii', newline='\n') as f:
        f.write('# YY Gem flare, the flaring star\'s V light over its quiescent light. ' + comment + '\n')
        f.write('# tools/castor/yy_gem_flares.py. Minutes from the time in castor.toml.\n')
        f.write('minutes,flux\n')
        for t, value in rows:
            f.write('%.2f,%.4f\n' % (t, value))


def share(excess, phase):
    """sin^2 of the patch's radius (the scene's peak_disk_share for a site fixed on the star)."""
    ratio = planck(9000.0, 550.0) / planck(3820.0, 550.0)
    return excess / ratio * (1.0 - LIMB / 3.0) / max(site_mu(phase), 0.2)


def main(directory):
    v_u = ((V_MAIN_PEAK - BASE_V) / (BASE_V - SKY_V)) / ((U_MAIN_PEAK - BASE_U) / (BASE_U - SKY_U_28))
    print('V excess / U excess at the main peak of March 28: %.3f' % v_u)

    raw = [(t, (c - SKY_U) / (base_march_24(t) - SKY_U)) for t, c in MARCH_24]
    peak = max(f for _, f in raw)
    scale = (10.0 ** (0.4 * PEAK_U_MAG) - 1.0) / (peak - 1.0)
    print('March 24: read U peak %.2f (%.2f mag), scaled to %.2f mag' % (peak, 2.5 * math.log10(peak), PEAK_U_MAG))
    rows = [(-1.2, 1.0)] + [(t, 1.0 + 2.0 * v_u * scale * (f - 1.0)) for t, f in raw]
    rows = [(t, max(f, 1.0)) for t, f in rows]
    write(directory + '/yy_gem_flare_1971_03_24.csv', 'Moffett & Bopp 1971, Fig. 1 (U), V from March 28.', rows)
    excess = max(f for _, f in rows) - 1.0
    t_peak = max(rows, key=lambda r: r[1])[0]
    jd_peak = 2441034.5 + (4.0 + 3.0 / 60.0 + t_peak / 60.0) / 24.0
    print('March 24: star V peak %.3f at +%.2f min, pair +%.3f mag; phase at the peak %.3f, at 04:09 UT %.3f (paper 0.154)'
          % (1.0 + excess, t_peak, 2.5 * math.log10(1.0 + excess / 2.0), phase_of(jd_peak),
             phase_of(2441034.5 + (4.0 + 9.0 / 60.0) / 24.0)))
    print('  peak_disk_share %.4f (site mu %.2f)' % (share(excess, phase_of(jd_peak)), site_mu(phase_of(jd_peak))))

    v = [(t, (c - BASE_V) / (BASE_V - SKY_V)) for t, c in MARCH_28_V]
    first_peak = (B_FIRST_PEAK - BASE_B) / (BASE_B - SKY_V) * ((V_MAIN_PEAK - BASE_V) / (B_MAIN_PEAK - BASE_B))
    rows = [(-1.0, 1.0), (0.6, 1.0 + 2.0 * first_peak)] + [(t, 1.0 + 2.0 * max(e, 0.0)) for t, e in v]
    write(directory + '/yy_gem_flare_1971_03_28.csv', 'Moffett & Bopp 1971, Fig. 3 (V; the first peak from B).', rows)
    excess = max(f for _, f in rows) - 1.0
    t_peak = max(rows, key=lambda r: r[1])[0]
    jd_peak = 2441038.5 + (4.0 + 23.0 / 60.0 + t_peak / 60.0) / 24.0
    print('March 28: star V peak %.3f at +%.2f min, pair +%.3f mag; phase at 04:43 UT %.3f (paper 0.095)'
          % (1.0 + excess, t_peak, 2.5 * math.log10(1.0 + excess / 2.0), phase_of(2441038.5 + (4.0 + 43.0 / 60.0) / 24.0)))
    print('  peak_disk_share %.4f (site mu %.2f)' % (share(excess, phase_of(jd_peak)), site_mu(phase_of(jd_peak))))


main(sys.argv[1])
