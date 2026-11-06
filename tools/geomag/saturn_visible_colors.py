"""Colors and layer ratio of Saturn's visible aurora (prints the constants of
src/scene/saturn_aurora.hpp).

    python -I saturn_visible_colors.py

Dyudina et al. 2015 (Icarus, arXiv:1506.00664) give the spectrum only as a plot.
FIG8 holds what was read from their Fig. 8 at 300 dpi, by locating the error bar
caps of each curve (+-10% or so; the axis is logarithmic, 431 px per decade).
Values are R/nm relative to the clear filter.

Color follows their Appendix A.2: the RED, GRN and BL1 brightnesses in I/F units
are taken as the true color. That is also what the renderer's linear RGB means
(1 = a white Lambertian surface in sunlight, aurora_colors.py), and it needs no
guess at the spectrum between the filters: building a spectrum from the filters
and integrating CIE functions moved the high aurora's hue from blue to magenta
depending on how the continuum was interpolated. Going from R/nm to I/F divides by
the sun's photons per nm at the filter center (5772 K black body, the white the
renderer uses; the real sun's lines change the blue by some 15%).

The layers' ratio is the boxes' clear-filter brightness (labels in Fig. 8, kR):
edge-on, the high layer at its peak height is 11 / 150 of the low layer at its
peak. The high box runs along dimmer longitudes than the bright patch under the
other two, so this is a lower limit on the high layer, by a factor of a few.

Standard library only.
"""

import math

H = 6.62607015e-34  # J s
C = 299792458.0     # m/s
K = 1.380649e-23    # J/K
T_SUN = 5772.0      # K, IAU 2015 B3 nominal

BANDS = {'RED': 650.0, 'GRN': 568.0, 'BL1': 451.0}  # filter centers, nm (Fig. 8 labels and curves)

FIG8 = {  # clear kR, RED, GRN, BL1
    'deep': (120.0, 1.03, 0.41, 0.42),
    'brightest': (150.0, 1.00, 0.48, 0.54),
    'high': (11.0, 0.93, 0.41, 1.08),
}

# Layers of aurora_model.cpp: Chapman peak and scale height (km), both assumed.
LOW = (450.0, 140.0)
HIGH = (1100.0, 180.0)
RADIUS_KM = 54900.0  # 1-bar radius near 70 deg latitude (radii_km of saturn.toml)


def sun_photons(nm):
    """Photons per nm of a black body, up to a constant."""
    lam = nm * 1e-9
    return 1.0 / lam ** 4 / (math.exp(H * C / (lam * K * T_SUN)) - 1.0)


def color_per_clear(box):
    """I/F in (RED, GRN, BL1) per unit clear-filter brightness, with RED's sun as 1."""
    red = sun_photons(BANDS['RED'])
    return [box[1], box[2] * red / sun_photons(BANDS['GRN']), box[3] * red / sun_photons(BANDS['BL1'])]


def chapman(h, peak, scale):
    z = (h - peak) / scale
    return math.exp(-z - math.exp(-z)) / scale


def tangent(height, layer):
    """Unit-column layer seen edge-on at this tangent height."""
    b = RADIUS_KM + height
    total = 0.0
    step = 5.0
    s = step / 2
    while s < 12000.0:
        total += 2.0 * chapman(math.hypot(b, s) - RADIUS_KM, *layer) * step
        s += step
    return total


def main():
    for name, box in FIG8.items():
        rgb = color_per_clear(box)
        print(f'{name:10} per clear kR {[round(v, 3) for v in rgb]}, chromaticity {[round(v / max(rgb), 3) for v in rgb]}')

    target = FIG8['high'][0] / FIG8['brightest'][0]
    low_at_low = max(tangent(h, LOW) for h in range(200, 700, 10))
    low_at_high = tangent(HIGH[0], LOW)
    high_at_low = tangent(LOW[0], HIGH)
    high_at_high = tangent(HIGH[0], HIGH)
    column = (target * low_at_low - low_at_high) / (high_at_high - target * high_at_low)
    print(f'edge-on per unit column: low {low_at_low:.3f} at its peak, {low_at_high:.4f} at {HIGH[0]:.0f} km; high {high_at_high:.3f}')
    print(f'high / low column for {target:.4f}: {column:.4f}')

    low = color_per_clear(FIG8['brightest'])
    high = color_per_clear(FIG8['high'])
    print('\n// saturn_aurora.hpp')
    print(f'kVisibleLowRgb = {{{low[0]:.3f}, {low[1]:.3f}, {low[2]:.3f}}}')
    print(f'kVisibleHighRgb = {{{high[0]:.3f}, {high[1]:.3f}, {high[2]:.3f}}}')
    print(f'kVisibleHighColumn = {column:.3f}')


if __name__ == '__main__':
    main()
