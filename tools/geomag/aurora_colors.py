"""Colors of the aurora's emission lines in the renderer's units (prints the
constants of src/scene/aurora.cpp).

    python aurora_colors.py <ciexyz31_1.csv>

The CIE 1931 2-degree color matching functions at 1 nm come from CVRL
(http://cvrl.ucl.ac.uk/database/data/cmfs/ciexyz31_1.csv). The renderer's
linear RGB is I/F: 1 is a white Lambertian surface lit by the sun at normal
incidence, and the sun is white (the eye adapts to it). So a line of 1 kR
(1e13 / (4 pi) photons s^-1 m^-2 sr^-1) is converted to XYZ with the CMFs,
divided by the luminance of that white surface, E_sun / pi, and white-balanced
channel by channel against sunlight (as the renderer does for a second star).
E_sun and the sun's color come from a 5772 K black body scaled to 1361 W m^-2
(IAU 2015 Resolution B3 nominal values, Prsa et al. 2016, AJ 152, 41).
Monochromatic lines lie outside the sRGB gamut: white is added until no
channel is negative, keeping the luminance.

Standard library only.
"""

import math
import sys

H = 6.62607015e-34   # J s (CODATA 2018, exact)
C = 299792458.0      # m/s
K = 1.380649e-23     # J/K
KM = 683.0           # lm/W, luminous efficacy at 540 THz (SI)
T_SUN = 5772.0       # K
S_SUN = 1361.0       # W m^-2
SIGMA = 5.670374419e-8

# Lines: wavelength (nm). 391.4 and 427.8 nm: N2+ first negative (0-0), (0-1);
# 557.7 nm: O(1S); 630.0 nm: O(1D) (with its weaker 636.4 nm companion).
LINES = {'violet': 391.4, 'blue': 427.8, 'green': 557.7, 'red': 630.0}

XYZ_TO_SRGB = [  # IEC 61966-2-1
    (3.2406, -1.5372, -0.4986),
    (-0.9689, 1.8758, 0.0415),
    (0.0557, -0.2040, 1.0570),
]
LUMA = (0.2126, 0.7152, 0.0722)


def load_cmf(path):
    cmf = {}
    with open(path, encoding='ascii') as f:
        for line in f:
            p = line.strip().split(',')
            if len(p) == 4:
                cmf[int(p[0])] = tuple(float(x) for x in p[1:])
    return cmf


def cmf_at(cmf, nm):
    lo = math.floor(nm)
    f = nm - lo
    return tuple((1 - f) * a + f * b for a, b in zip(cmf[lo], cmf[lo + 1]))


def planck(nm, t):
    lam = nm * 1e-9
    return 2 * H * C * C / lam ** 5 / (math.exp(H * C / (lam * K * t)) - 1)  # W m^-2 sr^-1 m^-1


def to_srgb(xyz):
    return [sum(m * v for m, v in zip(row, xyz)) for row in XYZ_TO_SRGB]


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    cmf = load_cmf(sys.argv[1])
    # Solar spectral irradiance: pi B_lambda scaled to the total S_SUN.
    scale = S_SUN / (SIGMA * T_SUN ** 4)
    sun = [0.0, 0.0, 0.0]
    for nm, bar in cmf.items():
        e = math.pi * planck(nm, T_SUN) * scale * 1e-9  # W m^-2 per nm
        for k in range(3):
            sun[k] += KM * bar[k] * e
    e_sun = sun[1]  # lux
    white = e_sun / math.pi  # cd m^-2 of the white surface
    sun_rgb = to_srgb([v / e_sun for v in sun])
    sun_y = sum(l * v for l, v in zip(LUMA, sun_rgb))
    sun_rgb = [v / sun_y for v in sun_rgb]
    print(f'// E_sun = {e_sun:.0f} lx, white surface {white:.0f} cd/m2, sunlight linear sRGB {sun_rgb}')

    for name, nm in LINES.items():
        radiance = 1e13 / (4 * math.pi) * H * C / (nm * 1e-9)  # W m^-2 sr^-1 per kR
        xyz = [KM * b * radiance / white for b in cmf_at(cmf, nm)]
        rgb = [v / s for v, s in zip(to_srgb(xyz), sun_rgb)]
        y = sum(l * v for l, v in zip(LUMA, rgb))
        lift = max(0.0, -min(rgb))
        rgb = [v + lift for v in rgb]
        rgb = [v * y / sum(l * u for l, u in zip(LUMA, rgb)) for v in rgb]
        print(f'{name:6} {nm} nm: {{{rgb[0]:.4e}, {rgb[1]:.4e}, {rgb[2]:.4e}}} per kR (luminance {y:.4e})')


if __name__ == '__main__':
    main()
