"""The height profile and the September 2017 timeline of Mars's diffuse aurora
(prints the constants of src/scene/mars_aurora.hpp).

    python -I mars_aurora_profile.py

Schneider et al. 2018 (GRL 45, 7391, doi:10.1029/2018GL077772) give both only as
plots, read here from 260-400 dpi renderings of the PDF by locating the curve and
the markers.

LIMB is the CO2+ UVD of their Fig. 2, orbit 5730: intensity seen edge-on at each
tangent height, so a line of sight at 25 km still crosses the layer at 60 km on the
way. It is not the emission at that height: comparing it directly with an emission
profile makes the layer look wider below than it is. The emission is fitted through
the line-of-sight integral: a Gaussian with its own sigma on each side.

MARKERS is their Fig. 3, CO Cameron bands: limb scans and the limb in apoapse images,
(day of September 2017, kR). The sites differ, and a few hours apart they differ by a
factor of 5, so the timeline takes the brightest of each 6 hours, keeps the first
peak (the limb scan of orbit 5730, the same orbit as Fig. 2), and raises anything
under the noise floor (about 50 R in the caption, 0.04 kR where the line is drawn) to it.

Standard library only.
"""

import math

RADIUS_KM = 3396.19

LIMB = ((23, 11.2), (30, 13.4), (35, 14.6), (40, 16.0), (45, 17.3), (50, 18.6), (55, 19.7),
        (57.5, 20.0), (60, 18.9), (65, 16.2), (70, 13.1), (75, 10.3), (80, 7.1), (85, 4.7),
        (90, 3.2), (95, 1.8), (100, 1.0), (110, 0.4), (120, 0.2))

MARKERS = ((7.01, 0.021), (8.43, 0.013), (10.40, 0.100), (11.15, 0.038), (11.21, 0.099), (11.27, 0.047),
           (11.58, 0.083), (11.76, 0.119), (11.88, 0.287), (11.95, 0.163), (12.01, 0.069), (12.13, 0.058),
           (12.50, 0.103), (12.87, 0.396), (12.99, 0.786), (13.06, 0.697), (13.12, 3.906), (13.24, 2.834),
           (13.43, 0.885), (13.61, 1.503), (13.73, 0.589), (13.86, 0.868), (13.98, 0.323), (14.17, 0.220),
           (14.35, 0.052), (14.47, 0.283), (14.60, 0.139), (14.72, 0.056), (14.91, 0.049), (15.21, 0.083),
           (15.34, 0.142), (16.69, 0.011), (17.06, 0.013), (17.18, 0.051), (17.43, 0.070), (17.80, 0.014),
           (17.87, 0.024), (18.05, 0.051), (18.17, 0.091), (18.54, 0.064), (18.67, 0.049), (18.92, 0.073),
           (19.28, 0.053), (20.39, 0.050), (20.89, 0.099), (21.14, 0.145), (21.20, 0.087), (21.26, 0.299),
           (21.39, 0.140), (21.50, 0.312), (21.57, 0.023), (21.63, 0.153), (22.00, 0.172), (22.24, 0.043),
           (22.38, 0.043), (22.62, 0.060), (22.74, 0.078), (23.48, 0.027), (23.73, 0.030), (24.22, 0.071),
           (24.47, 0.019), (25.21, 0.057), (25.58, 0.016))

NOISE_FLOOR_KR = 0.04


def layer(h, peak, low, high):
    return math.exp(-0.5 * ((h - peak) / (low if h < peak else high)) ** 2)


def edge_on(b, peak, low, high, top=160.0, n=400):
    end = math.sqrt((RADIUS_KM + top) ** 2 - (RADIUS_KM + b) ** 2)
    total = 0.0
    for j in range(n):
        u = (j + 0.5) / n
        total += layer(math.hypot(RADIUS_KM + b, end * u * u) - RADIUS_KM, peak, low, high) * 2.0 * end * u / n
    return 2.0 * total


def fit():
    best = None
    for peak in range(50, 76):
        for low in range(4, 27):
            for high in range(8, 35):
                model = [edge_on(h, peak, low, high) for h, _ in LIMB]
                scale = sum(m * o for m, (_, o) in zip(model, LIMB)) / sum(m * m for m in model)
                err = sum((scale * m - o) ** 2 for m, (_, o) in zip(model, LIMB))
                if best is None or err < best[0]:
                    best = (err, peak, low, high, scale)
    return best


def main():
    err, peak, low, high, scale = fit()
    print(f'peak {peak} km, sigma {low} km below, {high} km above; rms {math.sqrt(err / len(LIMB)):.2f} kR')
    for h, observed in LIMB:
        print(f'  {h:5.1f} km: observed {observed:5.1f} kR, model {scale * edge_on(h, peak, low, high):5.1f} kR,'
              f' old Gaussian {edge_on(h, 60, 15, 15) * scale:5.1f}')
    print(f'edge-on peak at {max((edge_on(h, peak, low, high), h) for h in range(20, 100))[1]} km (observed 57.5);'
          f' old Gaussian at {max((edge_on(h, 60, 15, 15), h) for h in range(20, 100))[1]} km')

    bins = {}
    for day, kr in MARKERS:
        key = math.floor(day * 4)
        if key not in bins or kr > bins[key][1]:
            bins[key] = (day, kr)
    print('\n// mars_aurora.hpp: day of September 2017, Cameron kR')
    knots = [bins[k] for k in sorted(bins)]
    print(', '.join(f'{{{day:.2f}, {max(kr, NOISE_FLOOR_KR):.3f}}}' for day, kr in knots))
    print(f'{len(knots)} knots')


if __name__ == '__main__':
    main()
