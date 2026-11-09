"""The light curve of Proxima's superflare of 2016 March 18, for alpha_centauri.toml.

Usage:
    python -E -P tools/proxima/superflare_curve.py <arXiv 1804.02001 PDF> <out.csv>

Needs PyMuPDF. Reads the Evryscope points from the vector Figure 1 of Howard et
al. 2018 (ApJL 860, L30): flux relative to quiescence, two-minute integrations.
Around the peak they are replaced by the flare template of Davenport et al.
2014 (ApJ 797, 122, Eqs. 1 and 4) with the paper's instantaneous peak of 68
times the quiescent flux, its scale time and peak time fitted so that the
template averaged over the two-minute integrations gives the measured points.
Time is from the peak time of the paper, 08:32:10 UT.

Also prints the area of the flare on the disk at its peak and the energy in
the g band, which the paper gives independently (10^32.8 erg).
"""

import math
import sys

import fitz

# Figure 1 axes (page 2, points): ticks at -20 and 80 min, and at 0 and 5 in flux.
X_0, X_PER_MIN = 139.12, (282.65 - 103.24) / 100.0
Y_0, Y_PER_UNIT = 238.20, (238.20 - 222.05) / 5.0
X_MIN, X_MAX = 71.95, 291.55     # the marks outside the axes are clipped in the figure

PEAK = 68.0                      # Howard et al. 2018, Sect. 3.1: flux at the peak over quiescence
EXPOSURE_MIN = 2.0               # Sect. 2: two-minute cadence
L0_ERG_S = 10.0 ** 28.6          # Sect. 3.2: quiescent luminosity in the g band
FG = 0.19                        # Sect. 3.2: share of a 9000 K flare's light in the g band
SIGMA_CGS = 5.670374419e-5
T_FLARE = 9000.0
R_STAR_CM = 107277.0e5           # scene: 0.1542 R_sun, Kervella et al. 2017
SITE_MU = math.cos(math.radians(30.0))  # the scene's assumed site, 30 deg from the sub-Earth point


def template(x):
    """Davenport et al. 2014, Eqs. 1 and 4; x in units of t_1/2 from the peak."""
    if x < -1.0:
        return 0.0
    if x <= 0.0:
        return max(1.0 + 1.941 * x - 0.175 * x ** 2 - 2.246 * x ** 3 - 1.125 * x ** 4, 0.0)
    return 0.6890 * math.exp(-1.600 * x) + 0.3030 * math.exp(-0.2783 * x)


def read_points(pdf):
    points = []
    for d in fitz.open(pdf)[1].get_drawings():
        if d['type'] == 'fs' and len(d['items']) == 8 and all(i[0] == 'c' for i in d['items']):
            r = d['rect']
            x, y = 0.5 * (r.x0 + r.x1), 0.5 * (r.y0 + r.y1)
            if X_MIN <= x <= X_MAX:
                points.append(((x - X_0) / X_PER_MIN, (Y_0 - y) / Y_PER_UNIT))
    return sorted(points)


def binned(t_center, peak_min, half_min):
    """The template's excess flux averaged over one integration."""
    n = 40
    total = 0.0
    for i in range(n):
        t = t_center + EXPOSURE_MIN * ((i + 0.5) / n - 0.5)
        total += template((t - peak_min) / half_min)
    return (PEAK - 1.0) * total / n


def fit_template(points):
    used = [p for p in points if -3.5 < p[0] < 4.5]
    best = (1e30, 0.0, 0.0)
    for i in range(400):
        peak_min = -2.0 + i * 0.01
        for j in range(1, 300):
            half_min = j * 0.01
            cost = sum((binned(t, peak_min, half_min) - (f - 1.0)) ** 2 for t, f in used)
            if cost < best[0]:
                best = (cost, peak_min, half_min)
    print('template fit: t_peak %.2f min, t_1/2 %.2f min, rms %.2f' % (best[1], best[2], math.sqrt(best[0] / len(used))))
    for t, f in used:
        print('  bin %+.2f min: flux %.2f, template %.2f' % (t, f, 1.0 + binned(t, best[1], best[2])))
    return best[1], best[2]


def main(pdf, out):
    points = read_points(pdf)
    print(len(points), 'points, %.1f .. %.1f min' % (points[0][0], points[-1][0]))
    peak_min, half_min = fit_template(points)
    start = peak_min - half_min
    # Back to the measured points where they pass the template (the third weaker flare is real).
    after = [p for p in points if p[0] > peak_min]
    stop = next(t for t, f in after if f - 1.0 > (PEAK - 1.0) * template((t - peak_min) / half_min))
    rows = [(t, f) for t, f in points if t < start]
    steps = int((stop - start) / (0.04 * half_min)) + 1
    rows += [(start + (stop - start) * k / steps,
              1.0 + (PEAK - 1.0) * template((start + (stop - start) * k / steps - peak_min) / half_min))
             for k in range(steps)]
    rows.append((peak_min, PEAK))
    rows += [(t, f) for t, f in points if t >= stop]
    rows.sort()
    print('template from %.2f to %.2f min, flux there %.2f' % (start, stop, 1.0 + (PEAK - 1.0) * template((stop - peak_min) / half_min)))
    with open(out, 'w', encoding='ascii', newline='\n') as f:
        f.write('# Proxima superflare of 2016 March 18: Howard et al. 2018, Fig. 1, around the peak\n')
        f.write('# Davenport et al. 2014; tools/proxima/superflare_curve.py. Minutes from 08:32:10 UT.\n')
        f.write('minutes,flux\n')
        for t, value in rows:
            f.write('%.3f,%.3f\n' % (t, max(value, 1.0)))

    energy = 0.0
    for (t0, f0), (t1, f1) in zip(rows, rows[1:]):
        energy += 0.5 * (f0 + f1 - 2.0) * (t1 - t0) * 60.0 * L0_ERG_S
    area = (PEAK - 1.0) * L0_ERG_S / (4.0 * FG * SIGMA_CGS * T_FLARE ** 4)
    share = area / (math.pi * R_STAR_CM ** 2)
    print('g band energy 10^%.2f erg (paper 10^32.8)' % math.log10(energy))
    print('peak disk share %.4f (flat 9000 K patch), cap radius at mu %.2f: %.1f deg'
          % (share, SITE_MU, math.degrees(math.asin(math.sqrt(share / SITE_MU)))))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
