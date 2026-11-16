"""The dipoles of YY Gem's Zeeman-Doppler maps, for castor.toml (PyMuPDF, NumPy).

Usage:
    python -I tools/castor/yy_gem_zdi_dipoles.py <arXiv 1902.04157 PDF>

Kochukhov & Shulyak 2019 (ApJ 873, 69; arXiv:1902.04157) map both stars from 36
Stokes V profiles of January 2012 (ESPaDOnS), up to degree 10, but publish the
maps only as Fig. 5 (page 8): Hammer-Aitoff projections, the central meridian at
longitude 180, longitude increasing to the right in the sense of the rotation,
filled contours in steps of about 0.08 kG between -0.48 and +0.48.

The figure is vector art. This renders the radial-field panels at 300 dpi, takes
each pixel's value from the nearest colour of the colour bar (pixels far from any,
the grid and labels, are dropped), inverts the projection (equal-area, so pixels
weigh alike) and fits Br = b . r_hat by least squares: on the sphere that is the
degree-1 part of the map, and |b| is Br at the positive magnetic pole.

Check: the mean |Br| comes out 168 and 100 G against the paper's 168 and 98 G
(Table 3). Their longitude faces Earth at phase phi when it is 360 (1 - phi), for
both stars, phase 0 at the primary (Ca's) eclipse; so a pole at longitude L faces
Earth at phase 1 - L / 360. The dipole holds about two thirds of the variance of
Br; degree 4 leaves a quarter of the rms.
"""

import math
import sys

import fitz
import numpy as np

DPI = 300
F = DPI / 150.0                       # the boxes below are in a 150-dpi render
COLOUR_BAR = (957, 387, 538, 384, 540)  # x, sampled y range, y of +0.48 and -0.48 kG
PANELS = {"Ca (YY Gem A)": (256, 569, 384, 542), "Cb (YY Gem B)": (598, 911, 384, 542)}


def render(path):
    page = fitz.open(path)[7]
    pm = page.get_pixmap(dpi=DPI)
    return np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n)[:, :, :3].astype(float)


def swatches(img):
    x, y0, y1, ytop, ybottom = COLOUR_BAR
    ys = np.arange(int(y0 * F), int(y1 * F))
    values = np.interp(ys, [ytop * F, ybottom * F], [0.48, -0.48])
    colours, groups = [], []
    for c, v in zip(img[ys, int(x * F)], values):
        if not colours or np.abs(c - colours[-1]).sum() > 15:
            colours.append(c)
            groups.append([v])
        else:
            groups[-1].append(v)
    return np.array(colours), np.array([np.mean(g) for g in groups])


def panel(img, colours, values, box):
    x0, x1, y0, y1 = (int(v * F) for v in box)
    sub = img[y0:y1, x0:x1]
    h, w = sub.shape[:2]
    yy, xx = np.mgrid[0:h, 0:w]
    X = ((xx + 0.5) / w * 2 - 1) * 2 * math.sqrt(2)
    Y = (1 - (yy + 0.5) / h * 2) * math.sqrt(2)
    inside = (X / (2 * math.sqrt(2))) ** 2 + (Y / math.sqrt(2)) ** 2 < 0.96  # clear of the outline
    z = np.sqrt(np.clip(1 - (X / 4) ** 2 - (Y / 2) ** 2, 0, 1))
    lon = np.radians(180.0) + 2 * np.arctan2(z * X, 2 * (2 * z * z - 1))
    lat = np.arcsin(np.clip(z * Y, -1, 1))
    dist = np.abs(sub[:, :, None, :] - colours[None, None, :, :]).sum(-1)
    ok = inside & (dist.min(-1) < 25)
    return lon[ok], lat[ok], values[dist.argmin(-1)[ok]]


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    img = render(sys.argv[1])
    colours, values = swatches(img)
    for name, box in PANELS.items():
        lon, lat, br = panel(img, colours, values, box)
        r_hat = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], 1)
        b, *_ = np.linalg.lstsq(r_hat, br, rcond=None)
        polar = np.linalg.norm(b)
        colat = math.degrees(math.acos(b[2] / polar))
        pole_lon = math.degrees(math.atan2(b[1], b[0])) % 360.0
        # The magnetic pole on the side of the visible (north, i < 90) rotation pole.
        positive = colat <= 90.0
        tilt = colat if positive else 180.0 - colat
        visible_lon = pole_lon if positive else (pole_lon + 180.0) % 360.0
        residual = br - r_hat @ b
        print(f"{name}: mean |Br| {1000 * np.mean(np.abs(br)):.0f} G over {len(br)} pixels")
        print(f"  dipole: polar_g = {1000 * polar:.0f}, tilt_deg = {tilt:.0f}, "
              f"visible_pole = \"{'positive' if positive else 'negative'}\" at longitude {visible_lon:.0f}, "
              f"phase = {(1.0 - visible_lon / 360.0) % 1.0:.3f}")
        print(f"  Br rms {1000 * np.sqrt(np.mean(br ** 2)):.0f} G, {1000 * np.sqrt(np.mean(residual ** 2)):.0f} G "
              f"left after the dipole")


if __name__ == "__main__":
    main()
