#!/usr/bin/env python3
"""Fits the radial profiles of beta Pictoris' scattered-light disks.

Golimowski et al. 2006 (AJ 131, 3109, arXiv:astro-ph/0602292) give the surface
brightness along the disk's spine as power laws of the projected distance. This
finds the column of dust of an optically thin, axisymmetric disk that shows
those laws when seen from Earth: the line-of-sight integral below, repeated with
a multiplicative correction at each radius. It prints the `radial_au` and
`column` lists of assets/scenes/beta_pictoris.toml and the widths of the result.

tests/stellar_system_tests.cpp measures the same laws with the program's own
integrator (scene/disk.cpp), so the two check each other.

Needs numpy: python -E -P tools/betapic/fit_disk.py
"""

import numpy as np

# Geometry, in the scene's convention (see its [disk] section).
MAIN_NODE, MAIN_I = 29.3, 88.6       # [LA12] via [MA19]; [MA19] Table 1, model 2
# The secondary disk's spine: [GO06] Table 2 and [RE24] Table 2 put it ~5.3 deg
# from the main disk's; 6 deg from edge-on, near side swapped [RE24] Sect. 5.2.
SECOND_NODE, SECOND_I = 34.6, 96.0

# Hedman & Stark 2015 Table 7, the G ring's three Henyey-Greenstein functions.
PHASE_G = (0.995, 0.665, 0.035)
PHASE_W = (0.643, 0.176, 0.181)

# Vertical shape: two Gaussians (sigmas in au, share of the second), constant
# height inside the flare radius and growing with r beyond it [GO06] Sect. 3.1.2.
MAIN_VERTICAL = (3.2, 16.0, 0.42)
SECOND_VERTICAL = (3.2, 3.2, 0.0)
FLARE_AU = 120.0

# [GO06] Table 3: deconvolved F606W indices, the mean of both extensions, in
# regions 1-3; region 4 from the images before deconvolution (the note says the
# deconvolved values there are suspect). Boundaries halfway across the gaps
# between the fitted regions. Inside region 1 [AP15] finds the slope continuous
# down to 0.5 arcsec; outside region 4 the last index simply continues.
COMPOSITE_INDEX = [(69.0, 1.265), (117.5, 1.885), (193.0, 4.47), (1e9, 3.815)]
# [GO06] Table 4: the secondary disk at 80-150 au, the mean of both extensions.
# It cannot be told from the main disk inside 80 au; there it is taken to rise
# like the composite. Sect. 3.2.2: 1.5 mag arcsec^-2 fainter at 100 au.
SECOND_INDEX = [(80.0, 1.265), (1e9, 4.175)]
SECOND_AT_100 = 10.0 ** (-0.4 * 1.5)

INNER_AU, OUTER_AU = 5.0, 2000.0
SECOND_OUTER_AU = 500.0


def plane(node_deg, i_deg):
    """The normal, in (north, east, away) coordinates (VisualOrbit's W)."""
    node, i = np.radians(node_deg), np.radians(i_deg)
    return np.array([-np.sin(node) * np.sin(i), np.cos(node) * np.sin(i), -np.cos(i)])


def phase(cos_theta):
    total = 0.0
    for g, w in zip(PHASE_G, PHASE_W):
        d = 1.0 + g * g - 2.0 * g * cos_theta
        total = total + w * (1.0 - g * g) / (d * np.sqrt(d))
    return total


def profile(index, radius):
    """A continuous broken power law, 1 at the radius of the first break."""
    radius = np.asarray(radius, dtype=float)
    out = np.zeros_like(radius)
    start, value = 0.0, 1.0
    first = True
    for edge, slope in index:
        inside = (radius >= start) & (radius < edge) if not first else radius < edge
        ref = start if not first else edge
        out[inside] = value * (radius[inside] / ref) ** (-slope)
        if first:
            first = False
        else:
            value = value * (edge / start) ** (-slope) if edge < 1e8 else value
        start = edge
    return out


def edge_taper(r, inner, outer):
    low = np.clip((r - inner) / inner, 0.0, 1.0)
    high = np.clip((outer - r) / (0.25 * outer), 0.0, 1.0)
    return low * low * (3.0 - 2.0 * low) * high * high * (3.0 - 2.0 * high)


class Disk:
    def __init__(self, node, i, vertical, outer, index):
        self.normal = plane(node, i)
        self.vertical = vertical
        self.outer = outer
        self.ln_r = np.linspace(np.log(INNER_AU), np.log(outer), 96)
        r = np.exp(self.ln_r)
        self.column = profile(index, r) * r  # a first guess: I(p) ~ column(p) / p
        self.column *= edge_taper(r, INNER_AU, outer)

    def radiance(self, points):
        """Emissivity at (north, east, away) points, per unit path."""
        distance = np.sqrt(np.sum(points * points, axis=-1))
        z = points @ self.normal
        r = np.sqrt(np.maximum(distance * distance - z * z, 1e-12))
        sigma1, sigma2, share = self.vertical
        stretch = (1.0 + (r / FLARE_AU) ** 4) ** 0.25
        s1, s2 = sigma1 * stretch, sigma2 * stretch
        vertical = (1.0 - share) * np.exp(-0.5 * (z / s1) ** 2) / s1 + share * np.exp(-0.5 * (z / s2) ** 2) / s2
        column = np.interp(np.log(np.maximum(r, 1e-6)), self.ln_r, self.column, left=0.0, right=0.0)
        cos_theta = -points[..., 2] / distance
        return column * vertical * phase(cos_theta) / (distance * distance)


def axes(pa_deg):
    pa = np.radians(pa_deg)
    major = np.array([np.cos(pa), np.sin(pa), 0.0])         # toward the NE extension
    across = np.array([np.sin(pa), -np.cos(pa), 0.0])       # toward the NW side
    return major, across


def line_integral(disks, major_offset, across_offset):
    """Surface brightness at sky offsets (arrays of the same shape)."""
    major, across = axes(MAIN_NODE)
    scale = np.maximum(0.3 * np.abs(major_offset), 3.0)
    u = np.linspace(-1.0, 1.0, 801)
    stretch = np.arcsinh(2500.0 / scale)[..., None]
    length = scale[..., None] * np.sinh(u * stretch)
    weight = scale[..., None] * np.cosh(u * stretch) * stretch * (u[1] - u[0])
    p = (major_offset[..., None, None] * major + across_offset[..., None, None] * across +
         length[..., None] * np.array([0.0, 0.0, 1.0]))
    total = 0.0
    for disk in disks:
        total = total + np.sum(disk.radiance(p) * weight, axis=-1)
    return total


def spine(disks, projected, side=1.0):
    """Brightest value across the disk at each projected distance."""
    offsets = np.linspace(-30.0, 30.0, 121)
    major = np.repeat(side * projected[:, None], offsets.size, axis=1)
    across = np.repeat(offsets[None, :], projected.size, axis=0)
    image = line_integral(disks, major, across)
    best = np.argmax(image, axis=1)
    return image[np.arange(projected.size), best], offsets[best]


def both_sides(disks, projected):
    ne, _ = spine(disks, projected, 1.0)
    sw, _ = spine(disks, projected, -1.0)
    return 0.5 * (ne + sw)


def fit(disk, others, target, projected, rounds):
    grid = np.exp(disk.ln_r)
    for k in range(rounds):
        model = both_sides([disk] + others, projected)
        correction = np.log(target / model)
        disk.column *= np.exp(0.8 * np.interp(np.log(grid), np.log(projected), correction))
        worst = np.max(np.abs(correction))
    return worst


def slope(disks, lo, hi):
    p = np.geomspace(lo, hi, 12)
    values = both_sides(disks, p)
    return -np.polyfit(np.log(p), np.log(values), 1)[0]


def widths(disks, projected):
    offsets = np.linspace(-80.0, 80.0, 641)
    major = np.full(offsets.shape, projected)
    image = line_integral(disks, major, offsets)
    peak = image.max()
    return (np.ptp(offsets[image >= 0.5 * peak]), np.ptp(offsets[image >= 0.1 * peak]))


def print_table(name, disk, count):
    r = np.geomspace(INNER_AU, disk.outer, count)
    values = np.exp(np.interp(np.log(r), disk.ln_r, np.log(np.maximum(disk.column, 1e-300))))
    values /= np.interp(np.log(115.0), disk.ln_r, disk.column)
    print(f"# {name}")
    print("radial_au = [" + ", ".join(f"{x:.4g}" for x in r) + "]")
    print("column = [" + ", ".join(f"{x:.4g}" for x in values) + "]")


def main():
    projected = np.geomspace(20.0, 1500.0, 40)
    main_disk = Disk(MAIN_NODE, MAIN_I, MAIN_VERTICAL, OUTER_AU, COMPOSITE_INDEX)
    second = Disk(SECOND_NODE, SECOND_I, SECOND_VERTICAL, SECOND_OUTER_AU, SECOND_INDEX)
    composite_target = profile(COMPOSITE_INDEX, projected)
    second_points = projected[projected < 0.7 * SECOND_OUTER_AU]
    second_target = profile(SECOND_INDEX, second_points)

    for _ in range(4):
        fit(second, [], second_target, second_points, 6)
        # The level relative to the main disk at 100 au [GO06].
        fit(main_disk, [second], composite_target, projected, 6)
        ratio = both_sides([second], np.array([100.0]))[0] / both_sides([main_disk], np.array([100.0]))[0]
        second.column *= SECOND_AT_100 / ratio
        second_target = second_target * (SECOND_AT_100 / ratio)
    worst = fit(main_disk, [second], composite_target, projected, 12)
    print(f"# worst misfit of the composite: {100.0 * worst:.2f}%")

    disks = [main_disk, second]
    for lo, hi, index in ((39, 67, 1.265), (71, 108, 1.885), (127, 193, 4.47), (193, 258, 3.815)):
        print(f"# slope {lo}-{hi} au: {slope(disks, lo, hi):.3f} (target {index})")
    for p in (60.0, 100.0, 200.0):
        fwhm, fw01 = widths(disks, p)
        print(f"# at {p:.0f} au: FWHM {fwhm:.1f} au, full width at 0.1 max {fw01:.1f} au")
    _, offset = spine([main_disk], np.array([40.0, 100.0, 250.0]))
    print("# main disk's spine, toward NW (au) at 40, 100, 250 au:", np.round(offset, 2))
    _, offset = spine([second], np.array([100.0, 130.0]))
    print("# secondary's spine, toward NW (au) at 100, 130 au, NE:", np.round(offset, 2))
    _, offset = spine([second], np.array([100.0, 130.0]), -1.0)
    print("# secondary's spine, toward NW (au) at 100, 130 au, SW:", np.round(offset, 2))

    print_table("Main disk", main_disk, 40)
    second.column /= np.interp(np.log(115.0), main_disk.ln_r, main_disk.column)
    # print_table normalizes each to its own value at 115 au; keep the ratio.
    ratio = np.interp(np.log(115.0), second.ln_r, second.column)
    print(f"# Secondary disk (column at 115 au relative to the main disk's: {ratio:.4f})")
    print_table("Secondary disk", second, 24)


if __name__ == "__main__":
    main()
