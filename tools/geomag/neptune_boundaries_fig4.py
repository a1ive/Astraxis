"""Reads the curves of Mejnertsen et al. 2016 (JGR Space Physics 121, 7497, CC BY 4.0), Figure 4,
from the article's PDF, and bakes them into assets/geomag/neptune_boundaries.txt: Neptune's
simulated bow shock and magnetopause in the y = 0 and z = 0 planes, with the dipole
"Earth-like" and "pole-on". See assets/geomag/SOURCES.md.

    python -I tools/geomag/neptune_boundaries_fig4.py <the article's PDF> [assets/geomag/neptune_boundaries.txt]

The figure is vector art: each curve is a filled outline of its stroke (0.6 pt wide, about
0.4 radii), so its points are within half that of the line. Curves are told apart by their
fill color and, of the two of a color in a panel, by which reaches farther sunward (the
shock). Needs PyMuPDF.

Each curve is cut at the Sun's line into two half profiles, the distance from the planet
against the angle from that line, so that a surface has four per state: toward +y (dusk), +z
(north), -y and -z. The output has a line "AXBOUNDARY1 <angles> <step in degrees>", then for
each surface (shock, magnetopause), state (earth_like, pole_on) and azimuth (0, 90, 180, 270
degrees from +y toward +z) a line "<surface> <state> <azimuth> <flaring>" followed by the
distances in planetary radii at 0, step, 2 step .. degrees. Beyond the last angle a profile
goes on as r (1 + cos theta)^flaring = constant, with the slope it ends with. The four
profiles of a surface and state start from the same distance on the Sun's line (the two
planes' mean).
"""

import math
import sys

import fitz

PAGE = 7  # journal page 7504
COLORS = {"earth_like": (1.0, 0.0, 0.0), "pole_on": (0.0, 0.58, 0.78)}
STEP_DEG = 2.5
# The frame ends at 60 radii across: how far round each surface has all its profiles.
LAST_DEG = {"shock": 80.0, "magnetopause": 115.0}
SMOOTH_DEG = 2.0  # the outlines are of dash-dotted and hand-smoothed curves


def half_profile(xy, sign):
    """Distance against angle from +x on one side (the sign of the other coordinate), in whole
    degrees from 0, from the outline's points: the mean of the two edges of the stroke, the
    gaps of a dash-dotted line filled along the curve."""
    bins = {}
    for x, v in xy:
        if v * sign < -0.3:
            continue
        angle = math.degrees(math.atan2(abs(v), x))
        # A degree and a half either way, so that both edges of the stroke are in each bin.
        for degree in range(max(round(angle) - 1, 0), round(angle) + 2):
            if abs(angle - degree) <= 1.5:
                bins.setdefault(degree, []).append(math.hypot(x, v))
    known = sorted(bins)
    last = known[-1]
    centers = {k: 0.5 * (min(bins[k]) + max(bins[k])) for k in known}
    out = []
    for degree in range(0, last + 1):
        if degree in centers:
            out.append(centers[degree])
            continue
        below = max((k for k in known if k < degree), default=None)
        above = min(k for k in known if k > degree)
        if below is None:
            out.append(centers[above])
        else:
            f = (degree - below) / (above - below)
            out.append(centers[below] + f * (centers[above] - centers[below]))
    return out


def smooth(values, width):
    out = []
    for i in range(len(values)):
        total = weight = 0.0
        for j in range(max(0, i - 3 * width), min(len(values), i + 3 * width + 1)):
            w = math.exp(-0.5 * ((j - i) / width) ** 2)
            total += w * values[j]
            weight += w
        out.append(total / weight)
    return out


def read_curves(pdf):
    """{surface: {state: {azimuth: [distance at each whole degree]}}}, and the outlines' points
    {(surface, state, plane): [(x, across)]}."""
    page = fitz.open(pdf)[PAGE]
    frames = set()
    outlines = []
    for drawing in page.get_drawings():
        points = []
        for item in drawing["items"]:
            if item[0] == "l":
                # The outline has few corners where the curve is smooth: points all along each side.
                steps = max(1, int(abs(item[2] - item[1]) / 0.25))
                points += [item[1] + (item[2] - item[1]) * (k / steps) for k in range(steps + 1)]
            elif item[0] == "c":
                # A cubic Bezier: the smooth stretches of the curves are these.
                a, b, c, d = item[1:5]
                steps = max(2, int((abs(b - a) + abs(c - b) + abs(d - c)) / 0.25))
                for k in range(steps + 1):
                    t = k / steps
                    points.append(a * (1 - t) ** 3 + b * (3 * t * (1 - t) ** 2) + c * (3 * t * t * (1 - t)) + d * t ** 3)
            elif item[0] == "re" and item[1].width > 100 and item[1].height > 150:
                frames.add((round(item[1].x0, 1), round(item[1].y0, 1), round(item[1].x1, 1), round(item[1].y1, 1)))
        fill = drawing.get("fill")
        for name, color in COLORS.items():
            if fill and len(points) > 100 and all(abs(a - b) < 0.02 for a, b in zip(fill, color)):
                outlines.append((name, points))
    frames = sorted(frames)
    if len(frames) != 2 or len(outlines) != 8:
        sys.exit(f"expected 2 panels and 8 curves, found {len(frames)} and {len(outlines)}")

    curves = {"shock": {}, "magnetopause": {}}
    outline_points = {}
    # The left panel is the y = 0 plane (across: z), the right one z = 0 (across: y).
    for (x0, y0, x1, y1), azimuths in zip(frames, ((90, 270), (0, 180))):
        pairs = {}
        for name, points in outlines:
            if not x0 - 1 <= points[0].x <= x1 + 1:
                continue
            # Axes: x from 40 at the left to -40, the other from 60 at the top to -60.
            xy = [(40.0 - 80.0 * (p.x - x0) / (x1 - x0), 60.0 - 120.0 * (p.y - y0) / (y1 - y0)) for p in points]
            pairs.setdefault(name, []).append(xy)
        for name, pair in pairs.items():
            pair.sort(key=max, reverse=True)
            for surface, xy in zip(("shock", "magnetopause"), pair):
                outline_points[(surface, name, "z" if azimuths[0] == 90 else "y")] = xy
                for azimuth, sign in zip(azimuths, (1.0, -1.0)):
                    curves[surface].setdefault(name, {})[azimuth] = half_profile(xy, sign)
    return curves, outline_points


def crossing(xy, across):
    """x where a curve crosses a level of the other coordinate: the middle of the stroke's two
    edges, straight from the outline (a check on the profiles, which are binned by angle)."""
    near = [x for x, v in xy if abs(v - across) < 0.15 and x > 0.0]
    return 0.5 * (min(near) + max(near)) if near else float("nan")


def main():
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    curves, outline_points = read_curves(sys.argv[1])
    for (surface, state, plane), xy in sorted(outline_points.items()):
        print(f"{surface:12s} {state:10s} x where {plane} = -20, -10, 10, 20: " +
              " ".join(f"{crossing(xy, v):5.2f}" for v in (-20.0, -10.0, 10.0, 20.0)))
    count = {surface: int(round(last / STEP_DEG)) + 1 for surface, last in LAST_DEG.items()}
    lines = [f"AXBOUNDARY1 {STEP_DEG}",
             "# Neptune's bow shock and magnetopause as simulated for Voyager 2's flyby: Mejnertsen et al. 2016,",
             "# JGR Space Physics 121, 7497, doi:10.1002/2015JA022272 (CC BY 4.0), Figure 4, read from the PDF's",
             "# vector curves by tools/geomag/neptune_boundaries_fig4.py; see SOURCES.md. Distances in Neptune",
             "# radii (24,765 km) against the angle from the Sun's line, every 2.5 deg from 0; azimuth from +y",
             "# (against the orbital motion) toward +z (north of the orbit)."]
    for surface in ("shock", "magnetopause"):
        for state in ("earth_like", "pole_on"):
            profiles = curves[surface][state]
            for azimuth, profile in profiles.items():
                if len(profile) <= LAST_DEG[surface] + SMOOTH_DEG:
                    sys.exit(f"{surface} {state} {azimuth}: the curve ends at {len(profile) - 1} deg")
            nose = sum(profile[0] for profile in profiles.values()) / 4.0
            summary = []
            for azimuth in (0, 90, 180, 270):
                smoothed = smooth(profiles[azimuth], int(SMOOTH_DEG))
                # Pin the Sun's line to the common nose, the correction fading over 20 deg.
                shift = nose - smoothed[0]
                smoothed = [r + shift * max(0.0, 1.0 - k / 20.0) for k, r in enumerate(smoothed)]
                at = lambda degrees: smoothed[int(degrees)] + (degrees - int(degrees)) * (smoothed[int(degrees) + 1] - smoothed[int(degrees)])
                values = [at(k * STEP_DEG) for k in range(count[surface])]
                # d ln r / d theta = flaring tan(theta / 2), over the last 10 degrees.
                end = LAST_DEG[surface]
                slope = (math.log(at(end)) - math.log(at(end - 10.0))) / math.radians(10.0)
                flaring = min(max(slope / math.tan(math.radians(end - 5.0) / 2.0), 0.3), 1.2)
                lines.append(f"{surface} {state} {azimuth} {flaring:.3f} " + " ".join(f"{r:.2f}" for r in values))
                summary.append(f"{azimuth:3d}: {at(30.0):5.1f} {at(60.0):5.1f} {values[-1]:5.1f} ({flaring:.2f})")
            print(f"{surface:12s} {state:10s} on the Sun's line {nose:5.2f}; at 30, 60 and {LAST_DEG[surface]:.0f} deg (flaring) by "
                  "azimuth " + "; ".join(summary))
    if len(sys.argv) == 3:
        with open(sys.argv[2], "w", encoding="ascii", newline="\n") as out:
            out.write("\n".join(lines) + "\n")
        print(f"{sys.argv[2]}: {len(lines) - 6} profiles")


if __name__ == "__main__":
    main()
