"""Reads the curves of Mejnertsen et al. 2016 (JGR Space Physics 121, 7497), Figure 4, from
the article's PDF: Neptune's simulated bow shock and magnetopause in the y = 0 and z = 0
planes, with the dipole "Earth-like" and "pole-on". See assets/geomag/SOURCES.md.

    python -I tools/geomag/neptune_boundaries_fig4.py <the article's PDF>

The figure is vector art: each curve is a filled outline of its stroke (0.6 pt wide, about
0.4 radii), so its points are within half that of the line. Curves are told apart by their
fill color and, of the two of a color in a panel, by which reaches farther sunward (the
shock). Needs PyMuPDF.
"""

import sys

import fitz

PAGE = 7  # journal page 7504
COLORS = {"Earth-like": (1.0, 0.0, 0.0), "pole-on": (0.0, 0.58, 0.78)}


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    page = fitz.open(sys.argv[1])[PAGE]
    frames = set()
    outlines = []
    for drawing in page.get_drawings():
        points = []
        for item in drawing["items"]:
            if item[0] == "l":
                points += [item[1], item[2]]
            elif item[0] == "re" and item[1].width > 100 and item[1].height > 150:
                frames.add((round(item[1].x0, 1), round(item[1].y0, 1), round(item[1].x1, 1), round(item[1].y1, 1)))
        fill = drawing.get("fill")
        for name, color in COLORS.items():
            if fill and len(points) > 100 and all(abs(a - b) < 0.02 for a, b in zip(fill, color)):
                outlines.append((name, points))
    frames = sorted(frames)
    if len(frames) != 2 or len(outlines) != 8:
        sys.exit(f"expected 2 panels and 8 curves, found {len(frames)} and {len(outlines)}")

    for panel, (x0, y0, x1, y1) in zip(("y = 0 (across: z)", "z = 0 (across: y)"), frames):
        print(panel)
        curves = {}
        for name, points in outlines:
            if not x0 - 1 <= points[0].x <= x1 + 1:
                continue
            # Axes: x from 40 at the left to -40, the other from 60 at the top to -60.
            xy = [(40.0 - 80.0 * (p.x - x0) / (x1 - x0), 60.0 - 120.0 * (p.y - y0) / (y1 - y0)) for p in points]
            curves.setdefault(name, []).append(xy)
        for name, pair in curves.items():
            pair.sort(key=max, reverse=True)
            for what, xy in zip(("bow shock", "magnetopause"), pair):
                def x_at(v):
                    # Both edges of the outline cross the level: their middle. The magnetopause
                    # is dash-dotted, so the nearest piece may be a little off the level.
                    off = min(abs(q[1] - v) for q in xy)
                    near = [q[0] for q in xy if abs(q[1] - v) < off + 0.6]
                    return 0.5 * (min(near) + max(near))
                sunward = max(xy)
                print(f"  {name:10s} {what:12s} on the axis {x_at(0.0):5.2f}; most sunward {sunward[0]:5.2f} at {sunward[1]:+5.1f}; "
                      f"x at -20, -10, +10, +20: {x_at(-20.0):5.1f} {x_at(-10.0):5.1f} {x_at(10.0):5.1f} {x_at(20.0):5.1f}")


if __name__ == "__main__":
    main()
