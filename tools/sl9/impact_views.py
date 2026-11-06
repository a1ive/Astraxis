"""Print what assets/bodies/jupiter.toml and the impact events of assets/scenes/jupiter.toml need,
from the tables tools/sl9/sl9_bake.py baked: where and when each fragment arrived, and camera offsets
for fragment G as the Earth and Galileo saw it.

Usage (needs the network; run outside the sandbox):
    python -I tools/sl9/impact_views.py

Site classes: [HAM95] Hammel et al. 1995, Science 267, 1288, Table 2. F, P2, T, U and V "were not
detected" (class 4); P1 faded from view before it arrived, so it is class 4 too.

Only the Python standard library is used.
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'tools', 'bake'))
import horizons_bake as hb  # noqa: E402

EPHEM = os.path.join(REPO, 'assets', 'ephem')
CLASSES = [('A', '2a'), ('B', '3'), ('C', '2a'), ('D', '3'), ('E', '2a'), ('F', '4'), ('G', '1'), ('H', '2a'),
           ('K', '1'), ('L', '1'), ('N', '3'), ('P1', '4'), ('P2', '4'), ('Q1', '2b'), ('Q2', '3'), ('R', '2b'),
           ('S', '2c'), ('T', '4'), ('U', '4'), ('V', '4'), ('W', '2c')]
# The peak flux at 945 nm, 1e-15 W/cm2/nm: Martin et al. 1995 Table 1, and Chapman et al. 1995 for K, N, W; L is 6% of Jupiter.
FLASH = {'G': 2.1, 'H': 1.8, 'L': 2.6, 'Q1': 1.1, 'K': 1.4, 'N': 0.57, 'W': 0.66}
RA, DEC = math.radians(268.056595), math.radians(64.495303)  # pck00011, as bodies/jupiter.toml
POLE = (math.cos(DEC) * math.cos(RA), math.cos(DEC) * math.sin(RA), math.sin(DEC))
DISTANCE_KM = 1.0e6       # far, so that the limb hides what it hid from far away, with a narrow view


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def unit(a):
    n = math.sqrt(dot(a, a))
    return tuple(x / n for x in a)


def body_axes(t):
    """The body frame's axes in ICRF: x through the prime meridian, z the pole (pck00011)."""
    w = math.radians(284.95 + 870.536 * t / hb.DAY)
    node = (-math.sin(RA), math.cos(RA), 0.0)
    east = cross(POLE, node)
    x = tuple(math.cos(w) * node[i] + math.sin(w) * east[i] for i in range(3))
    return x, cross(POLE, x), POLE


def vector(command, t):
    jd = t / hb.DAY + hb.J2000_JD
    return hb.horizons_vectors(command, '500@599', jd, jd + 2.0 / 1440.0, 1)[0][1]


def fmt(v):
    return '[' + ', '.join(f'{x:.1f}' for x in v) + ']'


def view(name, site, direction):
    """Offsets for a camera far along `direction`, with the site a little above the middle of the view."""
    d = unit(direction)
    up = unit(tuple(POLE[i] - dot(POLE, d) * d[i] for i in range(3)))
    aim = tuple(site[i] - 5000.0 * up[i] for i in range(3))
    camera = tuple(site[i] + DISTANCE_KM * d[i] for i in range(3))
    print(f'# {name}: angle from the site {math.degrees(math.acos(dot(unit(site), d))):.2f} deg')
    print(f'observer_offset_km = {fmt(camera)}')
    print(f'aim_offset_km = {fmt(aim)}')


def main():
    print('# --- for assets/bodies/jupiter.toml ---')
    ends = {}
    for name, klass in CLASSES:
        end = hb.read_eph(os.path.join(EPHEM, f'sl9_{name.lower()}.eph')).knots[-1]
        ends[name] = end
        x, y, z = body_axes(end[0])
        fixed = (dot(end[1], x), dot(end[1], y), dot(end[1], z))
        lat = math.degrees(math.asin(fixed[2] / math.sqrt(dot(fixed, fixed))))
        west = (-math.degrees(math.atan2(fixed[1], fixed[0]))) % 360.0
        print(f'[[impacts]]\nname = "{name}"\nclass = "{klass}"\ntime_tdb = "{hb.calendar(end[0])[:19]}"\n'
              f'lat_deg = {lat:.2f}\nwest_lon_deg = {west:.2f}' +
              (f'\nflash_share = {0.06 * FLASH[name] / 2.6:.4f}' if name in FLASH else '') + '\n')

    t, site, _ = ends['G']
    print('\n# --- the views of fragment G, ICRF ---')
    view('from the Earth', site, vector('399', t))
    view('from Galileo', site, vector('-77', t))

    # From overhead, turning with Jupiter.
    x, y, z = body_axes(t)
    fixed = (dot(site, x), dot(site, y), dot(site, z))
    normal = unit((fixed[0] / 71492.0 ** 2, fixed[1] / 71492.0 ** 2, fixed[2] / 66854.0 ** 2))
    north = unit(tuple((0.0, 0.0, 1.0)[i] - normal[2] * normal[i] for i in range(3)))
    print('\n# --- the site from overhead, body frame ---')
    print(f'observer_offset_km = {fmt(tuple(fixed[i] + DISTANCE_KM * normal[i] for i in range(3)))}')
    print(f'aim_offset_km = {fmt(tuple(fixed[i] - 8000.0 * north[i] for i in range(3)))}')


if __name__ == '__main__':
    main()
