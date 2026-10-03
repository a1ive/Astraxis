"""Phase corrections for "mean_elements" moons that also have a baked
ephemeris (ephemeris orbit with a mean_elements fallback in a scene file).

Usage:
    python fit_moon_phase.py <scene.toml> [body ...]

For each such body, compares the mean-element model (the same model as
src/ephem/mean_element_orbit.cpp: a precessing ellipse in the Laplace plane,
whose x axis is the plane's ascending node on the ICRF equator) with the baked
table every 2 days over the table's span. It fits three angles at the epoch,
keeping a, e, i and all rates: the mean longitude, the longitude of periapsis
and the node (the latter two only matter for eccentric or inclined orbits),
minimizing the mean squared angle between model and table (coarse grid, then
coordinate descent). It prints the corrected node_deg, arg_peri_deg and
mean_anomaly_deg. Fitting over the whole span rather than at one instant
averages out librations (e.g. Mimas' 70-year libration in its resonance with
Tethys).

Used for Saturn's moons, whose [ELEM] mean anomalies do not match SAT441 at
their epoch (see assets/scenes/saturn.toml). Only the standard library is used.
"""

import math
import os
import sys
import tomllib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import horizons_bake as hb  # noqa: E402

DAY = 86400.0
J2000_JD = 2451545.0


def unit_from_ra_dec(ra_deg, dec_deg):
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    return (math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def normalize(a):
    n = math.sqrt(dot(a, a))
    return tuple(x / n for x in a)


def solve_kepler(m, e):
    m = math.remainder(m, 2 * math.pi)
    big_e = m if e < 0.8 else math.copysign(math.pi, m)
    for _ in range(50):
        d = (big_e - e * math.sin(big_e) - m) / (1 - e * math.cos(big_e))
        big_e -= d
        if abs(d) < 1e-14:
            break
    return big_e


class MeanElements:
    def __init__(self, t):
        self.epoch = (t.get('epoch_jd_tdb', J2000_JD) - J2000_JD) * DAY
        self.a, self.e = t['a_km'], t['e']
        self.w0, self.m0 = math.radians(t['arg_peri_deg']), math.radians(t['mean_anomaly_deg'])
        self.i, self.node0 = math.radians(t['i_deg']), math.radians(t['node_deg'])

        def rate(period_s):
            return 2 * math.pi / period_s if period_s else 0.0
        year = 365.25 * DAY
        self.n = rate(t['period_days'] * DAY)
        self.w_rate = rate(t.get('apsis_period_years', 0.0) * year)
        self.node_rate = rate(t.get('node_period_years', 0.0) * year)
        self.pole = unit_from_ra_dec(*t['laplace_pole_ra_dec_deg'])
        x = cross((0.0, 0.0, 1.0), self.pole)
        self.x = normalize(x)
        self.y = cross(self.pole, self.x)

    def position(self, t):
        dt = t - self.epoch
        node = self.node0 + self.node_rate * dt
        w = self.w0 + self.w_rate * dt
        big_e = solve_kepler(self.m0 + self.n * dt, self.e)
        px, py = self.a * (math.cos(big_e) - self.e), self.a * math.sqrt(1 - self.e ** 2) * math.sin(big_e)
        # Rz(node) Rx(i) Rz(w) applied to the perifocal position.
        cw, sw, ci, si, cn, sn = math.cos(w), math.sin(w), math.cos(self.i), math.sin(self.i), math.cos(node), math.sin(node)
        x1, y1 = cw * px - sw * py, sw * px + cw * py
        y2, z2 = ci * y1, si * y1
        x3, y3 = cn * x1 - sn * y2, sn * x1 + cn * y2
        return tuple(x3 * self.x[k] + y3 * self.y[k] + z2 * self.pole[k] for k in range(3))


def fit(table, base):
    """Corrections (d_lambda, d_varpi, d_node) in radians and the rms / worst angle (deg)."""
    times = []
    t = table.start() + DAY
    while t < table.end() - DAY:
        times.append(t)
        t += 2.0 * DAY
    actual = [normalize(hb.eval_eph(table, t)) for t in times]

    def model_with(d):
        m = MeanElements.__new__(MeanElements)
        m.__dict__.update(base.__dict__)
        d_lambda, d_varpi, d_node = d
        m.node0 = base.node0 + d_node
        m.w0 = base.w0 + d_varpi - d_node
        m.m0 = base.m0 + d_lambda - d_varpi
        return m

    def cost(d, stats=False):
        m = model_with(d)
        angles = [math.acos(max(-1.0, min(1.0, dot(normalize(m.position(t)), a)))) for t, a in zip(times, actual)]
        rms = math.sqrt(sum(x * x for x in angles) / len(angles))
        return (rms, max(angles)) if stats else rms

    def mean_lambda_offset(d_varpi, d_node):
        m = model_with((0.0, d_varpi, d_node))
        offs = [math.atan2(dot(m.pole, cross(p, a)), dot(p, a)) for p, a in
                ((normalize(m.position(t)), a) for t, a in zip(times, actual))]
        return sum(offs) / len(offs)

    # Coarse grid over the periapsis and node (only where they matter), with the
    # mean longitude set to the average offset; then coordinate descent.
    varpis = [math.radians(k * 15.0) for k in range(24)] if base.e > 0.01 else [0.0]
    nodes = [math.radians(k * 15.0) for k in range(24)] if base.i > math.radians(1.0) else [0.0]
    best = None
    for dv in varpis:
        for dn in nodes:
            d = (mean_lambda_offset(dv, dn), dv, dn)
            c = cost(d)
            if best is None or c < best[0]:
                best = (c, d)
    d = list(best[1])
    step = math.radians(4.0)
    while step > math.radians(0.01):
        improved = False
        for k in range(3):
            if (k == 1 and len(varpis) == 1) or (k == 2 and len(nodes) == 1):
                continue
            for sign in (1.0, -1.0):
                trial = list(d)
                trial[k] += sign * step
                if cost(trial) < cost(d):
                    d = trial
                    improved = True
        if not improved:
            step *= 0.5
    rms, worst = cost(d, stats=True)
    return d, math.degrees(rms), math.degrees(worst)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    scene_path = sys.argv[1]
    only = set(sys.argv[2:])
    with open(scene_path, 'rb') as f:
        scene = tomllib.load(f)
    assets = os.path.dirname(os.path.dirname(os.path.abspath(scene_path)))
    for body in scene.get('bodies', []):
        orbit = body.get('orbit', {})
        fallback = orbit.get('fallback', {})
        if orbit.get('type') != 'ephemeris' or fallback.get('type') != 'mean_elements':
            continue
        if only and body['name'] not in only:
            continue
        table = hb.read_eph(os.path.join(assets, orbit['file']))
        base = MeanElements(fallback)
        (d_lambda, d_varpi, d_node), rms, worst = fit(table, base)
        node = (fallback['node_deg'] + math.degrees(d_node)) % 360.0
        w = (fallback['arg_peri_deg'] + math.degrees(d_varpi - d_node)) % 360.0
        m = (fallback['mean_anomaly_deg'] + math.degrees(d_lambda - d_varpi)) % 360.0
        print(f'{body["name"]}: d_lambda {math.degrees(d_lambda):+.2f}, d_varpi {math.degrees(d_varpi):+.2f}, '
              f'd_node {math.degrees(d_node):+.2f} deg; rms {rms:.2f}, worst {worst:.2f} deg -> '
              f'node_deg = {node:.2f}, arg_peri_deg = {w:.2f}, mean_anomaly_deg = {m:.2f}')


if __name__ == '__main__':
    main()
