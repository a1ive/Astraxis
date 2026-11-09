"""Numbers for assets/scenes/castor.toml (standard library only).

Usage:
    python -I tools/castor/scene_numbers.py

Prints the masses as GM, the barycentric scales, Castor C's orbit around AB,
and how well the published orbits reproduce what they were fitted to (the CHARA
positions and light-time corrections of Torres et al. 2022, Table 2).

Castor C's orbit is not known (Torres et al. 2022, Sect. 4.3). The scene uses
the one orbit that has their measured position and sky motion, their estimated
period, no line-of-sight velocity relative to AB, and C behind AB.
"""

import math

GM_SUN = 1.32712440041279419e11  # km^3/s^2, JPL DE440
AU_KM = 149597870.7
PC_AU = 648000.0 / math.pi
C_AU_DAY = 299792.458 * 86400.0 / AU_KM
JD_J2000 = 2451545.0
YEAR_DAYS = 365.25

# Torres et al. 2022, ApJ 941, 8, Tables 5 and 6.
DISTANCE_PC = 15.0703
M_AA, M_AB_STAR, M_BA, M_BB = 2.371, 0.3859, 1.789, 0.3865
M_A, M_B = M_AA + M_AB_STAR, M_BA + M_BB
OUTER = dict(P=459.1 * YEAR_DAYS, e=0.3382, i=115.107, node=41.304, w=251.84,
             T=JD_J2000 + (1959.59 - 2000.0) * YEAR_DAYS, M=M_A + M_B)
# omega_Aa is the primary's: the companion's is 180 deg away. Castor B is circular
# and T_B is Ba's maximum velocity: the companion is then at u = 180 deg.
INNER_A = dict(P=9.2127496, e=0.48769, i=35.00, node=95.100, w=264.968 + 180.0, T=2455817.7868, M=M_A)
INNER_B = dict(P=2.92835083, e=0.0, i=110.50, node=106.47, w=180.0, T=2456705.4942, M=M_B)

# Torres & Ribas 2002, ApJ 567, 1140, Sect. 5.
M_CA, M_CB = 0.5975, 0.6009
R_C, T_C = 0.6191, 3820.0
M_C = M_CA + M_CB

# Torres et al. 2022, Table 2: BJD - 2400000, light-time correction (d), rho (mas), theta (deg).
CHARA_A = [(54423.9754, -0.1121, 10.100, 37.56), (54426.9056, -0.1121, 3.482, 197.97),
           (54427.9267, -0.1121, 6.791, 280.64), (59275.6936, -0.0778, 9.475, 329.54),
           (59279.6593, -0.0777, 8.441, 49.01), (59538.0184, -0.0758, 7.934, 59.93),
           (59557.0147, -0.0756, 6.788, 78.93), (59568.8842, -0.0755, 7.419, 290.63)]
CHARA_B = [(59275.7346, 0.0988, 1.221, 24.43), (59279.6873, 0.0987, 2.986, 274.54),
           (59537.9664, 0.0962, 1.391, 160.78), (59557.0657, 0.0960, 1.786, 324.09),
           (59568.9552, 0.0959, 2.549, 306.31), (59569.0202, 0.0959, 2.843, 300.95)]

# SIMBAD, ICRS at J2000: Castor A (Hipparcos) and YY Gem (Gaia DR3).
RA_A, DEC_A = 113.64942845, 31.888276192
RA_C, DEC_C = 113.65603096981, 31.86949434486
# Torres et al. 2022, Sect. 4.2: C relative to AB, mas/yr.
DMU_E, DMU_N = -26.21, 2.88
PERIOD_C_YEARS = 14500.0     # their estimate


def semi_major_au(period_days, mass):
    return (mass * (period_days / YEAR_DAYS) ** 2) ** (1.0 / 3.0)


def kepler(o, jd):
    """(north, east, depth) of the companion relative to the primary, au."""
    a = semi_major_au(o['P'], o['M'])
    M = 2.0 * math.pi * (jd - o['T']) / o['P']
    M = math.remainder(M, 2.0 * math.pi)
    E = M + o['e'] * math.sin(M)
    for _ in range(50):
        E -= (E - o['e'] * math.sin(E) - M) / (1.0 - o['e'] * math.cos(E))
    nu = 2.0 * math.atan2(math.sqrt(1.0 + o['e']) * math.sin(E / 2.0), math.sqrt(1.0 - o['e']) * math.cos(E / 2.0))
    r = a * (1.0 - o['e'] * math.cos(E))
    u = nu + math.radians(o['w'])
    node, inc = math.radians(o['node']), math.radians(o['i'])
    return (r * (math.cos(u) * math.cos(node) - math.sin(u) * math.sin(node) * math.cos(inc)),
            r * (math.cos(u) * math.sin(node) + math.sin(u) * math.cos(node) * math.cos(inc)),
            r * math.sin(u) * math.sin(inc))


def sky(p):
    """Separation (arcsec) and position angle (deg) of (north, east, depth) in au."""
    rho = math.hypot(p[0], p[1]) / DISTANCE_PC
    return rho, math.degrees(math.atan2(p[1], p[0])) % 360.0


def compare_chara(name, orbit, table):
    print('CHARA ' + name + ': d(sep) mas, d(PA) deg, tangential mas')
    total = 0.0
    for bjd, dt, rho, theta in table:
        p = kepler(orbit, 2400000.0 + bjd + dt)
        sep, pa = sky(p)
        sep *= 1000.0
        d_pa = math.remainder(pa - theta, 360.0)
        tangential = rho * math.radians(d_pa)
        total += (sep - rho) ** 2 + tangential ** 2
        print('  %.4f  %+.3f  %+.2f  %+.3f' % (bjd, sep - rho, d_pa, tangential))
    print('  rms %.3f mas' % math.sqrt(total / len(table)))


def orbit_from_state(p, v, mass, jd):
    """Elements of a state (north, east, depth), au and au/day, in VisualOrbit's convention."""
    mu = mass * (2.0 * math.pi / YEAR_DAYS) ** 2
    r = math.sqrt(sum(x * x for x in p))
    h = (p[1] * v[2] - p[2] * v[1], p[2] * v[0] - p[0] * v[2], p[0] * v[1] - p[1] * v[0])
    h_len = math.sqrt(sum(x * x for x in h))
    inc = math.acos(h[2] / h_len)
    node = math.atan2(h[0], -h[1])
    energy = sum(x * x for x in v) / 2.0 - mu / r
    a = -mu / (2.0 * energy)
    ev = [((sum(x * x for x in v) - mu / r) * p[k] - sum(p[j] * v[j] for j in range(3)) * v[k]) / mu for k in range(3)]
    e = math.sqrt(sum(x * x for x in ev))
    u0 = (math.cos(node), math.sin(node), 0.0)
    u90 = (-math.sin(node) * math.cos(inc), math.cos(node) * math.cos(inc), math.sin(inc))
    u = math.atan2(sum(p[k] * u90[k] for k in range(3)), sum(p[k] * u0[k] for k in range(3)))
    w = math.atan2(sum(ev[k] * u90[k] for k in range(3)), sum(ev[k] * u0[k] for k in range(3)))
    nu = u - w
    E = 2.0 * math.atan2(math.sqrt(1.0 - e) * math.sin(nu / 2.0), math.sqrt(1.0 + e) * math.cos(nu / 2.0))
    M = E - e * math.sin(E)
    n = math.sqrt(mu / a ** 3)
    return dict(a=a, e=e, i=math.degrees(inc), node=math.degrees(node) % 360.0, w=math.degrees(w) % 360.0,
                T=jd - M / n, P=2.0 * math.pi / n)


def castor_c():
    b_a = kepler(OUTER, JD_J2000)
    share = M_B / OUTER['M']
    cos_dec = math.cos(math.radians(DEC_A))
    east = (RA_C - RA_A) * cos_dec * 3600.0 - share * b_a[1] / DISTANCE_PC
    north = (DEC_C - DEC_A) * 3600.0 - share * b_a[0] / DISTANCE_PC
    print('C from AB at J2000: %.2f" PA %.2f deg' % (math.hypot(east, north),
                                                     math.degrees(math.atan2(east, north)) % 360.0))
    d = DISTANCE_PC
    sky_p = (north * d, east * d)
    v = (DMU_N * 1e-3 * d / YEAR_DAYS, DMU_E * 1e-3 * d / YEAR_DAYS, 0.0)
    print('C sky speed %.2f km/s' % (math.hypot(v[0], v[1]) * AU_KM / 86400.0))
    mass = OUTER['M'] + M_C
    # Bisect the depth for the published period (C behind AB).
    lo, hi = 0.0, 3000.0
    for _ in range(80):
        z = 0.5 * (lo + hi)
        o = orbit_from_state((sky_p[0], sky_p[1], z), v, mass, JD_J2000)
        if o['P'] / YEAR_DAYS < PERIOD_C_YEARS:
            lo = z
        else:
            hi = z
    print('C depth %.1f au, a %.2f au, e %.4f, i %.3f, node %.3f, w %.3f, T year %.2f, P %.1f yr'
          % (z, o['a'], o['e'], o['i'], o['node'], o['w'], 2000.0 + (o['T'] - JD_J2000) / YEAR_DAYS,
             o['P'] / YEAR_DAYS))
    return mass


def main():
    print('Kepler: a_A %.4f mas (table 8.002), a_B %.4f mas (table 3.4442), a_AB %.4f" (table 6.722)'
          % (1000.0 * semi_major_au(INNER_A['P'], M_A) / DISTANCE_PC,
             1000.0 * semi_major_au(INNER_B['P'], M_B) / DISTANCE_PC,
             semi_major_au(OUTER['P'], OUTER['M']) / DISTANCE_PC))
    compare_chara('A', INNER_A, CHARA_A)
    compare_chara('B', INNER_B, CHARA_B)
    print('light time (d): epoch, A model / table, B model / table')
    for bjd, table_a, table_b in [(54423.9754, -0.1121, None), (59275.7, -0.0778, 0.0988), (59569.0, -0.0755, 0.0959)]:
        p = kepler(OUTER, 2400000.0 + bjd)
        z_a = -M_B / OUTER['M'] * p[2]
        z_b = M_A / OUTER['M'] * p[2]
        print('  %.1f  A %+.4f / %+.4f  B %+.4f / %s' % (bjd, -z_a / C_AU_DAY, table_a, -z_b / C_AU_DAY, table_b))
    mass = castor_c()
    print()
    masses = [('A/B outer', OUTER['M']), ('A pair', M_A), ('B pair', M_B), ('C pair', M_C), ('AB+C', mass)]
    for name, m in masses:
        print('GM %-10s %.9e  (%.4f M_sun)' % (name, m * GM_SUN, m))
    print('scale A in AB  %.9f   B in AB  %.9f' % (-M_B / OUTER['M'], M_A / OUTER['M']))
    print('scale Aa in A  %.9f   Ab in A  %.9f' % (-M_AB_STAR / M_A, M_AA / M_A))
    print('scale Ba in B  %.9f   Bb in B  %.9f' % (-M_BB / M_B, M_BA / M_B))
    print('scale AB       %.9f   C        %.9f' % (-M_C / mass, OUTER['M'] / mass))
    print('scale Ca in C  %.9f   Cb in C  %.9f' % (-M_CB / M_C, M_CA / M_C))
    print('L_C  %.4f' % (R_C ** 2 * (T_C / 5772.0) ** 4))
    print('sky of B from A: year, sep ", PA, 3D au')
    for year in (1718.0, 1778.0, 1803.0, 1826.0, 1959.59, 1970.0, 2000.0, 2026.85):
        p = kepler(OUTER, JD_J2000 + (year - 2000.0) * YEAR_DAYS)
        rho, pa = sky(p)
        print('  %.2f  %.2f  %.1f  %.1f' % (year, rho, pa, math.sqrt(sum(x * x for x in p))))


main()
