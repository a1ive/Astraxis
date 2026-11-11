"""Bake the lines from Io, Europa and Ganymede to their footprints; see assets/geomag/SOURCES.md.

    build/jovian_footprints/venv/Scripts/python tools/geomag/jupiter_flux_tube_bake.py

Requires JupiterMag==1.5.1 (offline only) and uses the same JRM33 + CON2020 as
jupiter_footprint_bake.py. For each moon and each of 72 System III east longitudes, at the
moon's nominal radius in the equatorial plane, one line in each hemisphere from the moon to the
900 km spheroid, resampled to 40 points that correspond from line to line (uniform in the
integral of ds / r, so denser near the planet). The points are turned by minus the seed's
longitude, so that lines interpolate without cutting the corner, and each has w, the fraction
of the Alfven wave's travel time spent by then.

The travel time goes as the integral of sqrt(rho) / B through the plasma sheet of Bagenal &
Delamere 2011 (JGR 116, A05209): electrons n0(rho) exp(-(z / H(rho))^2) around the centrifugal
equator (z as scene/jovian_plasma_sheet.hpp has it, after Phipps & Bagenal 2021), rho the
distance from the spin axis. Two numbers are fitted to Hue 2023's lead angles: one factor on H
for the three moons together, and the mass per electron of each moon.

Output: a line "AXFLUXTUBE2 72 40", then by moon: little-endian float32 x, y, z (in 71492 km),
w by longitude, hemisphere (north, south) and point; then the field at the moon's nominal place
at each longitude, float32 x, y, z in nT, turned like the points.
"""

import hashlib
import sys
from pathlib import Path

import JupiterMag as jm
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jupiter_footprint_bake import A, B, CON, HEIGHT, MOONS, RADII  # noqa: E402

LONGITUDES = 72
POINTS = 40
AXES = np.array(((A + HEIGHT) / A, (A + HEIGHT) / A, (B + HEIGHT) / A))
# Hue 2023 equation 1 and section 4.
SYNODIC_HOURS = (12.89, 11.22, 10.53)
MU0 = 4e-7 * np.pi
AMU = 1.66053907e-27
# VIP4's dipole longitude, as scene/jovian_plasma_sheet.hpp.
DIPOLE_EAST = np.radians(159.2)


# How far the moons go from the nominal radius (fraction) and the equatorial plane (deg): their
# eccentricities and inclinations (JUP365 mean elements) with some margin.
ORBIT_RANGE = ((0.006, 0.1), (0.012, 0.6), (0.003, 0.3))


def hue_lead_deg(moon, south, west_deg):
    """Hue 2023 table 1, as JovianAurora::lead_angle without its blend across 0 degrees."""
    x = np.radians(west_deg)
    c, s = np.cos, np.sin
    if moon == 0:
        return (4.14 - 1.89*c(x) - .70*s(x) + .22*c(2*x) - .12*s(2*x)) if south else \
               (4.26 + 2.64*c(x) + .50*s(x) + .20*c(2*x) + .126*s(2*x))
    if moon == 1:
        return (7.06 - 3.37*c(x) + .15*s(x)) if south else (5.78 + 3.29*c(.99*x) + 1.41*s(.99*x))
    return (12.97 - 5.39*c(.99*x) - 4.50*s(.99*x)) if south else (12.45 + 6.16*c(.98*x) + 3.64*s(.98*x))


def crossing(p, k):
    """The point where the segment k .. k + 1 of a trace meets the spheroid."""
    level = lambda v: np.sum((v / AXES) ** 2) - 1
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if level(p[k] + mid * (p[k+1] - p[k])) * level(p[k]) > 0:
            lo = mid
        else:
            hi = mid
    return hi


def trace(r, lat_deg, lon_deg, step=0.02):
    """For each seed: positions, |B| from the north foot to the south foot, and the seed's index."""
    r, lat, lon = [np.asarray(v, dtype=float).ravel() for v in np.broadcast_arrays(r, lat_deg, lon_deg)]
    lat, lon = np.radians(lat), np.radians(lon)
    seeds = np.column_stack((r * np.cos(lat) * np.cos(lon), r * np.cos(lat) * np.sin(lon), r * np.sin(lat)))
    max_len = int(8 * r.max() / step) + 1000
    t = jm.TraceField(seeds[:, 0], seeds[:, 1], seeds[:, 2], IntModel="jrm33", ExtModel="Con2020", MaxLen=max_len,
                      MaxStep=step, InitStep=step, MinStep=1e-6, ErrMax=1e-8, ionosphereShape="spheroid",
                      ai=AXES[0], bi=AXES[2])
    out = []
    for i, count in enumerate(t.nstep):
        if count >= max_len:
            raise RuntimeError("incomplete trace")
        p = np.column_stack((t.x[i, :count], t.y[i, :count], t.z[i, :count]))
        b = np.sqrt(t.Bx[i, :count] ** 2 + t.By[i, :count] ** 2 + t.Bz[i, :count] ** 2)
        if p[0, 2] < p[-1, 2]:
            p, b = p[::-1], b[::-1]
        level = np.sum((p / AXES) ** 2, axis=1) - 1
        inside = np.flatnonzero(level[:-1] * level[1:] <= 0)
        if len(inside) != 2:
            raise RuntimeError(f"expected two shell crossings, got {len(inside)}")
        first, last = inside
        u, v = crossing(p, first), crossing(p, last)
        line = np.vstack((p[first] + u * (p[first+1] - p[first]), p[first+1:last+1],
                          p[last] + v * (p[last+1] - p[last])))
        strength = np.concatenate(((b[first] + u * (b[first+1] - b[first]),), b[first+1:last+1],
                                   (b[last] + v * (b[last+1] - b[last]),)))
        seed = int(np.argmin(np.sum((line - seeds[i]) ** 2, axis=1)))
        line[seed] = seeds[i]
        out.append((line, strength, seed))
    return out


def arc_lengths(line):
    return np.concatenate(((0.0,), np.cumsum(np.linalg.norm(np.diff(line, axis=0), axis=1))))


def centrifugal_height(points):
    """Height above the centrifugal equator in RJ: centrifugal_equator_height_rj of
    scene/jovian_plasma_sheet.hpp (the tangent approximation after Phipps & Bagenal 2021)."""
    r = np.linalg.norm(points, axis=1)
    east = np.arctan2(points[:, 1], points[:, 0])
    offset = np.radians((1.66 * np.tanh(0.131 * r - 1.62) + 7.76) * np.cos(east - DIPOLE_EAST))
    return r * np.tan(np.arcsin(np.clip(points[:, 2] / r, -1.0, 1.0)) + offset)


def sheet_electrons(rho):
    """Electrons per cm3 on the centrifugal equator, Bagenal & Delamere 2011 (rho in RJ)."""
    x = rho / 6.0
    return 1987.0 * x ** -8.2 + 14.0 * x ** -3.2 + 0.05 * x ** -0.65


def sheet_scale_height(rho):
    """Their scale height H in RJ, of n = n0 exp(-(z / H)^2): 0.74 at Io, 1.7 at Europa, 2.8 at
    Ganymede."""
    x = np.log10(rho / 6.0)
    return 10.0 ** (-0.116 + 2.14 * x - 2.05 * x ** 2 + 0.491 * x ** 3 + 0.126 * x ** 4)


def travel_integrals(line, strength, seed, thickness):
    """Integral of sqrt(n) / B ds (sqrt(cm^-3) m / T) from the seed to each point; n the sheet's
    electrons with its scale height times `thickness`."""
    s = arc_lengths(line)
    rho = np.maximum(np.hypot(line[:, 0], line[:, 1]), 1.0)
    z = centrifugal_height(line)
    rate = np.sqrt(sheet_electrons(rho)) * np.exp(-0.5 * (z / (thickness * sheet_scale_height(rho))) ** 2) / (strength * 1e-9)
    total = np.concatenate(((0.0,), np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(s) * A * 1e3)))
    return total - total[seed]


def lead_fit(moon, lines, thickness):
    """Hue's lead angles at the lines' longitudes (north, then south), the model's with the mass
    per electron that fits them best, and that mass in amu."""
    west = -(np.arange(LONGITUDES) * 360.0 / LONGITUDES) % 360.0
    lead = np.concatenate((hue_lead_deg(moon, False, west), hue_lead_deg(moon, True, west)))
    integrals = np.empty(2 * LONGITUDES)
    for k, (line, strength, seed) in enumerate(lines):
        total = travel_integrals(line, strength, seed, thickness)
        integrals[k] = -total[0]
        integrals[LONGITUDES + k] = total[-1]
    scale = np.dot(lead, integrals) / np.dot(integrals, integrals)
    # lead / 360 * P_syn = sqrt(mu0 m n) * integral, n per cm3.
    mass = (scale * SYNODIC_HOURS[moon] * 3600.0 / 360.0) ** 2 / MU0 / AMU * 1e-6
    return lead, scale * integrals, mass


def hemisphere_points(line, integral, seed, south):
    """40 points from the seed to the foot, xyz in RJ and w, the fraction of the travel time."""
    if south:
        path, time = line[seed:], integral[seed:]
    else:
        path, time = line[seed::-1], -integral[seed::-1]
    steps = np.linalg.norm(np.diff(path, axis=0), axis=1)
    middle = 0.5 * (np.linalg.norm(path[1:], axis=1) + np.linalg.norm(path[:-1], axis=1))
    sigma = np.concatenate(((0.0,), np.cumsum(steps / middle)))
    grid = np.linspace(0.0, sigma[-1], POINTS)
    out = np.empty((POINTS, 4))
    for c in range(3):
        out[:, c] = np.interp(grid, sigma, path[:, c])
    out[:, 3] = np.interp(grid, sigma, time) / time[-1]
    out[0, :3], out[-1, :3] = path[0], path[-1]
    return out, time[-1]


def rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array(((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0)))


def table_line(table, moon_table_lon_deg, south):
    """Points of the table's line at a longitude, not yet turned back: as JovianAurora."""
    x = (moon_table_lon_deg % 360.0) * LONGITUDES / 360.0
    i = int(np.floor(x)) % LONGITUDES
    t = x - np.floor(x)
    return (1.0 - t) * table[i, int(south)] + t * table[(i + 1) % LONGITUDES, int(south)]


def miss_along(points, polyline):
    """The largest distance from the points to a polyline."""
    a, b = polyline[:-1], polyline[1:]
    d = b - a
    worst = 0.0
    for p in points:
        t = np.clip(np.sum((p - a) * d, axis=1) / np.maximum(np.sum(d * d, axis=1), 1e-30), 0.0, 1.0)
        worst = max(worst, np.sqrt(np.min(np.sum((a + t[:, None] * d - p) ** 2, axis=1))))
    return worst


def warped_line(table, seed, foot, south):
    """The renderer's line through `seed`: the nominal one, moved at the moon and at the foot.
    Keep in step with JovianAurora::flux_tube."""
    lon = np.degrees(np.arctan2(seed[1], seed[0]))
    turn = rotation(np.radians(lon))
    points = table_line(table, lon, south)[:, :3] @ turn.T
    r = np.linalg.norm(points, axis=1)
    near = (r - r[-1]) / (r[0] - r[-1])
    return points + near[:, None] * (seed - points[0]) + ((1.0 - near) ** 8)[:, None] * (foot - points[-1])


def main():
    if jm.__version__ != "1.5.1":
        raise RuntimeError("this bake requires JupiterMag 1.5.1")
    jm.Internal.Config(Model="jrm33", Degree=13)
    jm.Con2020.Config(**CON)
    lons = np.arange(LONGITUDES) * 360.0 / LONGITUDES
    payload = [f"AXFLUXTUBE2 {LONGITUDES} {POINTS}\n".encode("ascii")]
    rng = np.random.default_rng(20261105)
    traced = [trace(np.full(LONGITUDES, nominal), 0.0, lons) for nominal in RADII]

    # One thickness for the three moons: the least sum of their squared relative misses.
    def cost(thickness):
        total = 0.0
        for moon in range(3):
            lead, modelled, _ = lead_fit(moon, traced[moon], thickness)
            total += np.mean((modelled - lead) ** 2) / np.mean(lead) ** 2
        return total
    candidates = np.arange(0.8, 2.51, 0.05)
    thickness = float(candidates[int(np.argmin([cost(c) for c in candidates]))])
    print(f"scale heights of Bagenal & Delamere 2011 times {thickness:.2f}", flush=True)

    for moon, (name, nominal) in enumerate(zip(MOONS, RADII)):
        lines = traced[moon]
        lead, modelled, mass = lead_fit(moon, lines, thickness)
        _, unscaled, unscaled_mass = lead_fit(moon, lines, 1.0)
        print(f"{name}: lead angle rms miss {np.sqrt(np.mean((modelled - lead) ** 2)):.2f} deg, correlation "
              f"{np.corrcoef(lead, modelled)[0, 1]:.2f}, Hue {lead.min():.1f} .. {lead.max():.1f} deg, modelled "
              f"{modelled.min():.1f} .. {modelled.max():.1f} deg; {mass:.1f} amu per electron, {sheet_electrons(nominal):.3g} "
              f"electrons/cm3 and H = {thickness * sheet_scale_height(nominal):.2f} RJ at the moon "
              f"(with their H unchanged: rms {np.sqrt(np.mean((unscaled - lead) ** 2)):.2f} deg, {unscaled_mass:.1f} amu)", flush=True)

        table = np.empty((LONGITUDES, 2, POINTS, 4))
        for k, (line, strength, seed) in enumerate(lines):
            total = travel_integrals(line, strength, seed, thickness)
            turn = rotation(-np.radians(lons[k]))
            for south in (False, True):
                points, _ = hemisphere_points(line, total, seed, south)
                points[:, :3] = points[:, :3] @ turn.T
                table[k, int(south)] = points
        payload.append(table.astype("<f4").tobytes())

        # The field the moon sits in (nT), for its own field's surroundings (Ganymede).
        seeds = nominal * np.column_stack((np.cos(np.radians(lons)), np.sin(np.radians(lons)), np.zeros(LONGITUDES)))
        field = np.column_stack(jm.Internal.Field(seeds[:, 0], seeds[:, 1], seeds[:, 2])) + \
            np.column_stack(jm.Con2020.Field(seeds[:, 0], seeds[:, 1], seeds[:, 2]))
        strength = np.linalg.norm(field, axis=1)
        radial = np.sum(field * seeds, axis=1) / nominal
        print(f"  field at the moon: {strength.min():.1f} .. {strength.max():.1f} nT, Bz {field[:, 2].min():.1f} .. "
              f"{field[:, 2].max():.1f}, radial {radial.min():.1f} .. {radial.max():.1f}", flush=True)
        # Along the traced line at the seed, as it must be.
        for k, (line, _, seed) in enumerate(lines):
            along = line[seed + 1] - line[seed - 1]
            cosine = abs(np.dot(along, field[k])) / (np.linalg.norm(along) * strength[k])
            if cosine < 0.9999:
                raise RuntimeError(f"field at the seed is off the traced line: {cosine}")
        payload.append(np.stack([field[k] @ rotation(-np.radians(lons[k])).T for k in range(LONGITUDES)]).astype("<f4").tobytes())

        # Withheld moons where the orbit takes them, against a direct trace.
        r = nominal * (1.0 + rng.uniform(-1.0, 1.0, 60) * ORBIT_RANGE[moon][0])
        lat = rng.uniform(-1.0, 1.0, 60) * ORBIT_RANGE[moon][1]
        lon = rng.uniform(0.0, 360.0, 60)
        lon[:2] = (0.001, 359.999)
        exact_lines = trace(r, lat, lon, step=0.01)
        worst = worst_nominal = 0.0
        for k, (line, strength, seed) in enumerate(exact_lines):
            for south in (False, True):
                exact = line[seed:] if south else line[seed::-1]
                warped = warped_line(table, exact[0], exact[-1], south)
                worst = max(worst, miss_along(warped, exact), miss_along(exact[::10], warped))
                nominal_line = table_line(table, lon[k], south)[:, :3] @ rotation(np.radians(lon[k])).T
                worst_nominal = max(worst_nominal, miss_along(nominal_line, exact))
        print(f"  withheld lines: largest miss {worst * A:.0f} km, {worst_nominal * A:.0f} km without moving the "
              f"line to the moon and the foot", flush=True)

    out = Path("assets/geomag/jupiter_flux_tubes.bin")
    data = b"".join(payload)
    out.write_bytes(data)
    print(f"{out}: {len(data)} bytes, sha256 {hashlib.sha256(data).hexdigest()}")


if __name__ == "__main__":
    main()
