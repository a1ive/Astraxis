"""Turn the kOS dump of the stock Kerbol system into scene data.

    python kerbol_scene.py [kerbol_bodies.txt] > out.txt

kerbol_bodies.txt is the output of kerbol_dump.ks (KSP 1.12.5, kOS 1.6.0.1),
kept next to this script. The tool

1. checks the conventions: Kepler positions from the dumped elements, in the
   frame below, must reproduce the dumped positions; periods follow from the
   parent's mu alone; rotation is linear in UT;
2. prints the [[bodies]] blocks for assets/scenes/kerbol.toml (orbits, sizes,
   GMs, rotation; colors and atmospheres are added by hand);
3. prints each body's rotation angle and position relative to its parent at its
   sample UT in the scene frame (km), for tests/planet_scene_tests.cpp;
4. prints the atmosphere numbers used for [bodies.atmosphere].

Frame. kOS vectors are in Unity's left-handed frame with y toward the system's
north. Swapping y and z gives a right-handed frame; the scene's ICRF axes are
x = the solar prime vector (KSP's reference direction, from which longitudes of
ascending nodes and rotation angles are measured) and z = north. All KSP orbits,
moons included, are referred to this one plane: there are no local Laplace
planes and no axial tilts.

Time. KSP's UT (seconds since the start of a game) is used as TDB seconds since
J2000, so the elements' epoch (UT 0) is J2000.

Standard library only.
"""

import math
import os
import sys

# Earth's air, the reference for the Rayleigh depths: scene solar_system.toml [ATM]
# (Bruneton: 13.56e-6 m^-1 at 550 nm x 8 km = 0.1085 at 1 atm), at standard gravity.
EARTH_RAYLEIGH_DEPTH = 0.1085
STANDARD_GRAVITY = 9.80665
AIR_MOLAR_MASS = 0.0289644  # kg/mol, the value KSP gives Kerbin's and Laythe's air

# Scene order and names (KSP calls Kerbol "Sun").
ORDER = ["Sun", "Moho", "Eve", "Gilly", "Kerbin", "Mun", "Minmus", "Duna", "Ike", "Dres",
         "Jool", "Laythe", "Vall", "Tylo", "Bop", "Pol", "Eeloo"]
SCENE_NAME = {"Sun": "Kerbol"}


def parse(path):
    bodies, cur = [], None
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if line == "[[body]]":
            cur = {"atm_sample": []}
            bodies.append(cur)
            continue
        if not line or line.startswith("#") or "=" not in line or cur is None:
            continue
        key, text = (s.strip() for s in line.split("=", 1))
        if text.startswith("["):
            value = [float(x) for x in text.strip("[]").split(",")]
        elif text.startswith('"'):
            value = text.strip('"')
        elif text in ("True", "False"):
            value = text == "True"
        else:
            value = float(text)
        if key == "atm_sample":
            cur[key].append(value)
        else:
            cur[key] = value
    return {b["name"]: b for b in bodies}


def zup(v):
    return (v[0], v[2], v[1])


def scene_vector(b, v):
    """A raw kOS vector in the scene frame (x = solar prime vector, z = north)."""
    p = zup(v)
    s = zup(b["solar_prime_vector_raw"])
    ex = (s[0], s[1])
    ey = (-s[1], s[0])
    return (p[0] * ex[0] + p[1] * ex[1], p[0] * ey[0] + p[1] * ey[1], p[2])


def kepler_position(b, mu, t):
    a, e = b["semi_major_axis_m"], b["eccentricity"]
    i, node, w = (math.radians(b[k]) for k in ("inclination_deg", "lan_deg", "arg_periapsis_deg"))
    m = math.radians(b["mean_anomaly_at_epoch_deg"]) + math.sqrt(mu / a**3) * (t - b["epoch_s"])
    ea = m
    for _ in range(50):
        ea -= (ea - e * math.sin(ea) - m) / (1 - e * math.cos(ea))
    xp, yp = a * (math.cos(ea) - e), a * math.sqrt(1 - e * e) * math.sin(ea)
    cn, sn, ci, si, cw, sw = math.cos(node), math.sin(node), math.cos(i), math.sin(i), math.cos(w), math.sin(w)
    return ((cn * cw - sn * sw * ci) * xp + (-cn * sw - sn * cw * ci) * yp,
            (sn * cw + cn * sw * ci) * xp + (-sn * sw + cn * cw * ci) * yp,
            sw * si * xp + cw * si * yp)


def initial_rotation(b):
    t = b["sample_ut_start"]
    return (b["rotation_angle_deg"] - 360.0 * t / b["rotation_period_s"]) % 360.0


def check(by):
    for name in ORDER:
        b = by[name]
        assert b["sample_ut_start"] == b["sample_ut_end"], name + ": sample spans two ticks"
        w0 = initial_rotation(b)
        assert abs(w0 - round(w0)) < 1e-6 or abs(w0 - 360.0) < 1e-6, name + ": rotation not linear in UT"
        if "parent" not in b:
            continue
        mu = by[b["parent"]]["mu_m3_s2"]
        period = 2 * math.pi * math.sqrt(b["semi_major_axis_m"] ** 3 / mu)
        assert abs(period - b["period_s"]) < 1e-6, name + ": period is not from the parent's mu"
        predicted = kepler_position(b, mu, b["sample_ut_start"])
        actual = scene_vector(b, b["position_rel_parent_raw_m"])
        error = math.dist(predicted, actual)
        assert error < 1.0, f"{name}: Kepler position off by {error} m"
    print("# conventions check: ok", file=sys.stderr)


def fmt(x):
    return repr(float(x))


def body_block(by, name):
    b = by[name]
    out = ["[[bodies]]", f'name = "{SCENE_NAME.get(name, name)}"']
    if "parent" in b:
        out.append(f'parent = "{SCENE_NAME.get(b["parent"], b["parent"])}"')
    else:
        out.append('kind = "star"')
    out.append(f"radii_km = [{fmt(b['radius_m'] / 1000.0)}]")
    out.append(f"gm_km3_s2 = {fmt(b['mu_m3_s2'] / 1e9)}")
    out.append("pole_ra_dec_deg = [270.0, 90.0]")
    w0 = initial_rotation(b)
    w0 = 0.0 if abs(w0 - 360.0) < 1e-6 else round(w0, 6)
    out.append(f"prime_meridian_deg = [{fmt(w0)}, {fmt(360.0 * 86400.0 / b['rotation_period_s'])}]")
    if "parent" in b:
        out += ["[bodies.orbit]", 'type = "mean_elements"',
                f"a_km = {fmt(b['semi_major_axis_m'] / 1000.0)}",
                f"e = {fmt(b['eccentricity'])}",
                f"i_deg = {fmt(b['inclination_deg'])}",
                f"node_deg = {fmt(b['lan_deg'])}",
                f"arg_peri_deg = {fmt(b['arg_periapsis_deg'])}",
                f"mean_anomaly_deg = {fmt(b['mean_anomaly_at_epoch_deg'])}",
                f"period_days = {fmt(b['period_s'] / 86400.0)}"]
    return "\n".join(out)


def atmosphere(by, name):
    b = by[name]
    samples = b["atm_sample"]
    p0, t0 = samples[0][1], samples[0][2]
    n0 = p0 / t0
    # e-folding height of the number density (p / T) above the surface
    h_n = next(h for h, p, t in samples if t > 0 and p / t < n0 / math.e)
    g = b["mu_m3_s2"] / b["radius_m"] ** 2
    column = p0 * (STANDARD_GRAVITY / g) * (AIR_MOLAR_MASS / b["atm_molar_mass_kg_mol"])
    return (f"# {name}: top {b['atm_height_m'] / 1000:g} km; p0 {p0:.6g} atm (curve; ASL value "
            f"{b['atm_sea_level_pressure_atm']:.6g}), T0 {t0:.4g} K, M {b['atm_molar_mass_kg_mol']:.6g} kg/mol, "
            f"g {g:.5g} m/s^2; number density 1/e at {h_n / 1000:.2f} km; column / Earth's = {column:.4g}, "
            f"Rayleigh depth with air's cross-section = {EARTH_RAYLEIGH_DEPTH * column:.4g}")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                               "kerbol_bodies.txt")
    by = parse(path)
    check(by)
    for name in ORDER:
        print(body_block(by, name))
        print()
    print("# Sample UT, rotation angle (deg) and position relative to the parent (scene frame, km)")
    for name in ORDER:
        b = by[name]
        p = scene_vector(b, b.get("position_rel_parent_raw_m", [0.0, 0.0, 0.0]))
        print(f'{{"{SCENE_NAME.get(name, name)}", {b["sample_ut_start"]!r}, {b["rotation_angle_deg"]!r}, '
              f'{{{p[0] / 1000:.6f}, {p[1] / 1000:.6f}, {p[2] / 1000:.6f}}}}},')
    print()
    k = by["Kerbin"]
    print(f"# Kerbin solar day: {1 / (1 / k['rotation_period_s'] - 1 / k['period_s'])!r} s")
    for name in ORDER:
        if by[name]["atm_sample"] and name != "Sun":
            print(atmosphere(by, name))


if __name__ == "__main__":
    main()
