// Tests for the planets and their moons in the scenes: Solar System, Jupiter,
// Earth-Moon, Saturn (with the Janus-Epimetheus swaps), Mars' rotation, and the
// fictional Kerbol system.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/element_table.hpp"
#include "ephem/ephemeris.hpp"
#include "ephem/kepler.hpp"
#include "ephem/mean_element_orbit.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace astraxis;

namespace {

// The mean-element orbit behind a motion's fallback chain (ephemeris -> ... -> mean elements).
const MeanElementOrbit* mean_elements_of(const MotionSource* motion)
{
    while (motion) {
        if (const auto* mean = dynamic_cast<const MeanElementOrbit*>(motion)) {
            return mean;
        }
        const auto* eph = dynamic_cast<const EphemerisMotion*>(motion);
        motion = eph ? eph->fallback() : nullptr;
    }
    return nullptr;
}

// Mean longitude rate of each Galilean moon vs. the IAU synchronous rotation rate.
void test_galilean_rates()
{
    Scene scene = load_jupiter();
    const double pck_rates[] = {203.4889538, 101.3747235, 50.3176081, 21.5710715}; // deg/day
    double lambda0[4] = {};

    for (int m = 0; m < 4; ++m) {
        const MeanElementOrbit* orbit = mean_elements_of(scene.bodies[m + 1].motion.get());
        check(orbit != nullptr, "moon falls back to MeanElementOrbit");
        if (!orbit) {
            continue;
        }
        const double span_days = 1000.0;
        const KeplerElements a = orbit->elements_at(0.0);
        const KeplerElements b = orbit->elements_at(span_days * kSecondsPerDay);
        lambda0[m] = a.node + a.arg_peri + a.mean_anomaly;

        // Unwrapped longitude advance over the span.
        const double revolutions = std::floor(pck_rates[m] * span_days / 360.0);
        const double dl = wrap_two_pi((b.node + b.arg_peri + b.mean_anomaly) - lambda0[m]) + revolutions * kTwoPi;
        const double rate = dl * kRadToDeg / span_days;
        check(std::abs(rate - pck_rates[m]) < 0.002, "mean longitude rate matches PCK", rate);

        // Distance from Jupiter stays within the ellipse bounds.
        const double r = glm::length(orbit->eval(12345.0 * kSecondsPerDay).position);
        check(r > a.a * (1.0 - a.e) - 1.0 && r < a.a * (1.0 + a.e) + 1.0, "moon distance within [q, Q]", r);
    }

    // Laplace resonance: lambda_Io - 3 lambda_Europa + 2 lambda_Ganymede = 180 deg.
    const double laplace = wrap_pi(lambda0[0] - 3.0 * lambda0[1] + 2.0 * lambda0[2] - kPi) * kRadToDeg;
    check(std::abs(laplace) < 1.0, "Laplace resonance at J2000", laplace);

    // Inner moons: with apsides advancing and nodes regressing (Jupiter's
    // oblateness), the mean longitude rate matches the PCK synchronous rotation.
    // The periods of [ELEM] have 3 significant digits for Amalthea, hence 0.01.
    const char* inner[] = {"Amalthea", "Thebe", "Adrastea", "Metis"};
    const double inner_pck[] = {722.6314560, 533.7004100, 1206.9986602, 1221.2547301};
    for (int m = 0; m < 4; ++m) {
        const int b = scene.find(inner[m]);
        const MeanElementOrbit* orbit = b > 0 ? mean_elements_of(scene.bodies[static_cast<size_t>(b)].motion.get()) : nullptr;
        check(orbit != nullptr, "inner moon uses MeanElementOrbit");
        if (!orbit) {
            continue;
        }
        const double span_days = 100.0;
        const KeplerElements a = orbit->elements_at(0.0);
        const KeplerElements c = orbit->elements_at(span_days * kSecondsPerDay);
        const double revolutions = std::floor(inner_pck[m] * span_days / 360.0);
        const double dl = wrap_two_pi((c.node + c.arg_peri + c.mean_anomaly) - (a.node + a.arg_peri + a.mean_anomaly)) +
                          revolutions * kTwoPi;
        check(std::abs(dl * kRadToDeg / span_days - inner_pck[m]) < 0.01, "inner moon longitude rate matches PCK",
              dl * kRadToDeg / span_days);
    }
}

void test_jupiter_heliocentric()
{
    Scene scene = load_jupiter();
    // Perihelion ~4.95 au, aphelion ~5.46 au.
    for (int year = -50; year <= 50; ++year) {
        const double t = year * kDaysPerJulianYear * kSecondsPerDay;
        scene.update(t);
        const double r_au = glm::length(scene.sun_position()) / kAuKm;
        check(r_au > 4.9 && r_au < 5.5, "Jupiter heliocentric distance", r_au);
    }
}

// Jupiter scene: the inner moons against Horizons, and Io's blend from the mean
// elements into the Galileo-span table.
void test_jupiter_moon_ephemerides()
{
    Scene scene = load_jupiter();
    const int io = scene.find("Io");
    auto at = [&](int body, double t) { return scene.icrf_state_at(body, t).position; };
    double t = 0.0;
    parse_utc("2026-10-02", &t);

    // Inner moons (mean elements, period from the PCK spin rate) against Horizons
    // (@ 500@599, 2026-Oct-02 00:01:09 TDB): within a few degrees after 26 years.
    const char* inner[] = {"Amalthea", "Thebe", "Adrastea", "Metis"};
    const glm::dvec3 inner_hzn[] = {{-1.627599848191866E+05, -6.953652662401732E+04, -3.708388618213439E+04},
                                    {-5.620351095369970E+04, 1.969718451580069E+05, 9.037534633891605E+04},
                                    {9.314251567946289E+04, -8.105204208221256E+04, -3.714128408138068E+04},
                                    {-7.473411994107516E+04, -9.321647784889850E+04, -4.566251549820685E+04}};
    const double t_inner = t + 69.0; // TDB - UTC
    double worst_inner = 0.0;
    for (int m = 0; m < 4; ++m) {
        const glm::dvec3 p = at(scene.find(inner[m]), t_inner);
        const double c = glm::dot(p, inner_hzn[m]) / (glm::length(p) * glm::length(inner_hzn[m]));
        worst_inner = std::max(worst_inner, std::acos(std::clamp(c, -1.0, 1.0)) * kRadToDeg);
    }
    check(worst_inner < 3.0, "inner moons near Horizons in 2026 (deg)", worst_inner);
    std::printf("info: inner moons vs Horizons in 2026: worst %.2f deg\n", worst_inner);

    // Io fades from the mean elements into the Galileo-span table over a day: no
    // jump between one-minute steps beyond Io's own motion (~1,000 km/min).
    const auto* eph = dynamic_cast<const EphemerisMotion*>(scene.bodies[static_cast<size_t>(io)].motion.get());
    const auto* gll = eph ? dynamic_cast<const EphemerisMotion*>(eph->fallback()) : nullptr;
    check(gll != nullptr, "Io: Juno-span table, then Galileo-span table");
    if (gll) {
        const double t_start = gll->table().start();
        double worst_jump = 0.0;
        for (double tk = t_start - 2.0 * kSecondsPerDay; tk < t_start + 2.0 * kSecondsPerDay; tk += 60.0) {
            const State a = scene.icrf_state_at(io, tk);
            const glm::dvec3 b = at(io, tk + 60.0);
            worst_jump = std::max(worst_jump, glm::length(b - a.position - a.velocity * 60.0));
        }
        check(worst_jump < 50.0, "Io blends into the table without a jump (km/min)", worst_jump);
    }
}

// Io's shadow falls on Jupiter once per orbit (Io's orbit is nearly in Jupiter's equator
// and Jupiter's obliquity is only ~3 deg, so every conjunction with the sun line is a transit).
void test_io_shadow_transits()
{
    Scene scene = load_jupiter();
    double t0 = 0.0;
    parse_utc("2026-10-02 00:00", &t0);

    const double radius = scene.bodies[0].equatorial_radius_km;
    int transits = 0;
    bool in_transit = false;
    for (double t = t0; t < t0 + 10.0 * kSecondsPerDay; t += 300.0) {
        scene.update(t);
        const glm::dvec3 io = scene.bodies[1].world_position;
        const glm::dvec3 sun = glm::normalize(scene.sun_position());
        const double along = glm::dot(io, sun);
        const double perp = glm::length(io - along * sun);
        const bool now = along > 0.0 && perp < radius;
        if (now && !in_transit) {
            ++transits;
            if (transits == 1) {
                std::printf("info: first Io shadow transit after 2026-10-02: %s\n", format_utc(t).c_str());
            }
        }
        in_transit = now;
    }
    check(transits >= 5 && transits <= 6, "Io shadow transits in 10 days", transits);
}

// The Galilean moons are tidally locked: the IAU prime meridian (body-fixed +x,
// longitude 0) must point at Jupiter. This cross-checks the PCK rotation models
// against the JPL mean-element orbits (and fixes the texture longitude origin).
void test_tidal_locking()
{
    Scene scene = load_jupiter();
    double t0 = 0.0;
    parse_utc("2026-10-02 00:00", &t0);
    double worst_deg = 0.0;
    for (int day = 0; day < 3650; day += 7) {
        scene.update(t0 + day * kSecondsPerDay);
        for (size_t i = 1; i < scene.bodies.size(); ++i) {
            const Body& moon = scene.bodies[i];
            if (moon.kind != BodyKind::Planet) {
                continue; // spacecraft
            }
            const glm::dvec3 to_jupiter = glm::normalize(scene.bodies[0].world_position - moon.world_position);
            const glm::dvec3 prime = moon.orientation[0];
            worst_deg = std::max(worst_deg, std::acos(std::clamp(glm::dot(prime, to_jupiter), -1.0, 1.0)) * kRadToDeg);
        }
    }
    check(worst_deg < 5.0, "moon prime meridians face Jupiter (deg)", worst_deg);
}

void test_solar_system_scene()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    const int v1 = scene.find("Voyager 1");
    const int v2 = scene.find("Voyager 2");
    const int jupiter = scene.find("Jupiter");
    const int neptune = scene.find("Neptune");
    check(v1 > 0 && v2 > 0 && jupiter > 0 && neptune > 0, "solar_system.toml bodies");

    // Not launched yet in 1977-07; visible (linearly extrapolated) in 2150.
    double t = 0.0;
    parse_utc("1977-07-01", &t);
    scene.update(t);
    check(!scene.bodies[static_cast<size_t>(v1)].visible, "Voyager 1 invisible before launch");
    parse_utc("2150-01-01", &t);
    scene.update(t);
    check(scene.bodies[static_cast<size_t>(v1)].visible, "Voyager 1 extrapolated after 2100");

    // Closest approaches agree with the bake tool (same data, independent code).
    auto min_distance = [&](int craft, int planet, double t_center) {
        double best = 1e300;
        for (double dt = -3600.0; dt <= 3600.0; dt += 10.0) {
            const glm::dvec3 a = scene.icrf_state_at(craft, t_center + dt).position;
            const glm::dvec3 b = scene.icrf_state_at(planet, t_center + dt).position;
            best = std::min(best, glm::length(a - b));
        }
        return best;
    };
    const double ca_jupiter = (jd_from_calendar({1979, 3, 5, 12, 5, 34}) - kJ2000Jd) * kSecondsPerDay;
    const double d_jupiter = min_distance(v1, jupiter, ca_jupiter);
    check(std::abs(d_jupiter - 348066.0) < 100.0, "V1 Jupiter closest approach (km)", d_jupiter);
    const double ca_neptune = (jd_from_calendar({1989, 8, 25, 4, 2, 51}) - kJ2000Jd) * kSecondsPerDay;
    const double d_neptune = min_distance(v2, neptune, ca_neptune);
    check(std::abs(d_neptune - 33661.0) < 100.0, "V2 Neptune closest approach (km)", d_neptune);

    // Planet ephemerides agree with the Standish fallback to its stated accuracy.
    parse_utc("2026-01-01", &t);
    scene.update(t);
    const auto* eph = dynamic_cast<const EphemerisMotion*>(scene.bodies[static_cast<size_t>(jupiter)].motion.get());
    check(eph != nullptr, "Jupiter uses an ephemeris");

    // History trail of Voyager 1 starts at its current position.
    std::vector<glm::dvec3> points;
    std::vector<float> fades;
    scene.trail(v1, t, 2048, points, fades);
    check(points.size() > 100 && glm::length(points[0] - scene.bodies[static_cast<size_t>(v1)].world_position) < 1.0,
          "Voyager 1 trail head", static_cast<double>(points.size()));
    // Orbit trail of Jupiter (osculating ellipse around the Sun) starts at Jupiter.
    scene.trail(jupiter, t, 360, points, fades);
    check(points.size() == 360 &&
              glm::length(points[0] - scene.bodies[static_cast<size_t>(jupiter)].world_position) < 1.0,
          "Jupiter osculating orbit head", points.empty() ? -1.0 : glm::length(points[0]));

    // Planet-centered frame: Jupiter at the origin.
    for (size_t f = 0; f < scene.frames.size(); ++f) {
        if (scene.frames[f].name == "Jupiter-centered") {
            scene.set_active_frame(static_cast<int>(f));
        }
    }
    scene.update(t);
    check(glm::length(scene.bodies[static_cast<size_t>(jupiter)].world_position) < 1e-6, "Jupiter-centered frame");
}

// Solar System scene: moons and dwarf planets against Horizons, the two pairs
// around their barycenters, and the shapes.
void test_solar_system_bodies()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    auto at = [&](const char* name, double t) {
        const int i = scene.find(name);
        check(i >= 0, name);
        return i >= 0 ? scene.icrf_state_at(i, t).position : glm::dvec3(0.0);
    };
    auto angle_deg = [](const glm::dvec3& a, const glm::dvec3& b) {
        return std::acos(std::clamp(glm::dot(glm::normalize(a), glm::normalize(b)), -1.0, 1.0)) / kDegToRad;
    };

    // Mean-element moons against Horizons, 2026-Jan-01 00:00 TDB (JD 2461041.5), ICRF km,
    // relative to the planet's center. Bounds: the worst case over 1972-2100.
    const double t = tdb_from_jd_tdb(2461041.5);
    struct Ref {
        const char* moon;
        const char* planet;
        glm::dvec3 hzn;
        double max_deg;
    };
    const Ref refs[] = {
        {"Phobos", "Mars", {-2628.223, -8355.790, -2928.529}, 0.4}, // 4.3 with [ELEM]'s rates
        {"Deimos", "Mars", {10863.689, 20169.954, 5021.957}, 0.5},
        {"Miranda", "Uranus", {113206.352, -45282.815, 44805.770}, 4.4},
        {"Ariel", "Uranus", {173540.311, -55002.747, 57468.461}, 1.0},
        {"Umbriel", "Uranus", {-52480.299, -58662.732, 253016.586}, 0.3},
        {"Titania", "Uranus", {362308.026, -14664.726, -241389.292}, 0.4},
        {"Oberon", "Uranus", {289669.879, -193050.060, 467214.375}, 0.4},
        {"Triton", "Neptune", {-284586.643, -111100.883, 180291.265}, 0.2},
    };
    for (const Ref& r : refs) {
        const double a = angle_deg(at(r.moon, t) - at(r.planet, t), r.hzn);
        check(a < r.max_deg, r.moon, a);
        std::printf("info: %s vs Horizons in 2026: %.2f deg\n", r.moon, a);
    }

    // Earth and Moon around their barycenter (Horizons @ 500@3), Pluto and Charon
    // around theirs (@ 500@9), from one scaled relative table each.
    const glm::dvec3 moon_hzn(142572.086, 286065.541, 158212.898);
    const glm::dvec3 earth_hzn(-1753.642, -3518.617, -1946.025);
    const glm::dvec3 pluto_hzn(-1105.447, -785.154, 1644.354);
    const glm::dvec3 charon_hzn(9057.326, 6432.789, -13474.591);
    const double e_moon = glm::length(at("Moon", t) - at("Earth-Moon barycenter", t) - moon_hzn);
    const double e_earth = glm::length(at("Earth", t) - at("Earth-Moon barycenter", t) - earth_hzn);
    const double e_pluto = glm::length(at("Pluto", t) - at("Pluto system", t) - pluto_hzn);
    const double e_charon = glm::length(at("Charon", t) - at("Pluto system", t) - charon_hzn);
    check(e_moon < 25.0, "Moon around the Earth-Moon barycenter (km)", e_moon);
    check(e_earth < 1.0, "Earth around the Earth-Moon barycenter (km)", e_earth);
    check(e_pluto < 2.0, "Pluto around its barycenter (km)", e_pluto);
    check(e_charon < 12.0, "Charon around the Pluto barycenter (km)", e_charon);

    // The same geocentric Moon as the Earth-Moon scene's (separately baked) table.
    Scene earth_moon = load_scene_or_die("earth_moon.toml");
    const double t_a2 = tdb_from_jd_tdb(2461136.5); // 2026-04-06, during Artemis II
    const glm::dvec3 geo = earth_moon.icrf_state_at(earth_moon.find("Moon"), t_a2).position -
                           earth_moon.icrf_state_at(0, t_a2).position;
    const double e_geo = glm::length(at("Moon", t_a2) - at("Earth", t_a2) - geo);
    check(e_geo < 25.0, "geocentric Moon agrees between scenes (km)", e_geo);

    // Fallbacks at the end of the baked tables (faded in over fallback_blend_days):
    // the Moon's mean elements are off by degrees, Charon's by little.
    const double t_end = tdb_from_jd_tdb(2488067.5); // 2099-12-30, within the tables
    struct Join {
        const char* body;
        double max_deg;
    };
    for (const Join& j : {Join{"Charon", 0.1}, Join{"Moon", 12.0}}) {
        const int i = scene.find(j.body);
        const int parent = scene.bodies[static_cast<size_t>(i)].parent;
        const auto* scaled = dynamic_cast<const ScaledMotion*>(scene.bodies[static_cast<size_t>(i)].motion.get());
        const MotionSource& m = scaled ? scaled->inner() : *scene.bodies[static_cast<size_t>(i)].motion;
        const auto* eph = dynamic_cast<const EphemerisMotion*>(&m);
        check(eph != nullptr && eph->fallback() != nullptr && parent >= 0, "baked table with a fallback");
        if (!eph || !eph->fallback()) {
            continue;
        }
        const double a = angle_deg(eph->table().eval(t_end).position, eph->fallback()->eval(t_end).position);
        check(a < j.max_deg, j.body, a);
        std::printf("info: %s fallback at the end of its table: %.3f deg\n", j.body, a);
    }

    // The Pluto system and the small bodies continue on the two-body orbit from the
    // nearer end of their tables: continuous there in position and velocity. How
    // fast it drifts from the real (perturbed) orbit shows by running it back over
    // the table, from the 2100 end to 2026.
    const double t_2026 = tdb_from_jd_tdb(2461200.5);
    for (const Join& j : {Join{"Pluto system", 0.2}, Join{"Ceres", 2.5}, Join{"Pallas", 6.0}, Join{"Vesta", 1.0},
                          Join{"Eris", 0.2}, Join{"Haumea", 0.4}, Join{"Makemake", 0.3}, Join{"Gonggong", 0.2},
                          Join{"Quaoar", 0.1}, Join{"Sedna", 0.2}, Join{"Arrokoth", 0.1}}) {
        const MotionSource* m = scene.bodies[static_cast<size_t>(scene.find(j.body))].motion.get();
        const auto* eph = dynamic_cast<const EphemerisMotion*>(m);
        if (eph && eph->fallback()) {
            eph = dynamic_cast<const EphemerisMotion*>(eph->fallback()); // Arrokoth: the NH table first
        }
        check(eph != nullptr && eph->fallback() == nullptr && eph->table().reference_gm() > 0.0,
              (std::string(j.body) + ": Kepler-relative table, no fallback").c_str());
        if (!eph) {
            continue;
        }
        const EphemerisTable& table = eph->table();
        double worst_jump = 0.0;
        for (const bool at_end : {false, true}) {
            const EphemerisTable::Knot& k = at_end ? table.knots().back() : table.knots().front();
            const double outward = at_end ? 1.0 : -1.0; // 1 s outside the table
            const State outside = eph->eval(k.t + outward);
            worst_jump = std::max({worst_jump, glm::length(outside.position - (k.position + outward * k.velocity)),
                                   1e2 * glm::length(outside.velocity - k.velocity)}); // km, ~10 m/s
        }
        check(worst_jump < 1e-3 && eph->valid_at(table.start() - 1e9) && eph->valid_at(table.end() + 1e9),
              (std::string(j.body) + ": extrapolation continuous at both ends and valid beyond").c_str(), worst_jump);
        const EphemerisTable::Knot& last = table.knots().back();
        const State back = propagate_kepler({last.position, last.velocity}, table.reference_gm(), t_2026 - last.t);
        const double drift = angle_deg(back.position, table.eval(t_2026).position);
        check(drift < j.max_deg, (std::string(j.body) + ": two-body drift over 74 years (deg)").c_str(), drift);
        std::printf("info: %s two-body orbit from the 2100 end, run back to 2026: %.3f deg\n", j.body, drift);
    }

    // Shapes: triaxial Haumea, Quaoar's oblate spheroid.
    const Body& haumea = scene.bodies[static_cast<size_t>(scene.find("Haumea"))];
    check(haumea.equatorial_radius_km == 1161.0 && haumea.equatorial_radius_b_km == 852.0 &&
              haumea.polar_radius_km == 513.0,
          "Haumea semi-axes");
    const Body& quaoar = scene.bodies[static_cast<size_t>(scene.find("Quaoar"))];
    check(quaoar.equatorial_radius_b_km == 566.1 && quaoar.polar_radius_km == 511.2, "Quaoar spheroid");

    // Rings from bands: each profile sample averages over its cell, so the narrow
    // Uranian rings (1.5-2 km, samples ~1.2 km) keep their equivalent width
    // (sum of tau x width) and none falls between samples.
    for (const char* name : {"Uranus", "Neptune", "Haumea", "Quaoar"}) {
        const RingSystem& rings = scene.bodies[static_cast<size_t>(scene.find(name))].rings;
        double expected = 0.0;
        for (const RingBand& b : rings.bands) {
            expected += b.optical_depth * (b.outer_km - b.inner_km);
        }
        const double cell = (rings.outer_km - rings.inner_km) / static_cast<double>(rings.profile.size() - 1);
        double integral = 0.0;
        for (const glm::vec2& p : rings.profile) {
            integral += p.x * cell;
        }
        check(!rings.bands.empty() && std::abs(integral / expected - 1.0) < 1e-3, name, integral / expected);
    }
    const RingSystem& uranus_rings = scene.bodies[static_cast<size_t>(scene.find("Uranus"))].rings;
    bool all_present = uranus_rings.bands.size() == 10;
    for (const RingBand& b : uranus_rings.bands) {
        all_present = all_present && uranus_rings.optical_depth(0.5 * (b.inner_km + b.outer_km)) > 0.05;
    }
    check(all_present, "all ten Uranian rings in the profile");

    // Display hierarchy: the Moon and Charon are satellites of their primaries.
    check(scene.satellite_host(scene.find("Moon")) == scene.find("Earth") &&
              scene.satellite_host(scene.find("Earth")) == 0 &&
              scene.satellite_host(scene.find("Charon")) == scene.find("Pluto"),
          "satellite hosts in the Solar System scene");
}

// Solar System scene: the asteroid and Kuiper belts (two-body orbits from SBDB).
void test_solar_system_belts()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    check(scene.belts.size() == 2 && scene.belts[0].size() > 60000 && scene.belts[1].size() > 7000,
          "main belt and Kuiper belt loaded", scene.belts.empty() ? 0.0 : static_cast<double>(scene.belts[0].size()));
    if (scene.belts.size() != 2) {
        return;
    }
    auto angle_deg = [](const glm::dvec3& a, const glm::dvec3& b) {
        return std::acos(std::clamp(glm::dot(glm::normalize(a), glm::normalize(b)), -1.0, 1.0)) / kDegToRad;
    };

    // Sorted by H: the brightest are Vesta and Eris. Against the scene's own bodies
    // (Horizons): exact at the epoch, then drifting as a two-body orbit does.
    struct Probe {
        const SceneBelt* belt;
        const char* body;
        double max_at_epoch;
        double max_2000;
    };
    for (const Probe& p : {Probe{&scene.belts[0], "Vesta", 0.05, 3.0}, Probe{&scene.belts[1], "Eris", 0.03, 0.2}}) {
        const int body = scene.find(p.body);
        for (const double t : {p.belt->epoch_tdb, 0.0}) {
            const double a = angle_deg(belt_position(*p.belt, 0, t), scene.icrf_state_at(body, t).position);
            check(a < (t == 0.0 ? p.max_2000 : p.max_at_epoch), p.body, a);
            std::printf("info: %s from the belt elements vs Horizons, %s: %.3f deg\n", p.body,
                        t == 0.0 ? "2000" : "epoch", a);
        }
    }

    // The 3:1 Kirkwood gap at 2.50 au (resonance with Jupiter) is well below its
    // neighbourhood.
    auto count = [&](double lo, double hi) {
        size_t n = 0;
        for (size_t k = 0; k < scene.belts[0].size(); ++k) {
            const double a = scene.belts[0].elements[k * SceneBelt::kStride];
            n += (a >= lo && a < hi) ? 1 : 0;
        }
        return static_cast<double>(n) / (hi - lo);
    };
    const double gap = count(2.49, 2.51);
    const double around = 0.5 * (count(2.40, 2.47) + count(2.53, 2.60));
    check(gap < 0.3 * around, "Kirkwood gap at the 3:1 resonance", gap / around);
}

// Uranus' eccentric, inclined rings [FR24 = French et al. 2024, arXiv:2401.04634].
void test_uranus_rings()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    const int uranus = scene.find("Uranus");
    const Body& body = scene.bodies[static_cast<size_t>(uranus)];
    const RingSystem& rings = body.rings;
    check(rings.shaped() && std::abs(rings.epoch_tdb - tdb_from_jd_tdb(2446450.0)) < 1.0, "Uranus rings shaped");
    auto band = [&](const char* name) -> const RingBand* {
        for (const RingBand& b : rings.bands) {
            if (b.name == name) {
                return &b;
            }
        }
        return nullptr;
    };

    // The apse and node rates follow from the gravity field fitted to them,
    // [FR24] Table 17 (GM, R, adopted J2, J4, J6): the satellites' share and
    // higher-order terms are below 0.1%. Catches mistyped rates.
    constexpr double kGm = 5793950.3;
    constexpr double kR = 25559.0;
    constexpr double kJ2 = 3509.291e-6;
    constexpr double kJ4 = -35.522e-6;
    constexpr double kJ6 = 0.5e-6;
    double worst = 0.0;
    for (const RingBand& b : rings.bands) {
        if (b.peri_rate_deg_per_day == 0.0) {
            continue;
        }
        const double a = 0.5 * (b.inner_km + b.outer_km);
        const double x = (kR / a) * (kR / a);
        const double base = kGm / (a * a * a);
        const double n = std::sqrt(base * (1.0 + 1.5 * kJ2 * x - 15.0 / 8.0 * kJ4 * x * x + 35.0 / 16.0 * kJ6 * x * x * x));
        const double kappa = std::sqrt(base * (1.0 - 1.5 * kJ2 * x + 45.0 / 8.0 * kJ4 * x * x - 175.0 / 16.0 * kJ6 * x * x * x));
        const double nu = std::sqrt(base * (1.0 + 4.5 * kJ2 * x - 75.0 / 8.0 * kJ4 * x * x + 245.0 / 16.0 * kJ6 * x * x * x));
        const double to_deg_per_day = kSecondsPerDay / kDegToRad;
        worst = std::max(worst, std::abs(b.peri_rate_deg_per_day / ((n - kappa) * to_deg_per_day) - 1.0));
        if (b.a_sin_i_km > 0.0) {
            worst = std::max(worst, std::abs(b.node_rate_deg_per_day / ((n - nu) * to_deg_per_day) - 1.0));
        }
    }
    check(worst < 1e-3, "Uranian ring precession rates vs J2, J4, J6", worst);

    // Widths at periapsis and apoapsis against the ranges of French et al. 1991
    // (PDS Ring-Moon Systems Node), which are rounded to whole km.
    struct Width {
        const char* ring;
        double min_km;
        double max_km;
    };
    for (const Width& w : {Width{"Alpha", 4.0, 10.0}, Width{"Beta", 5.0, 11.0}, Width{"Epsilon", 20.0, 96.0}}) {
        const RingBand* b = band(w.ring);
        if (!b) {
            check(false, w.ring);
            continue;
        }
        const double at_peri = ring_edge_radius(b->outer_km, b->outer_ae_km, 0.0) - ring_edge_radius(b->inner_km, b->inner_ae_km, 0.0);
        const double at_apo = ring_edge_radius(b->outer_km, b->outer_ae_km, kPi) - ring_edge_radius(b->inner_km, b->inner_ae_km, kPi);
        check(std::abs(at_peri - w.min_km) < 1.5 && std::abs(at_apo - w.max_km) < 1.5, w.ring, at_apo);
    }
    // tau x width is kept: the epsilon ring's optical depth at both ends, 0.5-2.3 there.
    if (const RingBand* b = band("Epsilon")) {
        const double equivalent = b->optical_depth * (b->outer_km - b->inner_km);
        const double at_peri = equivalent / (ring_edge_radius(b->outer_km, b->outer_ae_km, 0.0) -
                                             ring_edge_radius(b->inner_km, b->inner_ae_km, 0.0));
        const double at_apo = equivalent / (ring_edge_radius(b->outer_km, b->outer_ae_km, kPi) -
                                            ring_edge_radius(b->inner_km, b->inner_ae_km, kPi));
        check(std::abs(at_peri - 2.3) < 0.15 && std::abs(at_apo - 0.5) < 0.05, "epsilon ring optical depth range", at_peri);
    }

    // Longitudes: from the ascending node of the equator on the ICRF equator, about
    // the pole of positive angular momentum, here taken from [FR24] Table 13 (the
    // scene uses the IAU pole, 0.002 deg away). At the epoch, the epsilon ring's
    // periapsis lies at 307.076 deg; ring 6 is highest 90 deg past its node.
    scene.update(rings.epoch_tdb);
    const glm::dmat3 body_to_icrf = scene.current_transform().axes * body.orientation;
    const glm::dvec3 pole = unit_toward(77.311327, 15.172795);
    const glm::dvec3 x = glm::normalize(glm::cross(glm::dvec3(0.0, 0.0, 1.0), pole));
    const glm::dvec3 y = glm::cross(pole, x);
    auto toward = [&](double longitude_deg) {
        return std::cos(longitude_deg * kDegToRad) * x + std::sin(longitude_deg * kDegToRad) * y;
    };
    const RingBand* epsilon = band("Epsilon");
    const RingBand* six = band("6");
    if (epsilon && six) {
        const RingBandShape shape = ring_band_shape(body, *epsilon, body_to_icrf, rings.epoch_tdb);
        const glm::dvec3 periapsis = body_to_icrf * glm::dvec3(shape.periapsis, 0.0);
        const double off = std::acos(std::clamp(glm::dot(periapsis, toward(307.076)), -1.0, 1.0)) / kDegToRad;
        check(off < 0.01, "epsilon ring periapsis direction", off);

        const RingBandShape tilt = ring_band_shape(body, *six, body_to_icrf, rings.epoch_tdb);
        const glm::dvec3 u = glm::transpose(body_to_icrf) * toward(89.727 + 90.0);
        const double height = (tilt.node.x * u.y - tilt.node.y * u.x) * glm::dot(body_to_icrf[2], pole);
        check(std::abs(height - 44.643) < 0.01, "ring 6 height 90 deg past its node", height);

        // A year later both have precessed at their rates.
        const double t = rings.epoch_tdb + 365.0 * kSecondsPerDay;
        scene.update(t);
        const glm::dmat3 later = scene.current_transform().axes * body.orientation;
        const glm::dvec3 moved = later * glm::dvec3(ring_band_shape(body, *epsilon, later, t).periapsis, 0.0);
        const double drift =
            std::acos(std::clamp(glm::dot(moved, toward(307.076 + 365.0 * 1.3632575)), -1.0, 1.0)) / kDegToRad;
        check(drift < 0.01, "epsilon ring apsidal precession", drift);
    }
}

// Ring arcs: Neptune's Adams arcs and Quaoar's dense Q1R arc.
void test_ring_arcs()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    auto find_band = [&](const Body& body, const char* name) -> const RingBand* {
        for (const RingBand& b : body.rings.bands) {
            if (b.name == name) {
                return &b;
            }
        }
        check(false, name);
        return nullptr;
    };

    // Carried from Voyager's 1989 position at [SOU22]'s mean motion, Fraternite's
    // middle must be where SPHERE saw it 27 years later: 217.43 +- 0.30 deg from
    // the ascending node of Neptune's equator on the ICRF equator, at the time the
    // light left Neptune, 28.959 au (Table 1) before the 05:33:51.958 UT frame.
    const Body& neptune = scene.bodies[static_cast<size_t>(scene.find("Neptune"))];
    if (const RingBand* fraternite = find_band(neptune, "Fraternite")) {
        double t = 0.0;
        parse_utc("2016-08-23 05:33:51.958", &t);
        t -= 28.959 * kAuKm / kSpeedOfLightKmS;
        scene.update(t);
        const glm::dmat3 body_to_icrf = scene.current_transform().axes * neptune.orientation;
        const double center = ring_band_shape(neptune, *fraternite, body_to_icrf, t).arc_center;
        const glm::dvec3 dir = body_to_icrf * glm::dvec3(std::cos(center), std::sin(center), 0.0);
        const glm::dvec3 pole = body_to_icrf[2];
        const glm::dvec3 x = glm::normalize(glm::cross(glm::dvec3(0.0, 0.0, 1.0), pole));
        const double longitude =
            wrap_two_pi(std::atan2(glm::dot(dir, glm::cross(pole, x)), glm::dot(dir, x))) / kDegToRad;
        // 0.3 deg each from that measurement and from the rate over 27 years; and
        // [SOU22] counted the interval in UT (13 leap seconds, 0.12 deg, short).
        check(std::abs(longitude - 217.43) < 0.6, "Fraternite in 2016 (SPHERE)", longitude);
        std::printf("info: Fraternite on 2016-08-23 at %.2f deg (SPHERE: 217.43 +- 0.30)\n", longitude);
    }
    if (const RingBand* liberte = find_band(neptune, "Liberte")) {
        const double y1989 = tdb_from_jd_tdb(2447757.0);
        const double y2009 = tdb_from_jd_tdb(2454832.5);
        check(liberte->fade(y1989 - kSecondsPerDay) == 1.0 &&
                  std::abs(liberte->fade(0.5 * (y1989 + y2009)) - 0.5) < 0.01 && liberte->fade(y2009 + kSecondsPerDay) == 0.0,
              "Liberte fades out by 2009");
    }

    // Quaoar's arc moves at the circular mean motion at its radius.
    const Body& quaoar = scene.bodies[static_cast<size_t>(scene.find("Quaoar"))];
    if (const RingBand* arc = find_band(quaoar, "Q1R dense arc")) {
        const double r = 0.5 * (arc->inner_km + arc->outer_km);
        const double n = std::sqrt(quaoar.body_gm_km3_s2 / (r * r * r)) * kSecondsPerDay / kDegToRad;
        check(std::abs(arc->arc_rate_deg_per_day / n - 1.0) < 1e-5 && quaoar.rings.shaped(), "Quaoar arc mean motion",
              arc->arc_rate_deg_per_day - n);
    }
}

void test_saturn_scene()
{
    Scene scene = load_scene_or_die("saturn.toml");
    auto tdb = [](const CalendarDateTime& c) { return (jd_from_calendar(c) - kJ2000Jd) * kSecondsPerDay; };
    auto at = [&](int body, double t) { return scene.icrf_state_at(body, t).position; };

    // Mean elements (phases fitted by tools/bake/fit_moon_phase.py) against the
    // Horizons tables wherever these have data. Mimas librates by +-44 deg; Hyperion,
    // in the 4:3 resonance with Titan, strays by up to 15 deg (rms 7).
    const char* moons[] = {"Mimas", "Enceladus", "Tethys", "Dione", "Rhea", "Titan", "Iapetus", "Hyperion"};
    for (const char* name : moons) {
        const int b = scene.find(name);
        const auto* eph = b > 0 ? dynamic_cast<const EphemerisMotion*>(scene.bodies[static_cast<size_t>(b)].motion.get()) : nullptr;
        if (!eph || !eph->fallback()) {
            check(false, "Saturn moon uses an ephemeris with mean-element fallback");
            continue;
        }
        double worst = 0.0;
        for (double t = eph->table().start(); t < eph->table().end(); t += 0.25 * kSecondsPerDay) {
            if (!eph->table().covers(t)) {
                continue;
            }
            const glm::dvec3 a = eph->table().eval(t).position;
            const glm::dvec3 m = eph->fallback()->eval(t).position;
            worst = std::max(worst, std::acos(std::clamp(glm::dot(a, m) / (glm::length(a) * glm::length(m)), -1.0, 1.0)) * kRadToDeg);
        }
        const std::string moon = name;
        const double limit = moon == "Mimas" ? 45.0 : moon == "Hyperion" ? 16.0 : 1.5;
        check(worst < limit, (moon + " mean elements near Horizons (deg)").c_str(), worst);
    }

    // Phoebe orbits backward: its orbital angular momentum points away from Saturn's pole.
    {
        const int phoebe = scene.find("Phoebe");
        const double t = tdb({2004, 6, 11, 0, 0, 0});
        const State s = scene.icrf_state_at(phoebe, t);
        const glm::dvec3 h = glm::normalize(glm::cross(s.position, s.velocity));
        const glm::dvec3 pole = unit_from_ra_dec(40.589 * kDegToRad, 83.537 * kDegToRad); // [PCK]
        const double tilt = std::acos(glm::dot(h, pole)) * kRadToDeg;
        check(tilt > 150.0, "Phoebe retrograde (deg from Saturn's pole)", tilt);
        // Outside the 1972-2100 table it keeps to a two-body orbit of the same size
        // (a = 12.93 million km, e = 0.164 [ELEM]).
        for (const double year : {1900.0, 2200.0}) {
            const double r = glm::length(at(phoebe, (year - 2000.0) * kDaysPerJulianYear * kSecondsPerDay));
            check(r > 10.0e6 && r < 16.0e6, "Phoebe beyond its table stays at its distance (km)", r);
        }
    }

    // Windowed table (Enceladus near Cassini only): gaps fall back to the mean
    // elements, and the position stays continuous across a window's edge.
    const int enceladus = scene.find("Enceladus");
    const auto* enc = dynamic_cast<const EphemerisMotion*>(scene.bodies[static_cast<size_t>(enceladus)].motion.get());
    check(enc && enc->table().max_gap() > 0.0 && !enc->table().covers(tdb({2006, 1, 1, 0, 0, 0})),
          "Enceladus table has gaps");
    if (enc) {
        double seg_start = 0.0;
        double seg_end = 0.0;
        enc->table().segment(tdb({2008, 10, 9, 19, 7, 44}), &seg_start, &seg_end);
        check(seg_end - seg_start > 2.0 * kSecondsPerDay && seg_end - seg_start < 30.0 * kSecondsPerDay,
              "Enceladus window around the E5 flyby (days)", (seg_end - seg_start) / kSecondsPerDay);
        double worst_jump = 0.0;
        for (double t = seg_start - kSecondsPerDay; t < seg_start + 2.0 * kSecondsPerDay; t += 60.0) {
            const State a = scene.icrf_state_at(enceladus, t);
            worst_jump = std::max(worst_jump, glm::length(at(enceladus, t + 60.0) - a.position - a.velocity * 60.0));
        }
        check(worst_jump < 50.0, "Enceladus enters its table window without a jump (km/min)", worst_jump);
    }

    // The ring profile (Cassini RSS, Rev 7): landmarks of the main rings.
    const RingSystem& rings = scene.bodies[0].rings;
    check(rings.profile.size() == 8192 && std::abs(rings.inner_km - 72771.0) < 5.0 && std::abs(rings.outer_km - 144989.0) < 5.0,
          "Saturn ring profile range", rings.inner_km);
    check(rings.optical_depth(110000.0) > 3.0 && rings.optical_depth(119000.0) < 0.3 && rings.optical_depth(80000.0) < 0.3 &&
              rings.optical_depth(130000.0) > 0.4 && rings.optical_depth(133580.0) < 0.05,
          "Saturn rings: opaque B ring, thin C ring and Cassini Division, A ring, Encke gap", rings.optical_depth(110000.0));

    // Synchronous moons: the PCK prime meridian faces Saturn (Horizons positions,
    // over one Iapetus orbit; Mimas within its 2010-02 window). The IAU meridians
    // are tied to surface features, which leaves offsets of a few degrees (Rhea
    // ~2.9); eccentricity swings the direction by up to 2e (Titan, Iapetus ~3.3,
    // Mimas 2.3 deg), and Iapetus' orbit is inclined ~7.6 deg to its Laplace plane.
    const char* locked[] = {"Rhea", "Titan", "Iapetus", "Mimas"};
    const double lock_limit[] = {4.0, 7.0, 10.0, 6.0};
    for (int m = 0; m < 4; ++m) {
        const bool mimas = m == 3;
        const double t_begin = mimas ? tdb({2010, 2, 13, 0, 0, 0}) : tdb({2008, 1, 1, 0, 0, 0});
        const double t_stop = mimas ? tdb({2010, 2, 14, 12, 0, 0}) : tdb({2008, 3, 21, 0, 0, 0});
        double worst_lock = 0.0;
        for (double t = t_begin; t < t_stop; t += (mimas ? 0.05 : 0.5) * kSecondsPerDay) {
            scene.update(t);
            const Body& moon = scene.bodies[static_cast<size_t>(scene.find(locked[m]))];
            const glm::dvec3 to_saturn = glm::normalize(scene.bodies[0].world_position - moon.world_position);
            worst_lock = std::max(worst_lock, std::acos(std::clamp(glm::dot(moon.orientation[0], to_saturn), -1.0, 1.0)) * kRadToDeg);
        }
        std::printf("info: %s prime meridian vs Saturn direction: worst %.2f deg\n", locked[m], worst_lock);
        check(worst_lock < lock_limit[m], "Saturn moon's prime meridian faces Saturn (deg)", worst_lock);
    }
}

// Janus and Epimetheus: the element tables against Horizons, and the swaps they hold.
void test_saturn_coorbitals()
{
    Scene scene = load_scene_or_die("saturn.toml");
    auto tdb = [](const CalendarDateTime& c) { return (jd_from_calendar(c) - kJ2000Jd) * kSecondsPerDay; };
    const int janus = scene.find("Janus");
    const int epimetheus = scene.find("Epimetheus");
    check(janus > 0 && epimetheus > 0, "Saturn scene has Janus and Epimetheus");
    if (janus <= 0 || epimetheus <= 0) {
        return;
    }
    const Body& jb = scene.bodies[static_cast<size_t>(janus)];
    const Body& eb = scene.bodies[static_cast<size_t>(epimetheus)];
    const auto* jt = dynamic_cast<const ElementTableOrbit*>(jb.motion.get());
    const auto* et = dynamic_cast<const ElementTableOrbit*>(eb.motion.get());
    check(jt && et, "Janus and Epimetheus use element tables");
    if (!jt || !et) {
        return;
    }

    // JPL Horizons (sat415_merged_DE437), geometric ICRF positions relative to Saturn's
    // center at 00:00 TDB; the tables are thinned to 10 km of the hourly positions.
    struct Sample {
        CalendarDateTime when;
        glm::dvec3 janus, epimetheus;
    };
    const Sample samples[] = {
        {{1981, 8, 25, 0, 0, 0}, {146309.3005, 39949.2094, -15146.0185}, {141984.8432, -47511.3035, -8922.8207}},
        {{2006, 1, 1, 0, 0, 0}, {88185.4289, 121435.7220, -16942.6047}, {101919.2617, 112568.3887, -17943.7833}},
        {{2047, 10, 17, 0, 0, 0}, {68940.2097, -135792.5591, 4497.5133}, {-109497.5635, 105791.0986, 818.7919}},
    };
    for (const Sample& s : samples) {
        const double t = tdb(s.when);
        const double dj = glm::length(scene.icrf_state_at(janus, t).position - s.janus);
        const double de = glm::length(scene.icrf_state_at(epimetheus, t).position - s.epimetheus);
        check(dj < 10.5 && de < 10.5, "Janus and Epimetheus tables match Horizons (km)", std::max(dj, de));
    }
    // The velocity (precessing ellipse) agrees with the change in position.
    {
        const double t = tdb({2006, 1, 1, 0, 0, 0});
        const State s = jt->eval(t);
        const glm::dvec3 drift = (jt->eval(t + 30.0).position - jt->eval(t - 30.0).position) / 60.0;
        check(glm::length(s.velocity - drift) < 1e-3, "Janus velocity matches its path (km/s)", glm::length(s.velocity - drift));
    }

    // The swaps: the semi-major axes cross every four years, the latest of Cassini's
    // first year on 2006-01-21 [DNC]. Each time the pair's orbital momentum is
    // conserved, GM_J da_J + GM_E da_E = 0: the a jumps (means over 100-200 days
    // before and after) are in the inverse ratio of the [SPHY] masses (0.2775).
    const double mass_ratio = eb.body_gm_km3_s2 / jb.body_gm_km3_s2;
    const ElementTable& ja = jt->table();
    const ElementTable& ea = et->table();
    auto a_diff = [&](double t) { return ja.at(t).a - ea.at(t).a; };
    int swaps = 0;
    double worst_ratio = 0.0;
    double t_2006 = 0.0;
    for (double t = ja.start(); t + kSecondsPerDay < ja.end(); t += kSecondsPerDay) {
        const double d0 = a_diff(t);
        const double d1 = a_diff(t + kSecondsPerDay);
        if ((d0 < 0.0) == (d1 < 0.0)) {
            continue;
        }
        const double tc = t + kSecondsPerDay * d0 / (d0 - d1);
        if (tc < ja.start() + 200.0 * kSecondsPerDay || tc > ja.end() - 200.0 * kSecondsPerDay) {
            continue;
        }
        ++swaps;
        auto mean_a = [&](const ElementTable& table, double from_days, double to_days) {
            double sum = 0.0;
            int n = 0;
            for (double d = from_days; d < to_days; d += 1.0, ++n) {
                sum += table.at(tc + d * kSecondsPerDay).a;
            }
            return sum / n;
        };
        const double da_j = mean_a(ja, 100.0, 200.0) - mean_a(ja, -200.0, -100.0);
        const double da_e = mean_a(ea, 100.0, 200.0) - mean_a(ea, -200.0, -100.0);
        worst_ratio = std::max(worst_ratio, std::abs(-da_j / da_e / mass_ratio - 1.0));
        check(std::abs(da_e) > 70.0 && std::abs(da_e) < 82.0 && (da_j < 0.0) == (da_e > 0.0),
              "Epimetheus' orbit jumps ~76 km at a swap, Janus' the other way (km)", da_e);
        if (std::abs(tc - tdb({2006, 1, 21, 0, 0, 0})) < 200.0 * kSecondsPerDay) {
            t_2006 = tc;
        }
        t += 300.0 * kSecondsPerDay;
    }
    std::printf("info: %d Janus-Epimetheus swaps, a jumps within %.2f%% of the mass ratio\n", swaps, 100.0 * worst_ratio);
    check(swaps == 19, "Janus and Epimetheus swap every four years 1974-2046", swaps);
    check(worst_ratio < 0.015, "swap a jumps in the inverse mass ratio (relative error)", worst_ratio);
    const double off_days = (t_2006 - tdb({2006, 1, 21, 12, 0, 0})) / kSecondsPerDay;
    check(std::abs(off_days) < 1.0, "2006 swap on January 21 (days)", off_days);

    // Never closer than ~10,000 km: Horizons' hourly samples give 10,235 km on 2006-01-20.
    double closest = 1e9;
    for (double t = tdb({2005, 12, 20, 0, 0, 0}); t < tdb({2006, 2, 20, 0, 0, 0}); t += 600.0) {
        closest = std::min(closest, glm::length(scene.icrf_state_at(janus, t).position - scene.icrf_state_at(epimetheus, t).position));
    }
    check(closest > 10100.0 && closest < 10400.0, "closest approach of the 2006 swap (km)", closest);

    // Outside 1972-2050 the elements run on: still on the shared orbit.
    for (const double year : {1900.0, 2100.0}) {
        const double r = glm::length(scene.icrf_state_at(janus, (year - 2000.0) * kDaysPerJulianYear * kSecondsPerDay).position);
        check(r > 150000.0 && r < 153000.0, "Janus beyond its table stays on its orbit (km)", r);
    }

    // Synchronous: the long axis (the shape models' x) toward Saturn.
    scene.update(tdb({2006, 1, 1, 0, 0, 0}));
    for (const Body* b : {&jb, &eb}) {
        const glm::dvec3 to_saturn = glm::normalize(scene.bodies[0].world_position - b->world_position);
        const double off = std::acos(std::clamp(glm::dot(b->orientation[0], to_saturn), -1.0, 1.0)) * kRadToDeg;
        check(off < 0.01, "Janus/Epimetheus long axis toward Saturn (deg)", off);
    }
}

// Earth-Moon scene: the Moon-centered rotating frame has the Moon at the origin
// and the Earth on -x.
void test_earth_moon_frame()
{
    Scene scene = load_scene_or_die("earth_moon.toml");
    const int moon = scene.find("Moon");
    for (size_t f = 0; f < scene.frames.size(); ++f) {
        if (scene.frames[f].name == "Earth-Moon rotating, Moon-centered") {
            scene.set_active_frame(static_cast<int>(f));
        }
    }
    scene.update(tdb_from_jd_tdb(jd_from_calendar({2023, 1, 1, 0, 0, 0})));
    const glm::dvec3 earth = scene.bodies[0].world_position;
    check(glm::length(scene.bodies[static_cast<size_t>(moon)].world_position) < 1e-6 && earth.x < -3.5e5 &&
              std::abs(earth.y) < 1e-3 && std::abs(earth.z) < 1e-3,
          "Moon-centered rotating frame", earth.x);
}

// Earth's prime meridian faces the Sun at 12:00 UTC (to within the equation of
// time, < 4.5 deg); the Moon's near side faces the Earth (libration < 10 deg).
void test_earth_moon_rotation()
{
    Scene scene = load_scene_or_die("space_telescopes.toml");
    scene.set_active_frame(0); // inertial
    const int moon = scene.find("Moon");
    for (const char* date : {"2023-02-11 12:00", "2026-06-21 12:00", "2029-11-03 12:00"}) {
        double t = 0.0;
        parse_utc(date, &t);
        scene.update(t);
        const glm::dvec3 pole = scene.bodies[0].orientation[2];
        const glm::dvec3 meridian = scene.bodies[0].orientation[0];
        glm::dvec3 sun = glm::normalize(scene.sun_position());
        sun = glm::normalize(sun - pole * glm::dot(sun, pole)); // project on the equator
        const double angle = std::acos(std::clamp(glm::dot(meridian, sun), -1.0, 1.0)) * kRadToDeg;
        check(angle < 4.5, "Greenwich faces the Sun at noon UTC (deg)", angle);

        const glm::dvec3 near_side = scene.bodies[static_cast<size_t>(moon)].orientation[0];
        const glm::dvec3 to_earth = glm::normalize(-scene.bodies[static_cast<size_t>(moon)].world_position);
        const double lib = std::acos(std::clamp(glm::dot(near_side, to_earth), -1.0, 1.0)) * kRadToDeg;
        check(lib < 10.0, "Moon near side faces Earth (deg)", lib);
    }
}

// Mars' axes in every scene that shows it: pck00011 with its long-period terms (the
// constant terms alone put the pole 1.5 deg off and the prime meridian 1-1.4 deg off).
// Checked against the IAU 2009 model (the "old values" in pck00011: pole 317.68143 -
// 0.1061 T, 52.88650 - 0.0609 T, W = 176.630 + 350.89198226 d), which has no such
// terms; the two agree to 0.003 deg.
void test_mars_rotation()
{
    for (const char* file : {"mars.toml", "solar_system.toml", "isee3.toml", "mercury.toml", "parker.toml",
                             "rosetta.toml"}) {
        Scene scene = load_scene_or_die(file);
        const int mars = scene.find("Mars");
        check(mars >= 0, "Mars in the scene");
        double worst = 0.0;
        for (const double jd : {2451545.0, 2459256.0, 2466000.0}) { // J2000, 2021-02-10, 2039-06
            const double d = jd - kJ2000Jd;
            const double centuries = d / kDaysPerJulianCentury;
            const glm::dmat3 iau2009 =
                iau_pole_frame((317.68143 - 0.1061 * centuries) * kDegToRad, (52.88650 - 0.0609 * centuries) * kDegToRad) *
                rotation_z(std::fmod(176.630 + 350.89198226 * d, 360.0) * kDegToRad);
            const glm::dmat3 axes = scene.body_to_icrf(scene.bodies[static_cast<size_t>(mars)], d * kSecondsPerDay);
            for (int c = 0; c < 3; ++c) {
                worst = std::max(worst, std::acos(std::min(1.0, glm::dot(axes[c], iau2009[c]))) * kRadToDeg);
            }
        }
        check(worst < 0.005, (std::string("Mars' axes match IAU 2009 (deg): ") + file).c_str(), worst);
    }
}

// Phobos and Deimos rotate synchronously in Cassini states: pck00011's periodic
// terms swing their spin poles around the Laplace poles with the orbits' nodes, so
// the poles stay close to the orbit normals (with the mean poles alone they are off
// by the inclinations, 1.1 and 1.8 deg). Phobos' long axis swings away from Mars by
// the optical libration (2e = 1.72 deg) plus its forced libration (1.143 deg), which
// is in antiphase: the orbital forcing is faster than its free libration. In phase it
// would be 0.6 deg, a quarter period off 2.1 deg; without the quadratic PM term the
// axis would lag by 12.7 deg in 2100. Deimos' orbit is nearly circular. Orbits from the
// baked tables and the mean elements alike (within 0.3 and 0.5 deg of Horizons), 1972-2100.
void test_mars_moon_rotation()
{
    for (const char* file : {"mars.toml", "solar_system.toml"}) {
        Scene scene = load_scene_or_die(file);
        const int mars = scene.find("Mars");
        struct Moon {
            const char* name;
            double max_pole_deg;
            double axis_deg; // expected worst angle of the long axis from Mars
            double axis_tolerance_deg;
        };
        for (const Moon& m : {Moon{"Phobos", 0.15, 1.72 + 1.143, 0.35}, Moon{"Deimos", 0.15, 0.0, 0.5}}) {
            const int moon = scene.find(m.name);
            check(moon >= 0 && mars >= 0, m.name);
            if (moon < 0 || mars < 0) {
                continue;
            }
            double worst_pole = 0.0;
            double worst_axis = 0.0;
            for (double jd = 2441317.5; jd < 2488069.5; jd += 9.7) { // 1972-01-01 .. 2100
                const double t = tdb_from_jd_tdb(jd);
                const State a = scene.icrf_state_at(mars, t);
                const State b = scene.icrf_state_at(moon, t);
                const glm::dvec3 r = b.position - a.position;
                const glm::dvec3 normal = glm::normalize(glm::cross(r, b.velocity - a.velocity));
                const glm::dmat3 axes = scene.body_to_icrf(scene.bodies[static_cast<size_t>(moon)], t);
                worst_pole = std::max(worst_pole, std::acos(std::min(1.0, glm::dot(axes[2], normal))) * kRadToDeg);
                const glm::dvec3 to_mars = glm::transpose(axes) * -r;
                worst_axis = std::max(worst_axis, std::abs(std::atan2(to_mars.y, to_mars.x)) * kRadToDeg);
            }
            const std::string where = std::string(m.name) + ", " + file;
            check(worst_pole < m.max_pole_deg, ("spin pole along the orbit normal (deg): " + where).c_str(), worst_pole);
            check(std::abs(worst_axis - m.axis_deg) < m.axis_tolerance_deg,
                  ("long axis toward Mars, librating (deg): " + where).c_str(), worst_axis);
            std::printf("info: %s: pole within %.3f deg of the orbit normal, long axis within %.2f deg of Mars\n",
                        where.c_str(), worst_pole, worst_axis);
        }
    }
}

} // namespace

// The Kerbol system against the running game (KSP 1.12.5, kOS dump, tools/kerbol/):
// each body's position relative to its parent and its rotation angle at the UT
// it was sampled, in the scene frame (x = KSP's solar prime vector, z = north),
// as printed by tools/kerbol/kerbol_scene.py.
void test_kerbol_scene()
{
    struct Sample {
        const char* body;
        double ut;
        double rotation_deg;
        glm::dvec3 position_km;
    };
    static const Sample kSamples[] = {
        {"Kerbol", 332.039999999988, 0.27669999999999, {0.0, 0.0, 0.0}},
        {"Moho", 332.319999999988, 190.098872066116, {-564669.323110, -6287327.876418, -198883.630460}},
        {"Eve", 332.339999999988, 1.48624099378876, {-9595693.991220, -2558836.137193, 436.442853}},
        {"Gilly", 332.899999999987, 9.24151477614565, {-25617.763405, -13663.859803, 4858.164644}},
        {"Kerbin", 332.159999999988, 95.5489925593853, {-13599827.571630, 18574.462196, 0.0}},
        {"Mun", 332.279999999988, 230.860678034095, {-1724.710820, 11875.410418, 0.0}},
        {"Minmus", 332.299999999988, 232.961089108911, {-45665.321098, 9977.135858, 4912.759825}},
        {"Duna", 332.459999999988, 91.8267629794643, {15516536.385679, -15288701.172025, 30.357090}},
        {"Ike", 332.579999999987, 1.82742226468935, {-702.093138, 3140.644127, 10.962961}},
        {"Dres", 332.939999999987, 28.4442068965516, {-45884866.640276, -8035236.234702, -4075491.109634}},
        {"Jool", 332.599999999987, 3.32599999999987, {34303640.676963, 55627044.143561, 164253.209534}},
        {"Laythe", 332.719999999987, 92.260800540243, {-27164.513797, -1029.099694, 0.0}},
        {"Vall", 332.839999999987, 1.13080443440434, {26151.405958, 34324.904522, 0.0}},
        {"Tylo", 332.879999999987, 0.565464348648838, {-68497.653779, -566.945094, -0.247376}},
        {"Bop", 332.859999999987, 230.220069724901, {-42041.149573, 103209.548204, 29190.902773}},
        {"Pol", 332.919999999987, 25.1328870732527, {11171.114416, 163784.243205, 12134.885248}},
        {"Eeloo", 332.959999999987, 31.1595889003081, {-72411834.780918, 86640570.684734, 11977933.186606}},
    };

    Scene scene = load_scene_or_die("kerbol.toml");
    scene.set_active_frame(0); // inertial
    check(scene.calendar == Calendar::Kerbin, "Kerbol scene uses the Kerbin calendar");
    double worst_km = 0.0;
    double worst_rotation_deg = 0.0;
    for (const Sample& s : kSamples) {
        const int b = scene.find(s.body);
        check(b >= 0, s.body);
        if (b < 0) {
            continue;
        }
        scene.update(s.ut);
        const Body& body = scene.bodies[static_cast<size_t>(b)];
        const glm::dvec3 parent =
            body.parent >= 0 ? scene.bodies[static_cast<size_t>(body.parent)].icrf_position : glm::dvec3(0.0);
        worst_km = std::max(worst_km, glm::length(body.icrf_position - parent - s.position_km));
        const glm::dvec3 meridian = body.orientation[0];
        const double angle = std::atan2(meridian.y, meridian.x) * kRadToDeg;
        const double off = std::abs(wrap_pi((angle - s.rotation_deg) * kDegToRad)) * kRadToDeg;
        worst_rotation_deg = std::max(worst_rotation_deg, off + std::abs(meridian.z));
    }
    check(worst_km < 0.01, "Kerbol bodies at the game's positions (km)", worst_km);
    check(worst_rotation_deg < 1e-6, "Kerbol bodies at the game's rotation angles (deg)", worst_rotation_deg);

    auto mean_motion = [&](const char* name) {
        const Body& body = scene.bodies[static_cast<size_t>(scene.find(name))];
        const MeanElementOrbit* orbit = mean_elements_of(body.motion.get());
        return orbit ? orbit->mean_motion() : 0.0; // rad/s
    };
    auto spin = [&](const char* name) { // rad/s
        return prime_meridian_rate(scene.bodies[static_cast<size_t>(scene.find(name))], 0.0) * kDegToRad /
               kSecondsPerDay;
    };

    // KSP's day: Kerbin's sidereal rotation and year make a solar day of 6 h exactly.
    const double solar_day = kTwoPi / (spin("Kerbin") - mean_motion("Kerbin"));
    check(std::abs(solar_day - kKerbinDaySeconds) < 1e-3, "Kerbin solar day = 6 h (s)", solar_day - kKerbinDaySeconds);

    // Tidally locked moons, and Duna locked to Ike (within the 2.8 ms of a float32 period).
    for (const char* moon : {"Mun", "Ike", "Laythe", "Vall", "Tylo", "Bop", "Pol"}) {
        const double n = mean_motion(moon);
        check(std::abs(spin(moon) - n) < 1e-12 * n, moon);
    }
    const double ike_period = kTwoPi / mean_motion("Ike");
    check(std::abs(kTwoPi / spin("Duna") - ike_period) < 0.01, "Duna turns with Ike's orbit (s)",
          kTwoPi / spin("Duna") - ike_period);

    // Laythe, Vall and Tylo: periods close to 1:2:4.
    const double vall = mean_motion("Laythe") / mean_motion("Vall");
    const double tylo = mean_motion("Laythe") / mean_motion("Tylo");
    check(std::abs(vall - 2.0) < 1e-5 && std::abs(tylo - 4.0) < 1e-4, "Laythe:Vall:Tylo = 1:2:4", tylo);

    // The schematic Mun transfer (tools/kerbol/kerbol_transfer.py): parking orbit 80 km
    // above Kerbin, the burns of the captions (across the 1 s step of each impulse), and
    // a circular orbit 30 km above the Mun after the capture.
    const int craft = scene.find("Mun transfer (schematic)");
    const int mun = scene.find("Mun");
    check(craft >= 0 && mun >= 0, "Kerbol transfer craft and the Mun");
    if (craft >= 0 && mun >= 0) {
        const MotionSource& c = *scene.bodies[static_cast<size_t>(craft)].motion;
        const MotionSource& m = *scene.bodies[static_cast<size_t>(mun)].motion;
        constexpr double kBurn = 7200.0;
        constexpr double kCapture = 31233.0;
        check(std::abs(glm::length(c.eval(kBurn - 60.0).position) - 680.0) < 0.1, "parking orbit 80 km up (km)",
              glm::length(c.eval(kBurn - 60.0).position) - 680.0);
        const double dv1 = glm::length(c.eval(kBurn + 1.5).velocity) - glm::length(c.eval(kBurn - 0.5).velocity);
        check(std::abs(dv1 * 1000.0 - 851.6) < 1.0, "transfer burn (m/s)", dv1 * 1000.0);
        auto mun_relative = [&](double t) {
            const State a = c.eval(t);
            const State b = m.eval(t);
            return State{a.position - b.position, a.velocity - b.velocity};
        };
        const double dv2 =
            glm::length(mun_relative(kCapture - 0.5).velocity) - glm::length(mun_relative(kCapture + 1.5).velocity);
        check(std::abs(dv2 * 1000.0 - 261.6) < 1.0, "capture burn (m/s)", dv2 * 1000.0);
        double worst = 0.0;
        for (double t = kCapture + 60.0; t < kCapture + 86000.0; t += 1234.5) {
            worst = std::max(worst, std::abs(glm::length(mun_relative(t).position) - 230.0));
        }
        check(worst < 0.1, "capture orbit 30 km above the Mun (km)", worst);
    }
}

void run_planet_scene_tests()
{
    test_galilean_rates();
    test_jupiter_heliocentric();
    test_io_shadow_transits();
    test_tidal_locking();
    test_jupiter_moon_ephemerides();
    test_solar_system_scene();
    test_solar_system_bodies();
    test_solar_system_belts();
    test_uranus_rings();
    test_ring_arcs();
    test_earth_moon_frame();
    test_earth_moon_rotation();
    test_mars_rotation();
    test_mars_moon_rotation();
    test_saturn_scene();
    test_saturn_coorbitals();
    test_kerbol_scene();
}
