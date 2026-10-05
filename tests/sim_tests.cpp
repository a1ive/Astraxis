// Minimal self-checking tests for core / ephem / scene (no framework).

#include "core/color.hpp"
#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/ephemeris.hpp"
#include "ephem/kepler.hpp"
#include "ephem/kerr_null.hpp"
#include "ephem/kerr_orbit.hpp"
#include "ephem/nbody.hpp"
#include "ephem/post_keplerian.hpp"
#include "ephem/visual_orbit.hpp"
#include "ephem/mean_element_orbit.hpp"
#include "scene/camera.hpp"
#include "scene/camera_director.hpp"
#include "scene/comet.hpp"
#include "scene/label_layout.hpp"
#include "scene/scene_loader.hpp"
#include "scene/simulation.hpp"
#include "scene/star_catalog.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace astraxis;

namespace {

int g_failures = 0;

void check(bool ok, const char* what, double value = 0.0)
{
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s (value = %.9g)\n", what, value);
    }
}

Scene load_jupiter()
{
    Scene scene;
    std::string error;
    if (!load_scene_file(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", scene, &error)) {
        std::printf("FAIL: cannot load jupiter.toml: %s\n", error.c_str());
        std::exit(1);
    }
    return scene;
}

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

void test_scene_loader()
{
    Scene scene = load_jupiter();
    check(scene.bodies.size() == 11, "jupiter.toml has 11 bodies", static_cast<double>(scene.bodies.size()));
    check(scene.bodies[0].rings.bands.size() == 5 && scene.bodies[0].rings.bands[1].name == "Main ring",
          "Jupiter's ring bands");
    check(scene.bodies[0].name == "Jupiter" && scene.bodies[0].parent == -1, "root is Jupiter");
    check(scene.bodies[4].name == "Callisto" && scene.bodies[4].parent == 0, "Callisto orbits Jupiter");
    check(std::abs(scene.bodies[0].polar_radius_km - 66854.0) < 1e-9, "Jupiter polar radius");
    check(scene.origin_heliocentric != nullptr, "heliocentric origin present");

    std::string error;
    Scene bad;
    check(!load_scene_string("name = 1", "bad1", bad, &error), "rejects non-string name");
    check(!load_scene_string("name = \"x\"\n[[bodies]]\nname = \"A\"\nradii_km = [1.0]\n"
                             "[[bodies]]\nname = \"B\"\nparent = \"C\"\nradii_km = [1.0]\n",
                             "bad2", bad, &error) &&
              error.find("parent 'C'") != std::string::npos,
          "rejects unknown parent");
    check(!load_scene_string("name = [", "bad3", bad, &error), "rejects TOML syntax error");
}

void test_star_catalog_and_color()
{
    std::vector<CatalogStar> stars;
    std::string error;
    check(load_star_catalog(ASTRAXIS_ASSET_DIR "/stars/bsc5.csv", stars, &error), "BSC5 loads");
    check(stars.size() > 9000, "BSC5 star count", static_cast<double>(stars.size()));
    double brightest = 99.0;
    for (const CatalogStar& s : stars) {
        brightest = std::min(brightest, s.vmag);
    }
    check(std::abs(brightest + 1.46) < 1e-9, "brightest star is Sirius (-1.46)", brightest);

    // The Sun (B-V = 0.65) is ~5800 K; blackbody colors go from red to blue.
    const double t_sun = temperature_from_bv(0.65);
    check(t_sun > 5600.0 && t_sun < 6000.0, "Sun temperature from B-V", t_sun);
    const glm::dvec3 cool = blackbody_linear_srgb(3000.0);
    const glm::dvec3 hot = blackbody_linear_srgb(15000.0);
    check(cool.r > cool.b && hot.b > hot.r, "blackbody color ordering");
    const glm::dvec3 d65ish = blackbody_linear_srgb(6500.0);
    check(d65ish.r > 0.9 && d65ish.g > 0.9 && d65ish.b > 0.85, "6500 K is near white", d65ish.b);
}

void test_kepler()
{
    const double eccentricities[] = {0.0, 0.01, 0.3, 0.7, 0.95};
    for (double e : eccentricities) {
        for (int k = -20; k <= 20; ++k) {
            const double m = k * 0.37;
            const double E = solve_kepler(m, e);
            const double residual = wrap_pi(E - e * std::sin(E) - m);
            check(std::abs(residual) < 1e-12, "Kepler equation residual", residual);
        }
    }
}

void test_calendar()
{
    CalendarDateTime c = calendar_from_jd(kJ2000Jd);
    check(c.year == 2000 && c.month == 1 && c.day == 1 && c.hour == 12 && c.minute == 0 && c.second == 0,
          "J2000 = 2000-01-01 12:00:00");

    c = calendar_from_jd(2299160.5);
    check(c.year == 1582 && c.month == 10 && c.day == 15, "Gregorian reform start 1582-10-15");

    c = calendar_from_jd(kUnixEpochJd);
    check(c.year == 1970 && c.month == 1 && c.day == 1 && c.hour == 0, "Unix epoch 1970-01-01");

    // 2026-10-02 00:00:00 UTC = Unix 1790899200.
    const double t = tdb_from_unix_utc(1790899200.0);
    check(format_utc(t) == "2026-10-02 00:00:00 UTC", "Unix -> UTC round trip");

    double parsed = 0.0;
    check(parse_utc("2026-10-02 00:00:00", &parsed) && std::abs(parsed - t) < 1e-3, "parse_utc", parsed - t);
    check(parse_utc("1979-03-05", &parsed) && format_utc(parsed) == "1979-03-05 00:00:00 UTC", "parse_utc date only");
    check(parse_utc("1066-10-14 09:30", &parsed) && format_utc(parsed) == "1066-10-14 09:30:00 UTC",
          "parse_utc Julian calendar");
    check(!parse_utc("2026-13-01", &parsed), "parse_utc rejects bad month");

    // Leap seconds: TDB - UTC = TAI-UTC + 32.184 s.
    check(std::abs(tdb_from_utc_jd(2451545.0) - (32.0 + 32.184)) < 1e-6, "TDB-UTC at J2000 (TAI-UTC 32 s)");
    parse_utc("1979-03-05 12:00:00", &parsed);
    const double jd_1979 = jd_from_calendar({1979, 3, 5, 12, 0, 0});
    check(std::abs(parsed - ((jd_1979 - kJ2000Jd) * kSecondsPerDay + 18.0 + 32.184)) < 1e-6,
          "TDB-UTC in 1979 (TAI-UTC 18 s)");
    check(format_utc(parsed) == "1979-03-05 12:00:00 UTC", "1979 UTC round trip");
    // Across the 2016-12-31 leap second the UTC clock advances 1 s less than TDB.
    double a = 0.0;
    double b = 0.0;
    parse_utc("2016-12-31 23:59:00", &a);
    parse_utc("2017-01-01 00:01:00", &b);
    check(std::abs((b - a) - 121.0) < 1e-4, "leap second at 2017-01-01", b - a);
    check(!parse_utc("hello", &parsed), "parse_utc rejects garbage");
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

Scene load_scene_or_die(const char* file)
{
    Scene scene;
    std::string error;
    if (!load_scene_file(std::string(ASTRAXIS_ASSET_DIR "/scenes/") + file, scene, &error)) {
        std::printf("FAIL: cannot load %s: %s\n", file, error.c_str());
        std::exit(1);
    }
    return scene;
}

double tdb_from_jd_tdb(double jd)
{
    return (jd - kJ2000Jd) * kSecondsPerDay;
}

void test_kepler_propagation()
{
    // Circular orbit (mu = 1, r = 1): a quarter period turns the state by 90 deg.
    const State circle{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    const State q = propagate_kepler(circle, 1.0, kPi / 2.0);
    check(glm::length(q.position - glm::dvec3(0.0, 1.0, 0.0)) < 1e-12 &&
              glm::length(q.velocity - glm::dvec3(-1.0, 0.0, 0.0)) < 1e-12,
          "Kepler quarter circle", glm::length(q.position - glm::dvec3(0.0, 1.0, 0.0)));

    // Eccentric ellipse and hyperbola: energy and angular momentum conserved,
    // and propagating back returns the start.
    for (const double speed : {1.2, 1.6}) { // escape speed is sqrt(2) ~ 1.414
        const State s0{{1.0, 0.0, 0.0}, {0.1, speed, 0.0}};
        const State s1 = propagate_kepler(s0, 1.0, 3.7);
        const State s2 = propagate_kepler(s1, 1.0, -3.7);
        auto energy = [](const State& s) { return 0.5 * glm::dot(s.velocity, s.velocity) - 1.0 / glm::length(s.position); };
        const double de = std::abs(energy(s1) - energy(s0));
        const double dh = glm::length(glm::cross(s1.position, s1.velocity) - glm::cross(s0.position, s0.velocity));
        check(de < 1e-12 && dh < 1e-12, "Kepler propagation conserves energy and momentum", std::max(de, dh));
        check(glm::length(s2.position - s0.position) < 1e-10, "Kepler propagation round trip",
              glm::length(s2.position - s0.position));
    }
}

void test_ephemeris()
{
    EphemerisTable v1;
    std::string error;
    check(v1.load(ASTRAXIS_ASSET_DIR "/ephem/voyager1.eph", &error), "voyager1.eph loads");
    check(v1.target_id() == -31 && v1.center_id() == 10, "voyager1.eph ids");

    // Horizons, Voyager 1 @ Sun, 1977-Sep-07 00:00 TDB (JD 2443393.5), ICRF km.
    const State s = v1.eval(tdb_from_jd_tdb(2443393.5));
    const glm::dvec3 horizons(1.457774659041131E+08, -3.572921412462115E+07, -1.539445114171495E+07);
    check(glm::length(s.position - horizons) < 20.0, "Voyager 1 matches Horizons (km)",
          glm::length(s.position - horizons));

    // Hermite velocity is the derivative of the interpolated position.
    const double t = tdb_from_jd_tdb(2444303.0); // 1979-03-05, near the Jupiter flyby (dense knots)
    const glm::dvec3 numeric = (v1.eval(t + 1.0).position - v1.eval(t - 1.0).position) / 2.0;
    check(glm::length(numeric - v1.eval(t).velocity) < 1e-3, "Hermite velocity consistent (km/s)",
          glm::length(numeric - v1.eval(t).velocity));

    EphemerisTable broken;
    check(!broken.load(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error), "rejects a non-ephemeris file");
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
        {"Phobos", "Mars", {-2628.223, -8355.790, -2928.529}, 4.3},
        {"Deimos", "Mars", {10863.689, 20169.954, 5021.957}, 0.7},
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

// Solar System scene: Pioneer 10/11 and New Horizons, and the planets baked
// back to 1972 for the Pioneer flybys.
void test_solar_system_missions()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    auto at = [&](const char* name, double t) { return scene.icrf_state_at(scene.find(name), t).position; };
    auto min_distance = [&](const char* a, const char* b, double t_center, double half, double step) {
        double best = 1e300;
        for (double dt = -half; dt <= half; dt += step) {
            best = std::min(best, glm::length(at(a, t_center + dt) - at(b, t_center + dt)));
        }
        return best;
    };

    // New Horizons' closest approaches against Horizons relative to the body centers
    // (@ 500@999 and @ 500@2486958, 1-minute samples, parabolic minimum).
    const double ca_pluto = tdb_from_jd_tdb(2457217.99279); // 2015-07-14 11:49:37 TDB
    const double d_pluto = min_distance("New Horizons", "Pluto", ca_pluto, 300.0, 1.0);
    check(std::abs(d_pluto - 13676.0) < 30.0, "New Horizons - Pluto closest approach (km)", d_pluto);
    const double ca_arrokoth = tdb_from_jd_tdb(2458484.73230); // 2019-01-01 05:34:31 TDB
    const double d_arrokoth = min_distance("New Horizons", "Arrokoth", ca_arrokoth, 300.0, 1.0);
    check(std::abs(d_arrokoth - 3538.0) < 60.0, "New Horizons - Arrokoth closest approach (km)", d_arrokoth);

    // Pioneer flybys agree with the bake tool's report (same data, independent code).
    struct Flyby {
        const char* craft;
        const char* planet;
        CalendarDateTime tdb;
        double km;
    };
    for (const Flyby& f : {Flyby{"Pioneer 10", "Jupiter", {1973, 12, 4, 2, 25, 46}, 202870.0},
                           Flyby{"Pioneer 11", "Jupiter", {1974, 12, 3, 5, 22, 11}, 113553.0},
                           Flyby{"Pioneer 11", "Saturn", {1979, 9, 1, 16, 30, 19}, 80008.0}}) {
        const double t = (jd_from_calendar(f.tdb) - kJ2000Jd) * kSecondsPerDay;
        const double d = min_distance(f.craft, f.planet, t, 600.0, 2.0);
        check(std::abs(d - f.km) < 100.0, f.craft, d);
    }

    // Lifetimes: not launched yet, then linearly extrapolated after the tables end (2050).
    double t = 0.0;
    parse_utc("1972-03-02", &t);
    scene.update(t);
    check(!scene.bodies[static_cast<size_t>(scene.find("Pioneer 10"))].visible, "Pioneer 10 before launch");
    parse_utc("2060-01-01", &t);
    scene.update(t);
    check(scene.bodies[static_cast<size_t>(scene.find("Pioneer 10"))].visible &&
              scene.bodies[static_cast<size_t>(scene.find("New Horizons"))].visible,
          "Pioneer 10 and New Horizons extrapolated after 2050");

    // The planets' tables now start in 1972, before the Pioneer 10 flyby.
    const auto* jup = dynamic_cast<const EphemerisMotion*>(scene.bodies[static_cast<size_t>(scene.find("Jupiter"))].motion.get());
    parse_utc("1973-12-01", &t);
    check(jup != nullptr && jup->fallback() != nullptr &&
              glm::length(jup->eval(t).position - jup->fallback()->eval(t).position) > 1000.0,
          "Jupiter in 1973 comes from the baked table, not the fallback");
}

// Signed volume of a closed mesh (positive: counter-clockwise seen from outside).
double signed_volume(const ShapeModel& shape)
{
    double volume = 0.0;
    for (size_t k = 0; k < shape.indices.size(); k += 3) {
        const glm::dvec3 a(shape.positions[shape.indices[k]]);
        const glm::dvec3 b(shape.positions[shape.indices[k + 1]]);
        const glm::dvec3 c(shape.positions[shape.indices[k + 2]]);
        volume += glm::dot(a, glm::cross(b, c)) / 6.0;
    }
    return volume;
}

// Largest difference (turns) between a mesh's map u and the east longitude of its
// vertex directions (seam vertices may be 0 or 1; poles skipped).
double map_u_error(const ShapeModel& shape)
{
    double worst = 0.0;
    for (size_t k = 0; k < shape.positions.size() && k < shape.map_u.size(); ++k) {
        const glm::vec3& p = shape.positions[k];
        if (std::hypot(p.x, p.y) < 1e-3 * glm::length(p)) {
            continue;
        }
        const double d = std::atan2(p.y, p.x) / (2.0 * kPi) - shape.map_u[k];
        worst = std::max(worst, std::abs(d - std::round(d)));
    }
    return worst;
}

// Arrokoth's shape model and pole against independent numbers in Porter et al. 2024,
// Table 2: the volume (equal-volume diameter 19.896 km) and the pole's inclination to
// the heliocentric orbit (100.39 deg) and to New Horizons' approach direction (41.096 deg).
void test_arrokoth_shape()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    const int arrokoth = scene.find("Arrokoth");
    const Body& body = scene.bodies[static_cast<size_t>(arrokoth)];
    check(body.shape != nullptr, "Arrokoth has a shape model");
    if (!body.shape) {
        return;
    }
    const ShapeModel& shape = *body.shape;
    glm::dvec3 lo(1e300), hi(-1e300);
    for (const glm::vec3& p : shape.positions) {
        lo = glm::min(lo, glm::dvec3(p));
        hi = glm::max(hi, glm::dvec3(p));
    }
    const glm::dvec3 size = hi - lo;
    check(std::abs(size.x - 34.546) < 0.01 && std::abs(size.y - 19.838) < 0.01 && std::abs(size.z - 13.822) < 0.01,
          "Arrokoth overall dimensions (km)", size.x);
    check(std::abs(size.x - 2.0 * body.equatorial_radius_km) < 0.01 &&
              std::abs(size.z - 2.0 * body.polar_radius_km) < 0.01,
          "Arrokoth radii are half its dimensions", body.equatorial_radius_km);

    const double d_equal = std::cbrt(6.0 * signed_volume(shape) / kPi);
    check(std::abs(d_equal - 19.896) < 0.01, "Arrokoth equal-volume diameter (km)", d_equal);
    const auto [min_albedo, max_albedo] = std::minmax_element(shape.albedo.begin(), shape.albedo.end());
    check(*min_albedo > 0.5f && *max_albedo < 2.5f, "Arrokoth relative albedo range", *max_albedo);
    check(shape.map_u.empty(), "Arrokoth's albedo is per vertex (no map coordinates)");

    auto angle_deg = [](const glm::dvec3& a, const glm::dvec3& b) {
        return std::acos(std::clamp(glm::dot(glm::normalize(a), glm::normalize(b)), -1.0, 1.0)) / kDegToRad;
    };
    const double ca = tdb_from_jd_tdb(2458484.73230); // 2019-01-01 05:34:31 TDB
    const State orbit = scene.icrf_state_at(arrokoth, ca);
    const double to_orbit = angle_deg(body.pole, glm::cross(orbit.position, orbit.velocity));
    check(std::abs(to_orbit - 100.39) < 0.05, "Arrokoth pole to its orbit normal (deg)", to_orbit);
    // Toward the approaching spacecraft, a day out (the approach is nearly straight).
    const glm::dvec3 approach =
        scene.icrf_state_at(scene.find("New Horizons"), ca - kSecondsPerDay).position -
        scene.icrf_state_at(arrokoth, ca - kSecondsPerDay).position;
    const double to_approach = angle_deg(body.pole, approach);
    check(std::abs(to_approach - (180.0 - 41.096)) < 0.1, "Arrokoth pole to the approach direction (deg)", to_approach);
    std::printf("info: Arrokoth pole to orbit normal %.2f deg, to the approaching spacecraft %.3f deg\n", to_orbit,
                to_approach);
}

// Phobos and Deimos: Thomas' shape models (west longitudes in the tables, east in the
// meshes) with map coordinates for their textures. Volumes against the mean radii in
// JPL SSD's satellite physical parameters (Archinal et al. 2018), and Stickney (49 W)
// as the check that the longitudes were turned around.
void test_mars_moon_shapes()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    struct Moon {
        const char* name;
        double mean_radius_km; // [SPHY]
        double tolerance_km;
    };
    for (const Moon& moon : {Moon{"Phobos", 11.08, 0.05}, Moon{"Deimos", 6.2, 0.1}}) {
        const Body& body = scene.bodies[static_cast<size_t>(scene.find(moon.name))];
        const std::string name = moon.name;
        check(body.shape != nullptr && !body.texture.empty(), (name + " has a shape model and a texture").c_str());
        if (!body.shape) {
            continue;
        }
        const ShapeModel& shape = *body.shape;
        check(shape.map_u.size() == shape.positions.size(), (name + " mesh has map coordinates").c_str());
        const double r_equal = std::cbrt(3.0 * signed_volume(shape) / (4.0 * kPi));
        check(std::abs(r_equal - moon.mean_radius_km) < moon.tolerance_km, (name + " equal-volume radius (km)").c_str(),
              r_equal);

        // Long axis toward Mars (the prime meridian), shortest along the pole, like the PCK ellipsoid.
        glm::dvec3 extent(0.0);
        for (const glm::vec3& p : shape.positions) {
            extent = glm::max(extent, glm::abs(glm::dvec3(p)));
        }
        check(extent.x > extent.y && extent.y > extent.z, (name + " axes ordered x > y > z").c_str(), extent.x);
        check(std::abs(extent.x - body.equatorial_radius_km) < 1.0, (name + " long axis near the PCK radius (km)").c_str(),
              extent.x);

        const double u_error = map_u_error(shape);
        check(u_error < 1e-6, (name + " map u matches the vertex longitude").c_str(), u_error);
    }

    // Stickney: the mean radius within 12 deg of its center is well below that of the
    // surrounding ring at 2 S 50 W, and not at the mirrored 2 S 50 E.
    const ShapeModel* phobos = scene.bodies[static_cast<size_t>(scene.find("Phobos"))].shape.get();
    if (!phobos) {
        return;
    }
    auto crater_depth = [&](double lon_east_deg) {
        const double lat = -2.0 * kDegToRad, lon = lon_east_deg * kDegToRad;
        const glm::dvec3 center(std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat));
        double inner = 0.0, ring = 0.0;
        int n_inner = 0, n_ring = 0;
        for (const glm::vec3& p : phobos->positions) {
            const double r = glm::length(glm::dvec3(p));
            const double angle = std::acos(std::clamp(glm::dot(glm::dvec3(p) / r, center), -1.0, 1.0)) / kDegToRad;
            if (angle < 12.0) {
                inner += r;
                ++n_inner;
            } else if (angle > 16.0 && angle < 24.0) {
                ring += r;
                ++n_ring;
            }
        }
        return ring / n_ring - inner / n_inner;
    };
    const double depth_west = crater_depth(-50.0);
    const double depth_east = crater_depth(50.0);
    check(depth_west > 0.5, "Stickney is a depression at 50 W (km)", depth_west);
    check(depth_west > depth_east + 0.5, "no Stickney at the mirrored 50 E (km)", depth_east);
    std::printf("info: Phobos mean radius within 12 deg of 2 S 50 W is %.2f km below the ring (50 E: %.2f km)\n",
                depth_west, depth_east);
}

// Vesta: the DLR Dawn HAMO DTM as a 1.5 deg grid, with the mosaic of the same release.
// Its volume against SBDB's equivalent diameter (522.77 +- 0.1 km, Park et al. 2025);
// the 2013 DTM and the grid's smoothing allow a few tenths of a km.
void test_vesta_shape()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    const Body& body = scene.bodies[static_cast<size_t>(scene.find("Vesta"))];
    check(body.shape != nullptr && !body.texture.empty(), "Vesta has a shape model and a texture");
    if (!body.shape) {
        return;
    }
    const ShapeModel& shape = *body.shape;
    const double d_equal = 2.0 * std::cbrt(3.0 * signed_volume(shape) / (4.0 * kPi));
    check(std::abs(d_equal - 522.77) < 0.6, "Vesta equal-volume diameter (km)", d_equal);
    glm::dvec3 lo(1e300), hi(-1e300);
    for (const glm::vec3& p : shape.positions) {
        lo = glm::min(lo, glm::dvec3(p));
        hi = glm::max(hi, glm::dvec3(p));
    }
    const glm::dvec3 size = hi - lo;
    check(size.x > size.y && size.y > size.z, "Vesta axes ordered x > y > z", size.x);
    const double u_error = map_u_error(shape);
    check(shape.map_u.size() == shape.positions.size() && u_error < 1e-6, "Vesta map u matches the vertex longitude",
          u_error);
}

// Pallas: the VLT/SPHERE MPCD shape model ([PALM] in solar_system.toml). Volume against
// Vernazza et al. 2021 (Table 1: D = 511 +- 4 km); the scene's IAU rotation against the
// spin solution of the ADAM model whose body frame the mesh shares (DAMIT model 4395,
// Marsset et al. 2020): lambda 42, beta -15 deg, P 7.81322 h, phi0 = 0 at JD 2433827.77154,
// body -> J2000 ecliptic = Rz(lambda) Ry(90 deg - beta) Rz(phi).
void test_pallas_shape()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    scene.set_active_frame(0); // inertial (ICRF axes)
    const int pallas = scene.find("Pallas");
    const Body& body = scene.bodies[static_cast<size_t>(pallas)];
    check(body.shape != nullptr && body.texture.empty(), "Pallas has a shape model (no texture)");
    if (!body.shape) {
        return;
    }
    const double d_equal = std::cbrt(6.0 * signed_volume(*body.shape) / kPi);
    check(std::abs(d_equal - 511.0) < 4.0, "Pallas equal-volume diameter (km)", d_equal);

    const double period_h = 360.0 * 24.0 / body.pm_rate_deg_per_day;
    check(std::abs(period_h - 7.81322) < 1e-5, "Pallas sidereal period (h)", period_h);
    // DAMIT's IAU conversion (from the pole rounded to RA 44, Dec 1): W0 = 41.7 deg.
    check(std::abs(body.pm_w0_deg - 41.7) < 0.15, "Pallas W0 near DAMIT's IAU value (deg)", body.pm_w0_deg);

    auto rotation_y = [](double a) {
        glm::dmat3 m(1.0); // columns
        m[0] = glm::dvec3(std::cos(a), 0.0, -std::sin(a));
        m[2] = glm::dvec3(std::sin(a), 0.0, std::cos(a));
        return m;
    };
    const double obliquity = 84381.448 / 3600.0 * kDegToRad; // IAU 1976
    const glm::dmat3 expected =
        rotation_x(obliquity) * rotation_z(42.0 * kDegToRad) * rotation_y((90.0 + 15.0) * kDegToRad);
    scene.update(tdb_from_jd_tdb(2433827.77154)); // phi = 0
    double worst = 0.0;
    for (int k : {0, 2}) {
        const double cos_angle = glm::dot(scene.bodies[static_cast<size_t>(pallas)].orientation[k], expected[k]);
        worst = std::max(worst, std::acos(std::clamp(cos_angle, -1.0, 1.0)) / kDegToRad);
    }
    check(worst < 0.05, "Pallas body axes match the DAMIT spin solution at t0 (deg)", worst);
}

// Amalthea and Thebe (Stooke, Jupiter scene) and Hyperion (Thomas, Cassini; Saturn
// scene): volumes against the mean radii in JPL SSD's satellite physical parameters
// (within their uncertainties). The Jovian moons' long axes point at Jupiter (+x).
void test_jupiter_saturn_small_moon_shapes()
{
    Scene jupiter = load_scene_or_die("jupiter.toml");
    Scene saturn = load_scene_or_die("saturn.toml");
    struct Moon {
        const Scene* scene;
        const char* name;
        double mean_radius_km; // [SPHY]
        double sigma_km;
        bool long_axis_x;
        bool textured;
    };
    for (const Moon& m : {Moon{&jupiter, "Amalthea", 83.5, 3.0, true, true}, Moon{&jupiter, "Thebe", 49.3, 4.0, true, false},
                          Moon{&saturn, "Hyperion", 135.0, 4.0, false, false}}) {
        const std::string name = m.name;
        const Body& body = m.scene->bodies[static_cast<size_t>(m.scene->find(m.name))];
        check(body.shape != nullptr, (name + " has a shape model").c_str());
        if (!body.shape) {
            continue;
        }
        const double r_equal = std::cbrt(3.0 * signed_volume(*body.shape) / (4.0 * kPi));
        check(std::abs(r_equal - m.mean_radius_km) < m.sigma_km, (name + " equal-volume radius (km)").c_str(), r_equal);
        check(body.texture.empty() != m.textured, (name + " texture as expected").c_str());
        if (m.textured) {
            check(map_u_error(*body.shape) < 1e-6, (name + " map u matches the vertex longitude").c_str());
        }
        if (m.long_axis_x) {
            glm::dvec3 extent(0.0);
            for (const glm::vec3& p : body.shape->positions) {
                extent = glm::max(extent, glm::abs(glm::dvec3(p)));
            }
            check(extent.x > extent.y && extent.x > extent.z, (name + " long axis toward Jupiter").c_str(), extent.x);
        }
    }

    // Hyperion spins at |omega|/n = 4.255 (Cassini, 2005-08-16 and 2005-09-25; Harbison et
    // al. 2011), n from [ELEM]'s mean-longitude period 21.276658 d: about 72 deg/day.
    const Body& hyperion = saturn.bodies[static_cast<size_t>(saturn.find("Hyperion"))];
    const double spin_over_n = hyperion.pm_rate_deg_per_day / (360.0 / 21.276658);
    check(std::abs(spin_over_n - 4.255) < 1e-4, "Hyperion spin rate / mean motion", spin_over_n);
}

// Halley (ISEE-3 scene): Stooke's shape model ([STKS]) against Szego 1991's independent
// Vega/Giotto model ([SZG]: a 15.3 x 7.2 x 7.22 km box, 365 km^3), and the long-axis mode
// of Belton et al. 1991 ([BLT]): from the long axis at 66.0 deg to M, 3.69 d around it and
// 7.1 d about itself, the abstract derives a total spin of 2.84 d at 21.4 deg from M.
void test_halley_shape()
{
    Scene scene = load_scene_or_die("isee3.toml");
    scene.set_active_frame(0); // inertial (ICRF axes)
    const int halley = scene.find("Halley");
    const Body& body = scene.bodies[static_cast<size_t>(halley)];
    check(body.shape != nullptr && body.texture.empty(), "Halley has a shape model (no texture)");
    check(body.free_precession.enabled && body.pm_rate_deg_per_day == 0.0, "Halley precesses freely");
    if (!body.shape) {
        return;
    }
    glm::dvec3 lo(1e300), hi(-1e300);
    for (const glm::vec3& p : body.shape->positions) {
        lo = glm::min(lo, glm::dvec3(p));
        hi = glm::max(hi, glm::dvec3(p));
    }
    const glm::dvec3 size = hi - lo;
    check(std::abs(size.z - 15.3) < 0.5, "Halley long axis along z (km)", size.z);
    check(std::abs(size.x - 7.2) < 0.5 && std::abs(size.y - 7.2) < 0.5, "Halley short axes (km)",
          std::max(size.x, size.y));
    const double volume = signed_volume(*body.shape);
    check(std::abs(volume / 365.0 - 1.0) < 0.15, "Halley volume near Szego's model (km^3)", volume);

    // Angular velocity from the orientation a minute apart: R' R^T = [omega]x.
    const double t0 = tdb_from_utc_jd(jd_from_calendar({1986, 3, 14, 0, 0, 0}));
    const double dt = 60.0;
    double worst_nutation = 0.0;
    for (double days : {0.0, 1.3, 2.9}) {
        const double t = t0 + days * kSecondsPerDay;
        scene.update(t);
        const glm::dmat3 r0 = scene.bodies[static_cast<size_t>(halley)].orientation;
        scene.update(t + dt);
        const glm::dmat3 r1 = scene.bodies[static_cast<size_t>(halley)].orientation;
        const glm::dmat3 w = (r1 - r0) * glm::transpose(r0) / dt;
        const glm::dvec3 omega(w[1][2], w[2][0], w[0][1]);
        const double period_d = kTwoPi / glm::length(omega) / kSecondsPerDay;
        const double spin_to_m = std::acos(glm::dot(glm::normalize(omega), body.pole)) / kDegToRad;
        const double nutation = std::acos(glm::dot(r0[2], body.pole)) / kDegToRad;
        worst_nutation = std::max(worst_nutation, std::abs(nutation - 66.0));
        check(std::abs(period_d - 2.84) < 0.01, "Halley total spin period (d)", period_d);
        check(std::abs(spin_to_m - 21.4) < 0.1, "Halley spin vector to M (deg)", spin_to_m);
        if (days == 0.0) {
            std::printf("info: Halley mesh %.2f x %.2f x %.2f km, %.1f km^3; total spin %.3f d at %.2f deg from M\n",
                        size.x, size.y, size.z, volume, period_d, spin_to_m);
        }
    }
    check(worst_nutation < 1e-6, "Halley long axis stays 66 deg from M (deg)", worst_nutation);
}

// Comets in the ISEE-3 scene ([SBDB], [SLV], [ESO] in isee3.toml). The ion tail points
// along the solar wind as the comet sees it: a few degrees off the anti-sun direction,
// toward where the comet came from. ICE flew through Giacobini-Zinner's tail 7,800 km
// behind the nucleus (Bame et al. 1986, Science 232, 356): at closest approach it should
// be downstream and within the tail's radius (1.2e4 km across, [SLV]) of the axis.
void test_comet_tails()
{
    Scene scene = load_scene_or_die("isee3.toml");
    const int gz = scene.find("Giacobini-Zinner");
    const int halley = scene.find("Halley");
    const int craft = scene.find("ISEE-3");
    const Body::Comet& gz_comet = scene.bodies[static_cast<size_t>(gz)].comet;
    const Body::Comet& halley_comet = scene.bodies[static_cast<size_t>(halley)].comet;
    check(gz_comet.enabled && halley_comet.enabled, "isee3.toml comets");
    check(gz_comet.m1 == 13.2 && gz_comet.k1 == 8.25 && halley_comet.m1 == 5.5 && halley_comet.k1 == 8.0,
          "comet magnitude laws [SBDB]");
    check(std::abs(comet_total_magnitude(5.5, 8.0, 0.1) - (5.5 - 8.0)) < 1e-12, "total magnitude law");

    auto utc = [](const CalendarDateTime& c) { return tdb_from_utc_jd(jd_from_calendar(c)); };
    const double t_ca = utc({1985, 9, 11, 11, 2, 43});
    const State sun = scene.sun_icrf_state_at(t_ca);
    const State comet = scene.icrf_state_at(gz, t_ca);
    const glm::dvec3 r = comet.position - sun.position;
    const glm::dvec3 v = comet.velocity - sun.velocity;
    const glm::dvec3 axis = ion_tail_direction(r, v, gz_comet.solar_wind_km_s);

    // Aberration: atan(v_perp / v_wind) away from the anti-sun direction, against the motion.
    const glm::dvec3 out = glm::normalize(r);
    const glm::dvec3 v_perp = v - glm::dot(v, out) * out;
    const double expected = std::atan2(glm::length(v_perp), gz_comet.solar_wind_km_s - glm::dot(v, out));
    const double aberration = std::acos(std::clamp(glm::dot(axis, out), -1.0, 1.0));
    check(std::abs(aberration - expected) < 1e-9, "ion tail aberration angle", aberration / kDegToRad);
    check(glm::dot(axis, v_perp) < 0.0, "ion tail lags behind the motion");

    const glm::dvec3 rel = scene.icrf_state_at(craft, t_ca).position - comet.position;
    const double downstream = glm::dot(rel, axis);
    const double off_axis = glm::length(rel - downstream * axis);
    check(downstream > 0.0, "ICE passed on the tail side (km)", downstream);
    check(off_axis < 0.5 * gz_comet.ion_tail_diameter_km, "ICE within the ion tail's radius (km)", off_axis);
    const double off_anti_sun = glm::length(rel - glm::dot(rel, out) * out);
    std::printf("info: G-Z tail %.2f deg from anti-sun; ICE %.0f km downstream, %.0f km off its axis (%.0f km off the "
                "anti-sun line)\n",
                aberration / kDegToRad, downstream, off_axis, off_anti_sun);
}

// Finson-Probstein dust (comet.hpp). A grain of radiation pressure ratio beta, released
// without ejection speed, drifts from the beta = 0 one by the radiation pressure
// acceleration: beta GM / r^2 outward, 1/2 a t^2 after a short time. Halley's tail on
// 1986-03-10: grains on the anti-sun side, curving back along the orbit (lagging the
// ion tail, which only aberrates by v / 400 km/s).
void test_dust_tail()
{
    const double gm = 1.32712440041279419e11; // [GM] in isee3.toml
    const State comet{glm::dvec3(0.6 * kAuKm, 0.0, 0.0), glm::dvec3(0.0, 50.0, 5.0)};
    const double dt = 6.0 * 3600.0;
    const glm::dvec3 drift = dust_grain_state(comet, glm::dvec3(0.0), 0.5, gm, dt).position -
                             dust_grain_state(comet, glm::dvec3(0.0), 0.0, gm, dt).position;
    const double expected = 0.5 * 0.5 * gm / (0.6 * kAuKm * 0.6 * kAuKm) * dt * dt;
    check(std::abs(glm::length(drift) / expected - 1.0) < 1e-3, "dust drift by radiation pressure (km)",
          glm::length(drift));
    check(glm::dot(glm::normalize(drift), glm::dvec3(1.0, 0.0, 0.0)) > 0.9999, "dust drift away from the sun");

    Scene scene = load_scene_or_die("isee3.toml");
    const int halley = scene.find("Halley");
    const double t = tdb_from_utc_jd(jd_from_calendar({1986, 3, 10, 0, 0, 0}));
    DustTail tail;
    tail.update(scene, halley, t);
    const auto& grains = tail.grains();
    const State sun = scene.sun_icrf_state_at(t);
    const State c = scene.icrf_state_at(halley, t);
    const glm::dvec3 out = glm::normalize(c.position - sun.position);
    const glm::dvec3 v = c.velocity - sun.velocity;
    const glm::dvec3 v_perp = glm::normalize(v - glm::dot(v, out) * out);
    glm::dvec3 mean(0.0);
    double total = 0.0;
    double far = 0.0;
    for (const DustGrain& g : grains) {
        mean += g.weight * (g.position - c.position);
        total += g.weight;
        far = std::max(far, glm::length(g.position - c.position));
    }
    mean /= total;
    const double anti_sun = glm::dot(mean, out);
    const double lag = glm::dot(mean, v_perp);
    const glm::dvec3 ion = ion_tail_direction(c.position - sun.position, v, 400.0);
    check(grains.size() > DustTail::kSlots * DustTail::kBetas * DustTail::kDirections * 9 / 10, "Halley dust grains",
          static_cast<double>(grains.size()));
    check(anti_sun > 0.0, "Halley dust tail on the anti-sun side (km)", anti_sun);
    check(lag / anti_sun < glm::dot(ion, v_perp) / glm::dot(ion, out), "Halley dust tail lags the ion tail");
    // Cost of a frame: the next update a minute later (cached slots, warm-started anomalies).
    const auto clock_start = std::chrono::steady_clock::now();
    for (int k = 1; k <= 20; ++k) {
        tail.update(scene, halley, t + 60.0 * k);
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - clock_start).count();
    std::printf("info: dust tail update %.2f ms per frame\n", ms / 20.0);
    // Giacobini-Zinner at ICE's pass: also a tail on the anti-sun side.
    {
        const int gz = scene.find("Giacobini-Zinner");
        const double t_gz = tdb_from_utc_jd(jd_from_calendar({1985, 9, 11, 11, 0, 0}));
        DustTail gz_tail;
        gz_tail.update(scene, gz, t_gz);
        const State gs = scene.sun_icrf_state_at(t_gz);
        const State gc = scene.icrf_state_at(gz, t_gz);
        glm::dvec3 gz_mean(0.0);
        double gz_total = 0.0;
        for (const DustGrain& g : gz_tail.grains()) {
            gz_mean += g.weight * (g.position - gc.position);
            gz_total += g.weight;
        }
        gz_mean /= gz_total;
        const double gz_anti_sun = glm::dot(gz_mean, glm::normalize(gc.position - gs.position));
        check(gz_anti_sun > 0.8 * glm::length(gz_mean), "Giacobini-Zinner dust tail on the anti-sun side (km)",
              gz_anti_sun);
    }
    std::printf("info: Halley dust 1986-03-10: %zu grains, weighted center %.3g km anti-sun, %.3g km along -v "
                "(%.1f deg behind the anti-sun line), farthest %.3g km\n",
                grains.size(), anti_sun, -lag, std::atan2(-lag, anti_sun) / kDegToRad, far);
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

void test_parker_scene()
{
    Scene scene = load_scene_or_die("parker.toml");
    const int parker = scene.find("Parker Solar Probe");
    const int venus = scene.find("Venus");
    check(parker > 0 && venus > 0, "parker.toml bodies");

    // Horizons, Parker Solar Probe (-96) and Venus (299) @ Sun, 2022-Jan-01 00:00 TDB, ICRF km.
    const double t_ref = tdb_from_jd_tdb(2459580.5);
    const glm::dvec3 parker_hzn(8.718605417334306E+07, -6.179560978012869E+07, -3.332843147602009E+07);
    const glm::dvec3 venus_hzn(-1.015245425299603E+07, 9.746038667700523E+07, 4.449519499322341E+07);
    const double e_parker = glm::length(scene.icrf_state_at(parker, t_ref).position - parker_hzn);
    const double e_venus = glm::length(scene.icrf_state_at(venus, t_ref).position - venus_hzn);
    check(e_parker < 5.0, "Parker matches Horizons (km)", e_parker);
    check(e_venus < 5.0, "Venus matches Horizons (km)", e_venus);

    // Each Venus flyby takes orbital energy away: the osculating period around the
    // Sun drops across every closest approach (bake tool times, TDB), ending at the
    // 88-day orbit of the mission design (Horizons spacecraft notes).
    const CalendarDateTime flybys[] = {{2018, 10, 3, 8, 45, 36},  {2019, 12, 26, 18, 15, 55}, {2020, 7, 11, 3, 24, 49},
                                   {2021, 2, 20, 20, 6, 53},  {2021, 10, 16, 9, 31, 58},  {2023, 8, 21, 12, 4, 3},
                                   {2024, 11, 6, 18, 44, 49}};
    auto period_days = [&](double t) {
        scene.update(t);
        double rp = 0.0;
        double ra = 0.0;
        double period = 0.0;
        scene.osculating_apsides(parker, &rp, &ra, &period);
        return period / kSecondsPerDay;
    };
    bool all_lower = true;
    double last_period = 0.0;
    for (const CalendarDateTime& d : flybys) {
        const double t = (jd_from_calendar(d) - kJ2000Jd) * kSecondsPerDay;
        const double before = period_days(t - 5.0 * kSecondsPerDay);
        last_period = period_days(t + 5.0 * kSecondsPerDay);
        all_lower = all_lower && last_period < before - 1.0;
    }
    check(all_lower, "every Venus flyby shortens the orbit", last_period);
    check(std::abs(last_period - 88.0) < 1.0, "final orbital period (days)", last_period);

    // The closest flyby (7) agrees with the bake tool and stays above the surface.
    const double ca7 = (jd_from_calendar(flybys[6]) - kJ2000Jd) * kSecondsPerDay;
    double best = 1e300;
    for (double dt = -600.0; dt <= 600.0; dt += 1.0) {
        best = std::min(best, glm::length(scene.icrf_state_at(parker, ca7 + dt).position -
                                          scene.icrf_state_at(venus, ca7 + dt).position));
    }
    check(std::abs(best - 6423.0) < 20.0, "Venus flyby 7 distance (km)", best);

    // Perihelion 22 at 9.86 solar radii, found by the scene's own periapsis search.
    const double p22 = (jd_from_calendar({2024, 12, 24, 11, 42, 5}) - kJ2000Jd) * kSecondsPerDay;
    const std::vector<double> peri = scene.periapsis_times(parker, p22 - 10.0 * kSecondsPerDay, p22 + 10.0 * kSecondsPerDay);
    const double r22 = peri.size() == 1 ? glm::length(scene.icrf_state_at(parker, peri[0]).position) : 0.0;
    check(peri.size() == 1 && std::abs(peri[0] - p22) < 60.0, "perihelion 22 time (s)",
          peri.empty() ? -1.0 : peri[0] - p22);
    check(std::abs(r22 / 695700.0 - 9.86) < 0.01, "perihelion 22 (solar radii)", r22 / 695700.0);

    // Sun-Venus rotating frame: Venus stays on the +x axis.
    for (size_t f = 0; f < scene.frames.size(); ++f) {
        if (scene.frames[f].name == "Sun-Venus rotating") {
            scene.set_active_frame(static_cast<int>(f));
        }
    }
    scene.update(t_ref);
    const glm::dvec3 w = scene.bodies[static_cast<size_t>(venus)].world_position;
    check(w.x > 1.0e8 && std::abs(w.y) < 1e-3 && std::abs(w.z) < 1e-3, "Venus fixed in rotating frame", w.y);

    // The history trail since launch, seen in the rotating frame, has no sharp
    // corners: the knots alone leave the aphelion loops as coarse polygons.
    double t_now = 0.0;
    parse_utc("2026-10-01", &t_now);
    scene.update(t_now);
    std::vector<glm::dvec3> points;
    std::vector<float> fades;
    scene.trail(parker, t_now, 2048, points, fades);
    double worst_turn_deg = 0.0;
    for (size_t k = 1; k + 1 < points.size(); ++k) {
        const glm::dvec3 a = points[k] - points[k - 1];
        const glm::dvec3 b = points[k + 1] - points[k];
        const double c = glm::dot(a, b) / (glm::length(a) * glm::length(b));
        worst_turn_deg = std::max(worst_turn_deg, std::acos(std::clamp(c, -1.0, 1.0)) * kRadToDeg);
    }
    check(points.size() > 2048 && worst_turn_deg < 20.0, "rotating-frame trail is smooth (deg)", worst_turn_deg);
    std::printf("info: Parker rotating-frame trail %zu points, sharpest turn %.1f deg\n", points.size(), worst_turn_deg);
}

// Galileo and Juno in the Jupiter scene: Kepler-relative tables against
// Horizons, moon flybys against the distances NASA/JPL published, and the
// blend from the Galilean moon tables into the mean elements.
void test_jupiter_missions()
{
    Scene scene = load_jupiter();
    const int galileo = scene.find("Galileo");
    const int juno = scene.find("Juno");
    const int io = scene.find("Io");
    const int europa = scene.find("Europa");
    const int ganymede = scene.find("Ganymede");
    check(galileo > 0 && juno > 0 && io > 0 && europa > 0 && ganymede > 0, "Jupiter mission bodies");
    auto at = [&](int body, double t) { return scene.icrf_state_at(body, t).position; };

    // Horizons, @ Jupiter center (500@599), 2021-Jun-07 16:57:18 TDB, ICRF km: Juno (-61)
    // and Ganymede (503) at the closest approach, between two knots of the table.
    const double t_g = tdb_from_jd_tdb(2459373.206458333);
    const glm::dvec3 juno_hzn(1.222322619851438E+05, -9.572422614074555E+05, -4.555370342133188E+05);
    const glm::dvec3 ganymede_hzn(1.197388830404732E+05, -9.586027392266861E+05, -4.578583386699203E+05);
    check(glm::length(at(juno, t_g) - juno_hzn) < 5.0, "Juno matches Horizons at the Ganymede flyby (km)",
          glm::length(at(juno, t_g) - juno_hzn));
    check(glm::length(at(ganymede, t_g) - ganymede_hzn) < 10.0, "Ganymede matches Horizons (km)",
          glm::length(at(ganymede, t_g) - ganymede_hzn));

    // Closest-approach altitudes (mean radii [SPHY]) against NASA/JPL news releases:
    // "Ride With Juno As It Flies Past the Solar System's Biggest Moon and Jupiter"
    // (2021-07-14): 1,038 km from Ganymede on 2021-06-07; "NASA's Juno Shares First
    // Image From Flyby of Jupiter's Moon Europa" (2022-09-29): Juno 352 km from Europa,
    // and Galileo 351 km on 2000-01-03. The releases give Earth-received times (here
    // ~39 min later than TDB at the spacecraft), so only the altitudes are compared.
    auto altitude = [&](int craft, int moon, const CalendarDateTime& tdb) {
        const double t0 = (jd_from_calendar(tdb) - kJ2000Jd) * kSecondsPerDay;
        double best = 1e300;
        for (double dt = -900.0; dt <= 900.0; dt += 1.0) {
            best = std::min(best, glm::length(at(craft, t0 + dt) - at(moon, t0 + dt)));
        }
        return best - scene.bodies[static_cast<size_t>(moon)].equatorial_radius_km;
    };
    const double alt_g = altitude(juno, ganymede, {2021, 6, 7, 16, 57, 17});
    const double alt_e = altitude(juno, europa, {2022, 9, 29, 9, 37, 37});
    const double alt_ge = altitude(galileo, europa, {2000, 1, 3, 18, 0, 46});
    check(std::abs(alt_g - 1038.0) < 15.0, "Juno-Ganymede flyby altitude (km)", alt_g);
    check(std::abs(alt_e - 352.0) < 15.0, "Juno-Europa flyby altitude (km)", alt_e);
    check(std::abs(alt_ge - 351.0) < 15.0, "Galileo-Europa flyby altitude (km)", alt_ge);

    // Galileo is there from 1995-11 until it enters Jupiter (2003-09-21), Juno from 2016-06.
    double t = 0.0;
    parse_utc("2003-09-22", &t);
    scene.update(t);
    check(!scene.bodies[static_cast<size_t>(galileo)].visible, "Galileo gone after its impact");
    parse_utc("2026-10-02", &t);
    scene.update(t);
    check(scene.bodies[static_cast<size_t>(juno)].visible, "Juno at Jupiter in 2026");

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

// Artemis I/II and CAPSTONE against NASA's published figures (distances only:
// mission blogs give times in local time zones, at the spacecraft or on Earth).
void test_earth_moon_scene()
{
    Scene scene = load_scene_or_die("earth_moon.toml");
    const int moon = scene.find("Moon");
    const int a1 = scene.find("Artemis I");
    const int a2 = scene.find("Artemis II");
    const int capstone = scene.find("CAPSTONE");
    check(moon > 0 && a1 > 0 && a2 > 0 && capstone > 0, "earth_moon.toml bodies");
    auto at = [&](int body, double t) { return scene.icrf_state_at(body, t).position; };
    auto tdb = [](const CalendarDateTime& c) { return (jd_from_calendar(c) - kJ2000Jd) * kSecondsPerDay; };
    const double moon_r = scene.bodies[static_cast<size_t>(moon)].equatorial_radius_km;
    const double earth_r = scene.bodies[0].equatorial_radius_km;

    // Extremum of f over [t0 - span, t0 + span] at one-second steps.
    auto extremum = [](auto f, double t0, double span, bool largest) {
        double best = largest ? -1e300 : 1e300;
        for (double dt = -span; dt <= span; dt += 1.0) {
            const double v = f(t0 + dt);
            best = largest ? std::max(best, v) : std::min(best, v);
        }
        return best;
    };
    auto above_moon = [&](int craft) { return [&, craft](double t) { return glm::length(at(craft, t) - at(moon, t)) - moon_r; }; };
    auto above_earth = [&](int craft) { return [&, craft](double t) { return glm::length(at(craft, t)) - earth_r; }; };
    constexpr double kMile = 1.609344;

    // "Artemis I - Flight Day Six: Orion Performs Lunar Flyby, Closest Outbound
    // Approach" (NASA blog, 2022-11-21): about 81 miles above the surface.
    const double a1_flyby = extremum(above_moon(a1), tdb({2022, 11, 21, 12, 58, 49}), 900.0, false);
    check(std::abs(a1_flyby - 81.0 * kMile) < 5.0, "Artemis I outbound flyby altitude (km)", a1_flyby);
    // "Artemis I - Flight Day 13: Orion Goes the (Max) Distance" (2022-11-28): 268,563 miles.
    const double a1_max = extremum(above_earth(a1), tdb({2022, 11, 28, 21, 7, 35}), 900.0, true);
    check(std::abs(a1_max - 268563.0 * kMile) < 15.0, "Artemis I farthest from the Earth (km)", a1_max);
    // "NASA's Artemis II Crew Eclipses Record for Farthest Human Spaceflight"
    // (2026-04-06): about 4,067 miles from the lunar surface, 252,756 miles from the
    // Earth at the farthest point. Both distances are from the surfaces (from the
    // Earth's center the maximum would be 6,378 km larger).
    const double a2_flyby = extremum(above_moon(a2), tdb({2026, 4, 6, 23, 1, 49}), 900.0, false);
    check(std::abs(a2_flyby - 4067.0 * kMile) < 15.0, "Artemis II lunar flyby altitude (km)", a2_flyby);
    const double a2_max = extremum(above_earth(a2), tdb({2026, 4, 6, 23, 6, 12}), 900.0, true);
    check(std::abs(a2_max - 252756.0 * kMile) < 15.0, "Artemis II farthest from the Earth (km)", a2_max);

    // CAPSTONE's near-rectilinear halo orbit is in 9:2 resonance with the synodic
    // month (29.530589 d, mean; Meeus, Astronomical Algorithms, ch. 49): the mean
    // time between perilunes over 2023-2025 is 2/9 of it. Perilunes lie over the
    // north pole (southern L2 family).
    const double t0 = tdb({2023, 1, 1, 0, 0, 0});
    const double t1 = tdb({2026, 1, 1, 0, 0, 0});
    const std::vector<double> peri = scene.periapsis_times(capstone, t0, t1);
    const double mean_days = peri.size() > 1 ? (peri.back() - peri.front()) / (peri.size() - 1) / kSecondsPerDay : 0.0;
    check(std::abs(mean_days - 2.0 / 9.0 * 29.530589) < 0.02, "CAPSTONE perilune interval (days)", mean_days);
    double lowest_lat = 90.0;
    for (const double tp : peri) {
        const glm::dvec3 r = scene.bodies[static_cast<size_t>(capstone)].motion->eval(tp).position;
        lowest_lat = std::min(lowest_lat, std::asin(glm::dot(r, scene.bodies[static_cast<size_t>(moon)].pole) / glm::length(r)) * kRadToDeg);
    }
    check(peri.size() > 150 && lowest_lat > 70.0, "CAPSTONE perilunes over the north pole (deg)", lowest_lat);

    // After splashdown the spacecraft is gone, but its trail lingers (fading) for
    // trail_linger_days and then disappears.
    std::vector<glm::dvec3> points;
    std::vector<float> fades;
    scene.update(tdb({2026, 4, 20, 0, 0, 0}));
    scene.trail(a2, tdb({2026, 4, 20, 0, 0, 0}), 2048, points, fades);
    const float min_fade = fades.empty() ? 0.0f : *std::min_element(fades.begin(), fades.end());
    check(!scene.bodies[static_cast<size_t>(a2)].visible && points.size() > 100 && min_fade > 0.25f,
          "Artemis II trail lingers, fading, after splashdown", static_cast<double>(points.size()));
    scene.update(tdb({2026, 5, 20, 0, 0, 0}));
    scene.trail(a2, tdb({2026, 5, 20, 0, 0, 0}), 2048, points, fades);
    check(points.empty(), "Artemis II trail gone 40 days after splashdown", static_cast<double>(points.size()));

    // Moon-centered rotating frame: the Moon at the origin, the Earth on -x.
    for (size_t f = 0; f < scene.frames.size(); ++f) {
        if (scene.frames[f].name == "Earth-Moon rotating, Moon-centered") {
            scene.set_active_frame(static_cast<int>(f));
        }
    }
    scene.update(t0);
    const glm::dvec3 earth = scene.bodies[0].world_position;
    check(glm::length(scene.bodies[static_cast<size_t>(moon)].world_position) < 1e-6 && earth.x < -3.5e5 &&
              std::abs(earth.y) < 1e-3 && std::abs(earth.z) < 1e-3,
          "Moon-centered rotating frame", earth.x);
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

    // Flyby altitudes (mean radii) against NASA: "Titan Flyby (T-70) - June 21, 2010"
    // (880 km), "Iapetus Flyby - Sept. 10, 2007" (about 1,640 km), "Enceladus Flyby -
    // October 09, 2008" (25 km), NASA Science Cassini mission pages.
    const int cassini = scene.find("Cassini");
    auto altitude = [&](int craft, const char* moon_name, const CalendarDateTime& c) {
        const int moon = scene.find(moon_name);
        const double t0 = tdb(c);
        double best = 1e300;
        for (double dt = -900.0; dt <= 900.0; dt += 1.0) {
            best = std::min(best, glm::length(at(craft, t0 + dt) - at(moon, t0 + dt)));
        }
        return best - scene.bodies[static_cast<size_t>(moon)].equatorial_radius_km;
    };
    const double t70 = altitude(cassini, "Titan", {2010, 6, 21, 1, 28, 48});
    const double iapetus = altitude(cassini, "Iapetus", {2007, 9, 10, 14, 16, 55});
    const double e5 = altitude(cassini, "Enceladus", {2008, 10, 9, 19, 7, 44});
    check(std::abs(t70 - 880.0) < 10.0, "Cassini-Titan T70 altitude (km)", t70);
    check(std::abs(iapetus - 1640.0) < 15.0, "Cassini-Iapetus altitude (km)", iapetus);
    check(std::abs(e5 - 25.0) < 5.0, "Cassini-Enceladus E5 altitude (km)", e5);

    // Huygens ends on Titan's surface.
    const int huygens = scene.find("Huygens");
    const double landed = altitude(huygens, "Titan", {2005, 1, 14, 11, 38, 49});
    check(std::abs(landed) < 5.0, "Huygens on Titan's surface (km)", landed);

    // Grand Finale: 22 dives between the rings and the planet (periapsis inside the
    // D ring's inner edge, ~66,900 km), then the plunge.
    int dives = 0;
    for (const double tp : scene.periapsis_times(cassini, tdb({2017, 4, 20, 0, 0, 0}), tdb({2017, 9, 15, 10, 0, 0}))) {
        dives += glm::length(at(cassini, tp)) < 66000.0 ? 1 : 0;
    }
    check(dives == 22, "Cassini's Grand Finale dives", dives);

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

// MESSENGER and BepiColombo: heliocentric cruise tables handed over to tables
// relative to Mercury (relative_to), checked against the Horizons notes of the
// two spacecraft and against NASA's MESSENGER Earth-flyby altitude.
void test_mercury_scene()
{
    Scene scene = load_scene_or_die("mercury.toml");
    const int mercury = scene.find("Mercury");
    const int venus = scene.find("Venus");
    const int earth = scene.find("Earth");
    const int messenger = scene.find("MESSENGER");
    const int bepi = scene.find("BepiColombo");
    check(mercury > 0 && venus > 0 && earth > 0 && messenger > 0 && bepi > 0, "mercury.toml bodies");
    auto at = [&](int body, double t) { return scene.icrf_state_at(body, t).position; };
    auto tdb = [](const CalendarDateTime& c) { return (jd_from_calendar(c) - kJ2000Jd) * kSecondsPerDay; };

    // Horizons, heliocentric (500@10), 2008-Jan-14 18:00 TDB, ICRF km: Mercury (199), Venus (299).
    const double t_ref = tdb_from_jd_tdb(2454480.25);
    const glm::dvec3 mercury_hzn(5.237913070701520E+07, 4.158221025197503E+06, -3.210073214907367E+06);
    const glm::dvec3 venus_hzn(-9.743411097601077E+07, -4.448008071855614E+07, -1.384639758080466E+07);
    check(glm::length(at(mercury, t_ref) - mercury_hzn) < 2.0, "Mercury matches Horizons (km)",
          glm::length(at(mercury, t_ref) - mercury_hzn));
    check(glm::length(at(venus, t_ref) - venus_hzn) < 5.0, "Venus matches Horizons (km)",
          glm::length(at(venus, t_ref) - venus_hzn));

    // Closest approach (+-1 h around the bake tool's time) minus a reference radius.
    auto altitude = [&](int craft, int planet, const CalendarDateTime& c, double radius) {
        const double t0 = tdb(c);
        double best = 1e300;
        for (double dt = -3600.0; dt <= 3600.0; dt += 1.0) {
            best = std::min(best, glm::length(at(craft, t0 + dt) - at(planet, t0 + dt)));
        }
        return best - radius;
    };
    const double r_mercury = scene.bodies[static_cast<size_t>(mercury)].equatorial_radius_km;
    const double r_earth = scene.bodies[static_cast<size_t>(earth)].equatorial_radius_km;

    // Horizons notes for MESSENGER (-236): first Mercury flyby 2008-01-14 19:05 UTC at
    // 200.6 km altitude, third 2009-09-29 21:56 UTC at 228 km; "2,347 kilometers over
    // central Mongolia" at the Earth flyby (2005-08-02). The Earth flyby needs the
    // Earth's center (scene: Earth around the Earth-Moon barycenter), not the barycenter.
    const double m1 = altitude(messenger, mercury, {2008, 1, 14, 19, 5, 44}, r_mercury);
    const double m3 = altitude(messenger, mercury, {2009, 9, 29, 21, 56, 0}, r_mercury);
    const double me = altitude(messenger, earth, {2005, 8, 2, 19, 4, 17}, r_earth);
    check(std::abs(m1 - 200.6) < 5.0, "MESSENGER first Mercury flyby altitude (km)", m1);
    check(std::abs(m3 - 228.0) < 5.0, "MESSENGER third Mercury flyby altitude (km)", m3);
    check(std::abs(me - 2347.0) < 30.0, "MESSENGER Earth flyby altitude (km)", me);
    std::printf("info: MESSENGER flybys: Mercury 1 %.1f km, Mercury 3 %.1f km, Earth %.1f km\n", m1, m3, me);

    // Horizons notes for BepiColombo (-121): Earth flyby 2020-04-10 04:24:58 UTC "@ 19066 km"
    // (from the Earth's center; ESA quoted ~12,700 km above the surface).
    const double be = altitude(bepi, earth, {2020, 4, 10, 4, 26, 7}, 0.0);
    check(std::abs(be - 19066.0) < 30.0, "BepiColombo Earth flyby distance (km)", be);
    std::printf("info: BepiColombo Earth flyby: %.0f km from the Earth's center\n", be);

    // Hand-over from the heliocentric cruise table to the Mercury-relative one: no jump.
    for (const auto& [craft, when] : {std::pair{messenger, CalendarDateTime{2011, 3, 12, 0, 0, 0}},
                                     std::pair{bepi, CalendarDateTime{2026, 10, 20, 0, 0, 0}}}) {
        const double t_switch = tdb(when);
        double worst = 0.0;
        for (double t = t_switch - 600.0; t < t_switch + 600.0; t += 60.0) {
            const State a = scene.icrf_state_at(craft, t);
            worst = std::max(worst, glm::length(at(craft, t + 60.0) - a.position - a.velocity * 60.0));
        }
        check(worst < 20.0, "spacecraft hand-over to the Mercury-relative table (km/min)", worst);
    }

    // MESSENGER ends at Mercury's surface (impact 2015-04-30 ~19:26 UTC).
    double t_end = tdb({2015, 4, 30, 19, 27, 0});
    scene.update(t_end);
    const double impact = glm::length(at(messenger, t_end) - at(mercury, t_end)) - r_mercury;
    check(scene.bodies[static_cast<size_t>(messenger)].visible && std::abs(impact) < 30.0, "MESSENGER impact altitude (km)", impact);
    scene.update(t_end + 600.0);
    check(!scene.bodies[static_cast<size_t>(messenger)].visible, "MESSENGER gone after the impact");
    scene.update(tdb({2004, 1, 1, 0, 0, 0}));
    check(!scene.bodies[static_cast<size_t>(messenger)].visible, "MESSENGER not there before its launch");

    // BepiColombo's science orbit (Horizons notes: MPO 480 x 1500 km, 2.3 h).
    double lowest = 1e300;
    double highest = 0.0;
    const double t_mpo = tdb({2027, 4, 9, 0, 0, 0});
    for (double t = t_mpo; t < t_mpo + kSecondsPerDay; t += 30.0) {
        const double h = glm::length(at(bepi, t) - at(mercury, t)) - r_mercury;
        lowest = std::min(lowest, h);
        highest = std::max(highest, h);
    }
    check(std::abs(lowest - 480.0) < 50.0 && std::abs(highest - 1500.0) < 100.0, "BepiColombo science orbit (km)", lowest);
    std::printf("info: BepiColombo in April 2027: %.0f x %.0f km\n", lowest, highest);

    // Trails: in the Sun-centered frame the 2013 trail reaches back into the cruise;
    // in the Mercury-centered frame it is capped at one day.
    const double t_orbit = tdb({2013, 1, 1, 0, 0, 0});
    std::vector<glm::dvec3> points;
    std::vector<float> fades;
    scene.update(t_orbit);
    scene.trail(messenger, t_orbit, 2048, points, fades);
    double farthest = 0.0;
    for (const glm::dvec3& q : points) {
        farthest = std::max(farthest, glm::length(q - scene.bodies[static_cast<size_t>(mercury)].world_position));
    }
    check(farthest > 1.0e7, "MESSENGER trail reaches back into the cruise (km from Mercury)", farthest);
    for (size_t f = 0; f < scene.frames.size(); ++f) {
        if (scene.frames[f].name == "Mercury-centered") {
            scene.set_active_frame(static_cast<int>(f));
        }
    }
    scene.update(t_orbit);
    scene.trail(messenger, t_orbit, 2048, points, fades);
    double widest = 0.0;
    for (const glm::dvec3& q : points) {
        widest = std::max(widest, glm::length(q));
    }
    check(points.size() > 50 && widest < 2.0e4, "MESSENGER trail in the Mercury frame: a day of orbits (km)", widest);
}

// ISEE-3 / ICE, reconstructed by tools/isee3/isee3_reconstruct.py from SSCWeb positions
// (1978-1983), the JPL navigation trajectory of the comet encounter and Horizons' 2014
// solution. Checked against figures published independently of those data.
void test_isee3_scene()
{
    Scene scene = load_scene_or_die("isee3.toml");
    const int sun = scene.find("Sun");
    const int earth = scene.find("Earth");
    const int moon = scene.find("Moon");
    const int craft = scene.find("ISEE-3");
    const int gz = scene.find("Giacobini-Zinner");
    const int halley = scene.find("Halley");
    check(sun >= 0 && earth > 0 && moon > 0 && craft > 0 && gz > 0 && halley > 0, "isee3.toml bodies");
    auto at = [&](int body, double t) { return scene.icrf_state_at(body, t).position; };
    auto utc = [](const CalendarDateTime& c) { return tdb_from_utc_jd(jd_from_calendar(c)); };
    auto closest = [&](int a, int b, double t0, double half, double step) {
        std::pair<double, double> best{1e300, t0};
        for (double t = t0 - half; t <= t0 + half; t += step) {
            const double d = glm::length(at(a, t) - at(b, t));
            if (d < best.first) {
                best = {d, t};
            }
        }
        return best;
    };

    // Launch 1978-08-12 15:12 UTC: the reconstruction begins at perigee an hour later.
    scene.update(utc({1978, 8, 12, 15, 0, 0}));
    check(!scene.bodies[static_cast<size_t>(craft)].visible, "ISEE-3 not there before the launch");
    const auto [perigee, t_perigee] = closest(craft, earth, utc({1978, 8, 12, 17, 0, 0}), 3600.0, 10.0);
    const double r_earth = scene.bodies[static_cast<size_t>(earth)].equatorial_radius_km;
    check(perigee - r_earth < 1000.0 && t_perigee - utc({1978, 8, 12, 15, 12, 0}) < 2 * 3600.0,
          "the transfer starts at a low perigee right after the launch (km)", perigee - r_earth);

    // Halo orbit "about 240 Earth radii upstream" (ESA/eoportal; L1 is 1.49 million km out).
    double sunward_sum = 0.0;
    int n = 0;
    for (double t = utc({1979, 1, 1, 0, 0, 0}); t < utc({1982, 1, 1, 0, 0, 0}); t += 5 * kSecondsPerDay, ++n) {
        const glm::dvec3 to_sun = glm::normalize(at(sun, t) - at(earth, t));
        sunward_sum += glm::dot(at(craft, t) - at(earth, t), to_sun);
    }
    const double sunward = sunward_sum / n / 6378.1366;
    check(std::abs(sunward - 240.0) < 20.0, "halo orbit around L1: mean sunward distance (Earth radii)", sunward);

    // Geotail passes: "distances as great as 237 Earth radii in the tail" (ICE summary in the
    // PDS Halley archive, DOCUMENT/ESA_1066/ICE.TXT).
    double tail = 0.0;
    for (double t = utc({1982, 6, 10, 0, 0, 0}); t < utc({1983, 12, 20, 0, 0, 0}); t += 3600.0) {
        const glm::dvec3 to_sun = glm::normalize(at(sun, t) - at(earth, t));
        tail = std::max(tail, -glm::dot(at(craft, t) - at(earth, t), to_sun));
    }
    tail /= 6378.1366;
    check(std::abs(tail - 237.0) < 5.0, "deepest geotail pass (Earth radii)", tail);

    // Lunar flybys: 1983-10-21 at 17,440 km (Wikipedia, "1983 in spaceflight"); the last one
    // 119.4 km above the surface (Horizons notes for -111; NASA/ESA mission summaries).
    const double r_moon = scene.bodies[static_cast<size_t>(moon)].equatorial_radius_km;
    const double flyby4 = closest(craft, moon, utc({1983, 10, 21, 16, 25, 0}), 3600.0, 2.0).first - r_moon;
    const auto [d5, t5] = closest(craft, moon, utc({1983, 12, 22, 18, 44, 0}), 3600.0, 1.0);
    check(std::abs(flyby4 - 17440.0) < 100.0, "lunar flyby 1983-10-21 (km)", flyby4);
    check(std::abs(d5 - r_moon - 119.4) < 5.0, "last lunar flyby 1983-12-22 (km above the surface)", d5 - r_moon);
    // "spent more than 30 minutes in its shadow" (TDA Progress Report 42-84, p. 178).
    double shadow_s = 0.0;
    for (double t = t5 - 7200.0; t < t5 + 7200.0; t += 10.0) {
        const glm::dvec3 p = at(craft, t);
        const glm::dvec3 to_moon = at(moon, t) - p;
        const glm::dvec3 to_sun = at(sun, t) - p;
        const double angle = std::acos(glm::dot(glm::normalize(to_moon), glm::normalize(to_sun)));
        if (angle < std::asin(r_moon / glm::length(to_moon)) + std::asin(695700.0 / glm::length(to_sun))) {
            shadow_s += 10.0;
        }
    }
    check(shadow_s > 25.0 * 60.0, "minutes in the Moon's shadow at the last flyby", shadow_s / 60.0);
    std::printf("info: ISEE-3 last lunar flyby %.1f km at %s, %.0f min in the Moon's shadow\n", d5 - r_moon,
                format_utc(t5).c_str(), shadow_s / 60.0);

    // Heliocentric table from 1983-12-31: no jump at the hand-over.
    {
        const double t_switch = utc({1983, 12, 31, 0, 0, 0});
        double worst = 0.0;
        for (double t = t_switch - 2 * kSecondsPerDay; t < t_switch + 2 * kSecondsPerDay; t += 600.0) {
            const State a = scene.icrf_state_at(craft, t);
            worst = std::max(worst, glm::length(at(craft, t + 600.0) - a.position - a.velocity * 600.0));
        }
        check(worst < 20.0, "ISEE-3 hand-over to the heliocentric table (km per 10 min)", worst);
    }

    // Giacobini-Zinner, 1985-09-11 11:02 UTC: 7,862 km from the nucleus at 20.7 km/s
    // (Horizons notes for -111; NASA).
    const auto [d_gz, t_gz] = closest(craft, gz, utc({1985, 9, 11, 11, 2, 0}), 600.0, 0.5);
    const double v_gz = glm::length(scene.icrf_state_at(craft, t_gz).velocity - scene.icrf_state_at(gz, t_gz).velocity);
    check(std::abs(d_gz - 7862.0) < 100.0, "Giacobini-Zinner closest approach (km)", d_gz);
    check(std::abs(t_gz - utc({1985, 9, 11, 11, 2, 0})) < 120.0, "Giacobini-Zinner closest approach time (s)",
          t_gz - utc({1985, 9, 11, 11, 2, 0}));
    check(std::abs(v_gz - 20.7) < 0.1, "Giacobini-Zinner flyby speed (km/s)", v_gz);

    // Halley: "flew between the Sun and comet Halley on 1986-Mar-28, 31 million km from the
    // comet" (Horizons notes for -111). Find when ICE is closest to the Sun-Halley line.
    double best_angle = 0.0;
    double t_align = 0.0;
    for (double t = utc({1986, 3, 1, 0, 0, 0}); t < utc({1986, 4, 30, 0, 0, 0}); t += 3600.0) {
        const glm::dvec3 p = at(craft, t);
        const double angle = std::acos(glm::dot(glm::normalize(at(sun, t) - p), glm::normalize(at(halley, t) - p)));
        if (angle > best_angle) {
            best_angle = angle;
            t_align = t;
        }
    }
    const double d_halley = glm::length(at(halley, t_align) - at(craft, t_align));
    check(std::abs(t_align - utc({1986, 3, 28, 12, 0, 0})) < 1.5 * kSecondsPerDay && std::abs(d_halley - 3.1e7) < 1.0e6,
          "between the Sun and Halley: distance (km)", d_halley);
    std::printf("info: ICE closest to the Sun-Halley line %s, %.1f million km from Halley\n",
                format_utc(t_align).c_str(), d_halley / 1e6);
    // Closest to Halley "28 x 10^6 km" (Wikipedia, International Cometary Explorer).
    const double halley_min =
        closest(craft, halley, utc({1986, 3, 25, 0, 0, 0}), 10 * kSecondsPerDay, 3600.0).first;
    check(std::abs(halley_min - 2.8e7) < 1.0e6, "closest to Halley (km)", halley_min);

    // 1987-2013, between the 1986-1987 maneuvers and the 2014 return: "an aphelion of 1.03 AU,
    // a perihelion of 0.93 AU and an inclination of 0.1 deg" (Wikipedia, International
    // Cometary Explorer) and a "355-day orbit" (Dunham, Farquhar et al., IAC-14.B6.3.4).
    {
        const double t95 = utc({1995, 1, 1, 0, 0, 0});
        const State s = scene.icrf_state_at(craft, t95);
        const double gm = scene.bodies[static_cast<size_t>(sun)].gm_km3_s2;
        const double r = glm::length(s.position);
        const double a = 1.0 / (2.0 / r - glm::dot(s.velocity, s.velocity) / gm);
        const glm::dvec3 h = glm::cross(s.position, s.velocity);
        const double e = glm::length(glm::cross(s.velocity, h) / gm - s.position / r);
        const double inc = std::acos(glm::dot(glm::normalize(h), ecliptic_pole_icrf())) * kRadToDeg;
        const double period = 2.0 * kPi * std::sqrt(a * a * a / gm) / kSecondsPerDay;
        check(std::abs(a * (1 - e) / kAuKm - 0.93) < 0.01 && std::abs(a * (1 + e) / kAuKm - 1.03) < 0.01,
              "ICE perihelion and aphelion in 1995 (au)", a * (1 - e) / kAuKm);
        check(inc < 0.3 && std::abs(period - 355.0) < 2.0, "ICE orbital period in 1995 (days)", period);
        std::printf("info: ICE in 1995: %.3f x %.3f au, i = %.2f deg, %.1f days\n", a * (1 - e) / kAuKm,
                    a * (1 + e) / kAuKm, inc, period);
    }

    // 2014 (Horizons -111, solution s45): Earth closest approach 2014-08-09.12355 TDB at
    // 0.001238 au; lunar flyby 2014-08-10 19:18:27 UTC, 18,401 km from the Moon's center
    // (Horizons notes).
    const auto [d_earth14, t_earth14] = closest(craft, earth, tdb_from_jd_tdb(2456878.62355), 3600.0, 10.0);
    check(std::abs(d_earth14 - 0.001238 * kAuKm) < 300.0, "Earth closest approach 2014 (km)", d_earth14);
    const auto [d_moon14, t_moon14] = closest(craft, moon, utc({2014, 8, 10, 19, 18, 27}), 3600.0, 2.0);
    check(std::abs(d_moon14 - 18401.0) < 50.0 && std::abs(t_moon14 - utc({2014, 8, 10, 19, 18, 27})) < 180.0,
          "lunar flyby 2014-08-10 (km from the center)", d_moon14);
    // Hand-over from the reconstruction to the Horizons table (2014-01-01): no jump.
    {
        const double t_switch = utc({2014, 1, 1, 0, 0, 0});
        double worst = 0.0;
        for (double t = t_switch - 2 * kSecondsPerDay; t < t_switch + 2 * kSecondsPerDay; t += 600.0) {
            const State a = scene.icrf_state_at(craft, t);
            worst = std::max(worst, glm::length(at(craft, t + 600.0) - a.position - a.velocity * 600.0));
        }
        check(worst < 20.0, "ISEE-3 hand-over to the 2014 table (km per 10 min)", worst);
    }
    // The last day of the 2014 table is not blended towards its fallback, which ended in
    // January (that would pull the end of the trail away). Horizons, 2014-12-31 18:00 TDB.
    const double t_last = tdb_from_jd_tdb(2457023.25);
    const glm::dvec3 last_hzn(2.452466292458899E+07, 1.580746653259222E+08, 7.169687407992832E+07);
    check(glm::length(at(craft, t_last) - at(sun, t_last) - last_hzn) < 5.0, "ISEE-3 at the end of its table (km)",
          glm::length(at(craft, t_last) - at(sun, t_last) - last_hzn));
    scene.update(utc({2015, 2, 1, 0, 0, 0}));
    check(!scene.bodies[static_cast<size_t>(craft)].visible, "ISEE-3 gone after the 2014 table ends");
}

// The Halley Armada in the ISEE-3 scene, reconstructed by tools/armada/armada_reconstruct.py
// from mission archives. The encounters of Giotto and Vega 1 are moved onto the published
// ones; Suisei's, Sakigake's and the Venus flyby are not, and are checked against figures
// published independently of the data.
void test_halley_armada()
{
    Scene scene = load_scene_or_die("isee3.toml");
    const int sun = scene.find("Sun");
    const int earth = scene.find("Earth");
    const int venus = scene.find("Venus");
    const int halley = scene.find("Halley");
    const int vega1 = scene.find("Vega 1");
    const int suisei = scene.find("Suisei");
    const int sakigake = scene.find("Sakigake");
    const int giotto = scene.find("Giotto");
    check(vega1 > 0 && suisei > 0 && sakigake > 0 && giotto > 0, "isee3.toml: the Armada");
    auto at = [&](int body, double t) { return scene.icrf_state_at(body, t).position; };
    auto utc = [](const CalendarDateTime& c) { return tdb_from_utc_jd(jd_from_calendar(c)); };
    auto closest = [&](int a, int b, double t0, double half, double step) {
        std::pair<double, double> best{1e300, t0};
        for (int pass = 0; pass < 3; ++pass, half = 2 * step, step /= 100.0) {
            const double center = best.first < 1e300 ? best.second : t0;
            for (double t = center - half; t <= center + half; t += step) {
                const double d = glm::length(at(a, t) - at(b, t));
                if (d < best.first) {
                    best = {d, t};
                }
            }
        }
        return best;
    };
    auto speed = [&](int a, int b, double t) {
        return glm::length(scene.icrf_state_at(a, t).velocity - scene.icrf_state_at(b, t).velocity);
    };
    auto visible = [&](int body, double t) {
        scene.update(t);
        return scene.bodies[static_cast<size_t>(body)].visible;
    };

    // Vega 1: launched 1984-12-15 [VEGA in isee3.toml]; its table starts near the Earth.
    check(!visible(vega1, utc({1984, 12, 15, 0, 0, 0})) && visible(vega1, utc({1984, 12, 15, 12, 0, 0})),
          "Vega 1 appears on its launch day");
    const auto [d_launch, t_launch] = closest(vega1, earth, utc({1984, 12, 15, 11, 0, 0}), 3600.0, 60.0);
    check(d_launch < 7000.0, "Vega 1 starts at the Earth (km from the center)", d_launch);
    // Venus: the lander arrived on 1985-06-11; the bus flew past at 39,000 km (Wikipedia,
    // Vega 1; NSSDC). The bridge across the data gap is not aimed at it.
    const auto [d_venus, t_venus] = closest(vega1, venus, utc({1985, 6, 11, 3, 0, 0}), 86400.0, 600.0);
    check(d_venus > 3.0e4 && d_venus < 5.0e4 && std::abs(t_venus - utc({1985, 6, 11, 3, 0, 0})) < 6 * 3600.0,
          "Vega 1 Venus flyby near the published 39,000 km (km)", d_venus);
    // Halley 1986-03-06 07:20:06 UT, 8,890 km, 79.2 km/s [VEGA].
    const auto [d_v1, t_v1] = closest(vega1, halley, utc({1986, 3, 6, 7, 20, 6}), 600.0, 10.0);
    check(std::abs(d_v1 - 8890.0) < 5.0, "Vega 1 closest approach (km)", d_v1);
    check(std::abs(t_v1 - utc({1986, 3, 6, 7, 20, 6})) < 1.0, "Vega 1 closest approach time (s)",
          t_v1 - utc({1986, 3, 6, 7, 20, 6}));
    check(std::abs(speed(vega1, halley, t_v1) - 79.2) < 0.1, "Vega 1 flyby speed (km/s)", speed(vega1, halley, t_v1));
    check(!visible(vega1, utc({1986, 7, 15, 0, 0, 0})), "Vega 1 gone after its positions end");

    // Suisei: 151,000 km "on the side facing the Sun" [JAXA], from the ISAS ephemeris unchanged.
    const auto [d_su, t_su] = closest(suisei, halley, utc({1986, 3, 8, 13, 6, 0}), 3 * 3600.0, 30.0);
    const glm::dvec3 to_sun = glm::normalize(at(sun, t_su) - at(halley, t_su));
    const double sunward = glm::dot(glm::normalize(at(suisei, t_su) - at(halley, t_su)), to_sun);
    check(std::abs(d_su / 151000.0 - 1.0) < 0.01, "Suisei closest approach (km)", d_su);
    check(sunward > 0.7, "Suisei passes on the sunward side (cosine)", sunward);
    check(visible(suisei, utc({1985, 8, 19, 6, 0, 0})) && !visible(suisei, utc({1985, 8, 18, 23, 0, 0})),
          "Suisei appears on its launch day");

    // Sakigake: "6.99 million km" on 1986-03-11 (Wikipedia); HelioWeb's rounded positions
    // give a few percent less.
    const auto [d_sk, t_sk] = closest(sakigake, halley, utc({1986, 3, 11, 0, 0, 0}), 2 * 86400.0, 600.0);
    check(std::abs(d_sk / 6.99e6 - 1.0) < 0.06 && std::abs(t_sk - utc({1986, 3, 11, 4, 0, 0})) < 12 * 3600.0,
          "Sakigake closest approach (km)", d_sk);

    // Giotto: 1986-03-14 00:03:01.84 UTC, 596 km [GIO]; "68 km/s" (ESA SP-1066). Only the
    // twelve days of ESOC's encounter file.
    const auto [d_g, t_g] = closest(giotto, halley, utc({1986, 3, 14, 0, 3, 2}), 600.0, 10.0);
    const double t_g_pub = utc({1986, 3, 14, 0, 3, 1}) + 0.84;
    check(std::abs(d_g - 596.0) < 2.0, "Giotto closest approach (km)", d_g);
    check(std::abs(t_g - t_g_pub) < 0.5, "Giotto closest approach time (s)", t_g - t_g_pub);
    check(std::abs(speed(giotto, halley, t_g) - 68.0) < 0.5, "Giotto flyby speed (km/s)", speed(giotto, halley, t_g));
    check(!visible(giotto, utc({1986, 3, 4, 0, 0, 0})) && visible(giotto, utc({1986, 3, 10, 0, 0, 0})) &&
              !visible(giotto, utc({1986, 3, 18, 0, 0, 0})),
          "Giotto shown only within ESOC's file");

    // Vega 2: no positions are published; reconstructed from its launch (1984-12-21), the
    // lander's entry at Venus (1985-06-15 02:06 UT) and the encounter: 1986-03-09 07:20:00 UT,
    // 8,030 km, 76.8 km/s [VEGA]. The speed is not imposed: it follows from the path.
    const int vega2 = scene.find("Vega 2");
    check(vega2 > 0, "isee3.toml: Vega 2");
    check(!visible(vega2, utc({1984, 12, 21, 9, 0, 0})) && visible(vega2, utc({1984, 12, 21, 18, 0, 0})),
          "Vega 2 appears on its launch day");
    const auto [d_v2, t_v2] = closest(vega2, halley, utc({1986, 3, 9, 7, 20, 0}), 600.0, 10.0);
    check(std::abs(d_v2 - 8030.0) < 5.0, "Vega 2 closest approach (km)", d_v2);
    check(std::abs(t_v2 - utc({1986, 3, 9, 7, 20, 0})) < 1.0, "Vega 2 closest approach time (s)",
          t_v2 - utc({1986, 3, 9, 7, 20, 0}));
    check(std::abs(speed(vega2, halley, t_v2) - 76.8) < 0.1, "Vega 2 flyby speed (km/s)", speed(vega2, halley, t_v2));
    const auto [d_venus2, t_venus2] = closest(vega2, venus, utc({1985, 6, 15, 2, 6, 0}), 86400.0, 600.0);
    check(d_venus2 > 7000.0 && d_venus2 < 6.0e4 && std::abs(t_venus2 - utc({1985, 6, 15, 2, 6, 0})) < 12 * 3600.0,
          "Vega 2 passes Venus at its lander's arrival (km)", d_venus2);
    check(!visible(vega2, utc({1986, 5, 15, 0, 0, 0})), "Vega 2 gone after the mission (1986-04)");

    // Pioneer Venus Orbiter around Venus: the periapses and an apoapsis of its VSO position
    // files (PDS PPI; distances = ALT + 6050 km): 1985-07-01 22:41:07 UTC 8,235.8 km,
    // 1986-04-27 23:25:26 UTC 8,300.3 km; 1985-07-02 10:41 UTC 70,682 km.
    const int pvo = scene.find("Pioneer Venus Orbiter");
    check(pvo > 0 && scene.bodies[static_cast<size_t>(pvo)].parent == venus, "isee3.toml: PVO around Venus");
    const auto [peri1, t_peri1] = closest(pvo, venus, utc({1985, 7, 1, 22, 41, 7}), 1800.0, 20.0);
    const auto [peri2, t_peri2] = closest(pvo, venus, utc({1986, 4, 27, 23, 25, 26}), 1800.0, 20.0);
    check(std::abs(peri1 - 8235.8) < 5.0 && std::abs(t_peri1 - utc({1985, 7, 1, 22, 41, 7})) < 30.0,
          "PVO periapsis 1985-07-01 (km)", peri1);
    check(std::abs(peri2 - 8300.3) < 5.0 && std::abs(t_peri2 - utc({1986, 4, 27, 23, 25, 26})) < 30.0,
          "PVO periapsis 1986-04-27 (km)", peri2);
    const double apo = glm::length(at(pvo, utc({1985, 7, 2, 10, 41, 6})) - at(venus, utc({1985, 7, 2, 10, 41, 6})));
    check(std::abs(apo - 70682.0) < 20.0, "PVO apoapsis 1985-07-02 (km)", apo);
    check(visible(pvo, utc({1986, 1, 20, 0, 0, 0})) && !visible(pvo, utc({1987, 2, 1, 0, 0, 0})),
          "PVO shown through 1986 (its data gaps bridged), not after");

    // Pioneer 7: integrated from its 1966 launch and published orbit, the period fitted to
    // NASA's time of closest approach, 23:36 UT on 1986-03-20; the distance is not imposed
    // (NASA: 12.1 million km). Downstream: on the comet's far side from the Sun [P7].
    const int pioneer7 = scene.find("Pioneer 7");
    check(pioneer7 > 0, "isee3.toml: Pioneer 7");
    const double t_p7_pub = utc({1986, 3, 20, 23, 36, 0});
    const auto [d_p7, t_p7] = closest(pioneer7, halley, t_p7_pub, 2 * 86400.0, 600.0);
    check(std::abs(t_p7 - t_p7_pub) < 120.0, "Pioneer 7 closest approach time (s)", t_p7 - t_p7_pub);
    check(std::abs(d_p7 / 12.1e6 - 1.0) < 0.02, "Pioneer 7 closest approach (km)", d_p7);
    const double p7_sunward = glm::dot(glm::normalize(at(pioneer7, t_p7) - at(halley, t_p7)),
                                       glm::normalize(at(sun, t_p7) - at(halley, t_p7)));
    check(p7_sunward < -0.3, "Pioneer 7 passes on the far side from the Sun (cosine)", p7_sunward);
    check(!visible(pioneer7, utc({1984, 12, 31, 0, 0, 0})) && visible(pioneer7, utc({1985, 1, 2, 0, 0, 0})) &&
              !visible(pioneer7, utc({1987, 1, 2, 0, 0, 0})),
          "Pioneer 7 shown in 1985-1986 only");
    std::printf("info: Pioneer 7 at Halley %.3f million km %s\n", d_p7 / 1e6, format_utc(t_p7).c_str());

    std::printf("info: Vega 2 at Halley %.0f km %s, at Venus %.0f km %s; PVO periapses %.1f, %.1f km\n", d_v2,
                format_utc(t_v2).c_str(), d_venus2, format_utc(t_venus2).c_str(), peri1, peri2);
    std::printf("info: Armada at Halley: Vega 1 %.0f km %s, Suisei %.0f km %s, Sakigake %.3g km %s, Giotto %.1f km "
                "%s; Vega 1 at Venus %.0f km %s\n",
                d_v1, format_utc(t_v1).c_str(), d_su, format_utc(t_su).c_str(), d_sk, format_utc(t_sk).c_str(), d_g,
                format_utc(t_g).c_str(), d_venus, format_utc(t_venus).c_str());
}

void test_jwst_scene()
{
    // Lagrange L2 of Sun / Earth-Moon: force balance in the rotating frame, ~1.5e6 km.
    const double mu = 3.04042340e-6;
    const double gamma = lagrange_gamma(mu, 2);
    const double x = 1.0 - mu + gamma;
    const double f = x - (1.0 - mu) / ((x + mu) * (x + mu)) - mu / ((x - 1.0 + mu) * (x - 1.0 + mu));
    check(std::abs(f) < 1e-12, "L2 force balance", f);
    check(gamma * kAuKm > 1.49e6 && gamma * kAuKm < 1.52e6, "L2 distance (km)", gamma * kAuKm);
    const double g1 = lagrange_gamma(mu, 1);
    const double x1 = 1.0 - mu - g1;
    const double f1 = x1 - (1.0 - mu) / ((x1 + mu) * (x1 + mu)) + mu / ((x1 - 1.0 + mu) * (x1 - 1.0 + mu));
    check(std::abs(f1) < 1e-12, "L1 force balance", f1);

    Scene scene = load_scene_or_die("jwst.toml");
    check(scene.frame_is_rotating(), "JWST scene starts in the rotating frame");
    const int jwst = scene.find("JWST");
    double worst = 0.0;
    for (int day = 200; day < 3500; day += 5) { // from mid-2022 into the predicted schedule
        double t = 0.0;
        parse_utc("2022-01-01", &t);
        t += day * kSecondsPerDay;
        scene.update(t);
        // The Sun stays on the -x axis of the rotating frame.
        const glm::dvec3 sun = glm::normalize(scene.sun_position());
        check(sun.x < -0.9999999, "Sun on -x in the rotating frame", sun.x);
        const glm::dvec3 l2 = scene.markers[0].world_position;
        check(l2.x > 1.4e6 && std::abs(l2.y) < 1.0 && std::abs(l2.z) < 1.0, "L2 marker on +x", l2.x);
        worst = std::max(worst, glm::length(scene.bodies[static_cast<size_t>(jwst)].world_position - l2));
    }
    // JWST's halo orbit amplitude is several hundred thousand km.
    check(worst < 1.2e6, "JWST stays near L2 (km)", worst);
}

void test_galactic_frame()
{
    // The Galactic centre direction (l = b = 0) in ICRS: RA 266.405, Dec -28.936 deg
    // (Hipparcos definition, ESA SP-1200 vol. 1, sec. 1.5.3).
    const glm::dvec3 gc = galactic_to_icrf() * glm::dvec3(1.0, 0.0, 0.0);
    double ra = std::atan2(gc.y, gc.x) * kRadToDeg;
    ra = ra < 0.0 ? ra + 360.0 : ra;
    const double dec = std::asin(gc.z) * kRadToDeg;
    check(std::abs(ra - 266.405) < 0.01 && std::abs(dec + 28.936) < 0.01, "Galactic centre direction", ra);
}

// Akeson et al. (2021) alpha Cen AB elements reproduce their own ALMA relative
// astrometry (Table 6: 2019.6505, PA 343.7729 deg, separation 5.27869").
void test_visual_orbit_convention()
{
    VisualOrbit ab;
    const double parallax_arcsec = 0.75081;
    ab.a_km = 17.4930 / parallax_arcsec * kAuKm;
    ab.e = 0.51947;
    ab.i_deg = 79.2430;
    ab.node_deg = 205.073;
    ab.arg_peri_deg = 231.519;
    ab.t_peri_tdb = tdb_from_julian_year(1955.564);
    ab.ra_deg = 219.85892215;
    ab.dec_deg = -60.83163195;
    const double gm = (1.0788 + 0.9092) * kSunGmKm3S2;

    const State s = visual_orbit_state(ab, gm, tdb_from_julian_year(2019.6505));
    const double ra = ab.ra_deg * kDegToRad;
    const double dec = ab.dec_deg * kDegToRad;
    const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
    const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
    const double distance_km = kAuKm / parallax_arcsec * 206264.80624709636; // 1 pc / parallax
    const double n = glm::dot(s.position, north) / distance_km * 206264.80624709636;
    const double e = glm::dot(s.position, east) / distance_km * 206264.80624709636;
    const double sep = std::sqrt(n * n + e * e);
    double pa = std::atan2(e, n) * kRadToDeg;
    pa = pa < 0.0 ? pa + 360.0 : pa;
    check(std::abs(sep - 5.27869) < 0.02, "alpha Cen AB separation (arcsec)", sep);
    check(std::abs(pa - 343.7729) < 0.3, "alpha Cen AB position angle (deg)", pa);

    // Kepler's third law with these elements gives the published total mass.
    const double period_s = kTwoPi * std::sqrt(ab.a_km * ab.a_km * ab.a_km / gm);
    check(std::abs(period_s / (kDaysPerJulianYear * kSecondsPerDay) - 79.762) < 0.05, "alpha Cen AB period", period_s);
}

void test_kerr_schwarzschild_precession()
{
    // S2-like orbit (GRAVITY 2020 elements) around a non-spinning hole.
    KerrOrbitSetup setup;
    setup.gm_km3_s2 = 4.261e6 * kSunGmKm3S2;
    setup.spin = 0.0;
    setup.elements.a_km = 0.125058 * 8246.7 * kAuKm;
    setup.elements.e = 0.884649;
    setup.elements.i_deg = 134.567;
    setup.elements.node_deg = 228.171;
    setup.elements.arg_peri_deg = 66.263;
    setup.elements.t_peri_tdb = tdb_from_julian_year(2018.379);
    setup.elements.ra_deg = 266.41683;
    setup.elements.dec_deg = -29.00782;
    const double period = kTwoPi * std::sqrt(std::pow(setup.elements.a_km, 3) / setup.gm_km3_s2);
    setup.elements_epoch_tdb = setup.elements.t_peri_tdb - 0.5 * period; // 2010 apocentre

    const double t_end = setup.elements_epoch_tdb + 4.2 * period;
    const KerrOrbitMotion orbit(setup, setup.elements_epoch_tdb, t_end);
    const auto& knots = orbit.table().knots();
    // A failed Newton solve ends the trajectory: the whole span must be covered.
    check(std::abs(knots.back().t - t_end) < 1.0, "S2 integration covers the span",
          (t_end - knots.back().t) / kSecondsPerDay);

    // Every Newton solve converges at the scene's step size (the residual floor
    // is ~1e-13 here; a tolerance below it would fail ~3% of steps).
    const KerrEquatorial& kerr = orbit.metric();
    const glm::dvec3 P = orbit.plane()[0];
    const glm::dvec3 Q = orbit.plane()[1];
    {
        const double m = kerr.mass();
        glm::dvec2 x(glm::dot(knots.front().position, P), glm::dot(knots.front().position, Q));
        glm::dvec2 u = kerr.u_from_velocity(
            x, glm::dvec2(glm::dot(knots.front().velocity, P), glm::dot(knots.front().velocity, Q)) /
                   kSpeedOfLightKmS);
        int failures = 0;
        for (double t = 0.0; t < period * kSpeedOfLightKmS;) {
            const double r = glm::length(x);
            const double dt = 0.003 * std::sqrt(r * r * r / m);
            failures += kerr.step(x, u, dt) ? 0 : 1;
            t += dt;
        }
        check(failures == 0, "S2 Kerr IMR Newton converges over one orbit", failures);
    }

    // Pericentre longitudes in the orbit plane: the knot with minimal r, refined
    // by a parabola through r(phi) of its neighbours (knots are ~14' apart in
    // angle at pericentre, comparable to the 12' precession being measured).
    auto phi_of = [&](const glm::dvec3& p) { return std::atan2(glm::dot(p, Q), glm::dot(p, P)); };
    std::vector<double> peri;
    for (size_t i = 1; i + 1 < knots.size(); ++i) {
        const double r0 = glm::length(knots[i - 1].position);
        const double r1 = glm::length(knots[i].position);
        const double r2 = glm::length(knots[i + 1].position);
        if (r1 < r0 && r1 <= r2) {
            const double p1 = phi_of(knots[i].position);
            const double a0 = wrap_pi(phi_of(knots[i - 1].position) - p1);
            const double a2 = wrap_pi(phi_of(knots[i + 1].position) - p1);
            // Vertex of the parabola through (a0, r0), (0, r1), (a2, r2).
            const double num = a0 * a0 * (r2 - r1) - a2 * a2 * (r0 - r1);
            const double den = a0 * (r2 - r1) - a2 * (r0 - r1);
            peri.push_back(p1 + 0.5 * num / den);
        }
    }
    check(peri.size() == 4, "S2: 4 pericentre passages", static_cast<double>(peri.size()));
    if (peri.size() >= 2) {
        const double advance = wrap_pi(peri.back() - peri.front()) / static_cast<double>(peri.size() - 1);
        const double m = setup.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
        const double expected =
            6.0 * kPi * m / (setup.elements.a_km * (1.0 - setup.elements.e * setup.elements.e));
        // ~12 arcmin per orbit; 1PN is accurate to O(M/r) ~ 1e-3 here.
        check(std::abs(advance / expected - 1.0) < 0.02, "S2 Schwarzschild precession / 1PN", advance / expected);
        std::printf("info: S2 periapsis advance %.2f arcmin/orbit (1PN %.2f)\n", advance * kRadToDeg * 60.0,
                    expected * kRadToDeg * 60.0);
    }
    // Radial period close to Keplerian (GR correction ~1e-3).
    check(knots.back().t > knots.front().t, "S2 integration span");
}

void test_kerr_conservation()
{
    // Strong field, rapidly spinning hole: energy and angular momentum stay constant.
    const double m = 10.0;
    const KerrEquatorial kerr(m, 0.9 * m);
    glm::dvec2 x(40.0 * m, 0.0);
    glm::dvec2 u = kerr.u_from_velocity(x, glm::dvec2(0.0, 0.17));
    const double e0 = kerr.energy(x, u);
    const double l0 = kerr.angular_momentum(x, u);
    // IMR is symplectic: the energy error is O(dt^2) and bounded (no drift).
    double worst_e[2] = {0.0, 0.0};
    double worst_l = 0.0;
    const int steps = 40000;
    for (int i = 0; i < steps; ++i) {
        const double r = glm::length(x);
        check(kerr.step(x, u, 0.01 * std::sqrt(r * r * r / m)), "Kerr IMR Newton converges");
        double& w = worst_e[i < steps / 2 ? 0 : 1];
        w = std::max(w, std::abs(kerr.energy(x, u) / e0 - 1.0));
        worst_l = std::max(worst_l, std::abs(kerr.angular_momentum(x, u) / l0 - 1.0));
    }
    check(worst_e[1] < 1e-6, "Kerr energy error small", worst_e[1]);
    check(worst_e[1] < 1.5 * worst_e[0] + 1e-12, "Kerr energy error bounded (no drift)", worst_e[1] / worst_e[0]);
    check(worst_l < 1e-6, "Kerr angular momentum conserved", worst_l);
    check(kerr.horizon_radius() > m && kerr.horizon_radius() < 2.0 * m, "Kerr horizon radius", kerr.horizon_radius());
}

void test_nbody()
{
    // Two bodies on a Keplerian orbit: one period later they are back.
    const double gm1 = kSunGmKm3S2;
    const double gm2 = 0.5 * kSunGmKm3S2;
    const double a = 10.0 * kAuKm;
    const double e = 0.5;
    const double r_peri = a * (1.0 - e);
    const double v_peri = std::sqrt((gm1 + gm2) * (1.0 + e) / r_peri);
    std::vector<NBodyParticle> p(2);
    p[0].gm = gm1;
    p[1].gm = gm2;
    p[1].position = glm::dvec3(r_peri, 0.0, 0.0);
    p[1].velocity = glm::dvec3(0.0, v_peri, 0.0);
    // Barycentric.
    const glm::dvec3 cm_v = (p[1].velocity * gm2) / (gm1 + gm2);
    const glm::dvec3 cm_r = (p[1].position * gm2) / (gm1 + gm2);
    for (NBodyParticle& q : p) {
        q.position -= cm_r;
        q.velocity -= cm_v;
    }
    const double period = kTwoPi * std::sqrt(a * a * a / (gm1 + gm2));
    const auto knots = integrate_nbody(p, 0.0, -period, period, period / 2000.0);

    auto relative_at = [&](double t) {
        EphemerisTable t1;
        EphemerisTable t2;
        t1.set_knots(knots[0]);
        t2.set_knots(knots[1]);
        return t2.eval(t).position - t1.eval(t).position;
    };
    const double err_fwd = glm::length(relative_at(period) - glm::dvec3(r_peri, 0.0, 0.0)) / a;
    const double err_bwd = glm::length(relative_at(-period) - glm::dvec3(r_peri, 0.0, 0.0)) / a;
    check(err_fwd < 1e-6 && err_bwd < 1e-6, "N-body returns after one period", std::max(err_fwd, err_bwd));

    std::vector<NBodyParticle> end = p;
    end[0].position = knots[0].back().position;
    end[0].velocity = knots[0].back().velocity;
    end[1].position = knots[1].back().position;
    end[1].velocity = knots[1].back().velocity;
    const double de = std::abs(nbody_energy(end) / nbody_energy(p) - 1.0);
    check(de < 1e-9, "N-body energy conserved", de);

    // Thinned to 24 knots per orbit, relative to body 0: Hermite interpolation
    // between knots stays close to the full-resolution path.
    NBodyOptions output;
    output.origin = {0};
    output.stride = {1, 2000 / 24};
    const auto thin = integrate_nbody(p, 0.0, -period, period, period / 2000.0, output);
    EphemerisTable full_table;
    EphemerisTable thin_table;
    std::vector<EphemerisTable::Knot> rel(knots[1].size());
    for (size_t k = 0; k < rel.size(); ++k) {
        rel[k] = {knots[1][k].t, knots[1][k].position - knots[0][k].position,
                  knots[1][k].velocity - knots[0][k].velocity};
    }
    full_table.set_knots(rel);
    thin_table.set_knots(thin[1]);
    check(thin[1].size() < 60, "thinned N-body knots", static_cast<double>(thin[1].size()));
    double worst = 0.0;
    for (int k = 0; k < 997; ++k) {
        const double t = -period + 2.0 * period * (k + 0.5) / 997.0;
        worst = std::max(worst, glm::length(thin_table.eval(t).position - full_table.eval(t).position) / r_peri);
    }
    check(worst < 2e-3, "thinned N-body interpolation error (/ r_peri)", worst);
}

// Earth's prime meridian faces the Sun at 12:00 UTC (to within the equation of
// time, < 4.5 deg); the Moon's near side faces the Earth (libration < 10 deg).
void test_earth_moon_rotation()
{
    Scene scene = load_scene_or_die("jwst.toml");
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

void test_alpha_centauri_scene()
{
    Scene scene = load_scene_or_die("alpha_centauri.toml");
    const int a = scene.find("Alpha Centauri A");
    const int b = scene.find("Alpha Centauri B");
    const int prox = scene.find("Proxima Centauri");
    const int sun = scene.find("Sun");

    // The integrated (three-body) AB orbit still reproduces the ALMA astrometry.
    const double t_alma = tdb_from_julian_year(2019.6505);
    const glm::dvec3 rel = scene.icrf_state_at(b, t_alma).position - scene.icrf_state_at(a, t_alma).position;
    const double ra = 219.85892215 * kDegToRad;
    const double dec = -60.83163195 * kDegToRad;
    const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
    const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
    const double arcsec_per_km = 0.75081 / kAuKm; // parallax / au
    const double n = glm::dot(rel, north) * arcsec_per_km;
    const double e = glm::dot(rel, east) * arcsec_per_km;
    check(std::abs(std::sqrt(n * n + e * e) - 5.27869) < 0.02, "alpha Cen N-body separation", std::sqrt(n * n + e * e));

    // Proxima ~12,947 au from the AB barycentre at the Hipparcos epoch [Kervella 2017].
    const double d = glm::length(scene.icrf_state_at(prox, tdb_from_julian_year(1991.25)).position) / kAuKm;
    check(std::abs(d - 12947.0) < 50.0, "Proxima distance (au)", d);

    // The Sun is ~1.33 pc away, in Cassiopeia.
    const double d_sun = glm::length(scene.icrf_state_at(sun, 0.0).position) / kParsecKm;
    check(std::abs(d_sun - 1.3319) < 1e-3, "Sun distance from alpha Cen (pc)", d_sun);

    // Proxima's planets [Suarez Mascareno et al. 2025, Table 3; Damasso et al. 2020,
    // Table 1]: circular orbits whose time of inferior conjunction T0 puts the planet
    // in front of the star, so Proxima's reflex RV is -K sin(2 pi (t - T0) / P)
    // (positive = receding) and the planet's own line-of-sight velocity relative
    // to Proxima is +v sin i sin(2 pi (t - T0) / P).
    struct ProximaPlanet {
        const char* name;
        double t0_bjd;
        double period_days;
        double a_au;
    };
    const ProximaPlanet planets[] = {
        {"Proxima Centauri d", 2460557.55, 5.12338, 0.02881},
        {"Proxima Centauri b", 2460548.59, 11.18465, 0.04848},
        {"Proxima Centauri c (unconfirmed)", 2455892.0, 1900.0, 1.48},
    };
    const glm::dvec3 away = unit_from_ra_dec(217.42894222 * kDegToRad, -62.67949019 * kDegToRad);
    const double sin_i = std::sin(47.0 * kDegToRad);
    scene.update(tdb_from_julian_year(2026.0)); // lighting goes by where the stars are
    for (const ProximaPlanet& p : planets) {
        const int index = scene.find(p.name);
        check(index >= 0 && scene.bodies[static_cast<size_t>(index)].parent == prox, p.name);
        check(scene.lighting_star(index) == prox, "lit by Proxima");
        const MotionSource& orbit = *scene.bodies[static_cast<size_t>(index)].motion;
        const double period = p.period_days * kSecondsPerDay;
        const double t0 = (p.t0_bjd - kJ2000Jd) * kSecondsPerDay;
        const double a_km = p.a_au * kAuKm;
        const double amplitude = kTwoPi / period * a_km * sin_i;
        double worst = 0.0;
        for (int k = 0; k < 12; ++k) {
            const double phase = k / 12.0;
            const double t = t0 + (7.0 + phase) * period;
            const double planet_los = glm::dot(orbit.eval(t).velocity, away);
            worst = std::max(worst, std::abs(planet_los - amplitude * std::sin(kTwoPi * phase)) / amplitude);
        }
        check(worst < 1e-6, "Proxima planet RV curve (relative error)", worst);
        // In front of the star at conjunction.
        const double depth = glm::dot(orbit.eval(t0).position, away);
        check(std::abs(depth / (a_km * sin_i) + 1.0) < 1e-6, "Proxima planet in front at conjunction", depth);
    }
}

// TransitOrbit: whatever the fit's convention, the planet is in front of the
// star at the given mid-transit time and keeps the given period.
void test_transit_orbit_convention()
{
    TransitOrbit t;
    t.period_days = 10.0;
    t.t_transit_tdb = 1.0e6;
    t.e = 0.2;
    t.arg_peri_deg = 40.0;
    t.i_deg = 88.0;
    t.node_deg = 30.0;
    t.ra_deg = 120.0;
    t.dec_deg = 35.0;
    const double gm = kSunGmKm3S2;
    const double ra = t.ra_deg * kDegToRad;
    const double dec = t.dec_deg * kDegToRad;
    const glm::dvec3 away(std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec));
    for (double u : {90.0, 270.0}) {
        t.transit_u_deg = u;
        const VisualOrbit o = visual_orbit_from_transit(t, gm);
        const glm::dvec3 r = visual_orbit_state(o, gm, t.t_transit_tdb).position;
        const double depth = glm::dot(glm::normalize(r), away);
        check(std::abs(depth + std::sin(t.i_deg * kDegToRad)) < 1e-9, "transit orbit: in front at mid-transit", depth);
        const double period_days = kTwoPi * std::sqrt(o.a_km * o.a_km * o.a_km / gm) / kSecondsPerDay;
        check(std::abs(period_days - t.period_days) < 1e-9, "transit orbit keeps the period", period_days);
    }
    // The two conventions differ only in where omega is measured from.
    t.transit_u_deg = 270.0;
    const VisualOrbit o270 = visual_orbit_from_transit(t, gm);
    t.transit_u_deg = 90.0;
    t.arg_peri_deg -= 180.0;
    const VisualOrbit o90 = visual_orbit_from_transit(t, gm);
    const double d = glm::length(visual_orbit_state(o270, gm, 0.0).position - visual_orbit_state(o90, gm, 0.0).position);
    check(d < 1e-3, "transit orbit conventions agree (km)", d);
}

// Mid-transit (minimum sky-projected separation, planet in front) of a planet
// around the root star, searched within +-window of t_guess; NaN if none.
double find_transit(const Scene& scene, int planet, const glm::dvec3& away, double t_guess, double window)
{
    const MotionSource& motion = *scene.bodies[static_cast<size_t>(planet)].motion;
    auto g = [&](double t) {
        const State s = motion.eval(t);
        const glm::dvec3 r = s.position - away * glm::dot(s.position, away);
        const glm::dvec3 v = s.velocity - away * glm::dot(s.velocity, away);
        return glm::dot(r, v);
    };
    const double step = window / 50.0;
    double t0 = t_guess - window;
    double g0 = g(t0);
    for (double t1 = t0 + step; t1 <= t_guess + window; t1 += step) {
        const double g1 = g(t1);
        if (g0 < 0.0 && g1 >= 0.0 && glm::dot(motion.eval(t1).position, away) < 0.0) {
            double lo = t0;
            double hi = t1;
            for (int k = 0; k < 60; ++k) {
                const double mid = 0.5 * (lo + hi);
                (g(mid) < 0.0 ? lo : hi) = mid;
            }
            return 0.5 * (lo + hi);
        }
        t0 = t1;
        g0 = g1;
    }
    return std::nan("");
}

glm::dvec3 unit_toward(double ra_deg, double dec_deg)
{
    const double ra = ra_deg * kDegToRad;
    const double dec = dec_deg * kDegToRad;
    return {std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec)};
}

// TRAPPIST-1 [Agol et al. 2021]: integrated from the Table 2 elements, the
// planets transit when the paper's posterior says they do (times_obs_and_posterior.txt
// in github.com/ericagol/TRAPPIST1_Spitzer: first, middle and last fitted transit,
// BJD_TDB - 2450000), with the Table 5 impact parameters.
void test_trappist1_scene()
{
    Scene scene = load_scene_or_die("trappist1.toml");
    const glm::dvec3 away = unit_toward(346.622368729, -5.041399251);
    const double r_star = scene.bodies[0].equatorial_radius_km;
    struct Planet {
        const char* name;
        double times[3];
        double impact; // b0 [Table 5]
    };
    const Planet planets[] = {
        {"TRAPPIST-1 b", {7322.517901, 7804.486285, 8769.942058}, 0.095},
        {"TRAPPIST-1 c", {7282.805834, 7796.232104, 8777.058537}, 0.109},
        {"TRAPPIST-1 d", {7560.801874, 7815.940948, 8374.794749}, 0.063},
        {"TRAPPIST-1 e", {7312.713901, 7831.152074, 8770.478717}, 0.191},
        {"TRAPPIST-1 f", {7321.522392, 7827.886782, 8785.389199}, 0.312},
        {"TRAPPIST-1 g", {7294.772224, 7961.825511, 8777.174343}, 0.379},
        {"TRAPPIST-1 h", {7662.550861, 7981.632640, 8769.837235}, 0.378},
    };
    double worst_min = 0.0;
    for (const Planet& p : planets) {
        const int index = scene.find(p.name);
        check(index > 0, p.name);
        if (index <= 0) {
            continue;
        }
        for (double t_post : p.times) {
            const double t_ref = tdb_from_jd_tdb(2450000.0 + t_post);
            const double t = find_transit(scene, index, away, t_ref, 0.1 * kSecondsPerDay);
            const double dt_min = std::isnan(t) ? 1e9 : (t - t_ref) / 60.0;
            worst_min = std::max(worst_min, std::abs(dt_min));
            if (std::abs(dt_min) > 3.0) {
                std::printf("  %s transit at %.6f: off by %.2f min\n", p.name, t_post, dt_min);
            }
        }
        const double t = find_transit(scene, index, away, tdb_from_jd_tdb(2450000.0 + p.times[0]), 0.1 * kSecondsPerDay);
        const glm::dvec3 r = scene.bodies[static_cast<size_t>(index)].motion->eval(t).position;
        const double impact = glm::length(r - away * glm::dot(r, away)) / r_star;
        check(std::abs(impact - p.impact) < 0.02, "TRAPPIST-1 impact parameter", impact);
    }
    std::printf("info: TRAPPIST-1 transit times vs Agol 2021 posterior: worst %.2f min\n", worst_min);
    check(worst_min < 3.0, "TRAPPIST-1 transit times match Agol 2021 (min)", worst_min);

    // The output is thinned to ~24 knots per orbit (span 2000-2050).
    const auto& b = static_cast<const EphemerisMotion&>(*scene.bodies[static_cast<size_t>(scene.find("TRAPPIST-1 b"))].motion);
    const double per_orbit = 1.510826 * static_cast<double>(b.table().knots().size()) / (50.0 * kDaysPerJulianYear);
    check(per_orbit > 20.0 && per_orbit < 30.0, "TRAPPIST-1 b knots per orbit", per_orbit);
}

// Mid-conjunction of `body` in front of `center` (minimum sky-projected
// separation), searched within +-window of t_guess; NaN if none.
double find_conjunction(const Scene& scene, int body, int center, const glm::dvec3& away, double t_guess,
                        double window)
{
    auto relative = [&](double t) {
        const State a = scene.icrf_state_at(body, t);
        const State b = scene.icrf_state_at(center, t);
        return State{a.position - b.position, a.velocity - b.velocity};
    };
    auto g = [&](double t) {
        const State s = relative(t);
        const glm::dvec3 r = s.position - away * glm::dot(s.position, away);
        const glm::dvec3 v = s.velocity - away * glm::dot(s.velocity, away);
        return glm::dot(r, v);
    };
    const double step = window / 50.0;
    double t0 = t_guess - window;
    double g0 = g(t0);
    for (double t1 = t0 + step; t1 <= t_guess + window; t1 += step) {
        const double g1 = g(t1);
        if (g0 < 0.0 && g1 >= 0.0 && glm::dot(relative(t1).position, away) < 0.0) {
            double lo = t0;
            double hi = t1;
            for (int k = 0; k < 60; ++k) {
                const double mid = 0.5 * (lo + hi);
                (g(mid) < 0.0 ? lo : hi) = mid;
            }
            return 0.5 * (lo + hi);
        }
        t0 = t1;
        g0 = g1;
    }
    return std::nan("");
}

// Kepler-47 [Orosz et al. 2019]: integrated from the Table 7 barycentric
// initial conditions (BJD - 2455000 below). In the Kepler years the planets
// transit the primary when the paper's model does (Table 2, "model time")
// to within the light-travel time across the orbits (2-8 min; ELC corrects
// for it, we do not). The model's predictions for 2019-2026 (Table 10) fall
// within a quarter of the predicted transit duration: c drifts ahead of ELC by
// about 0.2 primary radii on the sky (cause unknown; 1-sigma 0.04-0.26 d), b
// and d stay within minutes. The binary eclipses drift 0.6 s per orbit from
// the linear ephemeris of Eq. (1), most likely the GR terms ELC includes
// (Sect. 5.1) and we leave out.
void test_kepler47_scene()
{
    Scene scene = load_scene_or_die("kepler47.toml");
    const glm::dvec3 away = unit_toward(295.297909668, 46.920474260);
    const int a = scene.find("Kepler-47 A");
    const int b = scene.find("Kepler-47 B");
    check(a > 0 && b > 0, "Kepler-47 stars");
    if (a <= 0 || b <= 0) {
        return;
    }
    struct Transit {
        const char* name;
        double time;
        double duration_h; // 0: Kepler-era model time
    };
    const Transit transits[] = {
        {"Kepler-47 b", -30.80876, 0.0},  {"Kepler-47 b", 352.25805, 0.0},  {"Kepler-47 b", 1408.96429, 0.0},
        {"Kepler-47 d", 604.44003, 0.0},  {"Kepler-47 d", 1350.36740, 0.0}, {"Kepler-47 c", 246.65033, 0.0},
        {"Kepler-47 c", 550.47748, 0.0},  {"Kepler-47 c", 850.98582, 0.0},  {"Kepler-47 c", 1154.78852, 0.0},
        {"Kepler-47 b", 3522.04361, 4.86}, {"Kepler-47 b", 4961.81627, 3.34}, {"Kepler-47 b", 6402.38804, 10.99},
        {"Kepler-47 d", 3590.46599, 5.02}, {"Kepler-47 d", 5456.14904, 34.28}, {"Kepler-47 d", 6390.01925, 6.46},
        {"Kepler-47 c", 3575.38963, 5.78}, {"Kepler-47 c", 4483.55306, 6.61}, {"Kepler-47 c", 5996.00286, 6.10},
    };
    double worst_kepler_min = 0.0;
    double worst_fraction = 0.0;
    for (const Transit& p : transits) {
        const int index = scene.find(p.name);
        check(index > 0, p.name);
        if (index <= 0) {
            continue;
        }
        const double t_ref = tdb_from_jd_tdb(2455000.0 + p.time);
        const double t = find_conjunction(scene, index, a, away, t_ref, 0.5 * kSecondsPerDay);
        const double dt_min = std::isnan(t) ? 1e9 : (t - t_ref) / 60.0;
        if (p.duration_h == 0.0) {
            worst_kepler_min = std::max(worst_kepler_min, std::abs(dt_min));
        } else {
            worst_fraction = std::max(worst_fraction, std::abs(dt_min) / (60.0 * p.duration_h));
        }
    }
    std::printf("info: Kepler-47 transits vs Orosz 2019: Kepler era worst %.2f min, 2019-2026 worst %.2f durations\n",
                worst_kepler_min, worst_fraction);
    check(worst_kepler_min < 8.0, "Kepler-47 transits match Orosz 2019 Table 2 (min)", worst_kepler_min);
    check(worst_fraction < 0.25, "Kepler-47 transits within the Table 10 predictions (durations)", worst_fraction);

    // Primary eclipses: B in front of A, P = 7.44837568 d, T0 = BJD 2454963.246137.
    double worst_eclipse_min = 0.0;
    for (int n : {1, 190, 400, 816}) {
        const double t_ref = tdb_from_jd_tdb(2454963.246137 + 7.44837568 * n);
        const double t = find_conjunction(scene, b, a, away, t_ref, 0.1 * kSecondsPerDay);
        const double dt_min = std::isnan(t) ? 1e9 : (t - t_ref) / 60.0;
        worst_eclipse_min = std::max(worst_eclipse_min, std::abs(dt_min));
    }
    std::printf("info: Kepler-47 primary eclipses vs linear ephemeris 2009-2025: worst %.2f min\n",
                worst_eclipse_min);
    check(worst_eclipse_min < 10.0, "Kepler-47 primary eclipses near the linear ephemeris (min)", worst_eclipse_min);

    // Both stars light the planets, the red dwarf with ~2% of the light; the Sun not at all.
    scene.update(tdb_from_jd_tdb(2461041.5));
    StarLight lights[kMaxStarLights];
    const int count = scene.lighting_stars(scene.find("Kepler-47 c"), lights, kMaxStarLights);
    check(count == 2 && lights[0].star == a && lights[1].star == b, "Kepler-47 c lit by A and B", count);
    if (count == 2) {
        check(lights[1].relative_flux > 0.01 && lights[1].relative_flux < 0.03, "Kepler-47 B adds ~2% light",
              lights[1].relative_flux);
    }
}

// Kepler-64 [Schwamb et al. 2013]: integrated from the Table 7 joint solution,
// the binary eclipses and the planet transits Aa when they were measured
// (Table 3, 1-sigma 0.9 and 7 min; Table 2, 6-9 min; BJD - 2455000), and Ba +
// Bb stand where Gaia sees them.
void test_kepler64_scene()
{
    Scene scene = load_scene_or_die("kepler64.toml");
    const glm::dvec3 away = unit_toward(298.215070011, 39.955103210);
    const int aa = scene.find("Kepler-64 Aa");
    const int ab = scene.find("Kepler-64 Ab");
    const int b = scene.find("Kepler-64 b");
    check(aa > 0 && ab > 0 && b > 0, "Kepler-64 bodies");
    if (aa <= 0 || ab <= 0 || b <= 0) {
        return;
    }
    struct Event {
        int body;
        int center;
        double time;
    };
    const Event primary[] = {{ab, aa, -32.18064}, {ab, aa, 447.82520}, {ab, aa, 927.83150}};
    const Event secondary[] = {{aa, ab, -24.32048}, {aa, ab, 455.67885}, {aa, ab, 915.68830}};
    const Event transits[] = {{b, aa, 70.80674}, {b, aa, 344.11218}, {b, aa, 613.17869}, {b, aa, 885.91042}};
    auto worst = [&](const auto& events) {
        double w = 0.0;
        for (const Event& e : events) {
            const double t_ref = tdb_from_jd_tdb(2455000.0 + e.time);
            const double t = find_conjunction(scene, e.body, e.center, away, t_ref, 0.5 * kSecondsPerDay);
            w = std::max(w, std::isnan(t) ? 1e9 : std::abs(t - t_ref) / 60.0);
        }
        return w;
    };
    const double w_primary = worst(primary);
    const double w_secondary = worst(secondary);
    const double w_transit = worst(transits);
    std::printf("info: Kepler-64 vs Schwamb 2013: primary eclipses %.2f, secondary %.2f, transits %.2f min\n",
                w_primary, w_secondary, w_transit);
    check(w_primary < 3.0, "Kepler-64 primary eclipses (min)", w_primary);
    check(w_secondary < 20.0, "Kepler-64 secondary eclipses (min)", w_secondary);
    check(w_transit < 20.0, "Kepler-64 planet transits (min)", w_transit);

    // Ba + Bb: 0.7043" at position angle 123.29 deg from the A system (J2016.0, 1906.58 pc).
    const int bb = scene.find("Kepler-64 B");
    check(bb > 0, "Kepler-64 B");
    if (bb > 0) {
        const glm::dvec3 r = scene.icrf_state_at(bb, tdb_from_julian_year(2016.0)).position;
        const double ra = 298.215070011 * kDegToRad;
        const double dec = 39.955103210 * kDegToRad;
        const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
        const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
        const double pa = wrap_two_pi(std::atan2(glm::dot(r, east), glm::dot(r, north))) * kRadToDeg;
        const double sep = std::hypot(glm::dot(r, east), glm::dot(r, north)) / (1906.58 * kAuKm);
        check(std::abs(pa - 123.29) < 0.05, "Kepler-64 B position angle (deg)", pa);
        check(std::abs(sep - 0.7043) < 0.001, "Kepler-64 B separation (arcsec)", sep);
    }
}

// PSR B1620-26 [Thorsett et al. 1999]: the pulsar's line-of-sight motion
// reproduces the Table 1 timing orbit (Roemer delay x (1 - e) sin omega at
// periastron, semi-amplitude K), the inner binary moves 6.4 lt-s around the
// triple's barycentre and recedes through its node at JD 2449104 (Table 2),
// the two orbits are 40 deg apart [Sigurdsson & Thorsett 2005], and the sky is
// M4 seen from inside.
void test_psr_b1620_scene()
{
    Scene scene = load_scene_or_die("psr_b1620.toml");
    const glm::dvec3 away = unit_toward(245.9092575, -26.5316025);
    const int pulsar = scene.find("PSR B1620-26 A");
    const int wd = scene.find("PSR B1620-26 B");
    const int ab = scene.find("PSR B1620-26 AB");
    const int planet = scene.find("PSR B1620-26 b");
    check(pulsar > 0 && wd > 0 && ab > 0 && planet > 0, "PSR B1620-26 bodies");
    if (pulsar <= 0 || wd <= 0 || ab <= 0 || planet <= 0) {
        return;
    }
    const double lt_s = kSpeedOfLightKmS; // km per light-second
    auto pulsar_depth = [&](double t) {
        const State p = scene.icrf_state_at(pulsar, t);
        const State c = scene.icrf_state_at(ab, t);
        return std::make_pair(glm::dot(p.position - c.position, away), glm::dot(p.velocity - c.velocity, away));
    };
    const double x = 64.809460;
    const double e = 0.02531545;
    const double w = 117.1291 * kDegToRad;
    const double t_peri = tdb_from_jd_tdb(2448728.76242);
    const double depth_peri = pulsar_depth(t_peri).first / lt_s;
    check(std::abs(depth_peri - x * (1.0 - e) * std::sin(w)) < 1e-4, "PSR B1620-26 Roemer delay at periastron (s)",
          depth_peri);
    double k_max = 0.0;
    for (int i = 0; i < 2000; ++i) {
        k_max = std::max(k_max, pulsar_depth(t_peri + 191.44281 * kSecondsPerDay * i / 2000.0).second);
    }
    // Maximum radial velocity K (1 + e cos omega), K = 2 pi x c / (P sqrt(1 - e^2)).
    const double k = kTwoPi * x * lt_s / (191.44281 * kSecondsPerDay * std::sqrt(1.0 - e * e));
    check(std::abs(k_max / (k * (1.0 + e * std::cos(w))) - 1.0) < 1e-4, "PSR B1620-26 pulsar velocity amplitude",
          k_max);

    // The binary around the triple's barycentre (the root): node at JD 2449104, 6.4 lt-s.
    const double t_node = tdb_from_jd_tdb(2449104.0);
    const double quarter = 61.8 * kDaysPerJulianYear * kSecondsPerDay / 4.0;
    const State at_node = scene.icrf_state_at(ab, t_node);
    const double depth_q = glm::dot(scene.icrf_state_at(ab, t_node + quarter).position, away) / lt_s;
    check(std::abs(glm::dot(at_node.position, away)) / lt_s < 1e-6 && glm::dot(at_node.velocity, away) > 0.0,
          "PSR B1620-26 binary recedes through its node in 1993");
    check(std::abs(depth_q - 6.4) < 0.01, "PSR B1620-26 binary outer amplitude (lt-s)", depth_q);

    scene.update(tdb_from_julian_year(2026.0));
    const double tilt = std::acos(glm::dot(scene.orbit_normal(wd), scene.orbit_normal(planet))) * kRadToDeg;
    check(std::abs(tilt - 40.0) < 0.05, "PSR B1620-26 relative inclination (deg)", tilt);

    // The sky: ~13,500 members, many brighter than the naked-eye limit, the
    // brightest far brighter than Sirius.
    int naked_eye = 0;
    double brightest = 99.0;
    for (const CatalogStar& s : scene.sky.stars) {
        naked_eye += s.vmag < 6.5 ? 1 : 0;
        brightest = std::min(brightest, s.vmag);
    }
    std::printf("info: M4 from PSR B1620-26: %zu stars, %d brighter than V = 6.5, brightest V = %.2f\n",
                scene.sky.stars.size(), naked_eye, brightest);
    check(scene.sky.stars.size() > 13000, "M4 sky star count", static_cast<double>(scene.sky.stars.size()));
    check(naked_eye > 5000 && brightest < -4.0, "M4 sky is bright", brightest);

    // The pulsar's spin axis is the inner orbit's normal.
    check(scene.bodies[static_cast<size_t>(pulsar)].pulsar.enabled &&
              scene.bodies[static_cast<size_t>(pulsar)].pulsar.spin_axis_orbit_of == wd,
          "PSR B1620-26 spin axis");
}

// TIC 168789840 [Powell et al. 2021]: nested Keplerian orbits. Each binary's
// secondary passes in front of its primary at the Table 2 primary eclipses,
// close enough on the sky to eclipse; B stands at the Table 4 speckle position
// relative to AC; A and C keep the Table 8 period and eccentricity.
void test_tic168789840_scene()
{
    Scene scene = load_scene_or_die("tic168789840.toml");
    const glm::dvec3 away = unit_toward(63.5201633, -31.9228995);
    struct Binary {
        const char* primary;
        const char* secondary;
        double period_days;
        double t0_bjd;
    };
    const Binary binaries[] = {
        {"TIC 168789840 A1", "TIC 168789840 A2", 1.570013, 2458412.3855},
        {"TIC 168789840 B1", "TIC 168789840 B2", 8.217111, 2458411.9008},
        {"TIC 168789840 C1", "TIC 168789840 C2", 1.305883, 2458413.6822},
    };
    for (const Binary& bin : binaries) {
        const int p = scene.find(bin.primary);
        const int s = scene.find(bin.secondary);
        check(p > 0 && s > 0, bin.secondary);
        if (p <= 0 || s <= 0) {
            continue;
        }
        for (int n : {0, 1000}) {
            const double t_ref = tdb_from_jd_tdb(bin.t0_bjd + bin.period_days * n);
            const double t = find_conjunction(scene, s, p, away, t_ref, 0.2 * kSecondsPerDay);
            const double dt_min = std::isnan(t) ? 1e9 : (t - t_ref) / 60.0;
            check(std::abs(dt_min) < 0.1, "TIC 168789840 primary eclipse on the ephemeris (min)", dt_min);
            if (!std::isnan(t)) {
                const glm::dvec3 r = scene.icrf_state_at(s, t).position - scene.icrf_state_at(p, t).position;
                const double sky = glm::length(r - away * glm::dot(r, away));
                const double radii = scene.bodies[static_cast<size_t>(p)].equatorial_radius_km +
                                     scene.bodies[static_cast<size_t>(s)].equatorial_radius_km;
                check(sky < radii, "TIC 168789840 binaries eclipse", sky / radii);
            }
        }
    }

    // B relative to AC at the speckle epoch: PA 257.74 deg, 0.4230" (at 584 pc).
    const int ac = scene.find("TIC 168789840 AC");
    const int b = scene.find("TIC 168789840 B");
    check(ac > 0 && b > 0, "TIC 168789840 AC and B");
    if (ac > 0 && b > 0) {
        const double t = tdb_from_julian_year(2020.8236);
        const glm::dvec3 r = scene.icrf_state_at(b, t).position - scene.icrf_state_at(ac, t).position;
        const double ra = 63.5201633 * kDegToRad;
        const double dec = -31.9228995 * kDegToRad;
        const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
        const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
        const double pa = wrap_two_pi(std::atan2(glm::dot(r, east), glm::dot(r, north))) * kRadToDeg;
        const double sep = std::hypot(glm::dot(r, east), glm::dot(r, north)) / (584.0 * kAuKm) * 1000.0;
        check(std::abs(pa - 257.74) < 0.05, "TIC 168789840 B position angle (deg)", pa);
        check(std::abs(sep - 423.0) < 0.5, "TIC 168789840 B separation (mas)", sep);
    }

    // A and C: closest at periastron (BJD 2457662), a (1 - e) apart, every 3.7 years.
    const int a = scene.find("TIC 168789840 A");
    const int c = scene.find("TIC 168789840 C");
    check(a > 0 && c > 0, "TIC 168789840 A and C");
    if (a > 0 && c > 0) {
        auto separation = [&](double jd) {
            const double t = tdb_from_jd_tdb(jd);
            return glm::length(scene.icrf_state_at(c, t).position - scene.icrf_state_at(a, t).position) / kAuKm;
        };
        const double a_au = 3.67647;
        for (double jd : {2457662.0, 2457662.0 + 1351.425 * 3}) {
            check(std::abs(separation(jd) - a_au * (1.0 - 0.28)) < 1e-3, "TIC 168789840 AC periastron (au)",
                  separation(jd));
        }
        check(std::abs(separation(2457662.0 + 0.5 * 1351.425) - a_au * 1.28) < 1e-3, "TIC 168789840 AC apastron (au)",
              separation(2457662.0 + 0.5 * 1351.425));
    }
}

// Kepler-223 [Mills et al. 2016]: integrated from the best-fit initial
// conditions (Extended Data Table 3), the planets reproduce the transit-timing
// variations measured in each Kepler quarter (Extended Data Table 1: mean time
// of the quarter's transits, then the TTV relative to the given linear
// ephemeris with its -3, -1, +1, +3 sigma offsets; days, BJD - 2454900).
void test_kepler223_scene()
{
    Scene scene = load_scene_or_die("kepler223.toml");
    const glm::dvec3 away = unit_toward(298.318417682, 47.279530121);
    struct Quarter {
        double t_mean, m3, m1, ttv, p1, p3;
    };
    struct Planet {
        const char* name;
        double period_days;
        double t0;
        Quarter quarters[11];
    };
    const Planet planets[] = {
        {"Kepler-223 b", 7.3840154, 70.49489,
         {{123.32662, -0.0354, -0.0058, -0.0006, 0.0059, 0.0316},
          {239.02516, -0.0517, -0.0103, 0.0137, 0.0101, 0.0423},
          {416.51724, -0.0200, -0.0061, -0.0010, 0.0062, 0.0210},
          {521.69775, -0.0628, -0.0123, 0.0068, 0.0113, 0.0342},
          {699.18988, -0.0470, -0.0088, -0.0010, 0.0084, 0.0370},
          {797.79657, -0.0417, -0.0097, -0.0123, 0.0088, 0.0243},
          {886.54260, -0.0343, -0.0072, -0.0187, 0.0074, 0.0617},
          {1073.89526, -0.0500, -0.0118, 0.0730, 0.0140, 0.1150},
          {1162.64136, -0.0542, -0.0071, 0.0692, 0.0071, 0.0708},
          {1251.38745, -0.0217, -0.0062, 0.0507, 0.0065, 0.0333},
          {1458.46155, -0.0379, -0.0129, 0.0659, 0.0090, 0.0241}}},
        {"Kepler-223 c", 9.8487130, 71.37624,
         {{116.10564, -0.0362, -0.0103, -0.0168, 0.0133, 0.0518},
          {242.86336, -0.0683, -0.0405, 0.0023, 0.0077, 0.0747},
          {420.32413, -0.0254, -0.0077, 0.0264, 0.0076, 0.0476},
          {509.05453, -0.0266, -0.0077, 0.0166, 0.0090, 0.0304},
          {701.30371, -0.0585, -0.0121, -0.0155, 0.0126, 0.0535},
          {790.03418, -0.0302, -0.0075, -0.0318, 0.0084, 0.0268},
          {886.15869, -0.0173, -0.0046, -0.0737, 0.0048, 0.0537},
          {1071.01367, -0.1404, -0.0078, -0.0766, 0.0066, 0.0226},
          {1148.65283, -0.0411, -0.0067, -0.0959, 0.0064, 0.0349},
          {1252.17163, -0.0392, -0.0067, -0.0418, 0.0068, 0.0548},
          {1470.30054, -0.0361, -0.0056, -0.0449, 0.0051, 0.0179}}},
        {"Kepler-223 d", 14.7883997, 109.76775,
         {{132.10997, -0.0416, -0.0058, 0.0376, 0.0054, 0.0134},
          {248.65308, -0.0221, -0.0062, -0.0169, 0.0063, 0.0229},
          {427.57138, -0.0351, -0.0084, -0.0099, 0.0070, 0.0169},
          {519.49268, -0.0285, -0.0070, 0.0035, 0.0066, 0.0245},
          {711.54260, -0.0260, -0.0086, 0.0240, 0.0094, 0.0420},
          {800.18097, -0.0226, -0.0060, 0.0256, 0.0057, 0.0194},
          {898.66815, -0.0192, -0.0057, 0.0212, 0.0055, 0.0238},
          {1077.35193, -0.0530, -0.0080, 0.0020, 0.0077, 0.0210},
          {1169.50781, -0.0354, -0.0085, -0.0236, 0.0093, 0.0286},
          {1271.27771, -0.0272, -0.0131, -0.0578, 0.0132, 0.0328},
          {1483.43542, -0.0298, -0.0061, -0.1612, 0.0057, 0.0302}}},
        {"Kepler-223 e", 19.7213435, 68.10686,
         {{135.47421, -0.0303, -0.0060, -0.0067, 0.0053, 0.0187},
          {238.21753, -0.0232, -0.0067, 0.0022, 0.0072, 0.0458},
          {433.78842, -0.0542, -0.0095, 0.0302, 0.0084, 0.0298},
          {524.27625, -0.0244, -0.0061, -0.0106, 0.0063, 0.0296},
          {709.21222, -0.0432, -0.0071, 0.0022, 0.0065, 0.0208},
          {797.82037, -0.0240, -0.0060, 0.0090, 0.0061, 0.0240},
          {893.81256, -0.0357, -0.0216, 0.0297, 0.0242, 0.0513},
          {1079.88989, -0.0662, -0.0083, 0.0602, 0.0078, 0.0308},
          {1170.71301, -0.1067, -0.0118, 0.1167, 0.0110, 0.0453},
          {1263.01343, -0.0252, -0.0049, 0.1352, 0.0049, 0.0188},
          {1469.48169, -0.0393, -0.0097, 0.2283, 0.0100, 0.0467}}},
    };
    double residual_sq = 0.0;
    double signal_sq = 0.0;
    int within_3sigma = 0;
    int count = 0;
    for (const Planet& p : planets) {
        const int index = scene.find(p.name);
        check(index > 0, p.name);
        if (index <= 0) {
            continue;
        }
        for (const Quarter& q : p.quarters) {
            // Mean modelled TTV of the transits within half a quarter of the mean time.
            double sum = 0.0;
            int n_transits = 0;
            const int first = static_cast<int>(std::ceil((q.t_mean - 45.0 - p.t0) / p.period_days));
            const int last = static_cast<int>(std::floor((q.t_mean + 45.0 - p.t0) / p.period_days));
            for (int n = first; n <= last; ++n) {
                const double linear = 2454900.0 + p.t0 + n * p.period_days;
                const double t = find_transit(scene, index, away, tdb_from_jd_tdb(linear), 0.5 * kSecondsPerDay);
                if (!std::isnan(t)) {
                    sum += t / kSecondsPerDay + kJ2000Jd - linear;
                    ++n_transits;
                }
            }
            check(n_transits > 0, "Kepler-223 transits found in a quarter", q.t_mean);
            if (n_transits == 0) {
                continue;
            }
            const double model = sum / n_transits;
            ++count;
            residual_sq += (model - q.ttv) * (model - q.ttv);
            signal_sq += q.ttv * q.ttv;
            if (model >= q.ttv + q.m3 && model <= q.ttv + q.p3) {
                ++within_3sigma;
            } else {
                std::printf("  %s quarter at %.1f: model TTV %+.4f d, measured %+.4f d\n", p.name, q.t_mean, model,
                            q.ttv);
            }
        }
    }
    // The binned TTVs scatter around the photodynamical model by 1-3 sigma
    // ([MI16] Fig. 1), but the model must explain the signal (up to 0.23 d).
    const double rms_residual = std::sqrt(residual_sq / std::max(count, 1));
    const double rms_signal = std::sqrt(signal_sq / std::max(count, 1));
    std::printf("info: Kepler-223 quarterly TTVs: %d of %d within 3 sigma, rms %.4f d (signal %.4f d)\n",
                within_3sigma, count, rms_residual, rms_signal);
    check(count == 44 && within_3sigma == count, "Kepler-223 TTVs within 3 sigma", within_3sigma);
    check(rms_residual < 0.25 * rms_signal, "Kepler-223 TTV residuals small against the signal", rms_residual);
}

// Periastron advance and orbital decay of general relativity at first
// post-Newtonian order for two point masses (M_sun), e.g. Lorimer & Kramer,
// Handbook of Pulsar Astronomy (2005), eqs. 8.53 and 8.56.
double gr_omega_dot_deg_per_year(double m1, double m2, double pb_days, double e)
{
    const double t_sun = kSunGmKm3S2 / (kSpeedOfLightKmS * kSpeedOfLightKmS * kSpeedOfLightKmS); // s
    const double n = kTwoPi / (pb_days * kSecondsPerDay);
    const double rate = 3.0 * std::pow(t_sun * (m1 + m2) * n, 2.0 / 3.0) * n / (1.0 - e * e);
    return rate * kRadToDeg * kDaysPerJulianYear * kSecondsPerDay;
}

double gr_pb_dot(double m1, double m2, double pb_days, double e)
{
    const double t_sun = kSunGmKm3S2 / (kSpeedOfLightKmS * kSpeedOfLightKmS * kSpeedOfLightKmS);
    const double n = kTwoPi / (pb_days * kSecondsPerDay);
    const double e2 = e * e;
    const double f = (1.0 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e2 * e2) / std::pow(1.0 - e2, 3.5);
    return -192.0 * kPi / 5.0 * std::pow(t_sun * n, 5.0 / 3.0) * f * m1 * m2 / std::cbrt(m1 + m2);
}

// The post-Keplerian (DD) orbit: Keplerian without the relativistic terms,
// velocity = d position / dt, the periastron turning by k 2 pi per orbit, the
// quadratic phase of a changing period, and the Newtonian continuation.
void test_post_keplerian()
{
    VisualOrbit el;
    el.e = 0.6;
    el.i_deg = 50.0;
    el.node_deg = 30.0;
    el.arg_peri_deg = 100.0;
    el.t_peri_tdb = 1.0e8;
    el.ra_deg = 100.0;
    el.dec_deg = 20.0;
    const double gm = 2.8 * kSunGmKm3S2;
    const double pb = 0.3 * kSecondsPerDay;

    const PostKeplerianOrbit kepler = make_post_keplerian(el, pb, 0.0, 0.0, gm);
    const double t = el.t_peri_tdb + 12345.6;
    const State a = post_keplerian_state(kepler, gm, t);
    const State b = visual_orbit_state(kepler.elements, gm, t);
    check(glm::length(a.position - b.position) < 1e-6 && glm::length(a.velocity - b.velocity) < 1e-9,
          "post-Keplerian without relativistic terms is Keplerian", glm::length(a.position - b.position));

    // Exaggerated: the periastron turns 0.05 rad per radian of true anomaly.
    const double k = 0.05;
    const PostKeplerianOrbit fast = make_post_keplerian(el, pb, k * kTwoPi / pb, 0.0, gm);
    double worst = 0.0;
    for (double dt : {0.0, 3.0, 1234.5, 0.5 * pb, 7.9 * pb, -2.3 * pb}) {
        const double h = 0.5;
        const State s = post_keplerian_state(fast, gm, el.t_peri_tdb + dt);
        const glm::dvec3 fd = (post_keplerian_state(fast, gm, el.t_peri_tdb + dt + h).position -
                               post_keplerian_state(fast, gm, el.t_peri_tdb + dt - h).position) /
                              (2.0 * h);
        worst = std::max(worst, glm::length(fd - s.velocity) / glm::length(s.velocity));
    }
    check(worst < 1e-6, "post-Keplerian velocity is the derivative of position", worst);

    VisualOrbit turned = fast.elements;
    turned.arg_peri_deg += k * kTwoPi * 7.0 * kRadToDeg;
    const glm::dvec3 peri = glm::normalize(post_keplerian_state(fast, gm, el.t_peri_tdb + 7.0 * pb).position);
    const double off = std::acos(std::min(1.0, glm::dot(peri, visual_orbit_frame(turned)[0])));
    check(off < 1e-7, "periastron turns by 2 pi k per orbit", off);

    // A shrinking period: the N-th periastron comes ~Pb_dot N^2 / 2 periods early
    // (solving N = o - Pb_dot o^2 / 2 for the orbits o to second order in Pb_dot).
    const double pb_dot = -1e-6;
    const PostKeplerianOrbit decaying = make_post_keplerian(el, pb, 0.0, pb_dot, gm);
    const double n = 100.0;
    const double t_n = el.t_peri_tdb + pb * (n + 0.5 * pb_dot * n * n + 0.5 * pb_dot * pb_dot * n * n * n);
    const PostKeplerianPhase at_n = post_keplerian_phase(decaying, gm, t_n);
    check(std::abs(wrap_pi(at_n.mean_anomaly)) < 1e-6, "period decay: quadratic orbital phase", at_n.mean_anomaly);
    const double shrink = std::pow(1.0 + pb_dot * (t_n - el.t_peri_tdb) / pb, 2.0 / 3.0);
    check(std::abs(at_n.a_km / decaying.elements.a_km - shrink) < 1e-12, "semi-major axis follows the period",
          at_n.a_km);

    // Newton from t_ref: through the same point, with that moment's periastron, forever.
    const double t_ref = el.t_peri_tdb + 3.3 * pb;
    const VisualOrbit newton = newtonian_continuation(fast, gm, t_ref);
    const double moved = glm::length(visual_orbit_state(newton, gm, t_ref).position -
                                     post_keplerian_state(fast, gm, t_ref).position);
    check(moved < 1e-5, "Newtonian continuation starts at the same point", moved);
    check(std::abs(newton.arg_peri_deg * kDegToRad - post_keplerian_phase(fast, gm, t_ref).arg_peri_rad) < 1e-12,
          "Newtonian continuation keeps that moment's periastron");

    // The drawn ellipse starts at the body and closes.
    const PostKeplerianMotion motion(fast, gm);
    std::vector<glm::dvec3> ring;
    motion.sample_orbit(t, 64, ring);
    check(ring.size() == 64 && glm::length(ring.front() - motion.eval(t).position) < 1e-6 &&
              glm::length(ring.back() - ring.front()) < 1e-6,
          "post-Keplerian orbit line starts at the body and closes");
}

// The post-Keplerian orbit inside a split pair (ScaledMotion), or a ghost's ellipse.
const PostKeplerianMotion* post_keplerian_of(const Scene& scene, int body)
{
    const auto* scaled = dynamic_cast<const ScaledMotion*>(scene.bodies[static_cast<size_t>(body)].motion.get());
    return scaled ? dynamic_cast<const PostKeplerianMotion*>(&scaled->inner()) : nullptr;
}

const VisualOrbitMotion* ghost_orbit_of(const Scene& scene, int body)
{
    const auto* scaled = dynamic_cast<const ScaledMotion*>(scene.bodies[static_cast<size_t>(body)].motion.get());
    return scaled ? dynamic_cast<const VisualOrbitMotion*>(&scaled->inner()) : nullptr;
}

// A pulsar's distance behind the system's barycentre (the root) along the line of sight, in lt-s.
double pulsar_depth_s(const Scene& scene, int pulsar, const glm::dvec3& away, double t)
{
    return glm::dot(scene.icrf_state_at(pulsar, t).position - scene.icrf_state_at(0, t).position, away) /
           kSpeedOfLightKmS;
}

void test_psr_b1913_scene()
{
    Scene scene = load_scene_or_die("psr_b1913.toml");
    const int pulsar = scene.find("PSR B1913+16");
    const int companion = scene.find("B1913+16 companion");
    const int ghost_p = scene.find("PSR B1913+16 (Newton)");
    const int ghost_c = scene.find("Companion (Newton)");
    check(pulsar > 0 && companion > 0 && ghost_p > 0 && ghost_c > 0, "PSR B1913+16 bodies");
    if (pulsar <= 0 || companion <= 0 || ghost_p <= 0 || ghost_c <= 0) {
        return;
    }
    const PostKeplerianMotion* orbit = post_keplerian_of(scene, pulsar);
    const VisualOrbitMotion* newton = ghost_orbit_of(scene, ghost_p);
    check(orbit && newton, "PSR B1913+16 post-Keplerian orbit and Newtonian ghost");
    if (!orbit || !newton) {
        return;
    }

    // Roemer delay at periastron T0: x (1 - e) sin(omega0) [WH16 Table 2].
    const glm::dvec3 away = unit_toward(288.86666425, 16.10760744);
    const double x = 2.341776;
    const double e = 0.6171340;
    const double w = 292.54450 * kDegToRad;
    const double t0 = tdb_from_jd_tdb(2452145.40097849);
    const double depth = pulsar_depth_s(scene, pulsar, away, t0);
    check(std::abs(depth - x * (1.0 - e) * std::sin(w)) < 2e-5, "PSR B1913+16 Roemer delay at T0 (s)", depth);
    // Over one orbit the line-of-sight depth spans 2 x sqrt(1 - e^2 cos^2 omega).
    double z_min = 1e30;
    double z_max = -1e30;
    for (int i = 0; i < 4000; ++i) {
        const double z = pulsar_depth_s(scene, pulsar, away, t0 + 0.322997448918 * kSecondsPerDay * i / 4000.0);
        z_min = std::min(z_min, z);
        z_max = std::max(z_max, z);
    }
    const double span = 2.0 * x * std::sqrt(1.0 - e * e * std::cos(w) * std::cos(w));
    check(std::abs((z_max - z_min) / span - 1.0) < 1e-4, "PSR B1913+16 projected semi-major axis", z_max - z_min);

    // The measured advance and decay against GR for the [WH16] masses.
    const double omdot_gr = gr_omega_dot_deg_per_year(1.438, 1.390, 0.322997448918, e);
    check(std::abs(omdot_gr / 4.226585 - 1.0) < 2e-4, "PSR B1913+16 omega_dot is GR's for its masses", omdot_gr);
    const double pbdot_gr = gr_pb_dot(1.438, 1.390, 0.322997448918, e);
    check(std::abs(pbdot_gr / -2.40263e-12 - 1.0) < 1e-3, "PSR B1913+16 Pb_dot(GR) [WH16 Sec. 7]", pbdot_gr);
    std::printf("info: PSR B1913+16 GR: omega_dot %.5f deg/yr, Pb_dot %.5e\n", omdot_gr, pbdot_gr);

    // The Newtonian ghost leaves the real orbit on the day of discovery; by
    // 2026-10 the real periastron has turned ~221 deg away from it.
    double t_disc = 0.0;
    check(parse_utc("1974-07-02", &t_disc), "parse the discovery date");
    const double split = glm::length(scene.icrf_state_at(pulsar, t_disc).position -
                                     scene.icrf_state_at(ghost_p, t_disc).position);
    check(split < 1e-3, "PSR B1913+16 ghost starts on the real orbit (km)", split);
    const double t_now = tdb_from_julian_year(2026.76);
    const double turned = (post_keplerian_phase(orbit->orbit(), 3.7531078044e11, t_now).arg_peri_rad -
                           newton->orbit().arg_peri_deg * kDegToRad) *
                          kRadToDeg;
    const double expected = 4.226585 * (2026.76 - (2000.0 + t_disc / (kDaysPerJulianYear * kSecondsPerDay)));
    check(std::abs(turned - expected) < 0.01, "PSR B1913+16 periastron advance since 1974 (deg)", turned);
    std::printf("info: PSR B1913+16 periastron advance since discovery: %.2f deg\n", turned);

    // Ghosts pair among themselves and carry no mass; the pulsar's beams turn
    // about the real orbit's normal.
    check(scene.satellite_host(companion) == pulsar && scene.satellite_host(ghost_c) == ghost_p &&
              scene.satellite_host(ghost_p) == -1,
          "PSR B1913+16 hosts (ghosts pair with ghosts)");
    check(scene.bodies[static_cast<size_t>(ghost_p)].body_gm_km3_s2 == 0.0 &&
              std::abs(scene.bodies[static_cast<size_t>(pulsar)].body_gm_km3_s2 / kSunGmKm3S2 - 1.438) < 1e-6,
          "PSR B1913+16 masses (ghosts massless)");
    check(scene.bodies[static_cast<size_t>(pulsar)].pulsar.spin_axis_orbit_of == companion,
          "PSR B1913+16 spin axis along the orbit normal");
}

void test_psr_j1141_scene()
{
    Scene scene = load_scene_or_die("psr_j1141.toml");
    const int pulsar = scene.find("PSR J1141-6545");
    const int wd = scene.find("J1141-6545 white dwarf");
    const int ghost_p = scene.find("PSR J1141-6545 (Newton)");
    check(pulsar > 0 && wd > 0 && ghost_p > 0, "PSR J1141-6545 bodies");
    if (pulsar <= 0 || wd <= 0 || ghost_p <= 0) {
        return;
    }
    const PostKeplerianMotion* orbit = post_keplerian_of(scene, pulsar);
    check(orbit != nullptr, "PSR J1141-6545 post-Keplerian orbit");
    if (!orbit) {
        return;
    }
    const double gm = 3.03867692589e11;
    const glm::dvec3 away = unit_toward(175.27919583, -65.75531667);
    const double x = 1.858915;
    const double e = 0.171876;
    const double t0 = tdb_from_jd_tdb(2454000.4960283);
    const double depth = pulsar_depth_s(scene, pulsar, away, t0);
    check(std::abs(depth - x * (1.0 - e) * std::sin(80.6911 * kDegToRad)) < 2e-5,
          "PSR J1141-6545 Roemer delay at T0 (s)", depth);

    // [BBV08]'s measured advance is GR's for the [VK20] masses...
    const double omdot_gr = gr_omega_dot_deg_per_year(2.28967 - 1.02, 1.02, 0.19765096149, e);
    check(std::abs(omdot_gr / 5.3096 - 1.0) < 3e-4, "PSR J1141-6545 omega_dot is GR's for its masses", omdot_gr);
    const double pbdot_gr = gr_pb_dot(2.28967 - 1.02, 1.02, 0.19765096149, e);
    check(std::abs(pbdot_gr - -0.403e-12) < 0.025e-12, "PSR J1141-6545 Pb_dot within 1 sigma of GR", pbdot_gr);
    std::printf("info: PSR J1141-6545 GR: omega_dot %.5f deg/yr, Pb_dot %.4e\n", omdot_gr, pbdot_gr);

    // ...and the model, run back seven years from [VK20]'s T0, finds [BBV08]'s 1999
    // periastron (an independent timing solution) and its omega there. [VK20] used
    // tempo2, whose default time scale is TCB (TCB - TDB ~ 11-15 s then), [BBV08]
    // TEMPO (TDB): seconds of difference are expected.
    const double t_bbv = tdb_from_jd_tdb(2451370.3545515);
    const PostKeplerianPhase then = post_keplerian_phase(orbit->orbit(), gm, t_bbv);
    const double late_s = wrap_pi(then.mean_anomaly) / kTwoPi * 0.19765096149 * kSecondsPerDay;
    check(std::abs(late_s) < 30.0, "PSR J1141-6545 periastron of 1999 [BBV08] (s off)", late_s);
    const double w_pulsar = wrap_two_pi(then.arg_peri_rad - kPi) * kRadToDeg;
    check(std::abs(w_pulsar - 42.4561) < 0.003, "PSR J1141-6545 omega in 1999 [BBV08] (deg)", w_pulsar);
    std::printf("info: PSR J1141-6545 in 1999: periastron %.2f s off, omega %.4f deg\n", late_s, w_pulsar);

    double t_ref = 0.0;
    check(parse_utc("1999-07-10 20:30", &t_ref), "parse the J1141 reference date");
    const double split =
        glm::length(scene.icrf_state_at(pulsar, t_ref).position - scene.icrf_state_at(ghost_p, t_ref).position);
    check(split < 1e-3, "PSR J1141-6545 ghost starts on the real orbit (km)", split);
    check(std::abs(scene.bodies[static_cast<size_t>(wd)].body_gm_km3_s2 / kSunGmKm3S2 - 1.02) < 1e-6,
          "J1141-6545 white dwarf mass");
}

// The Newtonian ghosts of S2 and the Death Star start where the integrations do.
void test_sgr_a_ghosts()
{
    Scene scene = load_scene_or_die("sgr_a.toml");
    const int s2 = scene.find("S2");
    const int ghost = scene.find("S2 (Newton)");
    const int ds = scene.find("Death Star");
    const int ds_ghost = scene.find("Death Star (Newton)");
    check(s2 > 0 && ghost > 0 && ds > 0 && ds_ghost > 0, "Sgr A* ghost bodies");
    if (s2 <= 0 || ghost <= 0 || ds <= 0 || ds_ghost <= 0) {
        return;
    }
    // The integration starts from the elements at their 2010 apocentre (the
    // loader's "previous_apocentre": half the elements' Keplerian period before t_peri).
    const double a = 125.058 / 1000.0 * 8246.7 * kAuKm;
    const double apo = tdb_from_julian_year(2018.37900) - kPi * std::sqrt(a * a * a / 5.65487707e17);
    const double met = glm::length(scene.icrf_state_at(s2, apo).position - scene.icrf_state_at(ghost, apo).position);
    check(met < 1e3, "S2 and its ghost meet at the 2010 apocentre (km)", met);
    const double later = tdb_from_julian_year(2034.4);
    const double apart =
        glm::length(scene.icrf_state_at(s2, later).position - scene.icrf_state_at(ghost, later).position) / kAuKm;
    check(apart > 1.0, "S2 leaves its ghost by the next pericentre (au)", apart);
    std::printf("info: S2 vs Newton at the 2034 pericentre: %.1f au apart\n", apart);

    // The rolling Death Star restarts from its elements at the requested time.
    const double t = tdb_from_julian_year(2026.0);
    const double start = glm::length(scene.icrf_state_at(ds, t).position - scene.icrf_state_at(ds_ghost, t).position);
    check(start < 1e3, "Death Star starts on its ghost's ellipse (km)", start);
}

void test_sgr_a_scene()
{
    Scene scene = load_scene_or_die("sgr_a.toml");
    const int s2 = scene.find("S2");
    check(scene.bodies[0].kind == BodyKind::BlackHole, "Sgr A* is a black hole");
    // Horizon of a non-spinning 4.261e6 M_sun hole: 2 GM / c^2 ~ 1.26e7 km.
    check(std::abs(scene.bodies[0].equatorial_radius_km / 1.2584e7 - 1.0) < 1e-3, "Sgr A* horizon",
          scene.bodies[0].equatorial_radius_km);

    // Pericentre of the GR-integrated orbit (started at the 2010 apocentre).
    double best_t = 0.0;
    double best_r = 1e300;
    for (double y = 2017.5; y < 2019.5; y += 1e-4) {
        const double r = glm::length(scene.icrf_state_at(s2, tdb_from_julian_year(y)).position);
        if (r < best_r) {
            best_r = r;
            best_t = y;
        }
    }
    check(std::abs(best_t - 2018.379) < 0.01, "S2 pericentre epoch", best_t);
    check(std::abs(best_r / kAuKm - 119.0) < 3.0, "S2 pericentre distance (au)", best_r / kAuKm);
    std::printf("info: S2 pericentre %.4f at %.1f au\n", best_t, best_r / kAuKm);

    // The Death Star's rolling window follows the clock, including big jumps.
    const int hyp = scene.find("Death Star");
    check(hyp >= 0 && scene.bodies[static_cast<size_t>(hyp)].style == SurfaceStyle::DeathStar &&
              scene.bodies[static_cast<size_t>(hyp)].mark_periapsides,
          "Death Star body");
    check(scene.sky.milky_way == "textures/milky_way.jpg" && scene.sky.milky_way_brightness > 0.0, "Milky Way sky");
    check(std::filesystem::exists(std::filesystem::path(ASTRAXIS_ASSET_DIR) / scene.sky.milky_way),
          "Milky Way map present");
    const double m = scene.bodies[0].gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    double t = tdb_from_julian_year(2026.5);
    double r_min = 1e300;
    for (int i = 0; i < 2000; ++i) {
        t += 3600.0;
        r_min = std::min(r_min, glm::length(scene.icrf_state_at(hyp, t).position));
    }
    check(std::abs(r_min / m - 40.0) < 2.0, "Death Star pericentre ~40 M", r_min / m);
    const glm::dvec3 far = scene.icrf_state_at(hyp, tdb_from_julian_year(1900.0)).position;
    check(glm::length(far) > 30.0 * m && glm::length(far) < 400.0 * m, "rolling orbit after a jump", glm::length(far) / m);
}

// Exact Schwarzschild periapsis advance per orbit for turning points r_p, r_a
// (units of M): (du/dphi)^2 = 2 (u1 - u)(u - u2)(u3 - u), u3 = 1/2 - u1 - u2;
// with u = u2 + (u1 - u2) sin^2(theta) the period in phi is
// 4 * integral_0^{pi/2} dtheta / sqrt(2 (u3 - u2 - (u1 - u2) sin^2 theta)).
double schwarzschild_advance(double r_p, double r_a)
{
    const double u1 = 1.0 / r_p;
    const double u2 = 1.0 / r_a;
    const double u3 = 0.5 - u1 - u2;
    const int n = 20000;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        const double th = (i + 0.5) / n * 0.5 * kPi;
        const double s = std::sin(th);
        sum += 1.0 / std::sqrt(2.0 * (u3 - u2 - (u1 - u2) * s * s));
    }
    return 4.0 * sum * (0.5 * kPi / n) - kTwoPi;
}

// The Death Star's periapsides (as marked in the app) advance by the exact
// Schwarzschild amount for its integrated turning points (~18 deg per orbit,
// vs 16.6 deg at 1PN), and the event "at periapsis" really is at one.
void test_death_star_precession()
{
    Scene scene = load_scene_or_die("sgr_a.toml");
    const int ds = scene.find("Death Star");
    const double m = scene.bodies[0].gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);

    const double t0 = tdb_from_julian_year(2026.5);
    scene.update(t0);
    const double span = 10.0 * kSecondsPerDay; // inside the 30-day rolling window
    const std::vector<double> peri = scene.periapsis_times(ds, t0, t0 + span);
    check(peri.size() == 4 || peri.size() == 5, "periapsides in 10 days", static_cast<double>(peri.size()));
    if (peri.size() < 3) {
        return;
    }

    double r_p = 0.0;
    double r_a = 0.0;
    for (double tp : peri) {
        r_p += glm::length(scene.icrf_state_at(ds, tp).position) / m / static_cast<double>(peri.size());
    }
    for (double t = peri[0]; t < peri[1]; t += 60.0) {
        r_a = std::max(r_a, glm::length(scene.icrf_state_at(ds, t).position) / m);
    }
    const double expected = schwarzschild_advance(r_p, r_a);

    double worst = 0.0;
    for (size_t k = 1; k < peri.size(); ++k) {
        const glm::dvec3 a = scene.icrf_state_at(ds, peri[k - 1]).position;
        const glm::dvec3 b = scene.icrf_state_at(ds, peri[k]).position;
        const double advance = std::acos(std::clamp(glm::dot(a, b) / (glm::length(a) * glm::length(b)), -1.0, 1.0));
        worst = std::max(worst, std::abs(advance / expected - 1.0));
        // Prograde: the periapsis turns in the sense of the orbital motion.
        const glm::dvec3 h = glm::cross(a, scene.icrf_state_at(ds, peri[k - 1]).velocity);
        check(glm::dot(glm::cross(a, b), h) > 0.0, "apsidal advance is prograde");
    }
    check(worst < 0.01, "Death Star apsidal advance = exact Schwarzschild (rel. error)", worst);
    std::printf("info: Death Star r_p %.2f M, r_a %.1f M, advance %.2f deg/orbit (exact), 1PN %.2f\n", r_p, r_a,
                expected / kDegToRad, 6.0 * kPi / (0.5 * (r_p + r_a) * (1.0 - std::pow((r_a - r_p) / (r_a + r_p), 2.0))) / kDegToRad);

    // Event "Death Star at periapsis": a periapsis within the following 3 hours.
    for (const SceneEvent& e : scene.events) {
        if (e.focus == ds && e.name.find("periapsis") != std::string::npos) {
            scene.update(e.t_tdb);
            const std::vector<double> next = scene.periapsis_times(ds, e.t_tdb, e.t_tdb + 3.0 * 3600.0);
            check(next.size() == 1, "periapsis event is 0-3 h before a periapsis",
                  next.empty() ? -1.0 : (next[0] - e.t_tdb) / 3600.0);
        }
        if (e.from_orbit_normal >= 0) {
            check(e.from_orbit_normal == ds, "face-on event looks along the Death Star's orbit normal");
        }
    }
}

void test_kerr_null_geodesics()
{
    // Analytic metric derivatives (also used by the shader) vs finite differences.
    const double a = 0.9;
    double worst = 0.0;
    for (const glm::dvec3 x : {glm::dvec3(3.1, -2.2, 1.7), glm::dvec3(-0.4, 5.0, -3.3), glm::dvec3(12.0, 1.0, 0.2)}) {
        const KerrSchild k = kerr_schild(x, a);
        for (int j = 0; j < 3; ++j) {
            glm::dvec3 d(0.0);
            d[j] = 1e-6;
            const KerrSchild kp = kerr_schild(x + d, a);
            const KerrSchild km = kerr_schild(x - d, a);
            worst = std::max(worst, std::abs((kp.f - km.f) / 2e-6 - k.df[j]));
            for (int i = 0; i < 3; ++i) {
                worst = std::max(worst, std::abs((kp.l[i] - km.l[i]) / 2e-6 - k.dl[i][j]));
            }
        }
    }
    check(worst < 1e-6, "Kerr-Schild analytic derivatives", worst);

    // The ray stays null (H = 0) through a close passage.
    Photon ph = make_photon(glm::dvec3(-50.0, 6.0, 3.0), glm::dvec3(1.0, 0.0, 0.0), a);
    double worst_h = 0.0;
    for (int i = 0; i < 3000; ++i) {
        rk4_step(ph, a, ray_step(ph.x, a));
        worst_h = std::max(worst_h, std::abs(hamiltonian(ph, a)));
        if (glm::length(ph.x) > 60.0) {
            break;
        }
    }
    check(worst_h < 1e-6, "null geodesic stays null", worst_h);

    // Schwarzschild: capture iff b < 3 sqrt(3) = 5.196; weak-field deflection 4/b + 15 pi / (4 b^2).
    auto shoot = [](double spin, double b) {
        return trace_ray(glm::dvec3(-2000.0, b, 0.0), glm::dvec3(1.0, 0.0, 0.0), spin, 2100.0, 0.0, 0.0, 20000);
    };
    check(shoot(0.0, 5.15).outcome == RayOutcome::Captured, "b = 5.15 captured");
    check(shoot(0.0, 5.25).outcome == RayOutcome::Escaped, "b = 5.25 escapes");
    const RayResult far = shoot(0.0, 100.0);
    const double deflection = std::acos(std::clamp(far.direction.x, -1.0, 1.0));
    const double expected = 4.0 / 100.0 + 15.0 * kPi / (4.0 * 100.0 * 100.0);
    check(std::abs(deflection - expected) < 2e-3 * expected * 10.0, "weak-field deflection", deflection / expected);

    // Spin: the traced ray is time-reversed, so the real photon moves along -x and
    // has L_z = y. Prograde rays (y > 0) get closer before capture than retrograde ones.
    check(shoot(0.9, 4.0).outcome == RayOutcome::Escaped && shoot(0.9, -4.0).outcome == RayOutcome::Captured,
          "Kerr prograde/retrograde capture asymmetry");

    check(std::abs(isco_radius(0.0) - 6.0) < 1e-12, "ISCO a = 0", isco_radius(0.0));
    check(std::abs(isco_radius(0.9) - 2.3209) < 1e-3, "ISCO a = 0.9", isco_radius(0.9));
    check(std::abs(horizon_radius(0.0) - 2.0) < 1e-12, "horizon a = 0");

    // Doppler: a disk viewed edge-on from -x; the side moving toward the observer is blueshifted.
    const Photon toward = make_photon(glm::dvec3(0.0, 10.0, 0.0), glm::dvec3(1.0, 0.0, 0.0), 0.0);
    const Photon away = make_photon(glm::dvec3(0.0, -10.0, 0.0), glm::dvec3(1.0, 0.0, 0.0), 0.0);
    const double g_toward = disk_redshift(toward, 0.0);
    const double g_away = disk_redshift(away, 0.0);
    check(g_toward > 1.0 && g_away < 1.0, "disk Doppler shift sides", g_toward);
    // Face-on emission: only gravitational + transverse redshift, g = sqrt(1 - 3/r).
    const Photon face_on = make_photon(glm::dvec3(10.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 1.0), 0.0);
    check(std::abs(disk_redshift(face_on, 0.0) - std::sqrt(1.0 - 3.0 / 10.0)) < 1e-9, "face-on redshift sqrt(1-3/r)",
          disk_redshift(face_on, 0.0));
}

// Auto tour: many shots at varying times keep the camera finite and outside
// bodies. The time warp follows the director's shot warp, as in the app.
void test_camera_director(const char* scene_file, double base_warp, double minutes, uint32_t seed,
                          bool* transit_seen = nullptr)
{
    Scene scene = load_scene_or_die(scene_file);
    double t = 0.0;
    parse_utc("2026-10-02 00:00", &t);
    scene.update(t);

    OrbitCamera camera;
    camera.set_up_axis(scene.up_axis());
    camera.focus(0, scene);

    CameraDirector director;
    director.set_time_warp(base_warp);
    director.start(scene, camera, seed);

    bool kinds_seen[static_cast<int>(CameraDirector::ShotKind::Count)] = {};
    const double dt = 1.0 / 30.0;
    for (int frame = 0; frame < static_cast<int>(30.0 * 60.0 * minutes); ++frame) {
        const double warp = director.shot_warp() > 0.0 ? director.shot_warp() : base_warp;
        director.set_time_warp(warp);
        t += dt * warp;
        scene.update(t);
        director.update(dt, scene, camera);
        camera.update(dt, scene);
        kinds_seen[static_cast<int>(director.current_shot())] = true;
        if (director.current_shot() == CameraDirector::ShotKind::ShadowTransit) {
            check(director.shot_warp() > 0.0 && director.shot_warp() <= 600.0, "transit shot slows time");
            if (transit_seen) {
                *transit_seen = true;
            }
        }

        const CameraView view = camera.view(1920, 1080);
        const bool finite = std::isfinite(view.position.x) && std::isfinite(view.position.y) &&
                            std::isfinite(view.position.z);
        check(finite, "camera position finite");
        if (!finite) {
            return;
        }
        for (const Body& body : scene.bodies) {
            const double d = glm::length(view.position - body.world_position);
            if (d < body.equatorial_radius_km) {
                check(false, "camera inside a body", d);
                return;
            }
        }
    }
    int distinct = 0;
    for (bool seen : kinds_seen) {
        distinct += seen ? 1 : 0;
    }
    check(distinct >= 3, "director uses at least 3 shot kinds", distinct);
}

// What the info panel describes a body as orbiting, and its osculating orbit.
void test_orbit_center()
{
    std::printf("orbit center\n");
    auto period_days = [](const Scene& scene, int body) {
        const Scene::OrbitCenter c = scene.orbit_center(body);
        const State s = scene.icrf_state_at(body, scene.time());
        const State h = scene.icrf_state_at(c.body, scene.time());
        OrbitElements el;
        osculating_elements({s.position - h.position, s.velocity - h.velocity}, c.gm, &el);
        return el.period_s / kSecondsPerDay;
    };

    // Around a barycenter: the Moon orbits Earth with the pair's GM, Earth the Sun.
    Scene solar = load_scene_or_die("solar_system.toml");
    solar.update(tdb_from_jd_tdb(2461000.5));
    const int earth = solar.find("Earth");
    const int moon = solar.find("Moon");
    check(solar.orbit_center(moon).body == earth, "the Moon orbits Earth");
    check(std::abs(solar.orbit_center(moon).gm - 403503.235625) < 1e-6, "Earth-Moon GM");
    check(solar.orbit_center(earth).body == solar.find("Sun"), "Earth orbits the Sun");
    const double moon_days = period_days(solar, moon);
    const double earth_days = period_days(solar, earth);
    check(moon_days > 26.5 && moon_days < 28.5, "osculating lunar period (days)", moon_days);
    check(std::abs(earth_days - 365.25) < 3.0, "osculating Earth period (days)", earth_days);

    // Circumbinary planets orbit the barycenter; the secondary star its primary.
    // Periods as in kepler47.toml [OR19]: binary 7.448 d, planet b 49.46 d.
    Scene k47 = load_scene_or_die("kepler47.toml");
    k47.update(tdb_from_jd_tdb(2460000.5));
    const int b = k47.find("Kepler-47 b");
    const int star_b = k47.find("Kepler-47 B");
    check(k47.orbit_center(b).body == k47.find("Kepler-47"), "Kepler-47 b orbits the barycenter");
    check(k47.orbit_center(star_b).body == k47.find("Kepler-47 A"), "Kepler-47 B orbits A");
    const double binary_days = period_days(k47, star_b);
    const double b_days = period_days(k47, b);
    check(std::abs(binary_days - 7.448) < 0.05, "Kepler-47 binary period (days)", binary_days);
    check(std::abs(b_days - 49.5) < 2.0, "Kepler-47 b period (days)", b_days);

    // Unbound: no apoapsis or period; the periapsis is where the state started.
    OrbitElements el;
    check(osculating_elements({glm::dvec3(7000.0, 0.0, 0.0), glm::dvec3(0.0, 12.0, 0.0)}, 398600.0, &el),
          "hyperbolic elements");
    check(el.eccentricity > 1.0 && el.period_s == 0.0 && el.apoapsis_km == 0.0, "hyperbolic orbit", el.eccentricity);
    check(std::abs(el.periapsis_km - 7000.0) < 1e-6, "hyperbolic periapsis (km)", el.periapsis_km);
}

// Masses shown in the info panel (body_gm_km3_s2).
// Atmospheres: per-channel optical depths from the 550 nm values, Earth's
// against Bruneton's Rayleigh formula, Pluto's haze phase function against
// the New Horizons phase ratios, the same parameters in every scene, and the
// loader's checks.
void test_atmospheres()
{
    const Scene scene = load_scene_or_die("solar_system.toml");
    auto atmosphere = [&](const char* name) -> const Atmosphere& {
        return scene.bodies[static_cast<size_t>(scene.find(name))].atmosphere;
    };
    for (const char* name : {"Earth", "Venus", "Mars", "Titan", "Pluto", "Triton"}) {
        check(atmosphere(name).enabled, name);
    }
    check(!atmosphere("Moon").enabled && !atmosphere("Jupiter").enabled, "no atmosphere on the Moon and Jupiter");

    // Bruneton 2017 (demo.cc): beta_R = 1.24062e-6 lambda^-4 m^-1 (lambda in um), H_R = 8 km.
    const Atmosphere& earth = atmosphere("Earth");
    double worst = 0.0;
    for (int c = 0; c < 3; ++c) {
        const double um = kChannelWavelengthsNm[c] / 1000.0;
        const double expected = 1.24062e-6 * std::pow(um, -4.0) * 8000.0;
        worst = std::max(worst, std::abs(earth.rayleigh_depth[c] / expected - 1.0));
    }
    check(worst < 1e-3, "Earth Rayleigh depths match Bruneton's formula", worst);
    check(std::abs(earth.haze_depth.r - earth.haze_depth.b) < 1e-12, "Earth aerosols are gray (alpha = 0)");

    // Titan: opacity ~ lambda^-2.34 above 80 km (Tomasko et al. 2008).
    const Atmosphere& titan = atmosphere("Titan");
    const double ratio = titan.haze_depth.b / titan.haze_depth.r;
    check(std::abs(ratio - std::pow(440.0 / 680.0, -2.34)) < 1e-9, "Titan haze Angstrom exponent", ratio);
    // Triton: scattering optical depth ~ lambda^-2 (Ohno et al. 2021), haze to ~30 km inside the shell.
    const Atmosphere& triton = atmosphere("Triton");
    check(std::abs(triton.haze_depth.b / triton.haze_depth.r - std::pow(440.0 / 680.0, -2.0)) < 1e-9 &&
              triton.height_km >= 30.0,
          "Triton haze");
    check(titan.deck_altitude_km == 80.0 && atmosphere("Venus").deck_altitude_km == 74.0, "deck altitudes");

    // Pluto: Henyey-Greenstein ratio between 167 and 20 deg phase (scattering angles
    // 13 and 160 deg) within the haze I/F ratios of Cheng et al. 2017, Table 4
    // (peak: 0.3 / 0.02 = 15; at 45 km: 0.15 / 0.004 = 37).
    const double g = atmosphere("Pluto").haze_g;
    auto hg = [g](double theta_deg) {
        const double d = 1.0 + g * g - 2.0 * g * std::cos(theta_deg * kDegToRad);
        return (1.0 - g * g) / (d * std::sqrt(d));
    };
    const double pluto_ratio = hg(13.0) / hg(160.0);
    check(pluto_ratio > 15.0 && pluto_ratio < 37.0, "Pluto haze forward/back ratio", pluto_ratio);

    // Scenes that repeat a body carry the same atmosphere.
    for (const char* file : {"earth_moon.toml", "jwst.toml", "mercury.toml", "parker.toml", "saturn.toml"}) {
        const Scene other = load_scene_or_die(file);
        for (const Body& b : other.bodies) {
            const int k = scene.find(b.name);
            if (k < 0 || !scene.bodies[static_cast<size_t>(k)].atmosphere.enabled) {
                continue;
            }
            const Atmosphere& a = b.atmosphere;
            const Atmosphere& ref = scene.bodies[static_cast<size_t>(k)].atmosphere;
            check(a.enabled && a.height_km == ref.height_km && a.deck_altitude_km == ref.deck_altitude_km &&
                      a.rayleigh_depth == ref.rayleigh_depth && a.haze_depth == ref.haze_depth &&
                      a.haze_albedo == ref.haze_albedo && a.haze_g == ref.haze_g &&
                      a.haze_scale_height_km == ref.haze_scale_height_km &&
                      b.texture == scene.bodies[static_cast<size_t>(k)].texture,
                  (std::string(file) + ": " + b.name + " atmosphere as in solar_system.toml").c_str());
        }
    }

    std::string error;
    Scene bad;
    const char* head = R"(name = 'x'
[[bodies]]
)";
    auto bad_scene = [&](const char* body) {
        return !load_scene_string(std::string(head) + body, "bad_atmosphere", bad, &error);
    };
    check(bad_scene(R"(name = 'S'
kind = 'star'
radii_km = [1.0]
[bodies.atmosphere]
height_km = 1.0
)"),
          "rejects an atmosphere on a star");
    check(bad_scene(R"(name = 'P'
radii_km = [1.0]
[bodies.atmosphere]
height_km = 1.0
haze_depth = 0.1
)") && error.find("haze_scale_height_km") != std::string::npos,
          "requires a haze scale height");
    check(bad_scene(R"(name = 'P'
radii_km = [1.0]
[bodies.atmosphere]
height_km = 1.0
haze_albedo = [0.5, 1.5, 0.5]
)"),
          "rejects a haze albedo above 1");
}

// Plumes: vents from the IAU Gazetteer (west longitudes converted to east),
// the tiger stripes 30-36 km apart and perpendicular to the line through their
// centers, surface points on the ellipsoid, the same plumes in every scene,
// and the loader's checks.
void test_plumes()
{
    const Scene scene = load_scene_or_die("solar_system.toml");
    auto body = [&](const char* name) -> const Body& { return scene.bodies[static_cast<size_t>(scene.find(name))]; };
    auto plume = [&](const Body& b, const char* name) -> const Plume* {
        for (const Plume& p : b.plumes) {
            if (p.name == name) {
                return &p;
            }
        }
        return nullptr;
    };

    const Plume* pele = plume(body("Io"), "Pele");
    const Plume* tvashtar = plume(body("Io"), "Tvashtar");
    check(pele && pele->type == Plume::Type::Umbrella && std::abs(pele->lat_lon_deg.y - (360.0 - 255.28)) < 1e-9 &&
              pele->height_km == 300.0 && pele->radius_km == 600.0,
          "Pele: 255.28 W, 300 km high, 1,200 km across");
    check(tvashtar && std::abs(tvashtar->lat_lon_deg.y - (360.0 - 123.53)) < 1e-9 &&
              tvashtar->radius_km == 2.0 * tvashtar->height_km,
          "Tvashtar: 123.53 W, ballistic reach 2 H");

    // Tiger stripes: midpoints near the Gazetteer centers, lengths, spacing and orientation.
    const Body& enceladus = body("Enceladus");
    const char* stripes[4] = {"Alexandria Sulcus", "Cairo Sulcus", "Baghdad Sulcus", "Damascus Sulcus"};
    const double centers[4][2] = {{-75.63, 360.0 - 137.56}, {-81.62, 360.0 - 154.48},
                                  {-86.91, 360.0 - 230.54}, {-80.59, 360.0 - 285.87}};
    const double lengths[4] = {111.0, 165.0, 176.0, 125.0};
    glm::dvec3 mid[4];
    glm::dvec3 along[4];
    bool ok = true;
    for (int k = 0; k < 4; ++k) {
        const Plume* p = plume(enceladus, stripes[k]);
        if (!p || p->type != Plume::Type::Jets) {
            ok = false;
            continue;
        }
        const glm::dvec3 a = body_surface_point(enceladus, p->lat_lon_deg.x, p->lat_lon_deg.y);
        const glm::dvec3 b = body_surface_point(enceladus, p->end_lat_lon_deg.x, p->end_lat_lon_deg.y);
        // The chord's midpoint, raised to the surface (the sagitta is up to 15 km).
        mid[k] = glm::normalize(a + b) * glm::length(body_surface_point(enceladus, centers[k][0], centers[k][1]));
        along[k] = glm::normalize(b - a);
        const glm::dvec3 c = body_surface_point(enceladus, centers[k][0], centers[k][1]);
        ok = ok && glm::length(mid[k] - c) < 3.0 && std::abs(glm::length(b - a) / lengths[k] - 1.0) < 0.1;
    }
    check(ok, "tiger stripes centered on the Gazetteer centers, with their lengths");
    double worst_spacing = 0.0;
    double worst_angle = 0.0;
    for (int k = 0; k + 1 < 4; ++k) {
        const glm::dvec3 step = mid[k + 1] - mid[k];
        worst_spacing = std::max(worst_spacing, std::abs(glm::length(step) - 33.0));
        worst_angle = std::max(worst_angle, std::abs(glm::dot(glm::normalize(step), along[k])));
    }
    check(worst_spacing < 5.0, "tiger stripes ~35 km apart (km off)", worst_spacing);
    check(worst_angle < 0.2, "tiger stripes perpendicular to the line of centers (cos)", worst_angle);

    const Body& triton = body("Triton");
    const Plume* hili = plume(triton, "Hili");
    const Plume* mahilani = plume(triton, "Mahilani");
    check(hili && mahilani && hili->type == Plume::Type::Geyser && hili->height_km == 8.0 &&
              hili->tail_azimuth_deg == 270.0 && mahilani->lat_lon_deg == glm::dvec2(-50.5, 359.5),
          "Triton's geysers: 8 km, tails to the west");

    // Surface points lie on the ellipsoid; normals are unit and outward.
    double worst_surface = 0.0;
    for (const Body* b : {&body("Io"), &enceladus, &triton}) {
        for (double lat = -80.0; lat <= 80.0; lat += 40.0) {
            for (double lon = 0.0; lon < 360.0; lon += 60.0) {
                const glm::dvec3 p = body_surface_point(*b, lat, lon);
                const glm::dvec3 q = p / glm::dvec3(b->equatorial_radius_km, b->equatorial_radius_b_km, b->polar_radius_km);
                const glm::dvec3 n = body_surface_normal(*b, p);
                worst_surface = std::max({worst_surface, std::abs(glm::dot(q, q) - 1.0),
                                          std::abs(glm::length(n) - 1.0), glm::dot(n, glm::normalize(p)) > 0.99 ? 0.0 : 1.0});
            }
        }
    }
    check(worst_surface < 1e-12, "surface points and normals", worst_surface);

    for (const char* file : {"jupiter.toml", "saturn.toml"}) {
        const Scene other = load_scene_or_die(file);
        for (const Body& b : other.bodies) {
            const int k = scene.find(b.name);
            if (k < 0 || scene.bodies[static_cast<size_t>(k)].plumes.empty()) {
                continue;
            }
            const std::vector<Plume>& ref = scene.bodies[static_cast<size_t>(k)].plumes;
            bool same = b.plumes.size() == ref.size();
            for (size_t i = 0; same && i < ref.size(); ++i) {
                const Plume& x = b.plumes[i];
                const Plume& y = ref[i];
                same = x.name == y.name && x.type == y.type && x.lat_lon_deg == y.lat_lon_deg &&
                       x.end_lat_lon_deg == y.end_lat_lon_deg && x.count == y.count && x.height_km == y.height_km &&
                       x.radius_km == y.radius_km && x.spread_deg == y.spread_deg &&
                       x.optical_depth == y.optical_depth && x.albedo == y.albedo && x.g == y.g;
            }
            check(same, (std::string(file) + ": " + b.name + " plumes as in solar_system.toml").c_str());
        }
    }

    std::string error;
    Scene bad;
    const char* head = R"(name = 'x'
[[bodies]]
name = 'P'
radii_km = [100.0]
[[bodies.plumes]]
name = 'v'
lat_lon_deg = [0.0, 0.0]
height_km = 10.0
optical_depth = 0.1
)";
    auto bad_scene = [&](const char* rest) {
        return !load_scene_string(std::string(head) + rest, "bad_plume", bad, &error);
    };
    check(bad_scene("type = 'fountain'\n"), "rejects an unknown plume type");
    check(bad_scene("type = 'jets'\nend_lat_lon_deg = [1.0, 1.0]\ncount = 13\nspread_deg = 5.0\n"),
          "rejects too many jets");
    check(bad_scene("type = 'geyser'\nradius_km = 0.5\ntail_km = 10.0\n") &&
              error.find("tail_azimuth_deg") != std::string::npos,
          "requires a geyser's tail azimuth");
    check(!bad_scene("type = 'umbrella'\n") && bad.bodies[0].plumes[0].radius_km == 20.0,
          "umbrella radius defaults to 2 H");
}

void test_body_masses()
{
    // The planets' own GM ([SPHY], from the satellite ephemerides) plus their
    // moons add up to the DE440 system GM: Jupiter and Saturn leave out small
    // moons (2-2.5 km^3/s^2), Pluto and Charon (PLU060) differ from DE440 within
    // their uncertainties (0.4, 0.3).
    Scene solar = load_scene_or_die("solar_system.toml");
    struct System {
        const char* holder; // body whose gm_km3_s2 is the system GM
        std::vector<const char*> members;
        double tolerance;
    };
    const System systems[] = {
        {"Earth-Moon barycenter", {"Earth", "Moon"}, 1e-6},
        {"Mars", {"Mars", "Phobos", "Deimos"}, 0.003},
        {"Jupiter", {"Jupiter", "Io", "Europa", "Ganymede", "Callisto"}, 3.0},
        {"Saturn", {"Saturn", "Mimas", "Enceladus", "Tethys", "Dione", "Rhea", "Titan", "Iapetus"}, 3.0},
        {"Uranus", {"Uranus", "Miranda", "Ariel", "Umbriel", "Titania", "Oberon"}, 5.0},
        {"Neptune", {"Neptune", "Triton"}, 10.0},
        {"Pluto system", {"Pluto", "Charon"}, 0.5},
    };
    for (const System& s : systems) {
        double sum = 0.0;
        for (const char* m : s.members) {
            sum += solar.bodies[static_cast<size_t>(solar.find(m))].body_gm_km3_s2;
        }
        const double system_gm = solar.bodies[static_cast<size_t>(solar.find(s.holder))].gm_km3_s2;
        check(std::abs(sum - system_gm) < s.tolerance, s.holder, sum - system_gm);
    }

    // Split pairs: the body's share of the orbit's total GM.
    Scene tic = load_scene_or_die("tic168789840.toml");
    const double a1 = tic.bodies[static_cast<size_t>(tic.find("TIC 168789840 A1"))].body_gm_km3_s2 / kSunGmKm3S2;
    check(std::abs(a1 - 1.25) < 1e-6, "TIC 168789840 A1 mass (M_sun) [PO21 Table 7]", a1);
    Scene psr = load_scene_or_die("psr_b1620.toml");
    const double planet = psr.bodies[static_cast<size_t>(psr.find("PSR B1620-26 b"))].body_gm_km3_s2 / kSunGmKm3S2;
    check(std::abs(planet - 2.3865e-3) < 1e-6, "PSR B1620-26 b mass (M_sun) [SI03]", planet);
    // N-body members take their N-body GM.
    Scene k47 = load_scene_or_die("kepler47.toml");
    const double k47a = k47.bodies[static_cast<size_t>(k47.find("Kepler-47 A"))].body_gm_km3_s2 / kSunGmKm3S2;
    check(std::abs(k47a - 0.9573912601) < 1e-6, "Kepler-47 A mass (M_sun) [OR19 Table 7]", k47a);
}

// Moons orbit their planet, and the primary of a barycenter stands for its
// system; satellites fade out as their orbits shrink on screen.
void test_satellite_fades()
{
    Scene scene;
    auto add = [&](const char* name, int parent, BodyKind kind, const glm::dvec3& pos) {
        Body b;
        b.name = name;
        b.parent = parent;
        b.kind = kind;
        b.world_position = pos;
        scene.bodies.push_back(std::move(b));
    };
    const glm::dvec3 planet(1e8, 0.0, 0.0);
    const glm::dvec3 bary(-5e9, 0.0, 0.0);
    add("Sun", -1, BodyKind::Star, glm::dvec3(0.0));
    add("Planet", 0, BodyKind::Planet, planet);
    add("Moon", 1, BodyKind::Planet, planet + glm::dvec3(1e6, 0.0, 0.0));
    add("System", 0, BodyKind::Barycenter, bary);
    add("Primary", 3, BodyKind::Planet, bary + glm::dvec3(2000.0, 0.0, 0.0));
    add("Secondary", 3, BodyKind::Planet, bary - glm::dvec3(17000.0, 0.0, 0.0));

    const int hosts[] = {-1, 0, 1, 0, 0, 4};
    bool hosts_ok = true;
    for (int i = 0; i < 6; ++i) {
        hosts_ok = hosts_ok && scene.satellite_host(i) == hosts[i];
    }
    check(hosts_ok, "satellite hosts (barycenter primary stands for its system)");

    std::vector<float> fades;
    // 1e7 km from the planet: the moon's orbit spans 100 px.
    satellite_fades(scene, planet + glm::dvec3(0.0, 0.0, 1e7), 1000.0, 10.0, 28.0, fades);
    check(fades[2] == 1.0f && fades[4] == 1.0f && fades[5] == 0.0f, "fades near the planet");
    // Moon orbit at 19 px, halfway between hide and show.
    satellite_fades(scene, planet + glm::dvec3(0.0, 0.0, 1e6 * 1000.0 / 19.0), 1000.0, 10.0, 28.0, fades);
    check(std::abs(fades[2] - 0.5f) < 1e-3f, "fade halfway", fades[2]);
    // Far out: the planet's orbit spans 1 px, so it and its moon are hidden; the
    // dwarf-planet system still shows as its primary.
    satellite_fades(scene, glm::dvec3(0.0, 0.0, 1e11), 1000.0, 10.0, 28.0, fades);
    check(fades[0] == 1.0f && fades[1] == 0.0f && fades[2] == 0.0f && fades[4] == 1.0f && fades[5] == 0.0f,
          "fades far out");
}

// Simulation: scene loading, events, the clock and the tour's warp, as a host drives them.
void test_simulation()
{
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error), "simulation loads jupiter.toml");
    const size_t bodies = sim.scene().bodies.size();
    check(bodies > 0 && !sim.scene().events.empty(), "jupiter scene has bodies and events");
    check(sim.clock().warp == sim.scene().view.warp, "load sets the scene's warp", sim.clock().warp);
    check(!sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/no_such_scene.toml", &error) && !error.empty() &&
              sim.scene().bodies.size() == bodies,
          "failed load keeps the scene");

    sim.jump_to_event(0);
    const SceneEvent& event = sim.scene().events[0];
    check(sim.clock().t_tdb == event.t_tdb, "jump_to_event sets the time", sim.clock().t_tdb - event.t_tdb);
    check(sim.scene().time() == event.t_tdb, "jump_to_event updates the scene");

    SimClock& clock = sim.clock();
    clock.paused = false;
    clock.reverse = false;
    clock.warp = 1000.0;
    const double t0 = clock.t_tdb;
    const double anim0 = sim.animation_time();
    sim.update(0.05);
    check(std::fabs(clock.t_tdb - t0 - 50.0) < 1e-6, "update advances by real dt x warp", clock.t_tdb - t0);
    check(sim.scene().time() == clock.t_tdb, "update moves the scene to the clock");
    check(std::fabs(sim.animation_time() - anim0 - 0.05) < 1e-9, "animation time runs with the clock");
    clock.paused = true;
    sim.update(0.05);
    check(sim.animation_time() == anim0 + 0.05, "animation time stops while paused");
    clock.paused = false;

    OutputView view;
    sim.compute_view(1920, 1080, 1.0, view);
    check(view.body_fades.size() == bodies, "a fade per body");
    check(view.body_fades[static_cast<size_t>(sim.camera().target())] == 1.0f, "the focus always shows");
    check(view.camera.viewport.x == 1920.0f && view.camera.viewport.y == 1080.0f, "view has the output size");
    // Each output has its own view: a narrow one sees the same camera with another projection.
    OutputView narrow;
    sim.compute_view(400, 1080, 2.0, narrow);
    check(narrow.camera.position == view.camera.position && narrow.camera.viewport.x == 400.0f &&
              narrow.content_scale == 2.0,
          "per-output views");

    // The tour may impose its own warp; stopping it restores the user's.
    sim.start_tour();
    check(sim.touring(), "tour starts");
    bool tour_warped = false;
    for (int k = 0; k < 20 * 60 * 10; ++k) {
        sim.update(0.1);
        tour_warped = tour_warped || clock.warp != 1000.0;
    }
    check(tour_warped, "the jupiter tour sets a shot warp in 20 minutes");
    check(sim.touring(), "tour keeps running");
    sim.stop_tour();
    check(!sim.touring() && clock.warp == 1000.0, "stopping the tour restores the warp", clock.warp);

    // A running tour goes on in the next scene; a focus change ends it.
    sim.start_tour();
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/saturn.toml", &error) && sim.touring(),
          "tour survives a scene change");
    sim.compute_view(1920, 1080, 1.0, view);
    check(view.body_fades.size() == sim.scene().bodies.size(), "fades resized with the scene");
    sim.set_focus(0);
    check(!sim.touring(), "set_focus stops the tour");
}

// Every scene file, loaded once (by file stem).
const std::vector<std::pair<std::string, std::shared_ptr<const Scene>>>& all_scenes()
{
    static std::vector<std::pair<std::string, std::shared_ptr<const Scene>>> scenes;
    if (scenes.empty()) {
        for (const auto& file : std::filesystem::directory_iterator(ASTRAXIS_ASSET_DIR "/scenes")) {
            scenes.emplace_back(file.path().stem().string(),
                                std::make_shared<const Scene>(load_scene_or_die(file.path().filename().string().c_str())));
        }
    }
    return scenes;
}

// A body defined in several scenes (orbits are copied between scene files) moves
// the same way in each, relative to the same parent: inside its tables and in
// the fallbacks outside them.
void test_cross_scene_orbits()
{
    struct Entry {
        std::string scene;
        std::shared_ptr<const Scene> data;
        int body;
    };
    std::map<std::string, std::vector<Entry>> by_name; // "name / parent" -> entries
    for (const auto& [name, scene] : all_scenes()) {
        for (size_t i = 0; i < scene->bodies.size(); ++i) {
            const Body& b = scene->bodies[i];
            if (!b.motion || b.parent < 0) {
                continue;
            }
            const std::string key = b.name + " / " + scene->bodies[static_cast<size_t>(b.parent)].name;
            by_name[key].push_back({name, scene, static_cast<int>(i)});
        }
    }
    for (const auto& [key, entries] : by_name) {
        for (size_t k = 1; k < entries.size(); ++k) {
            const MotionSource& a = *entries[0].data->bodies[static_cast<size_t>(entries[0].body)].motion;
            const MotionSource& b = *entries[k].data->bodies[static_cast<size_t>(entries[k].body)].motion;
            double worst = 0.0;
            for (int year = 1900; year <= 2150; ++year) {
                for (int month = 0; month < 12; month += 5) {
                    const double t = tdb_from_unix_utc(0.0) + ((year - 1970) * 365.25 + month * 30.44) * kSecondsPerDay;
                    if (!a.valid_at(t) || !b.valid_at(t)) {
                        continue;
                    }
                    const glm::dvec3 pa = a.eval(t).position;
                    worst = std::max(worst, glm::length(b.eval(t).position - pa) / glm::length(pa));
                }
            }
            if (worst > 1e-3) {
                std::printf("  %s: %s vs %s differ by %.3g of the distance\n", key.c_str(), entries[0].scene.c_str(),
                            entries[k].scene.c_str(), worst);
            }
            check(worst <= 1e-3, ("same orbit in every scene: " + key).c_str(), worst);
        }
    }
}

// Where a table hands over to its fallback, the position fades between the two
// over fallback_blend_days, and the velocity leaves out the fade's own rate
// (EphemerisMotion::eval). Both are only good while the fallback stays close to
// the table there: the left-out rate, |dP/dt - V|, must stay a small part of
// the body's speed.
void test_ephemeris_blends()
{
    constexpr double kMaxShare = 0.1;
    for (const auto& [name, scene] : all_scenes()) {
        for (const Body& body : scene->bodies) {
            const MotionSource* m = body.motion.get();
            if (const auto* scaled = dynamic_cast<const ScaledMotion*>(m)) {
                m = &scaled->inner();
            }
            for (const auto* eph = dynamic_cast<const EphemerisMotion*>(m); eph && eph->fallback();
                 eph = dynamic_cast<const EphemerisMotion*>(eph->fallback())) {
                const EphemerisTable& table = eph->table();
                std::vector<double> edges;
                for (const EphemerisTable::Knot& k : table.knots()) {
                    double a = 0.0;
                    double b = 0.0;
                    table.segment(k.t, &a, &b);
                    if (edges.empty() || edges.back() != b) {
                        edges.push_back(a);
                        edges.push_back(b);
                    }
                }
                double worst = 0.0;
                for (double edge : edges) {
                    for (double t = edge - 31.0 * kSecondsPerDay; t <= edge + 31.0 * kSecondsPerDay; t += 3600.0) {
                        double a0 = 0.0, b0 = 0.0, a1 = 0.0, b1 = 0.0;
                        if (!table.covers(t - 60.0) || !table.covers(t + 60.0)) {
                            continue;
                        }
                        table.segment(t - 60.0, &a0, &b0);
                        table.segment(t + 60.0, &a1, &b1);
                        if (a0 != a1 || b0 != b1) {
                            continue;
                        }
                        const State s = eph->eval(t);
                        const glm::dvec3 rate = (eph->eval(t + 60.0).position - eph->eval(t - 60.0).position) / 120.0;
                        worst = std::max(worst, glm::length(rate - s.velocity) / glm::length(s.velocity));
                    }
                }
                if (worst > 0.02) {
                    std::printf("info: %s / %s: fallback blend leaves out %.1f%% of the speed\n", name.c_str(),
                                body.name.c_str(), 100.0 * worst);
                }
                check(worst < kMaxShare, (name + " / " + body.name + ": fallback close to the table at its edges").c_str(),
                      worst);
            }
        }
    }
}

// Label layout: priority order, easing, no overlapping labels once settled, picking.
void test_label_layout()
{
    Simulation sim;
    std::string error;
    if (!sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error)) {
        check(false, "label test loads jupiter.toml");
        return;
    }
    OutputView view;
    sim.compute_view(1280, 720, 1.0, view);
    const glm::vec2 screen(1280.0f, 720.0f);
    const float font = 13.0f;
    auto measure = [&](const std::string& text) { return glm::vec2(7.0f * static_cast<float>(text.size()), font); };
    const int focus = sim.camera().target();

    LabelLayout labels;
    labels.update(sim.scene(), view, focus, screen, font, measure, 0.0f);
    check(!labels.bodies().empty() && labels.bodies().front().body == focus, "the focus is placed first");
    bool all_hidden = true;
    for (const LabelLayout::BodyMark& mark : labels.bodies()) {
        all_hidden = all_hidden && mark.text_alpha == 0.0f;
    }
    check(all_hidden, "labels start transparent and ease in");

    for (int k = 0; k < 30; ++k) {
        labels.update(sim.scene(), view, focus, screen, font, measure, 1.0f / 30.0f);
    }
    const LabelLayout::BodyMark& first = labels.bodies().front();
    check(std::fabs(first.text_alpha - first.fade) < 1e-5f, "the focus label is fully shown after 1 s",
          first.text_alpha);
    std::vector<glm::vec4> shown;
    bool overlap = false;
    for (const LabelLayout::BodyMark& mark : labels.bodies()) {
        if (mark.text_alpha < 0.99f * mark.fade || mark.fade <= 0.05f) {
            continue;
        }
        const glm::vec2 size = measure(sim.scene().bodies[static_cast<size_t>(mark.body)].name);
        const glm::vec4 r(mark.text_position, mark.text_position + size);
        for (const glm::vec4& o : shown) {
            overlap = overlap || (r.x < o.z && o.x < r.z && r.y < o.w && o.y < r.w);
        }
        shown.push_back(r);
    }
    check(shown.size() > 1 && !overlap, "settled labels do not overlap", static_cast<double>(shown.size()));

    // (A moon over the focus's disk may win: the last mark under the point does.)
    check(labels.pick(first.position, 10.0f) >= 0, "picking at the focus finds a body");
    check(labels.pick(glm::vec2(-1000.0f, -1000.0f), 10.0f) == -1, "picking off screen finds nothing");

    labels.reset();
    labels.update(sim.scene(), view, focus, screen, font, measure, 0.0f);
    check(labels.bodies().front().text_alpha == 0.0f, "reset forgets the eased opacities");
}

// Auto tour around Sgr A*: black hole close-ups happen, and no shot or transition
// takes the camera into the hole's strong-field region.
void test_camera_director_black_hole()
{
    Scene scene = load_scene_or_die("sgr_a.toml");
    double t = 0.0;
    parse_utc("2026-10-02 00:00", &t);
    scene.update(t);
    const Body& hole = scene.bodies[0];
    const double m_km = hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    const int ds = scene.find("Death Star");

    OrbitCamera camera;
    camera.set_up_axis(hole.pole);
    camera.focus(0, scene);
    CameraDirector director;
    // 1 day/s: too fast to follow the Death Star (2.35 d orbit), so the face-on
    // shot must slow time down (as the app does with shot_warp()). The scene's
    // own 1 month/s would integrate ~100 years of its orbit here.
    const double base_warp = kSecondsPerDay;
    director.set_time_warp(base_warp);
    director.start(scene, camera, 5u);

    bool close_up = false;
    bool face_on = false;
    bool slowed = true;
    double min_distance = 1e300;
    const double dt = 1.0 / 30.0;
    for (int frame = 0; frame < 30 * 60 * 20; ++frame) {
        const double warp = director.shot_warp() > 0.0 ? director.shot_warp() : base_warp;
        director.set_time_warp(warp);
        t += dt * warp;
        scene.update(t);
        director.update(dt, scene, camera);
        camera.update(dt, scene);
        if (director.current_shot() == CameraDirector::ShotKind::BlackHoleCloseUp) {
            close_up = true;
        }
        if (director.current_shot() == CameraDirector::ShotKind::OrbitTopDown) {
            face_on = true;
            double r_p = 0.0;
            double r_a = 0.0;
            double period = 0.0;
            if (scene.osculating_apsides(ds, &r_p, &r_a, &period) && director.shot_warp() * 10.0 > period) {
                slowed = false;
            }
        }
        min_distance = std::min(min_distance, glm::length(camera.view(1920, 1080).position - hole.world_position));
    }
    check(close_up, "director picks a black hole close-up");
    check(face_on, "director picks a face-on orbit shot");
    check(slowed, "face-on orbit shot slows time to >= 10 s per orbit");
    check(min_distance > 11.0 * m_km, "camera stays beyond 11 M of the hole", min_distance / m_km);
}

} // namespace

int main()
{
    test_scene_loader();
    test_star_catalog_and_color();
    test_kepler();
    test_calendar();
    test_galilean_rates();
    test_jupiter_heliocentric();
    test_io_shadow_transits();
    test_tidal_locking();
    test_kepler_propagation();
    test_ephemeris();
    test_solar_system_scene();
    test_solar_system_bodies();
    test_solar_system_missions();
    test_solar_system_belts();
    test_arrokoth_shape();
    test_mars_moon_shapes();
    test_vesta_shape();
    test_pallas_shape();
    test_jupiter_saturn_small_moon_shapes();
    test_halley_shape();
    test_comet_tails();
    test_dust_tail();
    test_parker_scene();
    test_jupiter_missions();
    test_earth_moon_scene();
    test_saturn_scene();
    test_mercury_scene();
    test_isee3_scene();
    test_halley_armada();
    test_jwst_scene();
    test_galactic_frame();
    test_visual_orbit_convention();
    test_kerr_conservation();
    test_kerr_schwarzschild_precession();
    test_nbody();
    test_alpha_centauri_scene();
    test_transit_orbit_convention();
    test_trappist1_scene();
    test_kepler223_scene();
    test_kepler47_scene();
    test_tic168789840_scene();
    test_kepler64_scene();
    test_psr_b1620_scene();
    test_post_keplerian();
    test_psr_b1913_scene();
    test_psr_j1141_scene();
    test_earth_moon_rotation();
    test_kerr_null_geodesics();
    test_sgr_a_scene();
    test_death_star_precession();
    test_sgr_a_ghosts();
    test_satellite_fades();
    test_cross_scene_orbits();
    test_ephemeris_blends();
    test_simulation();
    test_label_layout();
    test_orbit_center();
    test_body_masses();
    test_atmospheres();
    test_plumes();
    for (uint32_t seed : {42u, 7u, 2026u, 99u, 12345u}) {
        test_camera_director("jupiter.toml", 3600.0, 20.0, seed);
    }
    for (uint32_t seed : {1u, 2u, 3u}) {
        test_camera_director("solar_system.toml", 2629800.0, 10.0, seed);
    }
    test_camera_director_black_hole();

    if (g_failures == 0) {
        std::printf("All tests passed.\n");
        return 0;
    }
    std::printf("%d failure(s).\n", g_failures);
    return 1;
}
