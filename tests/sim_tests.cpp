// Minimal self-checking tests for core / ephem / scene (no framework).

#include "core/color.hpp"
#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/ephemeris.hpp"
#include "ephem/kepler.hpp"
#include "ephem/kerr_null.hpp"
#include "ephem/kerr_orbit.hpp"
#include "ephem/nbody.hpp"
#include "ephem/visual_orbit.hpp"
#include "ephem/mean_element_orbit.hpp"
#include "scene/camera.hpp"
#include "scene/camera_director.hpp"
#include "scene/scene_loader.hpp"
#include "scene/star_catalog.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
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
    // osculating elements from 2026 (Pluto: 2016) drift, most for Ceres, which Jupiter
    // perturbs; the Moon's mean elements are off by degrees.
    const double t_end = tdb_from_jd_tdb(2488067.5); // 2099-12-30, within the tables
    struct Join {
        const char* body;
        double max_deg;
    };
    for (const Join& j : {Join{"Ceres", 6.0}, Join{"Eris", 0.1}, Join{"Sedna", 0.2}, Join{"Pluto system", 0.6},
                          Join{"Charon", 0.1}, Join{"Moon", 12.0}}) {
        const int i = scene.find(j.body);
        const int parent = scene.bodies[static_cast<size_t>(i)].parent;
        const auto* scaled = dynamic_cast<const ScaledMotion*>(scene.bodies[static_cast<size_t>(i)].motion.get());
        const MotionSource& m = scaled ? scaled->inner() : *scene.bodies[static_cast<size_t>(i)].motion;
        const auto* eph = dynamic_cast<const EphemerisMotion*>(&m);
        check(eph != nullptr && eph->fallback() != nullptr && parent >= 0, "baked table with a fallback");
        if (!eph || !eph->fallback()) {
            continue;
        }
        const double a = angle_deg(m.eval(t_end).position, eph->fallback()->eval(t_end).position);
        check(a < j.max_deg, j.body, a);
        std::printf("info: %s fallback at the end of its table: %.3f deg\n", j.body, a);
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

// Solar System scene: the asteroid and Kuiper belts (two-body orbits from SBDB).
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

    // Signed volume (positive: counter-clockwise seen from outside, a closed surface).
    double volume = 0.0;
    for (size_t k = 0; k < shape.indices.size(); k += 3) {
        const glm::dvec3 a(shape.positions[shape.indices[k]]);
        const glm::dvec3 b(shape.positions[shape.indices[k + 1]]);
        const glm::dvec3 c(shape.positions[shape.indices[k + 2]]);
        volume += glm::dot(a, glm::cross(b, c)) / 6.0;
    }
    const double d_equal = std::cbrt(6.0 * volume / kPi);
    check(std::abs(d_equal - 19.896) < 0.01, "Arrokoth equal-volume diameter (km)", d_equal);
    const auto [min_albedo, max_albedo] = std::minmax_element(shape.albedo.begin(), shape.albedo.end());
    check(*min_albedo > 0.5f && *max_albedo < 2.5f, "Arrokoth relative albedo range", *max_albedo);

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
    // Horizons tables wherever these have data. Mimas librates by +-44 deg.
    const char* moons[] = {"Mimas", "Enceladus", "Tethys", "Dione", "Rhea", "Titan", "Iapetus"};
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
        const bool mimas = std::string(name) == "Mimas";
        check(worst < (mimas ? 45.0 : 1.5), "Saturn moon mean elements near Horizons (deg)", worst);
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
    test_parker_scene();
    test_jupiter_missions();
    test_earth_moon_scene();
    test_saturn_scene();
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
    test_earth_moon_rotation();
    test_kerr_null_geodesics();
    test_sgr_a_scene();
    test_death_star_precession();
    test_satellite_fades();
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
