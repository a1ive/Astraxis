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

void test_voyager_scene()
{
    Scene scene = load_scene_or_die("voyager.toml");
    const int v1 = scene.find("Voyager 1");
    const int v2 = scene.find("Voyager 2");
    const int jupiter = scene.find("Jupiter");
    const int neptune = scene.find("Neptune");
    check(v1 > 0 && v2 > 0 && jupiter > 0 && neptune > 0, "voyager.toml bodies");

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

// Auto tour: many shots at varying times keep the camera finite and outside bodies.
void test_camera_director(uint32_t seed)
{
    Scene scene = load_jupiter();
    double t = 0.0;
    parse_utc("2026-10-02 00:00", &t);
    scene.update(t);

    OrbitCamera camera;
    camera.set_up_axis(scene.bodies[0].pole);
    camera.focus(0, scene);

    CameraDirector director;
    director.set_time_warp(3600.0);
    director.start(scene, camera, seed);

    bool kinds_seen[static_cast<int>(CameraDirector::ShotKind::Count)] = {};
    const double dt = 1.0 / 30.0;
    for (int frame = 0; frame < 30 * 60 * 20; ++frame) { // 20 minutes of real time
        t += dt * 3600.0;
        scene.update(t);
        director.update(dt, scene, camera);
        camera.update(dt, scene);
        kinds_seen[static_cast<int>(director.current_shot())] = true;

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
    test_voyager_scene();
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
    test_earth_moon_rotation();
    test_kerr_null_geodesics();
    test_sgr_a_scene();
    test_death_star_precession();
    for (uint32_t seed : {42u, 7u, 2026u, 99u, 12345u}) {
        test_camera_director(seed);
    }
    test_camera_director_black_hole();

    if (g_failures == 0) {
        std::printf("All tests passed.\n");
        return 0;
    }
    std::printf("%d failure(s).\n", g_failures);
    return 1;
}
