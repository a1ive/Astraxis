// Tests for the spacecraft trajectories in the scenes.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/ephemeris.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

using namespace astraxis;

namespace {

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
// Horizons, moon flybys against the distances NASA/JPL published.
void test_jupiter_missions()
{
    Scene scene = load_jupiter();
    const int galileo = scene.find("Galileo");
    const int juno = scene.find("Juno");
    const int europa = scene.find("Europa");
    const int ganymede = scene.find("Ganymede");
    check(galileo > 0 && juno > 0 && europa > 0 && ganymede > 0, "Jupiter mission bodies");
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
}

// Artemis I/II and CAPSTONE against NASA's published figures (distances only:
// mission blogs give times in local time zones, at the spacecraft or on Earth).
void test_earth_moon_missions()
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

    // ARTEMIS [ART]: P1 orbited L2 beyond the far side (from 2010-08-25), then moved to
    // L1 in January 2011 [ARTH]; P2 orbited L1 between the Earth and the Moon (sampled
    // daily, a month after arrival until shortly before leaving). P1's first lunar orbit:
    // "roughly 2,200 x 17,000 miles", which matches the distances from the Moon's center
    // (not altitudes).
    const int p1 = scene.find("ARTEMIS-P1");
    const int p2 = scene.find("ARTEMIS-P2");
    const int lcross = scene.find("LCROSS");
    const int ladee = scene.find("LADEE");
    check(p1 > 0 && p2 > 0 && lcross > 0 && ladee > 0, "earth_moon.toml: ARTEMIS, LCROSS, LADEE");
    // Signed distance beyond the Moon along the Earth-Moon line (> 0: far side).
    auto beyond_moon = [&](int craft, double t) {
        const glm::dvec3 m = at(moon, t);
        return glm::dot(at(craft, t) - m, glm::normalize(m));
    };
    auto least_beyond = [&](int craft, const CalendarDateTime& from, const CalendarDateTime& to, double sign) {
        double least = 1e300;
        for (double t = tdb(from); t < tdb(to); t += kSecondsPerDay) {
            least = std::min(least, sign * beyond_moon(craft, t));
        }
        return least;
    };
    check(least_beyond(p1, {2010, 9, 25, 0, 0, 0}, {2010, 12, 25, 0, 0, 0}, 1.0) > 0.0,
          "ARTEMIS-P1 stays beyond the Moon around L2");
    check(least_beyond(p1, {2011, 2, 1, 0, 0, 0}, {2011, 6, 15, 0, 0, 0}, -1.0) > 0.0,
          "ARTEMIS-P1 then stays on the Earth side around L1");
    check(least_beyond(p2, {2010, 11, 22, 0, 0, 0}, {2011, 6, 30, 0, 0, 0}, -1.0) > 0.0,
          "ARTEMIS-P2 stays on the Earth side around L1");
    auto from_moon = [&](int craft) { return [&, craft](double t) { return glm::length(at(craft, t) - at(moon, t)); }; };
    const double p1_peri = extremum(from_moon(p1), tdb({2011, 6, 27, 15, 16, 34}), 3600.0, false);
    const double p1_apo = extremum(from_moon(p1), tdb({2011, 6, 28, 15, 25, 3}), 3600.0, true);
    check(std::abs(p1_peri - 2200.0 * kMile) < 50.0 * kMile, "ARTEMIS-P1 first perilune (km)", p1_peri);
    check(std::abs(p1_apo - 17000.0 * kMile) < 500.0 * kMile, "ARTEMIS-P1 first apolune (km)", p1_apo);

    // LCROSS [LCR]: the table ends at the Shepherding Spacecraft's impact, 11:35:36.1 UTC
    // (+66.18 s to TDB), 3,809 m below the mean radius at 84.719 S. The scene's lunar
    // pole leaves out the IAU terms in E1 (up to ~4 deg), so the latitude is loose.
    const double t_lcross = tdb({2009, 10, 9, 11, 36, 42});
    const glm::dvec3 r_lcross = at(lcross, t_lcross - 0.5) - at(moon, t_lcross - 0.5);
    const glm::dvec3& moon_pole = scene.bodies[static_cast<size_t>(moon)].pole;
    const double lcross_lat = std::asin(glm::dot(glm::normalize(r_lcross), moon_pole)) * kRadToDeg;
    check(std::abs(glm::length(r_lcross) - (moon_r - 3.809)) < 3.0, "LCROSS impact radius (km)", glm::length(r_lcross));
    check(std::abs(lcross_lat + 84.719) < 4.5, "LCROSS impact latitude (deg)", lcross_lat);
    scene.update(t_lcross - 60.0);
    const bool lcross_before = scene.bodies[static_cast<size_t>(lcross)].visible;
    scene.update(t_lcross + 60.0);
    check(lcross_before && !scene.bodies[static_cast<size_t>(lcross)].visible, "LCROSS gone after its impact");

    // LADEE [LRO]: the data end at the impact, "approximately 3,800 miles per hour
    // (1,699 meters per second)", a few km above the mean radius.
    const double t_ladee = tdb({2014, 4, 18, 4, 33, 0});
    const State ladee_end = scene.icrf_state_at(ladee, t_ladee);
    const State moon_end = scene.icrf_state_at(moon, t_ladee);
    const double ladee_speed = glm::length(ladee_end.velocity - moon_end.velocity);
    const double ladee_alt = glm::length(ladee_end.position - moon_end.position) - moon_r;
    check(std::abs(ladee_speed - 1.699) < 0.02, "LADEE impact speed (km/s)", ladee_speed);
    check(ladee_alt > 0.0 && ladee_alt < 5.0, "LADEE altitude at the end of the data (km)", ladee_alt);

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
}

// Cassini and Huygens in the Saturn scene: flyby altitudes against NASA, the
// Huygens landing and the Grand Finale dives.
void test_saturn_missions()
{
    Scene scene = load_scene_or_die("saturn.toml");
    auto tdb = [](const CalendarDateTime& c) { return (jd_from_calendar(c) - kJ2000Jd) * kSecondsPerDay; };
    auto at = [&](int body, double t) { return scene.icrf_state_at(body, t).position; };

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
    // In the Moon's shadow 18:19-18:47 UTC, 28 minutes (Farquhar 2001, J. Astronaut. Sci. 49, 23,
    // Fig. 23; TDA Progress Report 42-84, p. 178, says "more than 30 minutes").
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
    check(std::abs(shadow_s - 28.0 * 60.0) < 2.0 * 60.0, "minutes in the Moon's shadow at the last flyby",
          shadow_s / 60.0);
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

// Rosetta and 67P: flybys against the Horizons notes for -226 and ESA's figures; the
// comet-relative phase against ESA's mission milestones.
void test_rosetta_scene()
{
    Scene scene = load_scene_or_die("rosetta.toml");
    const int earth = scene.find("Earth");
    const int lutetia = scene.find("Lutetia");
    const int comet = scene.find("67P");
    const int craft = scene.find("Rosetta");
    check(earth > 0 && lutetia > 0 && comet > 0 && craft > 0, "rosetta.toml bodies");
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

    // Earth flybys (Horizons notes): 2005-03-04 22:10 UTC 8,340 km from the center;
    // 2007-11-13 20:57 UTC 11,678 km; 2009-11-13 07:45:40 UTC ~2,490 km above the surface.
    // The Earth is the scene's (around the Earth-Moon barycenter, 50 km table tolerance).
    const double r_earth = scene.bodies[static_cast<size_t>(earth)].equatorial_radius_km;
    const auto [e1, te1] = closest(craft, earth, utc({2005, 3, 4, 22, 10, 0}), 3600.0, 5.0);
    const auto [e2, te2] = closest(craft, earth, utc({2007, 11, 13, 20, 57, 0}), 3600.0, 5.0);
    const auto [e3, te3] = closest(craft, earth, utc({2009, 11, 13, 7, 45, 40}), 3600.0, 5.0);
    check(std::abs(e1 - 8340.0) < 60.0 && std::abs(te1 - utc({2005, 3, 4, 22, 10, 0})) < 120.0,
          "Rosetta Earth flyby 2005 (km from the center)", e1);
    check(std::abs(e2 - 11678.0) < 60.0 && std::abs(te2 - utc({2007, 11, 13, 20, 57, 0})) < 120.0,
          "Rosetta Earth flyby 2007 (km from the center)", e2);
    check(std::abs(e3 - r_earth - 2490.0) < 60.0 && std::abs(te3 - utc({2009, 11, 13, 7, 45, 40})) < 120.0,
          "Rosetta Earth flyby 2009 (km above the surface)", e3 - r_earth);
    std::printf("info: Rosetta Earth flybys %.0f, %.0f, %.0f km from the center\n", e1, e2, e3);

    // 21 Lutetia, 2010-07-10 15:44:45 UTC, 3,235 km (Horizons notes).
    const auto [dl, tl] = closest(craft, lutetia, utc({2010, 7, 10, 15, 44, 45}), 1800.0, 0.5);
    check(std::abs(dl - 3235.0) < 100.0 && std::abs(tl - utc({2010, 7, 10, 15, 44, 45})) < 120.0,
          "Rosetta at Lutetia (km)", dl);

    // At the comet (ESA): arrival 2014-08-06 at ~100 km; Philae released 2014-11-12 08:35 UTC
    // from 22.5 km of the center; the 2015-02-14 flyby 6 km above the surface (~8 km from the
    // center, radii 2.4 x 1.55 x 1.2 km); touchdown 2016-09-30 10:39:28 UTC.
    const double arrival = glm::length(at(craft, utc({2014, 8, 6, 9, 6, 0})) - at(comet, utc({2014, 8, 6, 9, 6, 0})));
    check(std::abs(arrival - 100.0) < 10.0, "Rosetta at arrival (km from 67P)", arrival);
    const double release =
        glm::length(at(craft, utc({2014, 11, 12, 8, 35, 0})) - at(comet, utc({2014, 11, 12, 8, 35, 0})));
    check(std::abs(release - 22.5) < 0.5, "Rosetta at Philae's release (km from 67P)", release);
    const auto [flyby, t_flyby] = closest(craft, comet, utc({2015, 2, 14, 12, 41, 0}), 7200.0, 10.0);
    check(flyby > 6.5 && flyby < 9.0 && std::abs(t_flyby - utc({2015, 2, 14, 12, 41, 0})) < 600.0,
          "Rosetta close flyby 2015-02-14 (km from the center)", flyby);
    const double t_end = utc({2016, 9, 30, 10, 38, 30});
    check(glm::length(at(craft, t_end) - at(comet, t_end)) < 3.0, "Rosetta at the surface (km from the center)",
          glm::length(at(craft, t_end) - at(comet, t_end)));
    scene.update(utc({2016, 9, 30, 11, 0, 0}));
    check(!scene.bodies[static_cast<size_t>(craft)].visible, "Rosetta gone after its touchdown");

    // Hand-over from the heliocentric table to the comet-relative one (2014-05-01): no jump.
    const double t_switch = utc({2014, 5, 1, 0, 0, 0});
    double worst = 0.0;
    for (double t = t_switch - 2 * kSecondsPerDay; t < t_switch + 2 * kSecondsPerDay; t += 600.0) {
        const State a = scene.icrf_state_at(craft, t);
        worst = std::max(worst, glm::length(at(craft, t + 600.0) - a.position - a.velocity * 600.0));
    }
    check(worst < 5.0, "Rosetta hand-over to the comet-relative table (km per 10 min)", worst);
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

} // namespace

void run_mission_tests()
{
    test_solar_system_missions();
    test_parker_scene();
    test_jupiter_missions();
    test_earth_moon_missions();
    test_saturn_missions();
    test_mercury_scene();
    test_isee3_scene();
    test_rosetta_scene();
    test_halley_armada();
    test_jwst_scene();
}
