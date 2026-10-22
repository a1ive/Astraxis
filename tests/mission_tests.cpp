// Tests for the spacecraft trajectories in the scenes.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/ephemeris.hpp"
#include "ephem/kepler.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace astraxis;

namespace {

// Distance from the origin to the outermost surface of a shape along the direction of p (km).
double shape_radius(const ShapeModel& shape, const glm::dvec3& p)
{
    const glm::dvec3 d = glm::normalize(p);
    double best = 0.0;
    for (size_t i = 0; i + 2 < shape.indices.size(); i += 3) {
        const glm::dvec3 a(shape.positions[shape.indices[i]]);
        const glm::dvec3 e1 = glm::dvec3(shape.positions[shape.indices[i + 1]]) - a;
        const glm::dvec3 e2 = glm::dvec3(shape.positions[shape.indices[i + 2]]) - a;
        // Moller-Trumbore, ray from the origin.
        const glm::dvec3 h = glm::cross(d, e2);
        const double det = glm::dot(e1, h);
        if (std::abs(det) < 1e-15) {
            continue;
        }
        const glm::dvec3 s = -a;
        const double u = glm::dot(s, h) / det;
        const glm::dvec3 q = glm::cross(s, e1);
        const double w = glm::dot(d, q) / det;
        if (u < 0.0 || w < 0.0 || u + w > 1.0) {
            continue;
        }
        best = std::max(best, glm::dot(e2, q) / det);
    }
    return best;
}

// An orbit insertion event's ghost (SceneEvent::Ghost): the spacecraft's path
// without the burn. The burn's delta-v as the path shows it, against the
// published value (none if 0) within `tolerance` (relative); the ghost and the
// readout only during the event; the burn not under way at the time given for
// its start; the ghost on the approach path at the event's start (it drifts off
// slowly: the gravity it coasts in is not the whole field, and between knots a
// table is only as close as its tolerance, a few km), apart an hour after
// the burn, and gone by the end of its window: unbound by then. An impulsive
// burn (a reconstruction's) has no halfway. Returns the event's index, -1 if none.
int check_insertion_ghost(Scene& scene, const char* craft, double published_m_s, double tolerance,
                          bool impulsive = false)
{
    const std::string name(craft);
    const int body = scene.find(craft);
    int event = -1;
    for (size_t i = 0; i < scene.events.size(); ++i) {
        if (body >= 0 && scene.events[i].ghost.body == body) {
            event = static_cast<int>(i);
        }
    }
    check(event >= 0, (name + "'s orbit insertion has a ghost").c_str());
    if (event < 0) {
        return -1;
    }
    const SceneEvent::Ghost& g = scene.events[static_cast<size_t>(event)].ghost;
    const double dv = g.delta_v_km_s.empty() ? 0.0 : g.delta_v_km_s.back() * 1000.0;
    if (published_m_s > 0.0) {
        check(std::abs(dv / published_m_s - 1.0) < tolerance, (name + " orbit insertion delta-v (m/s)").c_str(), dv);
    }

    const int ghost = g.ghost;
    const Body& gb = scene.bodies[static_cast<size_t>(ghost)];
    const double burn_end = g.t_tdb + g.delta_v_span_s;
    const double mid = 0.5 * (g.burn_tdb + burn_end);
    scene.set_active_event(-1);
    scene.update(mid);
    double readout = 0.0;
    check(!gb.visible && !scene.burn_delta_v(body, &readout), (name + "'s ghost and readout only during its event").c_str());
    scene.set_active_event(event);
    scene.update(mid);
    const double low = impulsive ? 0.999 : 0.3;
    const double high = impulsive ? 1.001 : 0.7;
    check(gb.visible && scene.burn_delta_v(body, &readout) && readout > low * dv * 1e-3 && readout < high * dv * 1e-3,
          (name + "'s delta-v readout halfway through the burn (km/s)").c_str(), readout);
    // The ghost leaves the real path a minute before the time given. By then, less than half
    // a minute of the burn's mean thrust: the readout grows linearly between knots a minute
    // apart, so a burn starting in the interval after that time shows a little before it.
    // (An impulse shows over one step of the readout, about 4 s; until then a coast that
    // leaves out, e.g., sunlight pressure leaves a trace.)
    scene.update(impulsive ? g.burn_tdb - 4.0 : g.burn_tdb);
    scene.burn_delta_v(body, &readout);
    const double half_minute = impulsive ? 0.01 : dv / (g.delta_v_span_s - 120.0) * 30.0;
    check(readout * 1000.0 < half_minute, (name + "'s burn not under way at its given start (m/s)").c_str(),
          readout * 1000.0);
    scene.set_active_event(-1);

    const Body& parent = scene.bodies[static_cast<size_t>(gb.parent)];
    const double gm = parent.gm_km3_s2 > 0.0 ? parent.gm_km3_s2 : parent.body_gm_km3_s2;
    auto apart = [&](double t) {
        return glm::length(scene.icrf_state_at(ghost, t).position - scene.icrf_state_at(body, t).position);
    };
    const double before = apart(g.start_tdb);
    check(before < 10.0, (name + "'s ghost on the approach path at the event's start (km)").c_str(), before);
    const double after = apart(burn_end + 3600.0);
    check(after > 1000.0, (name + "'s ghost apart an hour after the burn (km)").c_str(), after);
    OrbitElements ghost_el;
    const State last = gb.motion->eval(g.end_tdb);
    check(!osculating_elements(last, gm, &ghost_el) || ghost_el.period_s == 0.0,
          (name + " would have flown past without the burn").c_str(), glm::length(last.position));
    std::printf("info: %s orbit insertion: %.1f m/s from the path (published %.0f); the ghost %.2f km off at the "
                "event's start (%.0f min before the burn), %.0f km an hour after, %.0f km from %s at the end\n",
                craft, dv, published_m_s, before, (g.burn_tdb - g.start_tdb) / 60.0, after,
                glm::length(last.position), parent.name.c_str());
    return event;
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

// Perihelion of a body in a scene with the Sun, and the osculating heliocentric orbit
// there (inclination to the J2000 ecliptic, IAU 1976 obliquity).
struct Perihelion {
    double t = 0.0;
    double q_km = 0.0;
    double e = 0.0;
    double incl_deg = 0.0;
    double v_km_s = 0.0;
    double v_inf_km_s = 0.0;
};

Perihelion find_perihelion(const Scene& scene, int body, double t0, double t1)
{
    const int sun = scene.find("Sun");
    auto rel = [&](double t) {
        const State b = scene.icrf_state_at(body, t);
        const State s = scene.icrf_state_at(sun, t);
        return State{b.position - s.position, b.velocity - s.velocity};
    };
    Perihelion p;
    p.q_km = 1e300;
    for (double t = t0; t <= t1; t += 60.0) {
        const double r = glm::length(rel(t).position);
        if (r < p.q_km) {
            p.q_km = r;
            p.t = t;
        }
    }
    const State s = rel(p.t);
    const double gm_sun = 1.32712440041279419e11; // [GM]
    OrbitElements el;
    osculating_elements(s, gm_sun, &el);
    p.e = el.eccentricity;
    const double obliquity = 84381.448 / 3600.0 * kDegToRad;
    const glm::dvec3 ecliptic_pole(0.0, -std::sin(obliquity), std::cos(obliquity));
    p.incl_deg = std::acos(glm::dot(glm::normalize(glm::cross(s.position, s.velocity)), ecliptic_pole)) / kDegToRad;
    p.v_km_s = glm::length(s.velocity);
    p.v_inf_km_s = std::sqrt(p.v_km_s * p.v_km_s - 2.0 * gm_sun / glm::length(s.position));
    return p;
}

// The interstellar objects in the Solar System scene: the heliocentric tables against
// Horizons' osculating elements at perihelion and against the papers.
void test_solar_system_interstellar()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    const int oumuamua = scene.find("1I/'Oumuamua");
    const int atlas = scene.find("3I/ATLAS");
    check(oumuamua > 0 && atlas > 0, "solar_system.toml has 1I/'Oumuamua and 3I/ATLAS");
    auto at = [&](const char* name, double t) { return scene.icrf_state_at(scene.find(name), t).position; };
    auto utc = [](const CalendarDateTime& c) { return tdb_from_utc_jd(jd_from_calendar(c)); };
    auto scan_min = [](double t0, double t1, double step, auto&& f) {
        std::pair<double, double> best{1e300, 0.0}; // value, time
        for (double t = t0; t <= t1; t += step) {
            best = std::min(best, std::pair{f(t), t});
        }
        return best;
    };
    auto shown_only_in = [&](int body, const char* what, std::initializer_list<std::pair<const char*, bool>> when) {
        for (const auto& [date, shown] : when) {
            double t = 0.0;
            parse_utc(date, &t);
            scene.update(t);
            check(scene.bodies[static_cast<size_t>(body)].visible == shown, what);
        }
    };

    // 1I/'Oumuamua. Horizons JPL#16 osculating elements at 2017-09-09 00:00 TDB: QR =
    // 38,304,506 km, Tp = JD 2458005.99498 TDB (11:53), e 1.199431, i 122.737 deg; 26.3 km/s
    // from afar and 87.3 km/s at perihelion, as the scene's notes and captions say. [JEW]'s
    // early orbit: q 0.254 au, e 1.197, i 122.6 deg.
    const Perihelion p1 = find_perihelion(scene, oumuamua, utc({2017, 9, 7, 0, 0, 0}), utc({2017, 9, 12, 0, 0, 0}));
    check(std::abs(p1.q_km - 38304506.0) < 100.0, "1I/'Oumuamua perihelion distance (km)", p1.q_km);
    check(std::abs(p1.t - tdb_from_jd_tdb(2458005.99498)) < 0.01 * kSecondsPerDay, "1I/'Oumuamua perihelion time (d)",
          (p1.t - tdb_from_jd_tdb(2458005.99498)) / kSecondsPerDay);
    check(std::abs(p1.e - 1.199431) < 0.0005 && std::abs(p1.incl_deg - 122.737) < 0.01,
          "1I/'Oumuamua eccentricity and inclination", p1.incl_deg);
    check(std::abs(p1.v_km_s - 87.3) < 0.1 && std::abs(p1.v_inf_km_s - 26.3) < 0.1, "1I/'Oumuamua speeds (km/s)",
          p1.v_km_s);
    check(std::abs(p1.q_km / kAuKm - 0.254) < 0.003 && std::abs(p1.e - 1.197) < 0.003 &&
              std::abs(p1.incl_deg - 122.6) < 0.2,
          "1I/'Oumuamua against the early orbit of Jewitt et al. 2017");
    // Closest to the Earth: Horizons geocentric vectors without light-time correction,
    // hourly, 2017-10-14 ~18:00 UT, 0.161754 au.
    const auto [d_earth, t_earth] = scan_min(utc({2017, 10, 13, 0, 0, 0}), utc({2017, 10, 16, 0, 0, 0}), 60.0,
                                             [&](double t) { return glm::length(at("1I/'Oumuamua", t) - at("Earth", t)); });
    check(std::abs(d_earth / kAuKm - 0.161754) < 0.00001, "1I/'Oumuamua closest to the Earth (au)", d_earth / kAuKm);
    check(std::abs(t_earth - utc({2017, 10, 14, 18, 0, 0})) < 3.0 * 3600.0,
          "1I/'Oumuamua closest to the Earth (h from 18:00 UT)", (t_earth - utc({2017, 10, 14, 18, 0, 0})) / 3600.0);
    // The tumbling disc ([MSH]): its normal keeps 46.5 deg from the angular momentum and
    // circles it every 10.81 h, so it is back after one period and 93 deg away after half.
    const Body& disc = scene.bodies[static_cast<size_t>(oumuamua)];
    const double t_model = utc({2017, 10, 27, 7, 31, 0});
    const double p_phi = 10.81 * 3600.0;
    const glm::dvec3 spin_axis = unit_toward(91.872, -17.016);
    auto normal = [&](double t) { return scene.body_to_icrf(disc, t)[2]; };
    auto angle = [](const glm::dvec3& a, const glm::dvec3& b) {
        return std::acos(std::clamp(glm::dot(a, b), -1.0, 1.0)) / kDegToRad;
    };
    const double nutation = angle(normal(t_model + 1000.0), spin_axis);
    check(std::abs(nutation - 46.5) < 0.01, "1I/'Oumuamua nutation (deg)", nutation);
    check(angle(normal(t_model), normal(t_model + p_phi)) < 0.01 &&
              std::abs(angle(normal(t_model), normal(t_model + 0.5 * p_phi)) - 93.0) < 0.01,
          "1I/'Oumuamua precesses in 10.81 h");
    std::printf("info: 1I/'Oumuamua: perihelion %.6f au, e %.6f, i %.3f deg, %.2f / %.2f km/s; Earth %.5f au\n",
                p1.q_km / kAuKm, p1.e, p1.incl_deg, p1.v_km_s, p1.v_inf_km_s, d_earth / kAuKm);
    shown_only_in(oumuamua, "1I/'Oumuamua shown only where the table is (2015-09-01 .. 2019-09-30)",
                  {{"2015-08-01", false}, {"2015-10-01", true}, {"2019-09-01", true}, {"2019-11-01", false}});

    // 2I/Borisov. Horizons JPL#54 osculating elements at 2019-12-08 00:00 TDB: QR =
    // 300,171,607 km, Tp = JD 2458826.05364 TDB (13:17), e 3.356479, i 44.0526 deg; 32.3 km/s
    // from afar ([GUZ]: ~32 km/s from the first orbit) and 43.9 km/s at perihelion.
    const int borisov = scene.find("2I/Borisov");
    check(borisov > 0, "solar_system.toml has 2I/Borisov");
    const Perihelion p2 = find_perihelion(scene, borisov, utc({2019, 12, 5, 0, 0, 0}), utc({2019, 12, 12, 0, 0, 0}));
    check(std::abs(p2.q_km - 300171607.0) < 100.0, "2I/Borisov perihelion distance (km)", p2.q_km);
    check(std::abs(p2.t - tdb_from_jd_tdb(2458826.05364)) < 0.01 * kSecondsPerDay, "2I/Borisov perihelion time (d)",
          (p2.t - tdb_from_jd_tdb(2458826.05364)) / kSecondsPerDay);
    check(std::abs(p2.e - 3.356479) < 0.0005 && std::abs(p2.incl_deg - 44.0526) < 0.01,
          "2I/Borisov eccentricity and inclination", p2.incl_deg);
    check(std::abs(p2.v_km_s - 43.9) < 0.1 && std::abs(p2.v_inf_km_s - 32.3) < 0.1, "2I/Borisov speeds (km/s)",
          p2.v_km_s);
    // Closest to the Earth: Horizons geocentric vectors without light-time correction,
    // hourly, 2019-12-28 ~04:00 UT, 1.936902 au.
    const auto [d_earth2, t_earth2] = scan_min(utc({2019, 12, 25, 0, 0, 0}), utc({2019, 12, 31, 0, 0, 0}), 600.0,
                                               [&](double t) { return glm::length(at("2I/Borisov", t) - at("Earth", t)); });
    check(std::abs(d_earth2 / kAuKm - 1.936902) < 0.00001, "2I/Borisov closest to the Earth (au)", d_earth2 / kAuKm);
    check(std::abs(t_earth2 - utc({2019, 12, 28, 4, 0, 0})) < 6.0 * 3600.0,
          "2I/Borisov closest to the Earth (h from 04:00 UT)", (t_earth2 - utc({2019, 12, 28, 4, 0, 0})) / 3600.0);
    std::printf("info: 2I/Borisov: perihelion %.6f au, e %.6f, i %.3f deg, %.2f / %.2f km/s; Earth %.5f au\n",
                p2.q_km / kAuKm, p2.e, p2.incl_deg, p2.v_km_s, p2.v_inf_km_s, d_earth2 / kAuKm);
    shown_only_in(borisov, "2I/Borisov shown only where the table is (2017-12-01 .. 2021-12-31)",
                  {{"2017-11-01", false}, {"2018-01-01", true}, {"2021-12-01", true}, {"2022-02-01", false}});

    // 3I/ATLAS. Horizons JPL#54 osculating elements at 2025-10-29 00:00 TDB: QR =
    // 202,921,707 km (1.356448 au), Tp = JD 2460977.98225 TDB (11:34). (The elements at the
    // solution epoch, 2026-02-19, give 1.356481 au and Tp 0.013 d later.) e 6.1393,
    // inclination 175.113 deg to the J2000 ecliptic; [REN]: 6.14 and 175 deg. 58.0 km/s
    // from afar and 68.3 km/s at perihelion, as the scene's notes and captions say.
    const Perihelion p3 = find_perihelion(scene, atlas, utc({2025, 10, 25, 0, 0, 0}), utc({2025, 11, 3, 0, 0, 0}));
    check(std::abs(p3.q_km - 202921707.0) < 100.0, "3I/ATLAS perihelion distance (km)", p3.q_km);
    check(std::abs(p3.t - tdb_from_jd_tdb(2460977.98225)) < 0.01 * kSecondsPerDay, "3I/ATLAS perihelion time (d)",
          (p3.t - tdb_from_jd_tdb(2460977.98225)) / kSecondsPerDay);
    check(std::abs(p3.e - 6.1393) < 0.002, "3I/ATLAS eccentricity", p3.e);
    check(std::abs(p3.incl_deg - 175.113) < 0.01, "3I/ATLAS inclination to the ecliptic (deg)", p3.incl_deg);
    check(std::abs(p3.v_km_s - 68.3) < 0.1 && std::abs(p3.v_inf_km_s - 58.0) < 0.1,
          "3I/ATLAS speed at perihelion (km/s)", p3.v_km_s);

    // [REN]: closest to Mars on 2025-10-03 04:38 UT at 0.194 au; Mars here is [APX], good
    // to far better than the 0.001 au tolerance.
    const auto [d_mars, t_mars] = scan_min(utc({2025, 10, 2, 0, 0, 0}), utc({2025, 10, 4, 12, 0, 0}), 60.0,
                                           [&](double t) { return glm::length(at("3I/ATLAS", t) - at("Mars", t)); });
    check(std::abs(d_mars / kAuKm - 0.194) < 0.001, "3I/ATLAS closest to Mars in the Solar System scene (au)",
          d_mars / kAuKm);
    check(std::abs(t_mars - utc({2025, 10, 3, 4, 38, 0})) < 2.0 * 3600.0, "3I/ATLAS closest to Mars (h from 04:38 UT)",
          (t_mars - utc({2025, 10, 3, 4, 38, 0})) / 3600.0);

    // [REN]: in solar conjunction, seen from the Earth, around 10-21 (the perihelion event).
    const auto [elong, t_conj] = scan_min(utc({2025, 10, 10, 0, 0, 0}), utc({2025, 11, 1, 0, 0, 0}), 600.0, [&](double t) {
        const glm::dvec3 e = at("Earth", t);
        return std::acos(glm::dot(glm::normalize(at("Sun", t) - e), glm::normalize(at("3I/ATLAS", t) - e)));
    });
    check(std::abs(t_conj - utc({2025, 10, 21, 12, 0, 0})) < 2.0 * kSecondsPerDay,
          "3I/ATLAS solar conjunction from the Earth (d from 10-21)",
          (t_conj - utc({2025, 10, 21, 12, 0, 0})) / kSecondsPerDay);
    std::printf("info: 3I/ATLAS: perihelion %.6f au, e %.4f, i %.3f deg, %.2f / %.2f km/s; Mars %.4f au; "
                "conjunction %.1f deg on day %.2f of October\n",
                p3.q_km / kAuKm, p3.e, p3.incl_deg, p3.v_km_s, p3.v_inf_km_s, d_mars / kAuKm, elong / kDegToRad,
                (t_conj - utc({2025, 10, 1, 0, 0, 0})) / kSecondsPerDay + 1.0);
    shown_only_in(atlas, "3I/ATLAS shown only where the table is (2024-01-01 .. 2027-12-31)",
                  {{"2023-12-01", false}, {"2024-02-01", true}, {"2027-12-01", true}, {"2028-02-01", false}});
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

    // Cassini's flyby. Horizons, -82 @ 500@599, 2000-12-31 00:00 TDB, ICRF km: a day
    // after closest approach, between two knots of the table.
    const int cassini = scene.find("Cassini");
    check(cassini > 0, "Cassini in the Jupiter scene");
    const double t_c = tdb_from_jd_tdb(2451909.5);
    const glm::dvec3 cassini_hzn(6.643951662465836E+06, -6.555686166043234E+06, -3.019377321932644E+06);
    check(glm::length(at(cassini, t_c) - cassini_hzn) < 5.0, "Cassini matches Horizons after its Jupiter flyby (km)",
          glm::length(at(cassini, t_c) - cassini_hzn));
    // NASA Science, "Cassini Celebrates 10 Years Since Jupiter Encounter" (2010-12-28):
    // within about 9.7 million km of the cloud tops, 2000-12-30 10:05 UTC.
    const double t_ca = (jd_from_calendar({2000, 12, 30, 10, 5, 0}) - kJ2000Jd) * kSecondsPerDay;
    double ca = 1e300;
    double t_min = 0.0;
    for (double dt = -6.0 * 3600.0; dt <= 6.0 * 3600.0; dt += 60.0) {
        const double d = glm::length(at(cassini, t_ca + dt));
        if (d < ca) {
            ca = d;
            t_min = dt;
        }
    }
    const double cloud_tops = ca - scene.bodies[static_cast<size_t>(scene.find("Jupiter"))].equatorial_radius_km;
    check(std::abs(cloud_tops - 9.7e6) < 0.05e6, "Cassini's closest approach to Jupiter's cloud tops (km)", cloud_tops);
    check(std::abs(t_min) < 600.0, "Cassini's closest approach time (s from NASA's)", t_min);
    // Jupiter turns the hyperbola by about 12 deg (jupiter.toml), near its equator.
    const glm::dvec3 v0 = scene.icrf_state_at(cassini, tdb_from_jd_tdb(2451818.5)).velocity;
    const glm::dvec3 v1 = scene.icrf_state_at(cassini, tdb_from_jd_tdb(2452000.5)).velocity;
    const double turn = glm::degrees(std::acos(glm::dot(glm::normalize(v0), glm::normalize(v1))));
    check(std::abs(turn - 12.0) < 1.0, "Jupiter turns Cassini's path (deg)", turn);
    const glm::dvec3 pole = scene.bodies[static_cast<size_t>(scene.find("Jupiter"))].pole;
    const State s_c = scene.icrf_state_at(cassini, t_c);
    const double tilt = glm::degrees(std::acos(glm::dot(glm::normalize(glm::cross(s_c.position, s_c.velocity)), pole)));
    check(tilt < 4.0, "Cassini's hyperbola near Jupiter's equator (deg)", tilt);
    parse_utc("2000-09-30", &t);
    scene.update(t);
    check(!scene.bodies[static_cast<size_t>(cassini)].visible, "Cassini not yet at Jupiter");
    parse_utc("2000-12-30", &t);
    scene.update(t);
    check(scene.bodies[static_cast<size_t>(cassini)].visible && scene.bodies[static_cast<size_t>(galileo)].visible,
          "Galileo and Cassini both at Jupiter");

    // Juno's orbit insertion: 542 m/s [JOI], read with Jupiter's J2 and J4 (without them,
    // the oblateness would count as hundreds of m/s of thrust).
    check_insertion_ghost(scene, "Juno", 542.0, 0.01);
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
    // The lunar orbit insertion: 737 mph planned [LADEEI].
    check_insertion_ghost(scene, "LADEE", 737.0 * 0.44704, 0.02);

    // Chang'e 5-T1's service module [L16] (reconstructed): the farthest distance of the
    // flown trajectory, 540,000 km (Table 4; the design's, flown "essentially identical");
    // around L2 beyond the Moon from the injection (2014-11-28) to the departure (01-04);
    // brake 1 at the "perilune arrival time" 19:04:33 UTC (+67.18 s to TDB); the orbit
    // after brake 3 about 200 km circular (Eq. 2).
    const int ce5 = scene.find("Chang'e 5-T1");
    check(ce5 > 0, "earth_moon.toml: Chang'e 5-T1");
    const double ce5_far = extremum([&](double t) { return glm::length(at(ce5, t)); },
                                    tdb({2014, 11, 9, 6, 43, 0}), 3600.0, true);
    check(std::abs(ce5_far - 540000.0) < 5000.0, "Chang'e 5-T1 farthest from the Earth (km)", ce5_far);
    check(least_beyond(ce5, {2014, 11, 29, 0, 0, 0}, {2015, 1, 4, 0, 0, 0}, 1.0) > 30000.0,
          "Chang'e 5-T1 stays beyond the Moon around L2");
    const double t_brake1 = tdb({2015, 1, 10, 19, 5, 40}) + 0.18;
    const double ce5_peri = extremum(above_moon(ce5), t_brake1, 600.0, false);
    check(std::abs(above_moon(ce5)(t_brake1) - ce5_peri) < 1.0, "Chang'e 5-T1 brake 1 at perilune (km above it)",
          above_moon(ce5)(t_brake1) - ce5_peri);
    double ce5_low = 1e300, ce5_high = 0.0;
    for (double t = tdb({2015, 1, 13, 0, 0, 0}); t < tdb({2015, 1, 13, 12, 0, 0}); t += 60.0) {
        ce5_low = std::min(ce5_low, above_moon(ce5)(t));
        ce5_high = std::max(ce5_high, above_moon(ce5)(t));
    }
    check(ce5_low > 150.0 && ce5_high < 250.0, "Chang'e 5-T1 final orbit about 200 km (km, highest)", ce5_high);
    // Brake 1's ghost: 219.7 m/s (Table 5), impulsive in the reconstruction. Around the Moon
    // alone the orbit before it is closed; with the Earth's tide the ghost leaves (unbound
    // 12 days later).
    const int ce5_brake = check_insertion_ghost(scene, "Chang'e 5-T1", 219.7, 0.001, true);
    if (ce5_brake >= 0) {
        const SceneEvent::Ghost& g = scene.events[static_cast<size_t>(ce5_brake)].ghost;
        OrbitElements before_brake;
        check(osculating_elements(scene.bodies[static_cast<size_t>(g.ghost)].motion->eval(g.t_tdb),
                                  scene.bodies[static_cast<size_t>(moon)].gm_km3_s2, &before_brake) &&
                  before_brake.period_s > 0.0,
              "Chang'e 5-T1's orbit around the Moon alone closed before brake 1");
    }

    // Chang'e 5-T1's rocket stage (WE0913A): close Earth flybys on 2022-01-20 and 02-08
    // [WE]; the table ends where it reaches the elevation of the impact site, 1,863 m above
    // the mean radius at 5.226 N, 234.486 E [LROC], ~15 deg off vertical [WE]. The site
    // checks the JPL fit against the crater found afterwards: 5.26 N, 234.89 E here,
    // although the scene's lunar rotation leaves out the IAU periodic terms (up to ~4 deg).
    const int stage = scene.find("Chang'e 5-T1 booster");
    check(stage > 0, "earth_moon.toml: Chang'e 5-T1 booster");
    auto from_earth = [&](double t) { return glm::length(at(stage, t)); };
    const double stage_peri1 = extremum(from_earth, tdb({2022, 1, 20, 23, 48, 0}), 1800.0, false);
    const double stage_peri2 = extremum(from_earth, tdb({2022, 2, 8, 21, 2, 0}), 1800.0, false);
    check(stage_peri1 < 40000.0 && stage_peri2 < 50000.0, "Chang'e 5-T1 booster Earth flybys (km, first)", stage_peri1);
    const double t_stage = tdb({2022, 3, 4, 12, 26, 36});
    const State stage_end = scene.icrf_state_at(stage, t_stage - 0.5);
    const State moon_stage = scene.icrf_state_at(moon, t_stage - 0.5);
    const glm::dvec3 r_stage = stage_end.position - moon_stage.position;
    const glm::dvec3 v_stage = stage_end.velocity - moon_stage.velocity;
    // Half a second before the end it is 1.3 km higher (2.5 km/s down); the Moon's own
    // table adds up to 2 km.
    const double stage_elev = glm::length(r_stage) - moon_r;
    check(stage_elev > 1.863 && stage_elev < 1.863 + 1.3 + 2.0, "Chang'e 5-T1 booster impact elevation (km)",
          stage_elev);
    const double stage_angle =
        std::asin(-glm::dot(glm::normalize(r_stage), glm::normalize(v_stage))) * kRadToDeg;
    check(std::abs(stage_angle - 75.0) < 3.0, "Chang'e 5-T1 booster impact angle (deg from horizontal)", stage_angle);
    const glm::dvec3 site = glm::transpose(body_rotation(scene.bodies[static_cast<size_t>(moon)], {}, t_stage)) *
                            glm::normalize(r_stage);
    const double site_lat = std::asin(site.z) * kRadToDeg;
    const double site_lon = std::atan2(site.y, site.x) * kRadToDeg + 360.0;
    check(std::abs(site_lat - 5.226) < 1.0 && std::abs(site_lon - 234.486) < 1.0,
          "Chang'e 5-T1 booster impact site (deg, longitude)", site_lon);
    scene.update(t_stage + 60.0);
    check(!scene.bodies[static_cast<size_t>(stage)].visible, "Chang'e 5-T1 booster gone after its impact");

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
    // Hyperion is irregular: the distance from its center, against Horizons (-82 @ 607,
    // 604.4 km at 2005-09-26 02:25:48 TDB; Thomas et al. 2007 give 618 km), within the
    // 10 km tolerance of Hyperion's own table.
    const double hyperion = altitude(cassini, "Hyperion", {2005, 9, 26, 2, 25, 48}) +
                            scene.bodies[static_cast<size_t>(scene.find("Hyperion"))].equatorial_radius_km;
    check(std::abs(hyperion - 604.4) < 10.0, "Cassini-Hyperion closest approach (km)", hyperion);
    // Phoebe likewise, against Horizons (-82 @ 609, 2,194.3 km at 2004-06-11 19:34:45 TDB;
    // NASA PIA18411 gives 2,068 km above the surface), within the tables' tolerances.
    const double phoebe = altitude(cassini, "Phoebe", {2004, 6, 11, 19, 34, 45}) +
                          scene.bodies[static_cast<size_t>(scene.find("Phoebe"))].equatorial_radius_km;
    check(std::abs(phoebe - 2194.3) < 15.0, "Cassini-Phoebe closest approach (km)", phoebe);

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

    // The orbit insertion: 626.17 m/s [SOI], read with Saturn's J2 and J4.
    check_insertion_ghost(scene, "Cassini", 626.17, 0.01);
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

    // The orbit insertion burn's ghost: 862 m/s [MESM], 1%.
    check_insertion_ghost(scene, "MESSENGER", 862.0, 0.01);

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

// Mars Express, MRO and MAVEN around Mars, against figures published by NASA and
// ESA independently of the Horizons tables.
void test_mars_scene()
{
    Scene scene = load_scene_or_die("mars.toml");
    const int mars = scene.find("Mars");
    const int phobos = scene.find("Phobos");
    const int mex = scene.find("Mars Express");
    const int mro = scene.find("MRO");
    const int maven = scene.find("MAVEN");
    check(mars == 0 && phobos > 0 && mex > 0 && mro > 0 && maven > 0, "mars.toml bodies");
    auto tdb = [](const CalendarDateTime& c) { return (jd_from_calendar(c) - kJ2000Jd) * kSecondsPerDay; };
    auto at = [&](int body, double t) { return scene.icrf_state_at(body, t).position; };
    const double gm = scene.bodies[static_cast<size_t>(mars)].gm_km3_s2;
    const double r_mean = scene.bodies[static_cast<size_t>(mars)].equatorial_radius_km;

    // Altitude above the IAU ellipsoid (3396.19 x 3376.20 km, the scene's pck00011 pole, which
    // moves 0.1 deg a century: taken in 2015; latitude geocentric, close enough at these heights).
    const glm::dvec3 pole = scene.body_to_icrf(scene.bodies[static_cast<size_t>(mars)], tdb({2015, 1, 1, 0, 0, 0}))[2];
    auto ellipsoid_altitude = [&](const glm::dvec3& p) {
        const double r = glm::length(p);
        const double s = glm::dot(p, pole) / r;
        const double a = 3396.19;
        const double b = 3376.20;
        return r - a * b / std::sqrt(b * b * (1.0 - s * s) + a * a * s * s);
    };

    // MRO's aerobraking ([MROS], NASA/JPL 2006-09-13): 426 passes through the upper
    // atmosphere from early April to 2006-08-30, the low point 98-105 km above the
    // surface, the high point lowered from 45,000 to 486 km.
    int passes = 0;
    double lowest = 1e300;
    for (int month = 3; month <= 9; ++month) {
        for (const double tp : scene.periapsis_times(mro, tdb({2006, month, month == 3 ? 12 : 1, 0, 0, 0}),
                                                     tdb({2006, month + 1, 1, 0, 0, 0}))) {
            const double h = ellipsoid_altitude(at(mro, tp));
            passes += h < 150.0 ? 1 : 0;
            lowest = std::min(lowest, h);
        }
    }
    check(std::abs(passes - 426) <= 3, "MRO aerobraking passes", passes);
    check(lowest > 95.0 && lowest < 108.0, "MRO aerobraking low point (km above the ellipsoid)", lowest);
    std::printf("info: MRO aerobraking: %d passes, lowest %.1f km\n", passes, lowest);
    OrbitElements el;
    osculating_elements(scene.icrf_state_at(mro, tdb({2006, 3, 12, 0, 0, 0})), gm, &el);
    check(std::abs(el.apoapsis_km - r_mean - 45000.0) < 500.0, "MRO capture orbit apoapsis (km)", el.apoapsis_km - r_mean);
    osculating_elements(scene.icrf_state_at(mro, tdb({2006, 8, 31, 0, 0, 0})), gm, &el);
    check(std::abs(el.apoapsis_km - r_mean - 486.0) < 20.0, "MRO apoapsis after aerobraking (km)", el.apoapsis_km - r_mean);

    // MRO's science orbit ([MROS]): 250 to 316 km above the surface from 2006-09-11.
    double low = 1e300;
    double high = 0.0;
    for (double t = tdb({2006, 9, 15, 0, 0, 0}); t < tdb({2006, 9, 16, 0, 0, 0}); t += 30.0) {
        const double h = ellipsoid_altitude(at(mro, t));
        low = std::min(low, h);
        high = std::max(high, h);
    }
    check(std::abs(low - 250.0) < 20.0 && std::abs(high - 316.0) < 20.0, "MRO science orbit (km)", low);
    std::printf("info: MRO in September 2006: %.0f x %.0f km\n", low, high);

    // MRO: two windows of data only.
    for (const auto& [when, shown] : {std::pair{CalendarDateTime{2006, 6, 1, 0, 0, 0}, true},
                                      std::pair{CalendarDateTime{2015, 1, 1, 0, 0, 0}, false},
                                      std::pair{CalendarDateTime{2025, 6, 1, 0, 0, 0}, true}}) {
        scene.update(tdb(when));
        check(scene.bodies[static_cast<size_t>(mro)].visible == shown, "MRO shown only inside its windows",
              static_cast<double>(when.year));
    }

    // Mars Express at Phobos ([MEXP], ESA 2013-12-23): 45 km above the surface on
    // 2013-12-29 (07:09 UTC), 67 km in March 2010. Distances from Phobos' center
    // minus its mean radius (11.08 km, [SPHY]); Phobos is irregular by a few km.
    auto phobos_pass = [&](const CalendarDateTime& c) {
        const double t0 = tdb(c);
        double best = 1e300;
        for (double dt = -1800.0; dt <= 1800.0; dt += 0.5) {
            best = std::min(best, glm::length(at(mex, t0 + dt) - at(phobos, t0 + dt)));
        }
        return best - 11.08;
    };
    const double p2013 = phobos_pass({2013, 12, 29, 7, 9, 0});
    const double p2010 = phobos_pass({2010, 3, 3, 21, 0, 0});
    check(std::abs(p2013 - 45.0) < 6.0, "Mars Express - Phobos 2013 (km above the mean radius)", p2013);
    check(std::abs(p2010 - 67.0) < 6.0, "Mars Express - Phobos 2010 (km above the mean radius)", p2010);
    std::printf("info: Mars Express at Phobos: %.1f km (2013), %.1f km (2010)\n", p2013, p2010);

    // Mars Express' period, "6h 43m" (Horizons notes), in 2004.
    osculating_elements(scene.icrf_state_at(mex, tdb({2004, 6, 1, 0, 0, 0})), gm, &el);
    check(std::abs(el.period_s / 60.0 - 403.0) < 5.0, "Mars Express period (min)", el.period_s / 60.0);

    // MAVEN's science orbit (Horizons notes: 150 x 6200 km, 75 deg inclination), then
    // gone at the end of the tracking data.
    low = 1e300;
    high = 0.0;
    glm::dvec3 normal(0.0);
    for (double t = tdb({2015, 3, 1, 0, 0, 0}); t < tdb({2015, 3, 2, 0, 0, 0}); t += 30.0) {
        const State s = scene.icrf_state_at(maven, t);
        const double h = ellipsoid_altitude(s.position);
        low = std::min(low, h);
        high = std::max(high, h);
        normal = glm::normalize(glm::cross(s.position, s.velocity));
    }
    const double incl = std::acos(glm::dot(normal, pole)) / kDegToRad;
    check(std::abs(low - 150.0) < 25.0 && std::abs(high - 6200.0) < 150.0, "MAVEN science orbit (km)", low);
    // The osculating inclination swings between ~73 and ~75.7 deg over the years.
    check(std::abs(incl - 75.0) < 2.5, "MAVEN inclination (deg)", incl);
    std::printf("info: MAVEN in March 2015: %.0f x %.0f km, %.1f deg\n", low, high, incl);
    scene.update(tdb({2026, 1, 15, 0, 0, 0}));
    check(!scene.bodies[static_cast<size_t>(maven)].visible, "MAVEN gone after the tracking data");

    // Hope at Deimos ([HOPE], Spaceflight Now 2023-04-25): about 100 km on 2023-03-10;
    // Deimos' mean radius 6.2 km [SPHY]. Its orbit (Horizons notes): 55 h, 25 deg.
    const int hope = scene.find("Hope");
    const int deimos = scene.find("Deimos");
    double hope_deimos = 1e300;
    for (double t = tdb({2023, 3, 10, 2, 0, 0}); t < tdb({2023, 3, 10, 2, 40, 0}); t += 0.5) {
        hope_deimos = std::min(hope_deimos, glm::length(at(hope, t) - at(deimos, t)));
    }
    check(std::abs(hope_deimos - 6.2 - 100.0) < 10.0, "Hope - Deimos 2023 (km above the mean radius)",
          hope_deimos - 6.2);
    const State hs = scene.icrf_state_at(hope, tdb({2022, 6, 1, 0, 0, 0}));
    osculating_elements(hs, gm, &el);
    const double hope_incl = std::acos(glm::dot(glm::normalize(glm::cross(hs.position, hs.velocity)), pole)) / kDegToRad;
    check(std::abs(el.period_s / 3600.0 - 55.0) < 1.0, "Hope period (h)", el.period_s / 3600.0);
    check(std::abs(hope_incl - 25.0) < 3.0, "Hope inclination (deg)", hope_incl);
    std::printf("info: Hope: %.0f km from Deimos (2023), %.1f h, %.1f deg\n", hope_deimos, el.period_s / 3600.0,
                hope_incl);

    // Comet Siding Spring ([CSS], NASA/JPL 2014): dust at about 56 km/s; MRO's HiRISE
    // imaged it from about 138,000 km at the closest; the greatest risk from 90 minutes
    // after the closest approach for about 20 minutes, when the orbiters "took shelter
    // behind Mars" (MRO and MAVEN here): on the far side of Mars from the dust, which
    // moves with the comet's velocity.
    const int comet = scene.find("Siding Spring");
    double t_ca = 0.0;
    double d_ca = 1e300;
    double mro_comet = 1e300;
    for (double t = tdb({2014, 10, 19, 17, 0, 0}); t < tdb({2014, 10, 19, 20, 0, 0}); t += 1.0) {
        const double d = glm::length(at(comet, t));
        if (d < d_ca) {
            d_ca = d;
            t_ca = t;
        }
        mro_comet = std::min(mro_comet, glm::length(at(comet, t) - at(mro, t)));
    }
    const glm::dvec3 v_comet = scene.icrf_state_at(comet, t_ca).velocity;
    check(std::abs(glm::length(v_comet) - 56.0) < 0.5, "Siding Spring speed past Mars (km/s)", glm::length(v_comet));
    check(std::abs(mro_comet - 138000.0) < 2000.0, "MRO - Siding Spring closest (km)", mro_comet);
    std::printf("info: Siding Spring: %.0f km from Mars' center, %.2f km/s; %.0f km from MRO\n", d_ca,
                glm::length(v_comet), mro_comet);
    const int odyssey = scene.find("Odyssey");
    const glm::dvec3 stream = glm::normalize(v_comet);
    for (const int craft : {mro, maven, odyssey}) {
        bool shielded = true;
        for (double t = t_ca + 90.0 * 60.0; t <= t_ca + 110.0 * 60.0; t += 60.0) {
            const glm::dvec3 p = at(craft, t);
            const double along = glm::dot(p, stream);
            shielded = shielded && along > 0.0 && glm::length(p - along * stream) < r_mean;
        }
        check(shielded, ("behind Mars during the comet's dust: " + scene.bodies[static_cast<size_t>(craft)].name).c_str());
    }

    // Odyssey ([ODY], JPL mission status 2001-11-13): the period "currently 15 hours"
    // during aerobraking.
    osculating_elements(scene.icrf_state_at(odyssey, tdb({2001, 11, 13, 0, 0, 0})), gm, &el);
    check(std::abs(el.period_s / 3600.0 - 15.0) < 1.0, "Odyssey period on 2001-11-13 (h)", el.period_s / 3600.0);

    // TGO ([TGOA], ESA 2018-02): from a four-day orbit of about 98,000 x 200 km; aerobraking
    // complete on 2018-02-20; a near-circular two-hour orbit about 400 km up by the end of
    // April 2018.
    const int tgo = scene.find("TGO");
    osculating_elements(scene.icrf_state_at(tgo, tdb({2016, 10, 21, 0, 0, 0})), gm, &el);
    check(std::abs(el.apoapsis_km - r_mean - 98000.0) < 1500.0 && std::abs(el.period_s / 86400.0 - 4.0) < 0.5,
          "TGO capture orbit (apoapsis altitude, km)", el.apoapsis_km - r_mean);
    int dips_before = 0;
    int dips_after = 0;
    for (const double tp : scene.periapsis_times(tgo, tdb({2018, 2, 10, 0, 0, 0}), tdb({2018, 3, 10, 0, 0, 0}))) {
        if (glm::length(at(tgo, tp)) - r_mean < 150.0) {
            (tp < tdb({2018, 2, 21, 0, 0, 0}) ? dips_before : dips_after) += 1;
        }
    }
    check(dips_before > 50 && dips_after == 0, "TGO aerobraking ends on 2018-02-20", dips_after);
    low = 1e300;
    high = 0.0;
    for (double t = tdb({2018, 4, 30, 0, 0, 0}); t < tdb({2018, 4, 30, 6, 0, 0}); t += 30.0) {
        const double h = glm::length(at(tgo, t)) - r_mean;
        low = std::min(low, h);
        high = std::max(high, h);
    }
    osculating_elements(scene.icrf_state_at(tgo, tdb({2018, 4, 30, 0, 0, 0})), gm, &el);
    check(std::abs(low - 400.0) < 60.0 && std::abs(high - 400.0) < 60.0 && std::abs(el.period_s / 3600.0 - 2.0) < 0.1,
          "TGO science orbit (km)", low);
    std::printf("info: TGO at the end of April 2018: %.0f x %.0f km, %.2f h\n", low, high, el.period_s / 3600.0);
    // ESA's long-term predict joins the data fit on 2025-06-21 with a ~30 km step (Horizons,
    // 16:16-16:17 TDB): bounded in the table, and TGO is shown today.
    double tgo_step = 0.0;
    for (double t = tdb({2025, 6, 21, 16, 0, 0}); t < tdb({2025, 6, 21, 16, 30, 0}); t += 60.0) {
        const State s = scene.icrf_state_at(tgo, t);
        tgo_step = std::max(tgo_step, glm::length(at(tgo, t + 60.0) - s.position - s.velocity * 60.0));
    }
    check(tgo_step < 45.0, "TGO data fit / predict seam (km per minute)", tgo_step);
    scene.update(tdb({2026, 10, 1, 0, 0, 0}));
    check(scene.bodies[static_cast<size_t>(tgo)].visible, "TGO shown in 2026 (ESA predicts)");

    // MOM ([MOMI], ISRO 2014-09-24): 421.7 x 76,993.6 km, inclination 150 deg, 72 h 51 m 51 s.
    const int mom = scene.find("MOM");
    const State ms = scene.icrf_state_at(mom, tdb({2014, 9, 24, 12, 0, 0}));
    osculating_elements(ms, gm, &el);
    const double mom_incl = std::acos(glm::dot(glm::normalize(glm::cross(ms.position, ms.velocity)), pole)) / kDegToRad;
    check(std::abs(el.periapsis_km - r_mean - 421.7) < 30.0 && std::abs(el.apoapsis_km - r_mean - 76993.6) < 300.0,
          "MOM capture orbit (apoapsis altitude, km)", el.apoapsis_km - r_mean);
    check(std::abs(el.period_s / 3600.0 - 72.864) < 0.3, "MOM period (h)", el.period_s / 3600.0);
    check(std::abs(mom_incl - 150.0) < 1.0, "MOM inclination (deg)", mom_incl);
    std::printf("info: MOM after arrival: %.0f x %.0f km, %.2f h, %.1f deg\n", el.periapsis_km - r_mean,
                el.apoapsis_km - r_mean, el.period_s / 3600.0, mom_incl);

    // Tianwen-1 (tools/tianwen1/). A state vector from the orbiter's telemetry (D. Estevez,
    // tianwen1_orbit_change_2024.script: 2024-12-18 23:24:53.892 UTC, GMAT "MarsInertial" =
    // Mars' equator of J2000 with the IAU 2000 pole RA 317.68143, Dec 52.88650 deg)
    // checks the MOMAG positions, which are independent of it.
    const int tw1 = scene.find("Tianwen-1");
    auto utc = [&](const CalendarDateTime& c) { return tdb_from_utc_jd(jd_from_calendar(c)); };
    {
        const double ra0 = 317.68143 * kDegToRad;
        const double de0 = 52.88650 * kDegToRad;
        const glm::dvec3 z(std::cos(de0) * std::cos(ra0), std::cos(de0) * std::sin(ra0), std::sin(de0));
        const glm::dvec3 x(-std::sin(ra0), std::cos(ra0), 0.0);
        const glm::dvec3 y = glm::cross(z, x);
        const glm::dvec3 telemetry = 4138.192797433121 * x + 1565.165195862828 * y + 3704.175680474208 * z;
        const double miss = glm::length(at(tw1, utc({2024, 12, 18, 23, 24, 53}) + 0.892) - telemetry);
        check(miss < 15.0, "Tianwen-1 MOMAG positions vs a telemetry state vector (km)", miss);
        std::printf("info: Tianwen-1 vs telemetry 2024-12-18: %.1f km\n", miss);
    }
    // [AMSAT] (Estevez et al., GRCon 2021): the plane change turned the orbit from 11 to 87
    // deg (about the IAU 2000 pole; the scene's pck00011 pole, with its long-period terms, is
    // within 0.001 deg of it and gives 11.0 and 87.3); the apoapsis lowering of 2021-02-20
    // gave 280 x 84,600 km; the 2-sol orbit (2 Mars sidereal days, 49.25 h) before the
    // landing; the deorbit burn aimed below the surface.
    auto tw1_incl = [&](double t) {
        const State s = scene.icrf_state_at(tw1, t);
        return std::acos(glm::dot(glm::normalize(glm::cross(s.position, s.velocity)), pole)) / kDegToRad;
    };
    // The capture orbit ([AMSAT]: 390 x 182,000 km, 11 deg).
    osculating_elements(scene.icrf_state_at(tw1, utc({2021, 2, 12, 0, 0, 0})), gm, &el);
    check(std::abs(el.periapsis_km - r_mean - 390.0) < 30.0 && std::abs(el.apoapsis_km - r_mean - 182000.0) < 3000.0,
          "Tianwen-1 capture orbit (apoapsis altitude, km)", el.apoapsis_km - r_mean);
    check(std::abs(tw1_incl(utc({2021, 2, 12, 0, 0, 0})) - 11.0) < 2.0, "Tianwen-1 capture orbit (deg)",
          tw1_incl(utc({2021, 2, 12, 0, 0, 0})));
    std::printf("info: Tianwen-1 capture orbit: %.0f x %.0f km\n", el.periapsis_km - r_mean, el.apoapsis_km - r_mean);
    const double i_before = tw1_incl(utc({2021, 2, 15, 0, 0, 0}));
    const double i_after = tw1_incl(utc({2021, 2, 15, 10, 0, 0}));
    check(std::abs(i_before - 11.0) < 0.5 && std::abs(i_after - 87.0) < 0.5, "Tianwen-1 plane change (deg)", i_after);
    std::printf("info: Tianwen-1 plane change: %.1f -> %.1f deg\n", i_before, i_after);
    osculating_elements(scene.icrf_state_at(tw1, utc({2021, 2, 20, 16, 0, 0})), gm, &el);
    check(std::abs(el.periapsis_km - r_mean - 280.0) < 20.0 && std::abs(el.apoapsis_km - r_mean - 84600.0) < 400.0,
          "Tianwen-1 after the phasing burn (apoapsis altitude, km)", el.apoapsis_km - r_mean);
    osculating_elements(scene.icrf_state_at(tw1, utc({2021, 5, 14, 16, 0, 0})), gm, &el);
    check(std::abs(el.period_s / 3600.0 - 49.25) < 0.1, "Tianwen-1 2-sol orbit (h)", el.period_s / 3600.0);
    // ... and in between (the 02-23 burn as Estevez predicted it, and the 05-14 vector's arc).
    for (const CalendarDateTime& when : {CalendarDateTime{2021, 3, 1, 0, 0, 0}, CalendarDateTime{2021, 4, 20, 0, 0, 0}}) {
        osculating_elements(scene.icrf_state_at(tw1, utc(when)), gm, &el);
        check(std::abs(el.period_s / 3600.0 - 49.25) < 1.5, "Tianwen-1 in the 2-sol parking orbit (h)",
              el.period_s / 3600.0);
    }
    osculating_elements(scene.icrf_state_at(tw1, utc({2021, 5, 14, 19, 0, 0})), gm, &el);
    check(el.periapsis_km < r_mean, "Tianwen-1 deorbited with the lander (periapsis below the surface, km)",
          el.periapsis_km - r_mean);
    osculating_elements(scene.icrf_state_at(tw1, utc({2021, 5, 14, 21, 30, 0})), gm, &el);
    check(el.periapsis_km - r_mean > 300.0, "Tianwen-1 back in a safe orbit (km)", el.periapsis_km - r_mean);
    // The relay orbit from 2021-05-17 ([AMSAT]: 1/3 Mars sidereal day, 8.21 h).
    for (const CalendarDateTime& when : {CalendarDateTime{2021, 5, 25, 0, 0, 0}, CalendarDateTime{2021, 8, 1, 0, 0, 0}}) {
        osculating_elements(scene.icrf_state_at(tw1, utc(when)), gm, &el);
        check(std::abs(el.period_s / 3600.0 - 8.21) < 0.15, "Tianwen-1 relay orbit (h)", el.period_s / 3600.0);
    }
    // The remote-sensing orbit (CNSA via Space.com, 2021-11: 265 x ~10,700 km, 7 h 05 min)
    // in the MOMAG positions.
    osculating_elements(scene.icrf_state_at(tw1, utc({2021, 11, 20, 0, 0, 0})), gm, &el);
    check(std::abs(el.period_s / 3600.0 - 7.083) < 0.05 && std::abs(el.apoapsis_km - r_mean - 10700.0) < 150.0,
          "Tianwen-1 remote-sensing orbit (h)", el.period_s / 3600.0);
    // Shown only where there are data: from the approach (2021-02-08) on, but not in the
    // MOMAG gap of 2022-11 .. 2023-04.
    for (const auto& [when, shown] : {std::pair{CalendarDateTime{2021, 2, 8, 0, 0, 0}, false},
                                      std::pair{CalendarDateTime{2021, 3, 15, 0, 0, 0}, true},
                                      std::pair{CalendarDateTime{2021, 9, 1, 0, 0, 0}, true},
                                      std::pair{CalendarDateTime{2023, 1, 15, 0, 0, 0}, false},
                                      std::pair{CalendarDateTime{2024, 6, 1, 0, 0, 0}, true}}) {
        scene.update(utc(when));
        check(scene.bodies[static_cast<size_t>(tw1)].visible == shown, "Tianwen-1 shown only where there are data",
              static_cast<double>(when.year));
    }

    // The orbit insertion events' ghosts: the approach hyperbolas without the burns. Every
    // 1-minute sample through each burn is a knot (ephem.toml keep_all, tianwen1_bake.py
    // KEEP_ALL); the delta-v is read knot to knot (unexplained_delta_v), so the cubic's swings
    // across the sudden start and end of a burn, or a seam between trajectory files where a
    // burn ends (Mars Express 0.9 km, TGO 99 km), do not count. Published values:
    // - Tianwen-1, Estevez's burn model [TW1]: 888 s at 3 kN, Isp 321.6 s, from 4,425 kg:
    //   321.6 g0 ln(4425 / (4425 - 888 * 3000 / (321.6 g0))) = 668 m/s; 1.5% (a model).
    // - Odyssey 1,427 m/s, planned [ODYK]: 1% (the readout +0.5%).
    // - Mars Express 806.762 m/s, planned and reached by burning longer [MEXN]: 1%.
    // - MRO "about 1,000 meters per second" [MROK]: 2%.
    // - MAVEN "about 2,750 mph" (1,229 m/s) [MAVM]: 1%.
    // - MOM 1,099 m/s [MOMI]: 1%.
    // - TGO "more than 1.5 km/s" [TGOM]: up to 1,600 m/s (the burn in the table runs a minute
    //   longer than ESA's 139; the readout stops before ESA's next file, 50 m/s off).
    // - Hope: none published.
    const double isp_g0 = 321.6 * 9.80665e-3; // km/s
    auto rocket_after = [&](double s) { return isp_g0 * std::log(4425.0 / (4425.0 - s * 3.0 / isp_g0)); };
    const int tw1_insertion = check_insertion_ghost(scene, "Tianwen-1", rocket_after(888.0) * 1000.0, 0.015);
    if (tw1_insertion >= 0) {
        const SceneEvent::Ghost& g = scene.events[static_cast<size_t>(tw1_insertion)].ghost;
        scene.set_active_event(tw1_insertion);
        scene.update(g.burn_tdb + 444.0);
        double readout = 0.0;
        check(scene.burn_delta_v(tw1, &readout) && std::abs(readout / rocket_after(444.0) - 1.0) < 0.03,
              "Tianwen-1's delta-v readout halfway through the burn (km/s)", readout);
        scene.set_active_event(-1);
    }
    check_insertion_ghost(scene, "Odyssey", 1427.0, 0.01);
    check_insertion_ghost(scene, "Mars Express", 806.762, 0.01);
    check_insertion_ghost(scene, "MRO", 1000.0, 0.02);
    check_insertion_ghost(scene, "MAVEN", 2750.0 * 0.44704, 0.01);
    check_insertion_ghost(scene, "MOM", 1099.0, 0.01);
    const int tgo_insertion = check_insertion_ghost(scene, "TGO", 1500.0, 1.0);
    if (tgo_insertion >= 0) {
        const std::vector<double>& dv = scene.events[static_cast<size_t>(tgo_insertion)].ghost.delta_v_km_s;
        const double total = dv.empty() ? 0.0 : dv.back();
        check(total > 1.5 && total < 1.6, "TGO orbit insertion delta-v (km/s)", total);
    }
    check_insertion_ghost(scene, "Hope", 0.0, 0.0);
    // Hope's capture orbit right after the burn: the orbit reached, "as close to Mars' planetary
    // surface as 1,063 kilometres" and "as far away as 42,461" [EMMO], not Horizons' arc through
    // the burn, which holds the planned burn (1,081 x 50,073 km; tools/bake/rejoin_burn.py). 30 km:
    // altitudes here are above the equatorial radius.
    {
        const double t_capture = tdb({2021, 2, 9, 20, 0, 0});
        const State a = scene.icrf_state_at(hope, t_capture);
        const State b = scene.icrf_state_at(mars, t_capture);
        OrbitElements capture;
        osculating_elements(State{a.position - b.position, a.velocity - b.velocity}, gm, &capture);
        check(std::abs(capture.periapsis_km - r_mean - 1063.0) < 30.0 &&
                  std::abs(capture.apoapsis_km - r_mean - 42461.0) < 30.0,
              "Hope's capture orbit (apoapsis altitude, km)", capture.apoapsis_km - r_mean);
    }

    // 3I/ATLAS ([REN], Ren et al. 2026): closest to Mars 2025-10-03 04:38 UT at 0.194 au,
    // ~86 km/s; 0.208 au from Tianwen-1 at the second imaging epoch, 2025-10-01 16:00 UTC
    // (Table 1).
    const int atlas = scene.find("3I/ATLAS");
    check(atlas > 0, "mars.toml has 3I/ATLAS");
    double t_atlas = 0.0;
    double d_atlas = 1e300;
    for (double t = utc({2025, 10, 3, 3, 0, 0}); t < utc({2025, 10, 3, 6, 0, 0}); t += 10.0) {
        const double d = glm::length(at(atlas, t));
        if (d < d_atlas) {
            d_atlas = d;
            t_atlas = t;
        }
    }
    const double v_atlas = glm::length(scene.icrf_state_at(atlas, t_atlas).velocity);
    check(std::abs(t_atlas - utc({2025, 10, 3, 4, 38, 0})) < 120.0, "3I/ATLAS closest to Mars (s from 04:38 UT)",
          t_atlas - utc({2025, 10, 3, 4, 38, 0}));
    check(std::abs(d_atlas / kAuKm - 0.194) < 0.001 && std::abs(v_atlas - 86.0) < 1.0,
          "3I/ATLAS distance from Mars (au)", d_atlas / kAuKm);
    const double t_epoch2 = utc({2025, 10, 1, 16, 0, 0});
    const double atlas_tw1 = glm::length(at(atlas, t_epoch2) - at(tw1, t_epoch2)) / kAuKm;
    check(std::abs(atlas_tw1 - 0.208) < 0.001, "3I/ATLAS seen from Tianwen-1 (au)", atlas_tw1);
    std::printf("info: 3I/ATLAS: %.4f au from Mars, %.2f km/s; %.4f au from Tianwen-1 on 2025-10-01\n",
                d_atlas / kAuKm, v_atlas, atlas_tw1);

    // Gravity-assist flybys of 2025, against the published closest approaches (UTC; heights
    // above Mars' mean radius, 3389.5 km [SPHY]). Europa Clipper ([ECM]): 2025-03-01 17:57,
    // 884 km. Hera ([HERA]): within 300 km of Deimos (ESA's simulation: 12:10 GMT) and about
    // 5000 km of Mars (12:50 GMT; Horizons notes: "~6000 km from Mars").
    const int clipper = scene.find("Europa Clipper");
    const int hera = scene.find("Hera");
    check(clipper > 0 && hera > 0, "mars.toml flyby spacecraft");
    auto closest = [&](int body, int other, double t0) {
        std::pair<double, double> best{1e300, 0.0}; // distance, time
        for (double t = t0 - 1800.0; t <= t0 + 1800.0; t += 0.5) {
            const glm::dvec3 d = at(body, t) - (other >= 0 ? at(other, t) : glm::dvec3(0.0));
            best = std::min(best, std::pair{glm::length(d), t});
        }
        return best;
    };
    const double r_mars = 3389.5;
    const double t_ecm = utc({2025, 3, 1, 17, 57, 0});
    const auto [ecm_d, ecm_t] = closest(clipper, -1, t_ecm);
    check(std::abs(ecm_d - r_mars - 884.0) < 3.0 && std::abs(ecm_t - t_ecm) < 60.0,
          "Europa Clipper at Mars (km above the mean radius)", ecm_d - r_mars);
    const double t_hd = utc({2025, 3, 12, 12, 10, 0});
    const auto [hd_d, hd_t] = closest(hera, deimos, t_hd);
    check(std::abs(hd_d - 300.0) < 10.0 && std::abs(hd_t - t_hd) < 90.0, "Hera - Deimos (km from the center)", hd_d);
    const double t_hm = utc({2025, 3, 12, 12, 50, 0});
    const auto [hm_d, hm_t] = closest(hera, -1, t_hm);
    check(hm_d - r_mars > 5000.0 && hm_d - r_mars < 6000.0 && std::abs(hm_t - t_hm) < 90.0,
          "Hera at Mars (km above the mean radius)", hm_d - r_mars);
    // Europa Clipper's speed around the Sun ([ECM]): 24.5 km/s arriving, 22.5 leaving.
    auto helio_speed = [&](int body, double t) {
        return glm::length(scene.icrf_state_at(body, t).velocity - scene.sun_icrf_state_at(t).velocity);
    };
    const double v_in = helio_speed(clipper, t_ecm - 2.0 * kSecondsPerDay);
    const double v_out = helio_speed(clipper, t_ecm + 2.0 * kSecondsPerDay);
    check(std::abs(v_in - 24.5) < 0.1 && std::abs(v_out - 22.5) < 0.1, "Europa Clipper's speed around the Sun (km/s)",
          v_out);
    std::printf("info: Europa Clipper %.1f km up, %+.0f s, %.2f -> %.2f km/s around the Sun; Hera %.1f km from "
                "Deimos (%+.0f s), %.0f km above Mars (%+.0f s)\n",
                ecm_d - r_mars, ecm_t - t_ecm, v_in, v_out, hd_d, hd_t - t_hd, hm_d - r_mars, hm_t - t_hm);
    // Each shown only around its flyby.
    for (const auto& [body, when, shown] : {std::tuple{clipper, CalendarDateTime{2025, 3, 1, 0, 0, 0}, true},
                                            std::tuple{clipper, CalendarDateTime{2025, 3, 12, 0, 0, 0}, false},
                                            std::tuple{hera, CalendarDateTime{2025, 3, 12, 0, 0, 0}, true},
                                            std::tuple{hera, CalendarDateTime{2025, 3, 1, 0, 0, 0}, false}}) {
        scene.update(utc(when));
        check(scene.bodies[static_cast<size_t>(body)].visible == shown, "flyby spacecraft shown only near Mars",
              static_cast<double>(when.day));
    }
}

// Osculating heliocentric elements in the J2000 ecliptic (obliquity 84381.448", as the
// SBDB and older orbit catalogs use): a (au), e, i, node, argument of perihelion, mean
// anomaly (deg).
struct EclipticElements {
    double a_au, e, i_deg, node_deg, peri_deg, mean_anomaly_deg;
};

EclipticElements ecliptic_elements(const State& s, double gm)
{
    const double eps = 84381.448 / 3600.0 * kDegToRad;
    auto to_ecliptic = [&](const glm::dvec3& p) {
        return glm::dvec3(p.x, std::cos(eps) * p.y + std::sin(eps) * p.z, -std::sin(eps) * p.y + std::cos(eps) * p.z);
    };
    const glm::dvec3 r = to_ecliptic(s.position);
    const glm::dvec3 v = to_ecliptic(s.velocity);
    const double rn = glm::length(r);
    const double a = 1.0 / (2.0 / rn - glm::dot(v, v) / gm);
    const glm::dvec3 h = glm::cross(r, v);
    const glm::dvec3 ev = glm::cross(v, h) / gm - r / rn;
    const double e = glm::length(ev);
    const glm::dvec3 node(-h.y, h.x, 0.0);
    auto wrap = [](double deg) { return deg < 0.0 ? deg + 360.0 : deg; };
    double peri = std::acos(glm::dot(node, ev) / (glm::length(node) * e)) * kRadToDeg;
    peri = ev.z < 0.0 ? 360.0 - peri : peri;
    double nu = std::acos(std::clamp(glm::dot(ev, r) / (e * rn), -1.0, 1.0));
    nu = glm::dot(r, v) < 0.0 ? 2.0 * kPi - nu : nu;
    const double ecc_anomaly = 2.0 * std::atan(std::sqrt((1.0 - e) / (1.0 + e)) * std::tan(nu / 2.0));
    return {a / kAuKm,
            e,
            std::acos(h.z / glm::length(h)) * kRadToDeg,
            wrap(std::atan2(node.y, node.x) * kRadToDeg),
            peri,
            wrap((ecc_anomaly - e * std::sin(ecc_anomaly)) * kRadToDeg)};
}

// Venus: Venus Express and Akatsuki against the orbits ESA and JAXA published, and
// Zoozve against the elements and the behavior of Mikkola et al. 2004.
void test_venus_scene()
{
    Scene scene = load_scene_or_die("venus.toml");
    const int venus = scene.find("Venus");
    const int vex = scene.find("Venus Express");
    const int akatsuki = scene.find("Akatsuki");
    const int zoozve = scene.find("Zoozve");
    check(venus > 0 && vex > 0 && akatsuki > 0 && zoozve > 0, "venus.toml bodies");
    auto tdb = [](const CalendarDateTime& c) { return (jd_from_calendar(c) - kJ2000Jd) * kSecondsPerDay; };
    auto rel = [&](int body, double t) {
        const State a = scene.icrf_state_at(body, t);
        const State b = scene.icrf_state_at(venus, t);
        return State{a.position - b.position, a.velocity - b.velocity};
    };
    const double gm = scene.bodies[static_cast<size_t>(venus)].gm_km3_s2;
    const double rv = 6051.8; // [PCK]
    const glm::dvec3 pole = unit_toward(272.76, 67.16);
    OrbitElements el;

    // Venus Express ([VOI], ESA): captured on 2006-04-11 into a 9-day polar orbit out to
    // 330,000-350,000 km; after the burns of 2006-04-15 .. 05-06, the 24-hour orbit
    // ([VEXE], ESA 2014-12-16) out to 66,000 km over the south pole, low over the north
    // pole (pericentre latitude ~80 deg, ESA's mission catalog).
    osculating_elements(rel(vex, tdb({2006, 4, 14, 0, 0, 0})), gm, &el);
    check(el.apoapsis_km - rv > 3.2e5 && el.apoapsis_km - rv < 3.5e5 && std::abs(el.period_s / 86400.0 - 9.0) < 0.3,
          "Venus Express capture orbit (apoapsis altitude, km)", el.apoapsis_km - rv);
    const State vs = rel(vex, tdb({2007, 1, 1, 0, 0, 0}));
    osculating_elements(vs, gm, &el);
    const double vex_incl = std::acos(glm::dot(glm::normalize(glm::cross(vs.position, vs.velocity)), pole)) * kRadToDeg;
    check(std::abs(el.period_s / 3600.0 - 24.0) < 0.1 && std::abs(el.apoapsis_km - rv - 66000.0) < 1000.0,
          "Venus Express science orbit (apoapsis altitude, km)", el.apoapsis_km - rv);
    check(std::abs(vex_incl - 90.0) < 1.5, "Venus Express polar orbit (deg)", vex_incl);
    // The capture burn's ghost; "approximately 1251 ms-1" [VOI], 1%. Venus Express' path
    // steps by 3.8 km where ESA's next file starts, at the burn's end: read knot to knot, the
    // delta-v does not count the cubic's swing there (it would add 13%).
    check_insertion_ghost(scene, "Venus Express", 1251.0, 0.01);
    // Akatsuki's second try ([VCO]: 1,228 s): the mission catalog's 134.8 m/s does not fit its
    // own burns on the same four thrusters (56.3 m/s in 303 s in July 2015 would make some 230
    // m/s); the table reads 199. Not compared.
    check_insertion_ghost(scene, "Akatsuki", 0.0, 0.0);
    const std::vector<double> peri06 = scene.periapsis_times(vex, tdb({2006, 6, 1, 0, 0, 0}), tdb({2006, 6, 2, 12, 0, 0}));
    const glm::dvec3 pp = peri06.empty() ? glm::dvec3(0.0) : rel(vex, peri06[0]).position;
    const double peri_lat = peri06.empty() ? 0.0 : std::asin(glm::dot(glm::normalize(pp), pole)) * kRadToDeg;
    check(std::abs(peri_lat - 80.0) < 3.0 && glm::length(pp) - rv > 200.0 && glm::length(pp) - rv < 350.0,
          "Venus Express pericentre over the north pole (deg)", peri_lat);

    // Aerobraking ([VEXE]): about 130-135 km from 18 June to 11 July 2014, then raised to
    // about 460 km by 26 July, the period just over 22 hours.
    int dips = 0;
    double lowest = 1e300;
    for (const double tp : scene.periapsis_times(vex, tdb({2014, 6, 1, 0, 0, 0}), tdb({2014, 7, 20, 0, 0, 0}))) {
        const double h = glm::length(rel(vex, tp).position) - rv;
        dips += h < 140.0 ? 1 : 0;
        lowest = std::min(lowest, h);
    }
    check(dips >= 23 && dips <= 32 && lowest > 125.0 && lowest < 136.0, "Venus Express aerobraking (lowest, km)", lowest);
    osculating_elements(rel(vex, tdb({2014, 7, 27, 0, 0, 0})), gm, &el);
    check(std::abs(el.periapsis_km - rv - 460.0) < 15.0 && el.period_s / 3600.0 > 22.0 && el.period_s / 3600.0 < 23.0,
          "Venus Express raised after aerobraking (km)", el.periapsis_km - rv);
    std::printf("info: Venus Express: pericentre at %.1f deg N, %d passes below 140 km in 2014 (lowest %.1f km)\n",
                peri_lat, dips, lowest);
    // Contact was lost on 2014-11-28: shown until then.
    scene.update(tdb({2014, 11, 27, 0, 0, 0}));
    check(scene.bodies[static_cast<size_t>(vex)].visible, "Venus Express shown until contact was lost");
    scene.update(tdb({2014, 12, 15, 0, 0, 0}));
    check(!scene.bodies[static_cast<size_t>(vex)].visible, "Venus Express gone after 2014-11-28");

    // Akatsuki's missed orbit insertion ([VCO]): VOI-1 on 2010-12-06 23:49 UTC; it passed
    // 550 km above Venus (Wikipedia, "Akatsuki (spacecraft)"; the bake tool: 549 km at
    // 2010-12-07 00:01:16 TDB), into a 203-day solar orbit of 0.61 x 0.74 au.
    double miss = 1e300;
    for (double t = tdb({2010, 12, 6, 23, 50, 0}); t < tdb({2010, 12, 7, 0, 15, 0}); t += 1.0) {
        miss = std::min(miss, glm::length(rel(akatsuki, t).position));
    }
    check(std::abs(miss - rv - 550.0) < 10.0, "Akatsuki misses Venus in 2010 (altitude, km)", miss - rv);
    osculating_elements(scene.icrf_state_at(akatsuki, tdb({2011, 2, 1, 0, 0, 0})), kSunGmKm3S2, &el);
    check(std::abs(el.period_s / 86400.0 - 203.0) < 1.0 && std::abs(el.periapsis_km / kAuKm - 0.61) < 0.01 &&
              std::abs(el.apoapsis_km / kAuKm - 0.74) < 0.01,
          "Akatsuki's orbit around the Sun (days)", el.period_s / 86400.0);

    // Captured on the second try ([AK15], JAXA 2015-12-09): about 400 x 440,000 km, 13 days
    // 14 hours, about 3 deg to Venus' orbital plane, moving the way Venus rotates (against
    // the IAU north pole: Venus spins retrograde).
    const State as = rel(akatsuki, tdb({2015, 12, 9, 0, 0, 0}));
    osculating_elements(as, gm, &el);
    const glm::dvec3 an = glm::normalize(glm::cross(as.position, as.velocity));
    const State venus_helio = scene.icrf_state_at(venus, tdb({2015, 12, 9, 0, 0, 0}));
    const glm::dvec3 venus_normal = glm::normalize(glm::cross(venus_helio.position, venus_helio.velocity));
    const double tilt = 180.0 - std::acos(glm::dot(an, venus_normal)) * kRadToDeg;
    check(std::abs(el.periapsis_km - rv - 400.0) < 50.0 && std::abs(el.apoapsis_km - rv - 440000.0) < 10000.0,
          "Akatsuki's capture orbit (apoapsis altitude, km)", el.apoapsis_km - rv);
    check(std::abs(el.period_s / 86400.0 - (13.0 + 14.0 / 24.0)) < 0.5, "Akatsuki's capture orbit (days)",
          el.period_s / 86400.0);
    check(glm::dot(an, pole) < 0.0 && tilt < 6.0, "Akatsuki orbits the way Venus spins, near its orbital plane (deg)",
          tilt);
    // [VCO]: ~360,000 km and 10.5 days in early 2016; after PC1 (2016-04-04) ~370,000 km and
    // 10.8 days.
    osculating_elements(rel(akatsuki, tdb({2016, 1, 15, 0, 0, 0})), gm, &el);
    check(std::abs(el.apoapsis_km - rv - 360000.0) < 10000.0 && std::abs(el.period_s / 86400.0 - 10.5) < 0.15,
          "Akatsuki in early 2016 (days)", el.period_s / 86400.0);
    osculating_elements(rel(akatsuki, tdb({2016, 4, 10, 0, 0, 0})), gm, &el);
    check(std::abs(el.apoapsis_km - rv - 370000.0) < 10000.0 && std::abs(el.period_s / 86400.0 - 10.8) < 0.15,
          "Akatsuki after PC1 (days)", el.period_s / 86400.0);
    // The heliocentric table hands over to the Venus-centered one 600,000 km out, with a
    // step within the tolerances of the two tables and Venus' (2 + 2 + 5 km; 4 km).
    const double handover = tdb({2015, 12, 4, 0, 0, 0});
    const State h0 = rel(akatsuki, handover - 30.0);
    const double seam = glm::length(rel(akatsuki, handover + 30.0).position - h0.position - h0.velocity * 60.0);
    check(seam < 9.0, "Akatsuki's tables join (km)", seam);
    // The final reconstruction ends on 2024-04-26 (contact lost late in April 2024 [AKE]).
    scene.update(tdb({2024, 4, 20, 0, 0, 0}));
    check(scene.bodies[static_cast<size_t>(akatsuki)].visible, "Akatsuki shown until April 2024");
    scene.update(tdb({2024, 6, 1, 0, 0, 0}));
    check(!scene.bodies[static_cast<size_t>(akatsuki)].visible, "Akatsuki gone after its data");

    // Zoozve against Mikkola et al. 2004 ([MIK], Table 1: elements at JD 2452900.5 from
    // Lowell's astorb.dat of 2003, J2000 ecliptic), on today's JPL orbit.
    const EclipticElements z = ecliptic_elements(scene.icrf_state_at(zoozve, tdb_from_jd_tdb(2452900.5)), kSunGmKm3S2);
    check(std::abs(z.a_au - 0.72368870) < 1e-4 && std::abs(z.e - 0.410427) < 2e-4, "Zoozve a, e (Mikkola 2004)", z.a_au);
    check(std::abs(z.i_deg - 8.9766) < 0.01 && std::abs(z.node_deg - 231.674) < 0.01 &&
              std::abs(z.peri_deg - 355.528) < 0.01 && std::abs(z.mean_anomaly_deg - 318.875) < 0.1,
          "Zoozve angles (Mikkola 2004; mean anomaly, deg)", z.mean_anomaly_deg);
    // It stays farther than ~0.2 au from Venus ([MIK] Fig. 3) ...
    double nearest = 1e300;
    for (double t = tdb({1900, 1, 10, 0, 0, 0}); t < tdb({2099, 12, 20, 0, 0, 0}); t += 2.0 * kSecondsPerDay) {
        nearest = std::min(nearest, glm::length(rel(zoozve, t).position));
    }
    check(nearest > 0.2 * kAuKm, "Zoozve never near Venus (au)", nearest / kAuKm);
    // ... and, in the frame rotating with Venus, goes around it backward once a Venus year
    // ([MIK] Fig. 1; 224.7 days, [APX]).
    double turned = 0.0;
    double last_angle = 0.0;
    const double t_loop = tdb({2026, 1, 1, 0, 0, 0});
    for (int k = 0; k <= 900; ++k) {
        const double t = t_loop + k * 224.7 / 900.0 * kSecondsPerDay;
        const State v = scene.icrf_state_at(venus, t);
        const glm::dvec3 x = glm::normalize(v.position);
        const glm::dvec3 zn = glm::normalize(glm::cross(v.position, v.velocity));
        const glm::dvec3 d = rel(zoozve, t).position;
        const double angle = std::atan2(glm::dot(d, glm::cross(zn, x)), glm::dot(d, x));
        if (k > 0) {
            turned += std::remainder(angle - last_angle, kTwoPi);
        }
        last_angle = angle;
    }
    check(std::abs(turned * kRadToDeg + 360.0) < 10.0, "Zoozve loops backward around Venus (deg per Venus year)",
          turned * kRadToDeg);
    // Closest to the Earth on 2026-11-02 19:10 TDB at 0.0402 au ([SBDB], close approach data;
    // the scene's Earth is the Earth-Moon barycenter, ~4,700 km off).
    const int earth = scene.find("Earth");
    double t_near = 0.0;
    double d_near = 1e300;
    for (double t = tdb({2026, 10, 28, 0, 0, 0}); t < tdb({2026, 11, 8, 0, 0, 0}); t += 60.0) {
        const double d = glm::length(scene.icrf_state_at(zoozve, t).position - scene.icrf_state_at(earth, t).position);
        if (d < d_near) {
            d_near = d;
            t_near = t;
        }
    }
    check(std::abs(d_near / kAuKm - 0.0402) < 2e-4 && std::abs(t_near - tdb({2026, 11, 2, 19, 10, 0})) < 3600.0,
          "Zoozve passes the Earth in 2026 (au)", d_near / kAuKm);
    // Gravity-assist flybys against the published closest approaches (UTC; distances from
    // the center, or altitudes above 6051.8 km): Galileo 1990-02-10 05:58:48, 16,123 km up
    // ([GLLM], PDS); Cassini 1998-04-26 about 284 km ([CAS1], JPL: 6:52 a.m. PDT received on
    // the Earth, 7.5 light minutes away, so about 13:44:30 UTC) and 1999-06-24 20:29:55,
    // 602.6 km ([CAS2], JPL); Juice 2025-08-31 05:28:51, 11,141.1 km from the center
    // (Horizons notes); MESSENGER 2007-06-05 23:08, 338 km ([MESS]); Solar Orbiter 2021-08-09
    // 04:42, 7995 km ([DVF], ESA's prediction a week before); BepiColombo 2021-08-10
    // 13:51:54, 552 km ([BCV2], ESA).
    struct Flyby {
        int body;
        CalendarDateTime utc;
        double km;     // from the center
        double tol_s;  // on the time (Cassini 1: from the Earth-received time)
    };
    const int galileo = scene.find("Galileo");
    const int cassini = scene.find("Cassini");
    const int juice = scene.find("Juice");
    const int messenger = scene.find("MESSENGER");
    const int solo = scene.find("Solar Orbiter");
    const int bepi = scene.find("BepiColombo");
    check(galileo > 0 && cassini > 0 && juice > 0 && messenger > 0 && solo > 0 && bepi > 0,
          "venus.toml flyby spacecraft");
    for (const Flyby& f : {Flyby{galileo, {1990, 2, 10, 5, 58, 48}, rv + 16123.0, 30.0},
                           Flyby{cassini, {1998, 4, 26, 13, 44, 30}, rv + 284.0, 120.0},
                           Flyby{cassini, {1999, 6, 24, 20, 29, 55}, rv + 602.6, 30.0},
                           Flyby{messenger, {2007, 6, 5, 23, 8, 0}, rv + 338.0, 60.0},
                           Flyby{solo, {2021, 8, 9, 4, 42, 0}, rv + 7995.0, 90.0},
                           Flyby{bepi, {2021, 8, 10, 13, 51, 54}, rv + 552.0, 30.0},
                           Flyby{juice, {2025, 8, 31, 5, 28, 51}, 11141.1, 30.0}}) {
        const double t0 = tdb_from_utc_jd(jd_from_calendar(f.utc));
        double best = 1e300;
        double t_best = 0.0;
        for (double t = t0 - 1800.0; t <= t0 + 1800.0; t += 0.5) {
            const double d = glm::length(rel(f.body, t).position);
            if (d < best) {
                best = d;
                t_best = t;
            }
        }
        const std::string what = scene.bodies[static_cast<size_t>(f.body)].name + " at Venus " +
                                 std::to_string(f.utc.year) + " (km from the center)";
        check(std::abs(best - f.km) < 3.0 && std::abs(t_best - t0) < f.tol_s, what.c_str(), best);
        std::printf("info: %s %d: %.1f km up, %+.0f s from the published time\n",
                    scene.bodies[static_cast<size_t>(f.body)].name.c_str(), f.utc.year, best - rv, t_best - t0);
    }
    // The boosts, two days before to two days after, as the releases count them. Galileo
    // ([GLLV]: "almost 5,000 mph", 2.24 km/s): the orbital energy gained around the Sun, as a
    // change of speed at Venus' distance. Cassini ([CAS1]: "about 7 kilometers per second"):
    // the turn of its velocity relative to Venus, |v_out - v_in| (its energy gain is 3.6 km/s).
    const double t_gll = tdb({1990, 2, 10, 6, 0, 0});
    auto energy = [&](int body, double t) {
        const State s = scene.icrf_state_at(body, t);
        return 0.5 * glm::dot(s.velocity, s.velocity) - kSunGmKm3S2 / glm::length(s.position);
    };
    const double r_venus = glm::length(scene.icrf_state_at(venus, t_gll).position);
    auto speed_at_venus = [&](double e) { return std::sqrt(2.0 * (e + kSunGmKm3S2 / r_venus)); };
    const double gll_gain = speed_at_venus(energy(galileo, t_gll + 2.0 * kSecondsPerDay)) -
                            speed_at_venus(energy(galileo, t_gll - 2.0 * kSecondsPerDay));
    const double t_cas = tdb({1998, 4, 26, 13, 46, 0});
    const double cas_gain = glm::length(rel(cassini, t_cas + 2.0 * kSecondsPerDay).velocity -
                                        rel(cassini, t_cas - 2.0 * kSecondsPerDay).velocity);
    check(std::abs(gll_gain - 2.2) < 0.15, "Galileo's speed gain at Venus (km/s)", gll_gain);
    check(std::abs(cas_gain - 7.0) < 0.5, "Cassini's turn at Venus 1998 (km/s)", cas_gain);
    // MESSENGER ([MESF]: "nearly 15,000 mph", 6.7 km/s), counted the same way as Galileo's.
    const double t_mes = tdb({2007, 6, 5, 23, 9, 0});
    const double r_venus_mes = glm::length(scene.icrf_state_at(venus, t_mes).position);
    auto speed_at_venus_mes = [&](double e) { return std::sqrt(2.0 * (e + kSunGmKm3S2 / r_venus_mes)); };
    const double mes_turn = speed_at_venus_mes(energy(messenger, t_mes - 2.0 * kSecondsPerDay)) -
                            speed_at_venus_mes(energy(messenger, t_mes + 2.0 * kSecondsPerDay));
    check(mes_turn > 6.0 && mes_turn < 6.71, "MESSENGER's speed loss at Venus 2007 (km/s)", mes_turn);
    // Venus Express is in orbit meanwhile, far out ([MESS]: its "upstream measurements").
    const double vex_out = glm::length(rel(vex, t_mes).position);
    check(vex_out > 3.0e4, "Venus Express out in its orbit (km)", vex_out);
    // Solar Orbiter and BepiColombo, "more than 575 000 km apart" at their closest ([DVF]).
    double closest_pair = 1e300;
    for (double t = tdb({2021, 8, 5, 0, 0, 0}); t < tdb({2021, 8, 14, 0, 0, 0}); t += 600.0) {
        closest_pair = std::min(closest_pair, glm::length(scene.icrf_state_at(solo, t).position -
                                                          scene.icrf_state_at(bepi, t).position));
    }
    check(closest_pair > 5.75e5 && closest_pair < 5.85e5, "Solar Orbiter and BepiColombo at their closest (km)",
          closest_pair);
    // The short tables against mercury.toml's heliocentric cruise tables, near the flybys.
    {
        Scene merc = load_scene_or_die("mercury.toml");
        const int m_venus = merc.find("Venus");
        double worst = 0.0;
        for (const auto& [name, mid] :
             {std::pair{"MESSENGER", t_mes}, std::pair{"BepiColombo", tdb({2021, 8, 10, 13, 53, 0})}}) {
            const int here = scene.find(name);
            const int there = merc.find(name);
            for (double t = mid - 86400.0; t <= mid + 86400.0; t += 600.0) {
                const glm::dvec3 d = merc.icrf_state_at(there, t).position - merc.icrf_state_at(m_venus, t).position;
                worst = std::max(worst, glm::length(rel(here, t).position - d));
            }
        }
        check(worst < 10.0, "Venus flybys match mercury.toml's cruise (km)", worst);
    }
    // Each shown only around its flybys.
    for (const auto& [body, when, shown] :
         {std::tuple{cassini, CalendarDateTime{1998, 4, 26, 0, 0, 0}, true},
          std::tuple{cassini, CalendarDateTime{1999, 1, 1, 0, 0, 0}, false},
          std::tuple{cassini, CalendarDateTime{1999, 6, 24, 0, 0, 0}, true},
          std::tuple{galileo, CalendarDateTime{1990, 3, 1, 0, 0, 0}, false},
          std::tuple{juice, CalendarDateTime{2025, 8, 31, 0, 0, 0}, true},
          std::tuple{juice, CalendarDateTime{2026, 1, 1, 0, 0, 0}, false},
          std::tuple{messenger, CalendarDateTime{2007, 6, 5, 0, 0, 0}, true},
          std::tuple{messenger, CalendarDateTime{2007, 7, 1, 0, 0, 0}, false},
          std::tuple{solo, CalendarDateTime{2021, 8, 10, 0, 0, 0}, true},
          std::tuple{bepi, CalendarDateTime{2021, 8, 10, 0, 0, 0}, true},
          std::tuple{bepi, CalendarDateTime{2021, 9, 1, 0, 0, 0}, false}}) {
        scene.update(tdb(when));
        check(scene.bodies[static_cast<size_t>(body)].visible == shown, "flyby spacecraft shown only near Venus",
              static_cast<double>(when.year));
    }
    std::printf("info: speed gained at Venus: Galileo %.2f km/s, Cassini %.2f km/s; MESSENGER lost %.2f km/s, "
                "Venus Express %.0f km out; Solar Orbiter and BepiColombo %.0f km apart at the closest\n",
                gll_gain, cas_gain, mes_turn, vex_out, closest_pair);

    std::printf("info: Akatsuki missed Venus by %.0f km (2010); Zoozve: Table 1 a %.6f e %.6f M %.3f, %.3f au from "
                "Venus at the closest (1900-2099), %.4f au from the Earth (2026)\n",
                miss - rv, z.a_au, z.e, z.mean_anomaly_deg, nearest / kAuKm, d_near / kAuKm);
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

    // Philae (tools/philae/philae_bake.py), in 67P's body-fixed frame: released from Rosetta
    // 2014-11-12 08:35 UTC; on the drawn shape at its first touchdown (Agilkia, 15:34:04 UTC)
    // and at rest (Abydos, 17:31:17 UTC), 1.27 km apart, and still there in 2016.
    const int philae = scene.find("Philae");
    check(philae > 0, "rosetta.toml has Philae");
    const Body& nucleus = scene.bodies[static_cast<size_t>(comet)];
    auto body_frame = [&](double t) {
        return glm::transpose(scene.body_to_icrf(nucleus, t)) * (at(philae, t) - at(comet, t));
    };
    auto altitude = [&](const glm::dvec3& p) { return glm::length(p) - shape_radius(*nucleus.shape, p); };
    scene.update(utc({2014, 11, 12, 8, 30, 0}));
    check(!scene.bodies[static_cast<size_t>(philae)].visible, "Philae not yet released");
    scene.update(utc({2014, 11, 12, 8, 40, 0}));
    check(scene.bodies[static_cast<size_t>(philae)].visible, "Philae released");
    const double t_sep = utc({2014, 11, 12, 8, 36, 0});
    check(glm::length(at(philae, t_sep) - at(craft, t_sep)) < 0.1, "Philae leaves from Rosetta (km)",
          glm::length(at(philae, t_sep) - at(craft, t_sep)));
    const glm::dvec3 agilkia = body_frame(utc({2014, 11, 12, 15, 34, 4}));
    const glm::dvec3 abydos = body_frame(utc({2014, 11, 12, 17, 32, 0}));
    const glm::dvec3 abydos_2016 = body_frame(utc({2016, 9, 30, 12, 0, 0}));
    check(std::abs(altitude(agilkia)) < 0.01, "Philae's first touchdown on the shape (km)", altitude(agilkia));
    check(std::abs(altitude(abydos)) < 0.01, "Philae at rest on the shape (km)", altitude(abydos));
    check(std::abs(glm::length(abydos - agilkia) - 1.27) < 0.02, "Agilkia to Abydos (km)",
          glm::length(abydos - agilkia));
    check(glm::length(abydos_2016 - abydos) < 1e-6, "Philae stays at Abydos", glm::length(abydos_2016 - abydos));
    std::printf("info: Philae on the 67P shape: Agilkia %+.1f m, Abydos %+.1f m, %.3f km apart\n",
                1000.0 * altitude(agilkia), 1000.0 * altitude(abydos), glm::length(abydos - agilkia));
    // At rest it turns with the comet: the velocity is the rotation's.
    const double t_rest = utc({2015, 3, 1, 0, 0, 0});
    const glm::dvec3 v = scene.icrf_state_at(philae, t_rest).velocity - scene.icrf_state_at(comet, t_rest).velocity;
    const glm::dvec3 dp = (at(philae, t_rest + 1.0) - at(comet, t_rest + 1.0)) -
                          (at(philae, t_rest - 1.0) - at(comet, t_rest - 1.0));
    check(glm::length(v - dp / 2.0) < 1e-7 && glm::length(v) > 1e-4, "Philae's velocity at rest (km/s)",
          glm::length(v));

    // 67P's rotation (spin_table, fitted to ESA's attitude kernels) against the periods of
    // Jorda et al. 2016 (Icarus 277, 257) as summed up by Kokotanekova et al. 2017 (MNRAS 471,
    // 2974, Sec. 2.1.6): 12.4041 h from July to October 2014, 12.4304 h on 2015-05-19,
    // 12.305 h just before perihelion on 2015-08-10, and ESA's 12.055 h at the end.
    auto period_h = [&](double t) { return 360.0 / prime_meridian_rate(nucleus, t) * 24.0; };
    const double p2014 = period_h(utc({2014, 9, 1, 0, 0, 0}));
    const double p_may = period_h(utc({2015, 5, 19, 0, 0, 0}));
    const double p_aug = period_h(utc({2015, 8, 10, 0, 0, 0}));
    const double p_end = period_h(utc({2016, 9, 30, 0, 0, 0}));
    check(std::abs(p2014 - 12.4041) < 0.0005, "67P period in 2014 (h)", p2014);
    check(std::abs(p_may - 12.4304) < 0.003, "67P period on 2015-05-19 (h)", p_may);
    check(std::abs(p_aug - 12.305) < 0.01, "67P period on 2015-08-10 (h)", p_aug);
    check(std::abs(p_end - 12.055) < 0.001, "67P period in 2016 (h)", p_end);
    std::printf("info: 67P period %.4f h (2014), %.4f h (2015-05-19), %.4f h (2015-08-10), %.4f h (2016)\n",
                p2014, p_may, p_aug, p_end);
    // In August-September 2014 it is the rotation of the Cheops frame (Scholten et al. 2015,
    // CHEOPS_REF_FRAME_V1.PDF: pole (69.54, 64.11) deg, W = 114.69 deg + 696.543884683 deg/d),
    // whose axes the shape model uses.
    double cheops = 0.0;
    for (double t = utc({2014, 8, 1, 0, 0, 0}); t < utc({2014, 10, 1, 0, 0, 0}); t += 0.37 * kSecondsPerDay) {
        const glm::dmat3 model = iau_pole_frame(69.54 * kDegToRad, 64.11 * kDegToRad) *
                                 rotation_z(wrap_two_pi((114.69 + 696.543884683 * t / kSecondsPerDay) * kDegToRad));
        const glm::dmat3 d = glm::transpose(model) * scene.body_to_icrf(nucleus, t);
        cheops = std::max(cheops, std::acos(std::clamp((d[0][0] + d[1][1] + d[2][2] - 1.0) / 2.0, -1.0, 1.0)));
    }
    check(cheops * kRadToDeg < 0.3, "67P against the Cheops frame in Aug-Sep 2014 (deg)", cheops * kRadToDeg);
    // Rosetta's last state (Horizons, relative to 67P, 10:38:52 UTC, ~36 s before it touched
    // down): just above the drawn surface, which it could only reach with the rotation right.
    const double t_last = utc({2016, 9, 30, 10, 38, 50});
    const double rosetta_altitude =
        altitude(glm::transpose(scene.body_to_icrf(nucleus, t_last)) * (at(craft, t_last) - at(comet, t_last)));
    check(rosetta_altitude > -0.03 && rosetta_altitude < 0.06, "Rosetta above the shape before touchdown (km)",
          rosetta_altitude);
    std::printf("info: 67P against the Cheops frame %.2f deg; Rosetta %.0f m above the shape at 10:38:50 UTC\n",
                cheops * kRadToDeg, 1000.0 * rosetta_altitude);
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

void test_space_telescopes_scene()
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

    Scene scene = load_scene_or_die("space_telescopes.toml");
    check(scene.frame_is_rotating(), "Space Telescopes scene starts in the rotating frame");
    auto at = [](const char* utc) {
        double t = 0.0;
        parse_utc(utc, &t);
        return t;
    };
    auto pos = [&scene](const char* name) { return scene.bodies[static_cast<size_t>(scene.find(name))].world_position; };

    // The L2 telescopes over their spans [HZN]: the halo orbits of JWST, Euclid and Roman
    // reach ~800,000 km from L2 ([EWG]: 400,000 to 800,000 km for Euclid's), Gaia's
    // Lissajous orbit ~350,000 km ([EWG]).
    struct Span {
        const char* name;
        const char* from;
        const char* to;
        double max_km;
    };
    for (const Span& s : {Span{"JWST", "2022-07-01", "2031-09-01", 9.0e5}, Span{"Euclid", "2023-09-01", "2031-10-01", 9.0e5},
                          Span{"Roman", "2026-10-01", "2028-10-01", 9.0e5}, Span{"Gaia", "2014-02-01", "2025-02-01", 4.5e5}}) {
        double worst = 0.0;
        for (double t = at(s.from); t < at(s.to); t += 5.0 * kSecondsPerDay) {
            scene.update(t);
            // The Sun stays on the -x axis of the rotating frame.
            const glm::dvec3 sun = glm::normalize(scene.sun_position());
            check(sun.x < -0.9999999, "Sun on -x in the rotating frame", sun.x);
            const glm::dvec3 l2 = scene.markers[0].world_position;
            check(l2.x > 1.4e6 && std::abs(l2.y) < 1.0 && std::abs(l2.z) < 1.0, "L2 marker on +x", l2.x);
            worst = std::max(worst, glm::length(pos(s.name) - l2));
        }
        check(worst < s.max_km, (std::string(s.name) + " stays near L2 (km)").c_str(), worst);
    }

    // Roman goes straight into its halo orbit: within a month of launch it is as far
    // from the Earth as L2 [HZN].
    scene.update(at("2026-09-30"));
    check(glm::length(pos("Roman")) > 1.3e6, "Roman out at L2 a month after launch (km)", glm::length(pos("Roman")));

    // Gaia left L2 in March 2025 and returns to the Earth's vicinity every 14 years, at
    // ~10 million km ([HZN] notes, [GF]). The first return, in 2039, is as described;
    // the later ones in the 100-year prediction come closer (5.7 million km in 2054).
    double gaia_min = 1e30;
    double gaia_min_t = 0.0;
    for (double t = at("2026-01-01"); t < at("2050-01-01"); t += kSecondsPerDay) {
        scene.update(t);
        if (glm::length(pos("Gaia")) < gaia_min) {
            gaia_min = glm::length(pos("Gaia"));
            gaia_min_t = t;
        }
    }
    check(gaia_min > 1.0e7 && gaia_min < 1.1e7, "Gaia comes back to ~10 million km (km)", gaia_min);
    const double gaia_years = (gaia_min_t - at("2025-03-28")) / kSecondsPerDay / 365.25;
    check(std::abs(gaia_years - 14.0) < 1.0, "Gaia's first return 14 years after it left L2 (years)", gaia_years);

    // Kepler's heliocentric orbit: 372.5 days [HZN], a little longer than the year, so it
    // falls behind the Earth (behind: -y in the rotating frame, the Earth moving to +y),
    // 170 million km by its retirement ([KR] caption), and comes back in 2099 [HZN].
    {
        scene.update(at("2012-01-01"));
        const double y0 = pos("Kepler").y;
        scene.update(at("2016-01-01"));
        const double y1 = pos("Kepler").y;
        check(y0 < -1.0e7 && y1 < y0 - 3.0e7, "Kepler trails the Earth and falls further behind (km)", y1);
        scene.update(at("2018-10-30"));
        check(std::abs(glm::length(pos("Kepler")) - 1.7e8) < 0.05e8, "Kepler's distance at retirement (km)",
              glm::length(pos("Kepler")));
        double back = 1e30;
        for (double t = at("2099-06-01"); t < at("2099-09-01"); t += kSecondsPerDay) {
            scene.update(t);
            back = std::min(back, glm::length(pos("Kepler")));
        }
        check(std::abs(back - 1.4e6) < 0.1e6, "Kepler comes back past the Earth in 2099 (km)", back);

        const double t = at("2014-01-01");
        const double dt = 60.0;
        const int rotating = scene.active_frame();
        scene.set_active_frame(0); // inertial, for the velocity
        scene.update(t - dt);
        const glm::dvec3 i0 = pos("Kepler") - scene.sun_position();
        scene.update(t + dt);
        const glm::dvec3 i1 = pos("Kepler") - scene.sun_position();
        const glm::dvec3 r = 0.5 * (i0 + i1);
        const glm::dvec3 v = (i1 - i0) / (2.0 * dt);
        const double gm_sun = 1.32712440041279419e11; // [GM]
        const double a = 1.0 / (2.0 / glm::length(r) - glm::dot(v, v) / gm_sun);
        const double period_days = kTwoPi * std::sqrt(a * a * a / gm_sun) / kSecondsPerDay;
        check(std::abs(period_days - 372.5) < 1.5, "Kepler's heliocentric period (days)", period_days);
        scene.set_active_frame(rotating);
    }

    // TESS's lunar flyby (the event caption): 10,000 km from the Moon's center on
    // 2018-05-17, lifting the perigee from under 10,000 km to over 100,000 km [HZN].
    {
        double closest = 1e30;
        double closest_t = 0.0;
        for (double t = at("2018-05-16"); t < at("2018-05-18"); t += 300.0) {
            scene.update(t);
            const double d = glm::length(pos("TESS") - pos("Moon"));
            if (d < closest) {
                closest = d;
                closest_t = t;
            }
        }
        check(std::abs(closest - 1.0e4) < 0.05e4, "TESS's lunar flyby distance (km)", closest);
        check(std::abs(closest_t - at("2018-05-17 06:40")) < 3600.0, "TESS's lunar flyby on 2018-05-17");
        double perigee_before = 1e30;
        double perigee_after = 1e30;
        for (double t = at("2018-04-20"); t < at("2018-08-01"); t += 1800.0) {
            scene.update(t);
            double& p = t < closest_t ? perigee_before : perigee_after;
            p = std::min(p, glm::length(pos("TESS")));
        }
        check(perigee_before < 1.0e4 && perigee_after > 1.0e5, "the flyby lifts TESS's perigee (km)", perigee_after);
    }

    // TESS: 13.7-day orbits in 2:1 resonance with the Moon, the Moon leading or lagging
    // the apogee by about 90 degrees [TR14]. Apogees (hourly samples) from 2019 to 2026.
    {
        std::vector<double> apogee_times;
        double worst_angle = 0.0;
        double prev2 = 0.0;
        double prev = 0.0;
        for (double t = at("2019-01-01"); t < at("2026-01-01"); t += 3600.0) {
            scene.update(t);
            const double r = glm::length(pos("TESS"));
            if (prev > prev2 && prev >= r && prev2 > 0.0) {
                apogee_times.push_back(t - 3600.0);
                const glm::dvec3 tess = glm::normalize(pos("TESS"));
                const glm::dvec3 moon = glm::normalize(pos("Moon"));
                const double angle = std::acos(std::clamp(glm::dot(tess, moon), -1.0, 1.0)) * kRadToDeg;
                worst_angle = std::max(worst_angle, std::abs(angle - 90.0));
            }
            prev2 = prev;
            prev = r;
        }
        // The Moon and the Sun pull single orbits to anywhere between 12.7 and 14.7 days.
        const double mean_period =
            (apogee_times.back() - apogee_times.front()) / static_cast<double>(apogee_times.size() - 1) / kSecondsPerDay;
        check(std::abs(mean_period - 13.7) < 0.15, "TESS's mean period (days)", mean_period);
        check(worst_angle < 45.0, "the Moon about 90 deg from TESS's apogee (deg off)", worst_angle);
        std::printf("info: TESS %zu apogees 2019-2025, mean period %.3f d, Moon-apogee angle within %.1f deg of 90\n",
                    apogee_times.size(), mean_period, worst_angle);
    }
}

// DART and Didymos: Dimorphos' orbit before and after the impact (solution JPL 547,
// Naidu et al. 2024 Table 1), the ghost that keeps the old orbit and the 2.6 mm/s
// readout; DART reaching Dimorphos at the impact time at the impact speed of Daly et
// al. 2023 (Table 1), on the side of the impact site; LICIACube's flyby (Horizons notes).
void test_didymos_scene()
{
    Scene scene = load_scene_or_die("didymos.toml");
    const int didymos = scene.find("Didymos");
    const int dimorphos = scene.find("Dimorphos");
    const int dart = scene.find("DART");
    const int licia = scene.find("LICIACube");
    const int earth = scene.find("Earth");
    check(didymos > 0 && dimorphos > 0 && dart > 0 && licia > 0 && earth > 0, "didymos.toml bodies");
    if (didymos <= 0 || dimorphos <= 0 || dart <= 0 || licia <= 0 || earth <= 0) {
        return;
    }
    const double impact = 717506133.3653613; // 2022-09-26 23:15:33.365 TDB [S547]
    auto rel = [&](int a, int b, double t) {
        return scene.icrf_state_at(a, t).position - scene.icrf_state_at(b, t).position;
    };

    // Dimorphos' angle around Didymos in its orbit plane (the plane at the impact).
    const State moon = scene.icrf_state_at(dimorphos, impact - 60.0);
    const State primary = scene.icrf_state_at(didymos, impact - 60.0);
    const glm::dvec3 r0 = moon.position - primary.position;
    const glm::dvec3 e1 = glm::normalize(r0);
    const glm::dvec3 normal = glm::normalize(glm::cross(r0, moon.velocity - primary.velocity));
    const glm::dvec3 e2 = glm::cross(normal, e1);
    auto angle = [&](const glm::dvec3& r) { return std::atan2(glm::dot(r, e2), glm::dot(r, e1)); };
    // Mean period over `orbits` from t0: the time to sweep orbits * 360 deg (hourly
    // samples, then bisection); 0 if it does not get there within twice the time.
    auto mean_period_h = [&](int body, double t0, int orbits) {
        double swept = 0.0;
        double prev = angle(rel(body, didymos, t0));
        double t = t0;
        const double target = 2.0 * kPi * orbits;
        for (int hours = 0; hours < 24 * orbits; ++hours) {
            const double a = angle(rel(body, didymos, t + 3600.0));
            const double step = std::remainder(a - prev, 2.0 * kPi);
            if (swept + step >= target) {
                double lo = t;
                double hi = t + 3600.0;
                for (int i = 0; i < 40; ++i) {
                    const double mid = 0.5 * (lo + hi);
                    const double s = std::remainder(angle(rel(body, didymos, mid)) - prev, 2.0 * kPi);
                    (swept + s >= target ? hi : lo) = mid;
                }
                return (lo - t0) / orbits / 3600.0;
            }
            swept += step;
            prev = a;
            t += 3600.0;
        }
        return 0.0;
    };
    const double radius = glm::length(r0);
    check(std::abs(radius - 1.18921) < 0.001, "Dimorphos' distance before the impact (km)", radius);
    const double before = mean_period_h(dimorphos, impact - 10.0 * kSecondsPerDay, 10);
    check(std::abs(before - 11.921493) < 2.0 / 3600.0, "Dimorphos' period before the impact (h)", before);
    // Settled at 11.3674 h. Single orbits swing with the eccentricity, the libration and
    // the apsidal precession (6.7 deg/d): the mean over two turns of the apsides (107 days)
    // from a month after the impact.
    const double settled = mean_period_h(dimorphos, impact + 30.0 * kSecondsPerDay, 226);
    check(std::abs(settled - 11.3674) < 3.0 / 3600.0, "Dimorphos' settled period (h)", settled);
    // Turning with Dimorphos (synchronous: x toward Didymos, z along the orbit pole)
    // matches Horizons' rotation for it before the impact (pole (69.70029, -72.69527),
    // W = 73.38929405 + 724.7244051 d + 1.0256e-6 d^2 deg) and its frame after.
    {
        double worst = 0.0;
        for (const double days : {-200.0, -10.0, -0.01}) {
            const double t = impact + days * kSecondsPerDay;
            const double d = t / kSecondsPerDay;
            const glm::dmat3 horizons =
                iau_pole_frame(69.70029 * kDegToRad, -72.69527 * kDegToRad) *
                rotation_z(std::fmod(73.38929405 + 724.7244051 * d + 1.0256e-6 * d * d, 360.0) * kDegToRad);
            const glm::dmat3 axes = scene.body_to_icrf(scene.bodies[static_cast<size_t>(dimorphos)], t);
            for (int c = 0; c < 3; ++c) {
                worst = std::max(worst, std::acos(std::min(1.0, glm::dot(axes[c], horizons[c]))) * kRadToDeg);
            }
        }
        check(worst < 0.05, "Dimorphos turns as in Horizons before the impact (deg)", worst);
    }

    // The ghost keeps the old orbit: within 0.2 deg of the circle at n0 after three days,
    // while the real moon is 100 deg ahead.
    int event = -1;
    for (size_t i = 0; i < scene.events.size(); ++i) {
        if (scene.events[i].ghost.body == dimorphos) {
            event = static_cast<int>(i);
        }
    }
    check(event >= 0, "Dimorphos has a ghost");
    if (event >= 0) {
        const SceneEvent::Ghost& g = scene.events[static_cast<size_t>(event)].ghost;
        const double n0 = 1.46401901971247768710e-04; // [S547]
        const double t_end = g.end_tdb - 60.0;
        const glm::dvec3 ghost = scene.icrf_state_at(g.ghost, t_end).position - scene.icrf_state_at(didymos, t_end).position;
        const double ghost_off =
            std::abs(std::remainder(angle(ghost) - n0 * (t_end - (impact - 60.0)), 2.0 * kPi)) * kRadToDeg;
        check(ghost_off < 0.2, "the ghost on the old orbit after three days (deg)", ghost_off);
        const double ahead = std::remainder(angle(rel(dimorphos, didymos, t_end)) - angle(ghost), 2.0 * kPi) * kRadToDeg;
        check(ahead > 90.0 && ahead < 115.0, "Dimorphos ahead of its ghost after three days (deg)", ahead);
        // The impulse: 2.630 mm/s along the path and 0.33 mm/s outward [S547].
        const double dv = g.delta_v_km_s.back() * 1e6;
        const double expected = std::hypot(2.630, 0.330);
        check(std::abs(dv - expected) < 0.03, "DART's push on Dimorphos (mm/s)", dv);
        scene.set_active_event(-1);
        scene.update(impact + 3600.0);
        double readout = 0.0;
        check(!scene.bodies[static_cast<size_t>(g.ghost)].visible && !scene.burn_delta_v(dimorphos, &readout),
              "Dimorphos' ghost and readout only during its event");
        scene.set_active_event(event);
        scene.update(impact - 5.0);
        double total = 0.0;
        check(scene.burn_delta_v(dimorphos, &readout, &total) && readout * 1e6 < 0.05 && std::abs(total * 1e6 - dv) < 1e-9,
              "no push before the impact (mm/s)", readout * 1e6);
        scene.update(impact + 3600.0);
        check(scene.bodies[static_cast<size_t>(g.ghost)].visible && scene.burn_delta_v(dimorphos, &readout) &&
                  std::abs(readout * 1e6 - dv) < 1e-6,
              "the whole push after the impact (mm/s)", readout * 1e6);
        scene.set_active_event(-1);
        std::printf("info: Dimorphos %.1f m from Didymos, %.6f h before the impact, %.4f h settled; push %.3f mm/s; "
                    "the ghost %.3f deg off the old orbit after three days, the moon %.1f deg ahead\n",
                    radius * 1e3, before, settled, dv, ghost_off, ahead);
    }

    // DART: at Dimorphos' center at the impact time, at 6.1449 km/s [DALY]; gone after.
    {
        const glm::dvec3 miss = rel(dart, dimorphos, impact - 1e-3);
        check(glm::length(miss) < 0.02, "DART at Dimorphos at the impact time (km)", glm::length(miss));
        const glm::dvec3 v =
            scene.icrf_state_at(dart, impact - 10.0).velocity - scene.icrf_state_at(dimorphos, impact - 10.0).velocity;
        check(std::abs(glm::length(v) - 6.1449) < 0.001, "DART's impact speed (km/s)", glm::length(v));
        check(!scene.bodies[static_cast<size_t>(dart)].motion->valid_at(impact + 1.0), "DART gone after the impact");
        // Daly et al. put the impact site at 8.84 S, 264.30 E, 25 m from the center of
        // figure: the point facing DART, which aims at the center here, is within about
        // asin(25 / 85) = 17 deg of it. 264 E is the leading side: the impact slowed
        // Dimorphos.
        const glm::dmat3 axes = scene.body_to_icrf(scene.bodies[static_cast<size_t>(dimorphos)], impact);
        const glm::dvec3 from = glm::transpose(axes) * glm::normalize(-v);
        const double lat = std::asin(from.z) * kRadToDeg;
        const double lon = std::fmod(std::atan2(from.y, from.x) * kRadToDeg + 360.0, 360.0);
        const double site_deg =
            std::acos(std::clamp(glm::dot(from, glm::dvec3(std::cos(-8.84 * kDegToRad) * std::cos(264.30 * kDegToRad),
                                                          std::cos(-8.84 * kDegToRad) * std::sin(264.30 * kDegToRad),
                                                          std::sin(-8.84 * kDegToRad))),
                                 -1.0, 1.0)) *
            kRadToDeg;
        check(site_deg < 20.0, "DART came from the impact site's side (deg off)", site_deg);
        const double earth_km = glm::length(rel(earth, dimorphos, impact));
        check(earth_km > 1.1e7 && earth_km < 1.15e7, "the impact 11 million km from Earth (km)", earth_km);
        // The heliocentric cruise hands over to the table relative to the system a day
        // before the impact without a visible step.
        double step = 0.0;
        for (double t = impact - 2.5 * kSecondsPerDay; t < impact - 0.5 * kSecondsPerDay; t += 600.0) {
            const glm::dvec3 a = scene.icrf_state_at(dart, t).position;
            const State s = scene.icrf_state_at(dart, t + 600.0);
            step = std::max(step, glm::length(s.position - a - 600.0 * s.velocity));
        }
        check(step < 5.0, "DART's handover to the relative table (km per 10 min off a straight line)", step);
        std::printf("info: DART reaches Dimorphos at %.4f km/s from %.1f S, %.1f E (Daly: 8.84 S, 264.30 E), %.0f km "
                    "from Earth\n",
                    glm::length(v), -lat, lon, earth_km);
    }

    // LICIACube: 57.7 km from Dimorphos at 23:17:11.5 UTC (Horizons notes).
    {
        double best = 1e300;
        double best_t = 0.0;
        for (double t = impact + 60.0; t < impact + 300.0; t += 0.1) {
            const double d = glm::length(rel(licia, dimorphos, t));
            if (d < best) {
                best = d;
                best_t = t;
            }
        }
        const double expected_t = impact + (3.0 * 60.0 + 11.5 - 24.183);
        check(std::abs(best - 57.7) < 0.2, "LICIACube's closest approach to Dimorphos (km)", best);
        check(std::abs(best_t - expected_t) < 2.0, "LICIACube's closest approach time (s off)", best_t - expected_t);
    }
}

} // namespace

void run_mission_tests()
{
    test_solar_system_missions();
    test_solar_system_interstellar();
    test_parker_scene();
    test_jupiter_missions();
    test_earth_moon_missions();
    test_saturn_missions();
    test_mercury_scene();
    test_mars_scene();
    test_venus_scene();
    test_isee3_scene();
    test_rosetta_scene();
    test_halley_armada();
    test_space_telescopes_scene();
    test_didymos_scene();
}
