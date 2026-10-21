// Tests for the scene framework: loader, simulation, labels, camera director.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/ephemeris.hpp"
#include "scene/camera.hpp"
#include "scene/camera_director.hpp"
#include "scene/label_layout.hpp"
#include "scene/scene_loader.hpp"
#include "scene/simulation.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace astraxis;

namespace {

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

void test_scene_loader()
{
    Scene scene = load_jupiter();
    check(scene.bodies.size() == 13, "jupiter.toml has 13 bodies (with Juno's ghost)",
          static_cast<double>(scene.bodies.size()));
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

// Event captions are placards: a few sentences, cleanly joined from the TOML's
// continued lines. Every event has one.
// A jump starts the caption's clock; loading a scene clears it.
void test_event_captions()
{
    for (const auto& [stem, scene] : all_scenes()) {
        for (const SceneEvent& e : scene->events) {
            if (e.caption.empty()) {
                check(false, (stem + " / " + e.name + " has a caption").c_str());
                continue;
            }
            const std::string what = stem + " / " + e.name + ": caption";
            size_t words = 0;
            bool in_word = false;
            for (const char c : e.caption) {
                words += (c != ' ' && !in_word) ? 1 : 0;
                in_word = c != ' ';
            }
            check(words >= 15 && words <= 80, (what + " has 15-80 words").c_str(), static_cast<double>(words));
            check(e.caption.front() != ' ' && e.caption.back() != ' ' &&
                      e.caption.find("  ") == std::string::npos && e.caption.find('\n') == std::string::npos,
                  (what + " is one clean paragraph").c_str());
        }
    }

    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/solar_system.toml", &error), "solar_system loads");
    check(sim.last_event() == -1, "no event before a jump");
    sim.jump_to_event(3);
    sim.update(0.5);
    check(sim.last_event() == 3 && std::abs(sim.event_age() - 0.5) < 1e-12, "a jump starts the caption clock",
          sim.event_age());
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error) && sim.last_event() == -1,
          "loading a scene clears the event");
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
    // Circular orbits: the fade radius is the distance from the host.
    for (int i = 1; i < 6; ++i) {
        scene.bodies[static_cast<size_t>(i)].fade_radius_km = glm::length(
            scene.bodies[static_cast<size_t>(i)].world_position - scene.bodies[static_cast<size_t>(hosts[i])].world_position);
    }

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

    // An eccentric orbit fades by its size, not by where the body is on it:
    // S2 (e = 0.88) is 15 times closer at pericentre than at apocentre, but its
    // fade radius stays the semi-major axis (125.058 mas at 8246.7 pc, [GR20]).
    // The Newtonian osculating a of the relativistic orbit is ~1% short at
    // pericentre (GM / r c^2 ~ 3.5e-4, amplified by 2 / (1 - e)).
    Scene sgr = load_scene_or_die("sgr_a.toml");
    const int s2 = sgr.find("S2");
    const double a_km = 125.058e-3 * 8246.7 * kAuKm;
    const double t_peri = tdb_from_julian_year(2018.379);
    for (const double t : {t_peri, t_peri - 0.5 * 16.0455 * kDaysPerJulianYear * kSecondsPerDay}) {
        sgr.update(t);
        const Body& b = sgr.bodies[static_cast<size_t>(s2)];
        check(std::abs(b.fade_radius_km / a_km - 1.0) < 0.02, "S2 fade radius is its semi-major axis",
              b.fade_radius_km / a_km);
    }
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

int find_frame(const Scene& scene, std::string_view name)
{
    for (size_t i = 0; i < scene.frames.size(); ++i) {
        if (scene.frames[i].name == name) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// East longitude (degrees, -180..180) of a display-frame position in a body-fixed frame.
double fixed_longitude(const glm::dvec3& p)
{
    return std::atan2(p.y, p.x) * kRadToDeg;
}

// Body-fixed frames turn with the origin's IAU axes: the body itself stands
// still, and its moons move across it at their rate relative to its rotation.
void test_body_fixed_frames()
{
    {
        Scene scene = load_scene_or_die("mars.toml");
        const int frame = find_frame(scene, "Mars-fixed");
        const int mars = scene.find("Mars");
        const int phobos = scene.find("Phobos");
        const int deimos = scene.find("Deimos");
        check(frame > 0 && scene.frames[static_cast<size_t>(frame)].type == DisplayFrame::Type::BodyFixed,
              "mars.toml has a Mars-fixed frame");
        scene.set_active_frame(frame);
        check(scene.frame_is_rotating(), "a body-fixed frame turns");
        const double t0 = tdb_from_jd_tdb(2460827.5); // 2025-06-01
        scene.update(t0);
        const Body& m = scene.bodies[static_cast<size_t>(mars)];
        double off_identity = 0.0;
        for (int c = 0; c < 3; ++c) {
            off_identity = std::max(off_identity, glm::length(m.orientation[c] - glm::dmat3(1.0)[c]));
        }
        check(off_identity < 1e-12, "Mars' own axes are the frame's", off_identity);
        check(glm::length(m.world_position) < 1e-6, "Mars at the origin of its fixed frame");
        check(glm::length(scene.up_axis() - glm::dvec3(0.0, 0.0, 1.0)) < 1e-12, "up is Mars' north pole");
        const glm::dmat3 axes = scene.body_to_icrf(m, t0);
        const Body& p = scene.bodies[static_cast<size_t>(phobos)];
        check(glm::length(p.world_position - glm::transpose(axes) * (p.icrf_position - m.icrf_position)) < 1e-6,
              "positions in the Mars-fixed frame");

        // Over an hour: Phobos gains on Mars' rotation (pole and prime meridian rates
        // of the scene's [PCK] values: 1128.845 - 350.892 = 777.95 deg/day eastward),
        // Deimos falls behind (285.162 - 350.892 = -65.73 deg/day).
        const glm::dvec3 phobos0 = p.world_position;
        const glm::dvec3 deimos0 = scene.bodies[static_cast<size_t>(deimos)].world_position;
        scene.update(t0 + 3600.0);
        const double phobos_rate =
            wrap_pi((fixed_longitude(p.world_position) - fixed_longitude(phobos0)) * kDegToRad) * kRadToDeg * 24.0;
        const double deimos_rate =
            wrap_pi((fixed_longitude(scene.bodies[static_cast<size_t>(deimos)].world_position) -
                     fixed_longitude(deimos0)) *
                    kDegToRad) *
            kRadToDeg * 24.0;
        check(std::abs(phobos_rate / 777.95 - 1.0) < 0.1, "Phobos moves east across Mars", phobos_rate);
        check(std::abs(deimos_rate / -65.73 - 1.0) < 0.1, "Deimos drifts west across Mars", deimos_rate);

        // Phobos' orbit becomes a history trail of two Mars days: about 4.4 laps of
        // the same near-circle, starting at Phobos.
        std::vector<glm::dvec3> points;
        std::vector<float> fades;
        scene.trail(phobos, scene.time(), 512, points, fades);
        const double a = 9375.0; // Phobos' a_km in mars.toml [ELEM]
        double worst = 0.0;
        double turned = 0.0;
        for (size_t k = 0; k < points.size(); ++k) {
            worst = std::max(worst, std::abs(glm::length(points[k]) / a - 1.0));
            if (k > 0) {
                turned += std::abs(wrap_pi((fixed_longitude(points[k]) - fixed_longitude(points[k - 1])) * kDegToRad));
            }
        }
        check(!points.empty() && glm::length(points[0] - p.world_position) < 1.0, "Phobos' trail starts at Phobos");
        check(worst < 0.03, "Phobos' trail circles Mars", worst);
        const double laps = turned / kTwoPi;
        check(std::abs(laps / (2.0 * 777.95 / 350.892) - 1.0) < 0.05, "two Mars days of Phobos' trail", laps);
    }
    {
        // Pluto and Charon are mutually locked: in the Pluto-fixed frame Charon stays put,
        // over the mean sub-Charon meridian (the zero of Pluto's longitudes; Charon's [PCK]
        // prime meridian is Pluto's plus 180 deg, facing Pluto). It holds from 1990 to 2050:
        // the PCK rotation and the Horizons orbit agree.
        Scene scene = load_scene_or_die("solar_system.toml");
        const int frame = find_frame(scene, "Pluto-fixed");
        const int charon = scene.find("Charon");
        check(frame > 0, "solar_system.toml has a Pluto-fixed frame");
        scene.set_active_frame(frame);
        for (double jd : {2447892.5, 2457218.0, 2469807.5}) { // 1990, 2015-07-14, 2050
            double lon_min = 1e9;
            double lon_max = -1e9;
            double lat_max = 0.0;
            for (int k = 0; k <= 8; ++k) { // one rotation
                scene.update(tdb_from_jd_tdb(jd) + k * 6.387 / 8.0 * kSecondsPerDay);
                const glm::dvec3 c = scene.bodies[static_cast<size_t>(charon)].world_position;
                const double lon = fixed_longitude(c);
                lon_min = std::min(lon_min, lon);
                lon_max = std::max(lon_max, lon);
                lat_max = std::max(lat_max, std::abs(std::asin(c.z / glm::length(c))) * kRadToDeg);
            }
            std::printf("info: Charon over Pluto from JD %.1f: east longitude %.2f .. %.2f, latitude < %.2f deg\n",
                        jd, lon_min, lon_max, lat_max);
            check(lon_max - lon_min < 0.2, "Charon stands still over Pluto", lon_max - lon_min);
            check(std::abs(lon_min) < 3.0 && std::abs(lon_max) < 3.0, "Charon over Pluto's prime meridian", lon_min);
            check(lat_max < 1.0, "Charon over Pluto's equator", lat_max);
        }

        // Only the Pluto system leaves trails: a planet's would sweep a circle tens of au
        // wide every 6.4 days, and Pluto's own stays at the origin.
        std::vector<glm::dvec3> points;
        std::vector<float> fades;
        scene.trail(charon, scene.time(), 256, points, fades);
        check(points.size() > 1, "Charon's trail in the Pluto-fixed frame");
        scene.trail(scene.find("Neptune"), scene.time(), 256, points, fades);
        check(points.empty(), "no Neptune trail in the Pluto-fixed frame");
        scene.trail(scene.find("Pluto"), scene.time(), 256, points, fades);
        check(points.empty(), "no trail for the frame's own body");
    }
}

// The sky turns with a rotating frame: it fades out where that would spin it by (Saturn
// turning with Janus at 2 days per second: over 1,000 deg/s), and stays where the frame
// turns slowly (the Earth-Moon frame at the same warp: 26 deg/s).
void test_sky_fades_in_fast_frames()
{
    auto settle = [](Simulation& sim) {
        for (int i = 0; i < 20; ++i) {
            sim.update(0.1);
        }
    };
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/saturn.toml", &error), "saturn.toml loads");
    settle(sim);
    check(sim.sky_visibility() == 1.0, "the sky shows in an inertial frame", sim.sky_visibility());
    const Scene& scene = sim.scene();
    const int janus_frame = find_frame(scene, "Saturn-Janus rotating");
    check(janus_frame > 0, "Saturn scene has the Janus frame");
    sim.set_frame(janus_frame);
    sim.clock().warp = 2.0 * kSecondsPerDay;
    // At a fixed date (the scene starts now), averaged over one orbit: the rate at any
    // one moment swings by about 2e = 1.4 percent around Janus' eccentric orbit. Mean
    // motion: the synchronous spin rate of pck00011 (BODY610_PM, 518.2359876 deg/day,
    // fitted near 1981); each orbit swap changes it by ~0.25 deg/day (Epimetheus'
    // BODY611_PM is 518.4907239).
    constexpr double kJanusDegPerDay = 518.2359876;
    double t0 = 0.0;
    check(parse_utc("2026-01-01", &t0), "parse the Janus frame date");
    sim.clock().t_tdb = t0;
    constexpr int kSamples = 96;
    const double period = 360.0 / kJanusDegPerDay * kSecondsPerDay;
    double rate = 0.0;
    double lo = 1e300, hi = 0.0;
    for (int i = 0; i < kSamples; ++i) {
        const double r = scene.frame_turn_rate(t0 + period * i / kSamples) * kRadToDeg * kSecondsPerDay;
        rate += r / kSamples;
        lo = std::min(lo, r);
        hi = std::max(hi, r);
    }
    check(std::abs(rate - kJanusDegPerDay) < 0.5, "the Janus frame turns at its mean motion (deg/day)", rate);
    std::printf("info: Janus frame over an orbit from 2026-01-01: mean %.2f, %.2f .. %.2f deg/day\n", rate, lo, hi);
    settle(sim);
    check(sim.sky_visibility() == 0.0, "the sky fades out while it would spin by", sim.sky_visibility());
    sim.set_frame(0);
    settle(sim);
    check(sim.sky_visibility() == 1.0, "and comes back in the inertial frame", sim.sky_visibility());

    Simulation earth_moon;
    check(earth_moon.load_scene(ASTRAXIS_ASSET_DIR "/scenes/earth_moon.toml", &error), "earth_moon.toml loads");
    earth_moon.set_frame(find_frame(earth_moon.scene(), "Earth-Moon rotating"));
    earth_moon.clock().warp = 2.0 * kSecondsPerDay;
    settle(earth_moon);
    check(earth_moon.sky_visibility() == 1.0, "the slowly turning Earth-Moon frame keeps its sky",
          earth_moon.sky_visibility());
}

// An event can shorten a trail while its frame is shown: in Parker's
// Venus-centered frame, the earlier flybys would otherwise cross Venus.
void test_event_trails()
{
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/parker.toml", &error), "parker.toml loads");
    const Scene& scene = sim.scene();
    const int parker = scene.find("Parker Solar Probe");
    auto find_event = [&](std::string_view prefix) {
        for (size_t i = 0; i < scene.events.size(); ++i) {
            if (scene.events[i].name.starts_with(prefix)) {
                return i;
            }
        }
        return scene.events.size();
    };
    const size_t flyby2 = find_event("Venus flyby 2");
    const size_t perihelion1 = find_event("Perihelion 1");
    const int venus_frame = find_frame(scene, "Venus-centered");
    const int rotating = find_frame(scene, "Sun-Venus rotating");
    check(parker > 0 && flyby2 < scene.events.size() && perihelion1 < scene.events.size() && venus_frame > 0 &&
              rotating > 0,
          "parker events and frames");
    check(scene.active_event() == -1 && scene.trail_history_days(parker) == 0.0, "no event trails after loading");

    // Passes within 1e6 km of Venus along the trail (Venus at the frame's origin).
    std::vector<glm::dvec3> points;
    std::vector<float> fades;
    auto venus_passes = [&] {
        scene.trail(parker, sim.clock().t_tdb, 2048, points, fades);
        int passes = 0;
        bool was_near = false;
        for (const glm::dvec3& p : points) {
            const bool now_near = glm::length(p) < 1.0e6;
            passes += now_near && !was_near ? 1 : 0;
            was_near = now_near;
        }
        return passes;
    };

    sim.jump_to_event(flyby2);
    check(scene.active_event() == static_cast<int>(flyby2) && scene.trail_history_days(parker) == 30.0 &&
              scene.trail_mode(parker) == TrailMode::History,
          "the flyby event sets Parker's trail", scene.trail_history_days(parker));
    const int venus = scene.find("Venus");
    check(scene.trail_mode(venus) == scene.bodies[static_cast<size_t>(venus)].trail, "other bodies keep their trails");
    const int with_event = venus_passes();
    check(with_event == 1, "flyby 2 trail passes Venus once", with_event);

    // Another frame ends the event; back in the Venus frame, flyby 1 shows up again.
    sim.set_frame(rotating);
    check(scene.active_event() == -1, "a frame change ends the event's trails");
    sim.set_frame(venus_frame);
    const int without_event = venus_passes();
    check(without_event == 2, "full trail passes Venus at both flybys", without_event);

    sim.jump_to_event(perihelion1);
    check(scene.active_event() == static_cast<int>(perihelion1) && scene.trail_history_days(parker) == 0.0,
          "an event without trails keeps the body's own");

    // An event's trail is not cut by its frame's cap (20 days in Rosetta's Sun-fixed
    // frame): the whole stay at 67P reaches back to the approach, millions of km out.
    // Picking the frame again ends the event's trails.
    Simulation rosetta_sim;
    check(rosetta_sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/rosetta.toml", &error), "rosetta.toml loads");
    const Scene& rosetta_scene = rosetta_sim.scene();
    const int rosetta = rosetta_scene.find("Rosetta");
    size_t whole_stay = rosetta_scene.events.size();
    for (size_t i = 0; i < rosetta_scene.events.size(); ++i) {
        if (rosetta_scene.events[i].name.starts_with("The whole stay")) {
            whole_stay = i;
        }
    }
    check(rosetta > 0 && whole_stay < rosetta_scene.events.size(), "rosetta's whole-stay event");
    auto farthest = [&] {
        rosetta_scene.trail(rosetta, rosetta_sim.clock().t_tdb, 2048, points, fades);
        double far = 0.0;
        for (const glm::dvec3& p : points) {
            far = std::max(far, glm::length(p)); // 67P at the frame's origin
        }
        return far;
    };
    rosetta_sim.jump_to_event(whole_stay);
    const double with_trail = farthest();
    check(rosetta_scene.frame().trail_history_days == 20.0 && with_trail > 1.0e6,
          "the event's 1000-day trail is not capped by the frame (km)", with_trail);
    rosetta_sim.set_frame(rosetta_scene.active_frame());
    const double capped = farthest();
    check(rosetta_scene.active_event() == -1 && capped < 2000.0, "the frame picked again caps the trail (km)",
          capped);
}

// from_body = "Sun" in a scene centered on a planet (no Sun body): the event's view is
// from the sun's direction, turned by phase_deg about the up axis.
void test_event_from_sun()
{
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/mars.toml", &error), "mars.toml loads");
    const Scene& scene = sim.scene();
    size_t index = scene.events.size();
    for (size_t i = 0; i < scene.events.size() && index == scene.events.size(); ++i) {
        index = scene.events[i].from_sun ? i : index;
    }
    check(scene.find("Sun") < 0 && index < scene.events.size(), "mars.toml has an event seen from the sun");
    if (index == scene.events.size()) {
        return;
    }
    const SceneEvent& e = scene.events[index];
    check(e.from_body < 0 && e.phase_deg != 0.0, "seen from the sun: no from_body, turned by a phase");
    sim.jump_to_event(index);
    sim.clock().paused = true;
    for (int k = 0; k < 40; ++k) {
        sim.update(0.1);
    }
    const glm::dvec3 focus = scene.bodies[static_cast<size_t>(e.focus)].world_position;
    double yaw = 0.0;
    double pitch = 0.0;
    double sun_yaw = 0.0;
    double sun_pitch = 0.0;
    sim.camera().angles_for_direction(sim.camera().view(1920, 1080).position - focus, &yaw, &pitch);
    sim.camera().angles_for_direction(scene.sun_position() - focus, &sun_yaw, &sun_pitch);
    const double turn = std::remainder(yaw - sun_yaw, kTwoPi) / kDegToRad;
    check(std::abs(pitch - sun_pitch) < 1e-3 && std::abs(turn - e.phase_deg) < 0.1,
          "the camera is turned from the sun's direction by phase_deg (deg)", turn);
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

    // Labels covered by a huge focus label come back only after their place
    // has stayed free for a moment, and never move relative to their body.
    glm::vec2 label_size(5000.0f);
    auto measure_sized = [&](const std::string&) { return label_size; };
    labels.reset();
    for (int k = 0; k < 30; ++k) {
        labels.update(sim.scene(), view, focus, screen, font, measure_sized, 1.0f / 30.0f);
    }
    std::vector<glm::vec2> offsets(sim.scene().bodies.size(), glm::vec2(0.0f));
    for (const LabelLayout::BodyMark& mark : labels.bodies()) {
        offsets[static_cast<size_t>(mark.body)] = mark.text_position - mark.position;
    }
    label_size = glm::vec2(1.0f);
    float waiting_max = 0.0f;
    for (int k = 0; k < 15; ++k) {
        labels.update(sim.scene(), view, focus, screen, font, measure_sized, 1.0f / 30.0f);
    }
    for (const LabelLayout::BodyMark& mark : labels.bodies()) {
        waiting_max = mark.body == focus ? waiting_max : std::max(waiting_max, mark.text_alpha);
    }
    check(waiting_max == 0.0f, "uncovered labels wait before they show again", waiting_max);
    for (int k = 0; k < 60; ++k) {
        labels.update(sim.scene(), view, focus, screen, font, measure_sized, 1.0f / 30.0f);
    }
    int reshown = 0;
    bool moved = false;
    for (const LabelLayout::BodyMark& mark : labels.bodies()) {
        reshown += mark.body != focus && mark.text_alpha > 0.99f * mark.fade && mark.fade > 0.05f ? 1 : 0;
        const glm::vec2 offset = mark.text_position - mark.position;
        moved = moved || glm::length(offset - offsets[static_cast<size_t>(mark.body)]) > 1e-3f;
    }
    check(reshown > 0, "uncovered labels show again after the wait", static_cast<double>(reshown));
    check(!moved, "labels keep their place relative to their body");
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

// Auto tour: many shots at varying times keep the camera finite and outside
// bodies. The time warp follows the director's shot warp, as in the app.
void test_camera_director(const char* scene_file, double base_warp, double minutes, uint32_t seed,
                          bool* transit_seen = nullptr, int min_kinds = 3)
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
            const double transit_warp = director.shot_warp() > 0.0 ? director.shot_warp() : base_warp;
            check(transit_warp <= 600.0, "transit shot slows time", transit_warp);
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
    check(distinct >= min_kinds, "director uses enough shot kinds", distinct);
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

// The flyby an event shows: Galileo passing Io on arrival day, 898 km up
// (the event's own name; Io's mean radius 1821.49 km [SPHY] as in jupiter.toml).
void test_find_flyby()
{
    Scene scene = load_scene_or_die("jupiter.toml");
    const int io = scene.find("Io");
    const int galileo = scene.find("Galileo");
    double t0 = 0.0;
    parse_utc("1995-12-07 16:47", &t0);
    const Flyby flyby = find_flyby(scene, io, galileo, t0, t0 + 7200.0);
    check(flyby.body == galileo, "Galileo flies past Io", flyby.body);
    const double altitude = flyby.distance_km - scene.bodies[static_cast<size_t>(io)].equatorial_radius_km;
    check(std::abs(altitude - 898.0) < 20.0, "Galileo-Io closest approach altitude (km)", altitude);
    check(flyby.t_tdb > t0 + 1800.0 && flyby.t_tdb < t0 + 5400.0, "about an hour after the event time (s)",
          flyby.t_tdb - t0);
    check(std::abs(glm::length(flyby.normal_icrf) - 1.0) < 1e-9, "flyby plane normal is a unit vector");
    // The other way round: the same pass.
    const Flyby reverse = find_flyby(scene, galileo, io, t0, t0 + 7200.0);
    check(reverse.body == io && std::abs(reverse.t_tdb - flyby.t_tdb) < 1.0, "from Galileo: Io, at the same time");
    // Before Galileo arrived: no flyby.
    double t1990 = 0.0;
    parse_utc("1990-01-01 00:00", &t1990);
    check(find_flyby(scene, io, galileo, t1990, t1990 + 7200.0).body < 0, "no Galileo in 1990");
    // Still closing in at the end of the interval: not a closest approach.
    check(find_flyby(scene, io, galileo, t0, t0 + 1800.0).body < 0, "approach cut off by the interval");
}

// Every event that names a flyby shows it: the body passes the focus, closer
// than the camera, within the search after the event's time. A jump to the
// event looks across the pass, as the tour does.
void test_event_flybys()
{
    int named = 0;
    for (const auto& [stem, scene_ptr] : all_scenes()) {
        const Scene& scene = *scene_ptr;
        if (std::none_of(scene.events.begin(), scene.events.end(), [](const SceneEvent& e) { return e.flyby >= 0; })) {
            continue;
        }
        Simulation sim;
        std::string error;
        check(sim.load_scene(std::string(ASTRAXIS_ASSET_DIR "/scenes/") + stem + ".toml", &error),
              (stem + " loads").c_str());
        for (size_t i = 0; i < sim.scene().events.size(); ++i) {
            const SceneEvent& e = sim.scene().events[i];
            if (e.flyby < 0) {
                continue;
            }
            ++named;
            const std::string what = stem + " / " + e.name;
            const double warp = e.warp > 0.0 ? e.warp : sim.scene().view.warp;
            const Flyby flyby = sim.event_flyby(i, warp);
            check(flyby.body == e.flyby, (what + ": the named flyby is found").c_str());
            if (flyby.body < 0) {
                continue;
            }
            const double distance = e.distance_km > 0.0 ? e.distance_km : 0.0;
            check(distance <= 0.0 || flyby.distance_km < distance, (what + ": passes within the shot (km)").c_str(),
                  flyby.distance_km);

            if (e.from_body >= 0 || e.from_sun || e.from_orbit_normal >= 0) {
                continue;
            }
            // The jump's view: the camera on the side of the pass's normal
            // (turned at most ~31 deg toward the sun).
            sim.jump_to_event(i);
            sim.clock().paused = true;
            for (int k = 0; k < 40; ++k) {
                sim.update(0.1);
            }
            const Scene& s = sim.scene();
            const glm::dvec3 normal = s.current_transform().direction_to_display(flyby.normal_icrf);
            const glm::dvec3 from_focus = glm::normalize(sim.camera().view(1920, 1080).position -
                                                         s.bodies[static_cast<size_t>(e.focus)].world_position);
            check(std::abs(glm::dot(from_focus, normal)) > 0.85, (what + ": seen across the pass").c_str(),
                  glm::dot(from_focus, normal));
            sim.clock().paused = false;
        }
    }
    check(named >= 80, "events name their flybys", named);
}

// The tour goes through the events: it fades out, moves the clock to before
// the next event, slows down to the event's warp as the event comes (a flyby:
// 20 s before closest approach), shows its caption, holds a flyby until after
// closest approach and then goes on with other shots at the user's warp.
void test_tour_events()
{
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error), "jupiter loads for the event tour");
    const Scene& scene = sim.scene();
    SimClock& clock = sim.clock();
    parse_utc("1995-06-01 00:00", &clock.t_tdb); // before all events
    clock.warp = 3600.0;
    clock.paused = false;
    clock.reverse = false;
    sim.update(0.0);
    sim.start_tour(11u);

    // The first event in time is the Galileo tour overview, then its Io flyby.
    const int io_event = 0;
    check(scene.events[static_cast<size_t>(io_event)].warp == 60.0, "Io flyby event is at 1 min/s");
    const Flyby io_flyby = find_flyby(scene, scene.find("Io"), scene.find("Galileo"), scene.events[0].t_tdb,
                                      scene.events[0].t_tdb + 6000.0);

    std::vector<int> shown;
    bool faded_out = false;
    bool slowed_down = false;
    bool arrived_ok = true;
    bool caption_ok = true;
    bool flyby_held = true;
    double io_arrival_lead = 0.0;
    bool warp_restored = false;
    int last = -1;
    bool was_arrived = false;
    double approach_max_warp = 0.0;
    const double dt = 0.05;
    for (int k = 0; k < static_cast<int>(12.0 * 60.0 / dt); ++k) {
        sim.update(dt);
        faded_out = faded_out || sim.fade() == 0.0;
        const int event = sim.tour_event();
        if (event >= 0 && event != last) {
            shown.push_back(event);
            approach_max_warp = 0.0;
        }
        if (event >= 0 && !sim.tour_event_arrived() && sim.fade() == 1.0) {
            approach_max_warp = std::max(approach_max_warp, clock.warp);
        }
        const bool arrived = event >= 0 && sim.tour_event_arrived();
        if (arrived && !was_arrived) {
            const SceneEvent& e = scene.events[static_cast<size_t>(event)];
            const double w = e.warp > 0.0 ? e.warp : 3600.0;
            arrived_ok = arrived_ok && clock.t_tdb >= e.t_tdb && clock.t_tdb < e.t_tdb + 100.0 * w;
            if (event == io_event) {
                io_arrival_lead = (io_flyby.t_tdb - clock.t_tdb) / w;
            }
            caption_ok = caption_ok && sim.last_event() == event && sim.event_age() == 0.0;
            if (approach_max_warp > 4.0 * w) {
                slowed_down = true;
            }
            check(clock.warp == w, "an event is held at its warp", clock.warp);
        }
        if (last == io_event && event != io_event && io_flyby.body >= 0) {
            flyby_held = flyby_held && clock.t_tdb > io_flyby.t_tdb;
        }
        if (event < 0 && last >= 0 && sim.director().shot_warp() == 0.0) {
            warp_restored = clock.warp == 3600.0;
        }
        last = event;
        was_arrived = arrived;
    }
    check(shown.size() >= 3, "the tour shows several events in 12 minutes", static_cast<double>(shown.size()));
    check(!shown.empty() && shown.front() == 1, "first the earliest event (Galileo's tour overview)",
          shown.empty() ? -1.0 : shown.front());
    check(shown.size() >= 2 && shown[1] == io_event, "then the next one (the Io flyby, 5 hours later)",
          shown.size() >= 2 ? shown[1] : -1.0);
    bool chronological = true;
    for (size_t i = 2; i < shown.size(); ++i) {
        chronological = chronological && scene.events[static_cast<size_t>(shown[i])].t_tdb >
                                             scene.events[static_cast<size_t>(shown[i - 1])].t_tdb;
    }
    check(chronological, "events in time order");
    check(faded_out, "moving the clock fades to black");
    check(arrived_ok, "the clock reaches each event's time");
    check(slowed_down, "the warp eases down to a flyby's own");
    check(caption_ok, "arriving at an event shows its caption");
    check(io_flyby.body >= 0 && flyby_held, "the Io flyby shot lasts until after closest approach");
    check(std::abs(io_arrival_lead - 20.0) < 0.2, "the Io approach ends 20 s before closest approach",
          io_arrival_lead);
    check(warp_restored, "after an event the user's warp returns");
    sim.stop_tour();
    check(sim.fade() == 1.0 && sim.tour_event() < 0 && clock.warp == 3600.0, "stopping the tour ends the event");

    // An event coming up anyway is reached without a jump: no fade, and the
    // clock only moves forward.
    parse_utc("1996-09-06 16:00", &clock.t_tdb); // 2 h before "Galileo at Ganymede"
    clock.warp = 600.0;
    sim.update(0.0);
    sim.start_tour(12u);
    const int ganymede_event = sim.tour_event();
    check(ganymede_event == 2, "an imminent event is taken at once", ganymede_event);
    bool no_jump = true;
    double t_prev = clock.t_tdb;
    for (int k = 0; k < static_cast<int>(60.0 / dt) && !sim.tour_event_arrived(); ++k) {
        sim.update(dt);
        no_jump = no_jump && sim.fade() == 1.0 && clock.t_tdb >= t_prev && clock.t_tdb - t_prev <= 600.0 * dt * 1.01;
        t_prev = clock.t_tdb;
    }
    check(no_jump, "approached without a jump, never faster than the user's warp");
    check(sim.tour_event_arrived() && sim.last_event() == ganymede_event, "and arrived at it");
    sim.stop_tour();
}

} // namespace

void run_scene_tests()
{
    test_scene_loader();
    test_satellite_fades();
    test_cross_scene_orbits();
    test_ephemeris_blends();
    test_simulation();
    test_event_trails();
    test_sky_fades_in_fast_frames();
    test_event_from_sun();
    test_body_fixed_frames();
    test_label_layout();
    test_orbit_center();
    test_body_masses();
    test_event_captions();
    test_find_flyby();
    test_event_flybys();
    test_tour_events();
    for (uint32_t seed : {42u, 7u, 2026u, 99u, 12345u}) {
        test_camera_director("jupiter.toml", 3600.0, 20.0, seed);
    }
    for (uint32_t seed : {1u, 2u, 3u}) {
        test_camera_director("solar_system.toml", 2629800.0, 10.0, seed);
    }
    for (uint32_t seed : {5u, 6u}) {
        test_camera_director("mars.toml", 600.0, 10.0, seed); // Phobos: 11 km, close to Mars
    }
    // Venus in 2026: no spacecraft and no moons, so only overviews and planet close-ups
    // (Zoozve, 0.1 km across, is one of the planets here).
    for (uint32_t seed : {7u, 8u}) {
        test_camera_director("venus.toml", 172800.0, 10.0, seed, nullptr, 2);
    }
    test_camera_director_black_hole();
}
