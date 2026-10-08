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

void run_scene_tests()
{
    test_scene_loader();
    test_satellite_fades();
    test_cross_scene_orbits();
    test_ephemeris_blends();
    test_simulation();
    test_label_layout();
    test_orbit_center();
    test_body_masses();
    test_event_captions();
    for (uint32_t seed : {42u, 7u, 2026u, 99u, 12345u}) {
        test_camera_director("jupiter.toml", 3600.0, 20.0, seed);
    }
    for (uint32_t seed : {1u, 2u, 3u}) {
        test_camera_director("solar_system.toml", 2629800.0, 10.0, seed);
    }
    test_camera_director_black_hole();
}
