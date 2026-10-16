#pragma once

#include "scene/camera.hpp"
#include "scene/scene.hpp"

#include <cstdint>
#include <random>

namespace astraxis {

// Auto tour ("idle mode"): picks a sequence of slow cinematic shots and flies
// the camera between them. Shots adapt to the scene: system overviews at the
// scale of a random orbit, lit-crescent planet close-ups, moons against their
// planet, spacecraft with their trail, moon shadow transits when one is
// happening, black holes seen from just above their disk, and precessing
// orbits seen face-on. Shots of fast motion (a spinning planet, a moon's
// orbit, a shadow transit) slow the time warp down while they last.
//
// Every few shots the tour shows one of the scene's events instead. Those
// move the clock, so Simulation runs them: the director only says when one is
// due (take_event_request) and holds it for as long as Simulation says.
class CameraDirector {
public:
    enum class ShotKind {
        Overview,
        PlanetCloseUp,
        MoonWithPlanet,
        ShadowTransit,
        SpacecraftFollow,
        BlackHoleCloseUp,
        OrbitTopDown,
        Event, // run by Simulation
        Count,
    };

    void start(const Scene& scene, OrbitCamera& camera, uint32_t seed);
    void stop(OrbitCamera& camera);
    bool active() const { return m_active; }

    void update(double real_dt, const Scene& scene, OrbitCamera& camera);

    // Simulated seconds per real second (0 when paused); used to skip shots
    // that would be over in a blink at high time warp.
    void set_time_warp(double warp) { m_time_warp = warp; }

    ShotKind current_shot() const { return m_kind; }
    static const char* shot_name(ShotKind kind);
    double time_left() const { return m_time_left; }
    // Time warp the current shot needs (e.g. to follow a fast orbit), or 0.
    double shot_warp() const { return m_active ? m_shot_warp : 0.0; }

    // Index of a moon whose shadow currently falls on its planet, or -1.
    static int find_shadow_transit(const Scene& scene);

    // Whether the scene has an event to show now, and whether the next one is
    // coming up soon enough to be shown without moving the clock; an imminent
    // event is taken at the next shot, others after one or two other shots.
    void set_events(bool available, bool imminent)
    {
        m_events_available = available;
        m_event_imminent = imminent;
    }
    // True (once) when the shot just begun is an event: the camera is left
    // to the caller, which either calls hold_event() or skip_event().
    bool take_event_request();
    // The event shot ends `seconds` from now and the tour goes on.
    void hold_event(double seconds) { m_time_left = seconds; }
    // No event can be shown after all: the next shot starts now.
    void skip_event(const Scene& scene, OrbitCamera& camera);

private:
    void next_shot(const Scene& scene, OrbitCamera& camera);
    double uniform(double lo, double hi);
    double random_sign();
    template <typename T>
    const T& pick(const std::vector<T>& v)
    {
        return v[static_cast<size_t>(uniform(0.0, static_cast<double>(v.size()))) % v.size()];
    }

    bool m_active = false;
    std::mt19937 m_rng{1u};
    ShotKind m_kind = ShotKind::Overview;
    int m_last_target = -1;
    double m_time_left = 0.0;
    double m_time_warp = 1.0;
    double m_shot_warp = 0.0;
    bool m_events_available = false;
    bool m_event_imminent = false;
    bool m_event_requested = false;
    int m_shots_until_event = 1;
};

// The closest approach of `other` to `body` during [t0, t1] (the closest of
// several). Only true minima count, not either end of the interval (where the
// distance is still falling or already rising), nor times when either body
// does not exist.
struct Flyby {
    int body = -1; // the body passing by; -1 = none
    double t_tdb = 0.0;
    double distance_km = 0.0;
    glm::dvec3 normal_icrf{0.0}; // unit normal of the relative path (r x v)
};
Flyby find_flyby(const Scene& scene, int body, int other, double t0, double t1);

} // namespace astraxis
