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
// orbits seen face-on.
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

    // Index of a moon whose shadow currently falls on the root body, or -1.
    static int find_shadow_transit(const Scene& scene);

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
};

} // namespace astraxis
