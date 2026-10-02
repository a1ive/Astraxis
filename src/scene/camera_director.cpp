#include "scene/camera_director.hpp"

#include "core/math.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace astraxis {

namespace {

constexpr double kTransitionSeconds = 7.0;
constexpr double kMaxTransitWarp = 600.0; // 10 min/s: a 2 h transit lasts ~12 s
// Face-on orbit shots: at most one orbit per kMinOrbitSeconds of real time,
// otherwise the warp is lowered to one orbit per kOrbitSeconds.
constexpr double kMinOrbitSeconds = 10.0;
constexpr double kOrbitSeconds = 25.0;

bool is_planet(const Body& b)
{
    return b.kind == BodyKind::Planet;
}

} // namespace

const char* CameraDirector::shot_name(ShotKind kind)
{
    switch (kind) {
    case ShotKind::Overview:
        return "Overview";
    case ShotKind::PlanetCloseUp:
        return "Planet close-up";
    case ShotKind::MoonWithPlanet:
        return "Moon";
    case ShotKind::ShadowTransit:
        return "Shadow transit";
    case ShotKind::SpacecraftFollow:
        return "Spacecraft";
    case ShotKind::BlackHoleCloseUp:
        return "Black hole";
    case ShotKind::OrbitTopDown:
        return "Orbit (face-on)";
    case ShotKind::Count:
        break;
    }
    return "";
}

double CameraDirector::uniform(double lo, double hi)
{
    return std::uniform_real_distribution<double>(lo, hi)(m_rng);
}

double CameraDirector::random_sign()
{
    return uniform(0.0, 1.0) < 0.5 ? -1.0 : 1.0;
}

void CameraDirector::start(const Scene& scene, OrbitCamera& camera, uint32_t seed)
{
    m_active = true;
    m_rng.seed(seed);
    m_last_target = -1;
    next_shot(scene, camera);
}

void CameraDirector::stop(OrbitCamera& camera)
{
    m_active = false;
    camera.set_drift(0.0);
}

int CameraDirector::find_shadow_transit(const Scene& scene)
{
    if (scene.bodies.empty() || !is_planet(scene.bodies[0])) {
        return -1;
    }
    const Body& root = scene.bodies[0];
    const glm::dvec3 sun = glm::normalize(scene.sun_position() - root.world_position);
    for (size_t i = 1; i < scene.bodies.size(); ++i) {
        const Body& b = scene.bodies[i];
        if (b.parent != 0 || !is_planet(b) || !b.visible) {
            continue;
        }
        const glm::dvec3 rel = b.world_position - root.world_position;
        const double along = glm::dot(rel, sun);
        const double perp = glm::length(rel - along * sun);
        // Keep a margin so the shadow is well inside the disk, not at the limb.
        if (along > 0.0 && perp < root.equatorial_radius_km * 0.8) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

void CameraDirector::update(double real_dt, const Scene& scene, OrbitCamera& camera)
{
    if (!m_active) {
        return;
    }
    m_time_left -= real_dt;
    if (m_time_left <= 0.0) {
        next_shot(scene, camera);
    }
}

void CameraDirector::next_shot(const Scene& scene, OrbitCamera& camera)
{
    if (scene.bodies.empty()) {
        return;
    }
    const Body& root = scene.bodies[0];

    // What this scene offers right now.
    std::vector<int> planets;     // visible planets/moons (close-up targets)
    std::vector<int> moons;       // planets orbiting a planet
    std::vector<int> spacecraft;  // visible spacecraft
    std::vector<int> black_holes; // visible black holes
    std::vector<int> precessing;  // bodies whose periapsides are marked
    std::vector<double> orbit_radii;
    for (size_t i = 0; i < scene.bodies.size(); ++i) {
        const Body& b = scene.bodies[i];
        if (!b.visible) {
            continue;
        }
        const int idx = static_cast<int>(i);
        if (is_planet(b) || b.kind == BodyKind::Star) {
            planets.push_back(idx); // close-up targets: planets, moons and stars
        }
        if (is_planet(b)) {
            if (b.parent >= 0 && is_planet(scene.bodies[static_cast<size_t>(b.parent)])) {
                moons.push_back(idx);
            }
        } else if (b.kind == BodyKind::Spacecraft) {
            spacecraft.push_back(idx);
        } else if (b.kind == BodyKind::BlackHole) {
            black_holes.push_back(idx);
        }
        if (b.mark_periapsides && b.motion) {
            precessing.push_back(idx);
        }
        if (b.parent == 0 && b.kind != BodyKind::Barycenter) {
            orbit_radii.push_back(glm::length(b.world_position - root.world_position));
        }
    }

    // Weighted choice among the available kinds; avoid repeating the last one
    // (except moons/spacecraft, which vary by target). Transits are rare, so
    // they are taken when available and watchable at the current warp.
    double weights[static_cast<int>(ShotKind::Count)] = {};
    weights[static_cast<int>(ShotKind::Overview)] = orbit_radii.empty() ? 0.0 : 0.25;
    weights[static_cast<int>(ShotKind::PlanetCloseUp)] = planets.empty() ? 0.0 : 0.30;
    weights[static_cast<int>(ShotKind::MoonWithPlanet)] = moons.empty() ? 0.0 : 0.45;
    weights[static_cast<int>(ShotKind::SpacecraftFollow)] = spacecraft.empty() ? 0.0 : 0.40;
    weights[static_cast<int>(ShotKind::BlackHoleCloseUp)] = black_holes.empty() ? 0.0 : 0.40;
    weights[static_cast<int>(ShotKind::OrbitTopDown)] = precessing.empty() ? 0.0 : 0.40;
    if (m_kind == ShotKind::Overview || m_kind == ShotKind::PlanetCloseUp || m_kind == ShotKind::BlackHoleCloseUp ||
        m_kind == ShotKind::OrbitTopDown) {
        weights[static_cast<int>(m_kind)] = 0.0;
    }

    ShotKind kind = ShotKind::Overview;
    const int transit = m_time_warp <= kMaxTransitWarp ? find_shadow_transit(scene) : -1;
    if (transit >= 0 && m_kind != ShotKind::ShadowTransit && uniform(0.0, 1.0) < 0.7) {
        kind = ShotKind::ShadowTransit;
    } else {
        double total = 0.0;
        for (double w : weights) {
            total += w;
        }
        if (total <= 0.0) {
            weights[static_cast<int>(ShotKind::Overview)] = 1.0;
            total = 1.0;
        }
        double r = uniform(0.0, total);
        for (int k = 0; k < static_cast<int>(ShotKind::Count); ++k) {
            if (weights[k] <= 0.0) {
                continue;
            }
            kind = static_cast<ShotKind>(k); // last positive kind absorbs rounding
            if (r < weights[k]) {
                break;
            }
            r -= weights[k];
        }
    }
    m_kind = kind;

    // Angles of the direction from `body` toward the sun (display frame).
    auto sun_angles = [&](int body, double* yaw, double* pitch) {
        camera.angles_for_direction(scene.sun_position() - scene.bodies[static_cast<size_t>(body)].world_position,
                                    yaw, pitch);
    };

    CameraShot shot;
    double duration = 40.0;
    double drift = 0.0;
    m_shot_warp = 0.0;

    switch (kind) {
    case ShotKind::Overview: {
        // Frame the orbit of a random body (inner or outer system).
        const double radius = orbit_radii.empty() ? scene.system_extent_km() : pick(orbit_radii);
        shot.target = 0;
        shot.yaw = uniform(-kPi, kPi);
        shot.pitch = uniform(6.0, 32.0) * kDegToRad;
        shot.distance = std::max(radius, root.equatorial_radius_km * 4.0) * uniform(1.8, 3.0);
        duration = uniform(40.0, 60.0);
        drift = random_sign() * uniform(0.006, 0.014);
        break;
    }

    case ShotKind::PlanetCloseUp: {
        // Partly lit: 30-120 deg away from the sun direction gives gibbous to crescent.
        std::vector<int> candidates;
        for (int p : planets) {
            if (p != m_last_target) {
                candidates.push_back(p);
            }
        }
        const int target = is_planet(root) && m_last_target != 0 ? 0 : pick(candidates.empty() ? planets : candidates);
        m_last_target = target;
        double sun_yaw = 0.0;
        double sun_pitch = 0.0;
        sun_angles(target, &sun_yaw, &sun_pitch);
        shot.target = target;
        shot.yaw = sun_yaw + random_sign() * uniform(30.0, 120.0) * kDegToRad;
        shot.pitch = uniform(-12.0, 25.0) * kDegToRad;
        shot.distance = scene.bodies[static_cast<size_t>(target)].equatorial_radius_km * uniform(3.0, 6.0);
        duration = uniform(35.0, 50.0);
        drift = random_sign() * uniform(0.006, 0.012);
        break;
    }

    case ShotKind::MoonWithPlanet: {
        // Prefer moons on the day side of their planet: the camera sits beyond the
        // moon looking back at the planet, so a night-side moon means staring
        // into the sun with both bodies backlit.
        std::vector<int> candidates;
        for (int i : moons) {
            const Body& m = scene.bodies[static_cast<size_t>(i)];
            const glm::dvec3 planet = scene.bodies[static_cast<size_t>(m.parent)].world_position;
            const glm::dvec3 radial = glm::normalize(m.world_position - planet);
            const glm::dvec3 sun = glm::normalize(scene.sun_position() - planet);
            if (i != m_last_target && glm::dot(radial, sun) > -0.2) {
                candidates.push_back(i);
            }
        }
        const int moon = pick(candidates.empty() ? moons : candidates);
        m_last_target = moon;
        const Body& m = scene.bodies[static_cast<size_t>(moon)];

        // Offset the camera from the radial direction toward the sun side.
        double radial_yaw = 0.0;
        double radial_pitch = 0.0;
        camera.angles_for_direction(m.world_position - scene.bodies[static_cast<size_t>(m.parent)].world_position,
                                    &radial_yaw, &radial_pitch);
        double sun_yaw = 0.0;
        double sun_pitch = 0.0;
        sun_angles(moon, &sun_yaw, &sun_pitch);
        const double toward_sun = wrap_pi(sun_yaw - radial_yaw) >= 0.0 ? 1.0 : -1.0;

        shot.target = moon;
        shot.relative_to_parent = true;
        shot.yaw = toward_sun * uniform(15.0, 40.0) * kDegToRad;
        shot.pitch = uniform(3.0, 15.0) * kDegToRad;
        shot.distance = m.equatorial_radius_km * uniform(5.0, 10.0);
        duration = uniform(30.0, 40.0);
        drift = random_sign() * uniform(0.003, 0.008);
        break;
    }

    case ShotKind::ShadowTransit: {
        // From near the sun direction, slightly off so the shadow is not hidden by the moon.
        double sun_yaw = 0.0;
        double sun_pitch = 0.0;
        sun_angles(0, &sun_yaw, &sun_pitch);
        shot.target = 0;
        shot.yaw = sun_yaw + random_sign() * uniform(12.0, 25.0) * kDegToRad;
        shot.pitch = sun_pitch + uniform(-5.0, 10.0) * kDegToRad;
        shot.distance = root.equatorial_radius_km * uniform(3.5, 5.5);
        duration = uniform(30.0, 40.0);
        drift = random_sign() * 0.004;
        break;
    }

    case ShotKind::SpacecraftFollow: {
        // Beyond the spacecraft, looking back toward where it came from (its parent),
        // far enough to see the trail.
        const int craft = pick(spacecraft);
        m_last_target = craft;
        const Body& c = scene.bodies[static_cast<size_t>(craft)];
        const double r = glm::length(c.world_position - scene.bodies[static_cast<size_t>(c.parent)].world_position);
        shot.target = craft;
        shot.relative_to_parent = true;
        shot.yaw = random_sign() * uniform(20.0, 60.0) * kDegToRad;
        shot.pitch = uniform(10.0, 30.0) * kDegToRad;
        shot.distance = std::max(r * uniform(0.25, 0.6), 1000.0);
        duration = uniform(35.0, 50.0);
        drift = random_sign() * uniform(0.004, 0.01);
        break;
    }

    case ShotKind::BlackHoleCloseUp: {
        // A few tens of M away, just above the disk plane (the scene's up axis is
        // the hole's pole): the lensed far side of the disk arches over the shadow.
        const int hole = pick(black_holes);
        m_last_target = hole;
        const double m_km = scene.bodies[static_cast<size_t>(hole)].gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
        shot.target = hole;
        shot.yaw = uniform(-kPi, kPi);
        shot.pitch = uniform(4.0, 18.0) * kDegToRad;
        shot.distance = m_km * uniform(30.0, 60.0);
        duration = uniform(40.0, 60.0);
        drift = random_sign() * uniform(0.004, 0.008);
        break;
    }

    case ShotKind::OrbitTopDown: {
        // Along the orbit normal, framing the apoapsis, so the trail and the
        // periapsis spokes show the rosette undistorted. No drift: yaw drift
        // turns about the scene's up axis, not about the orbit normal.
        const int body = pick(precessing);
        m_last_target = body;
        const Body& b = scene.bodies[static_cast<size_t>(body)];
        double r_peri = 0.0;
        double r_apo = 0.0;
        double period = 0.0;
        if (!scene.osculating_apsides(body, &r_peri, &r_apo, &period)) {
            r_apo = 2.0 * glm::length(b.world_position - scene.bodies[static_cast<size_t>(b.parent)].world_position);
        } else if (m_time_warp * kMinOrbitSeconds > period) {
            m_shot_warp = period / kOrbitSeconds;
        }
        double yaw = 0.0;
        double pitch = 0.0;
        camera.angles_for_direction(scene.orbit_normal(body), &yaw, &pitch);
        shot.target = b.parent;
        shot.yaw = yaw;
        shot.pitch = pitch;
        // The parent sits at a focus, so the rosette reaches r_apo on every side;
        // fit that within the vertical half field of view (22.5 deg).
        shot.distance = r_apo * uniform(2.8, 3.3);
        duration = uniform(40.0, 60.0);
        break;
    }

    case ShotKind::Count:
        break;
    }

    camera.set_drift(drift);
    camera.fly_to(shot, kTransitionSeconds, scene);
    m_time_left = duration + kTransitionSeconds;
}

} // namespace astraxis
