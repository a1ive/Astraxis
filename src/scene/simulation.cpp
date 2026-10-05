#include "scene/simulation.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "scene/scene_loader.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <utility>

namespace astraxis {

namespace {

// A moon's label, marker and orbit fade out as its orbit shrinks on screen
// from kSatelliteShowPx to kSatelliteHidePx (points, not pixels).
constexpr double kSatelliteHidePx = 10.0;
constexpr double kSatelliteShowPx = 28.0;
constexpr double kDiskMaxWarp = 60.0; // cap on the disk animation speed (ISCO period ~92 M)

} // namespace

bool Simulation::load_scene(const std::filesystem::path& file, std::string* error)
{
    Scene scene;
    if (!load_scene_file(file, scene, error)) {
        return false;
    }
    m_scene = std::move(scene);

    const bool touring = m_director.active();
    stop_tour();
    if (m_scene.view.start_now) {
        reset_to_now();
    } else {
        m_clock.t_tdb = m_scene.view.start_tdb;
    }
    m_clock.warp = m_scene.view.warp;
    m_scene.update(m_clock.t_tdb);

    m_camera = OrbitCamera{};
    m_camera.set_up_axis(m_scene.up_axis());
    const int focus = m_scene.view.focus;
    m_camera.focus(focus, m_scene);
    if (m_scene.view.distance_km > 0.0) {
        m_camera.fly_to({focus, m_scene.view.yaw_rad, m_scene.view.pitch_rad, m_scene.view.distance_km, false}, 0.0, m_scene);
    }
    if (touring) {
        start_tour();
    }
    return true;
}

void Simulation::update(double real_dt)
{
    m_clock.advance(real_dt);
    m_scene.update(m_clock.t_tdb);

    // Accretion-disk animation: follows the clock, but capped so that the inner
    // disk does not blur at high time warp.
    for (const Body& body : m_scene.bodies) {
        if (body.kind == BodyKind::BlackHole && body.gm_km3_s2 > 0.0) {
            const double t_m = body.gm_km3_s2 / std::pow(kSpeedOfLightKmS, 3.0); // seconds per M
            const double rate = m_clock.paused ? 0.0 : std::min(m_clock.warp, kDiskMaxWarp);
            m_disk_time = std::fmod(m_disk_time + real_dt * rate * (m_clock.reverse ? -1.0 : 1.0) / t_m, 1.0e6);
            break;
        }
    }

    if (!m_clock.paused) {
        m_animation_time = std::fmod(m_animation_time + real_dt, 1.0e6);
    }

    m_director.set_time_warp(m_clock.paused ? 0.0 : m_clock.warp);
    m_director.update(real_dt, m_scene, m_camera);
    apply_tour_warp();
    m_camera.update(real_dt, m_scene);
}

void Simulation::compute_view(int width, int height, double content_scale, OutputView& out) const
{
    CameraView& view = out.camera;
    view = m_camera.view(width, height);
    // Stars are ICRF directions: rotate them into the display frame.
    view.sky_view_proj = view.view_proj * glm::mat4(glm::mat3(glm::transpose(m_scene.current_transform().axes)));
    out.content_scale = content_scale;

    const double px_per_radian = 0.5 * view.viewport.y * view.focal_y;
    satellite_fades(m_scene, view.position, px_per_radian, kSatelliteHidePx * content_scale,
                    kSatelliteShowPx * content_scale, out.body_fades);
    const int target = m_camera.target();
    if (target >= 0 && target < static_cast<int>(out.body_fades.size())) {
        out.body_fades[static_cast<size_t>(target)] = 1.0f; // the focus always shows
    }
}

void Simulation::reset_to_now()
{
    const auto since_epoch = std::chrono::system_clock::now().time_since_epoch();
    m_clock.t_tdb = tdb_from_unix_utc(std::chrono::duration<double>(since_epoch).count());
}

void Simulation::set_focus(int body)
{
    if (body >= 0 && body < static_cast<int>(m_scene.bodies.size())) {
        stop_tour();
        if (body != m_camera.target()) {
            m_camera.focus(body, m_scene);
        }
    }
}

void Simulation::set_frame(int frame)
{
    if (frame == m_scene.active_frame()) {
        return;
    }
    m_scene.set_active_frame(frame);
    m_scene.update(m_clock.t_tdb);
    m_camera.set_up_axis(m_scene.up_axis());
}

void Simulation::jump_to_event(size_t index)
{
    const SceneEvent& e = m_scene.events[index];
    stop_tour();
    m_clock.t_tdb = e.t_tdb;
    if (e.warp > 0.0) {
        m_clock.warp = e.warp;
    }
    set_frame(e.frame >= 0 ? e.frame : 0); // events without a frame use the root frame
    m_scene.update(m_clock.t_tdb);
    if (e.focus >= 0) {
        const double distance =
            e.distance_km > 0.0 ? e.distance_km : m_camera.default_distance(e.focus, m_scene);
        double yaw = -1.2;
        double pitch = 0.3;
        if (e.from_orbit_normal >= 0) {
            m_camera.angles_for_direction(m_scene.orbit_normal(e.from_orbit_normal), &yaw, &pitch);
        } else if (e.from_body >= 0) {
            const glm::dvec3 to_body = m_scene.bodies[static_cast<size_t>(e.from_body)].world_position -
                                       m_scene.bodies[static_cast<size_t>(e.focus)].world_position;
            m_camera.angles_for_direction(to_body, &yaw, &pitch);
            yaw += e.phase_deg * kDegToRad;
        }
        m_camera.fly_to({e.focus, yaw, pitch, distance, false}, 2.0, m_scene);
    }
}

void Simulation::look_from_sun_side()
{
    stop_tour();
    // Slightly off the sun line, so moon shadows on the planet are not hidden behind the moons.
    const glm::dvec3 target = m_scene.bodies[static_cast<size_t>(m_camera.target())].world_position;
    m_camera.look_from(m_scene.sun_position() - target);
    m_camera.rotate(0.35, 0.1);
}

void Simulation::start_tour()
{
    const auto ticks = std::chrono::steady_clock::now().time_since_epoch().count();
    m_director.set_time_warp(m_clock.paused ? 0.0 : m_clock.warp);
    m_director.start(m_scene, m_camera, static_cast<uint32_t>(ticks));
    apply_tour_warp();
}

void Simulation::stop_tour()
{
    if (m_director.active()) {
        m_director.stop(m_camera);
        apply_tour_warp();
    }
}

void Simulation::apply_tour_warp()
{
    const double want = m_director.shot_warp();
    if (want == m_tour_warp) {
        return;
    }
    if (want > 0.0) {
        if (m_tour_warp <= 0.0) {
            m_saved_warp = m_clock.warp;
        }
        m_clock.warp = want;
    } else {
        m_clock.warp = m_saved_warp;
    }
    m_tour_warp = want;
}

} // namespace astraxis
