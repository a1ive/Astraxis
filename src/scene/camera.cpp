#include "scene/camera.hpp"

#include "core/math.hpp"

#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr double kFocusTransitionSeconds = 1.5;
constexpr double kMaxDistanceKm = 1e14; // ~670,000 au
constexpr double kPitchLimit = 89.0 * kDegToRad;
// The camera never gets closer than this many radii to any body's center
// (transitions between distant targets could otherwise cut through a planet).
constexpr double kBodyClearance = 1.1;
// Black holes (radius = horizon): stay outside ~12 M, where the view is still
// recognisable (shadow plus disk) rather than filled by the shadow.
constexpr double kBlackHoleClearance = 6.0;

double clearance_radius(const Body& b)
{
    return b.equatorial_radius_km * (b.kind == BodyKind::BlackHole ? kBlackHoleClearance : kBodyClearance);
}

double smoothstep(double x)
{
    x = std::clamp(x, 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

} // namespace

void OrbitCamera::set_up_axis(const glm::dvec3& up)
{
    m_basis = frame_from_pole(glm::normalize(up));
}

double OrbitCamera::default_distance(int body, const Scene& scene) const
{
    const Body& b = scene.bodies[static_cast<size_t>(body)];
    if (b.parent < 0) {
        // Root: frame the whole system.
        return scene.system_extent_km() * 2.6;
    }
    return b.equatorial_radius_km * 12.0;
}

double OrbitCamera::min_distance(int body, const Scene& scene) const
{
    const Body& b = scene.bodies[static_cast<size_t>(body)];
    return b.kind == BodyKind::BlackHole ? clearance_radius(b) : b.equatorial_radius_km * 1.15;
}

void OrbitCamera::focus(int body, const Scene& scene)
{
    if (body < 0 || body >= static_cast<int>(scene.bodies.size())) {
        return;
    }
    // Keep the current absolute viewing angles.
    fly_to({body, m_yaw + m_base_yaw, m_pitch + m_base_pitch, default_distance(body, scene), false},
           m_distance > 0.0 ? kFocusTransitionSeconds : 0.0, scene);
}

bool OrbitCamera::radial_angles(int body, const Scene& scene, double* yaw, double* pitch) const
{
    const Body& b = scene.bodies[static_cast<size_t>(body)];
    if (b.parent < 0) {
        return false;
    }
    const glm::dvec3 dir = b.world_position - scene.bodies[static_cast<size_t>(b.parent)].world_position;
    angles_for_direction(dir, yaw, pitch);
    return true;
}

void OrbitCamera::fly_to(const CameraShot& shot, double duration, const Scene& scene)
{
    // Current absolute angles, re-expressed in the new shot's reference.
    const double abs_yaw = m_yaw + m_base_yaw;
    const double abs_pitch = m_pitch + m_base_pitch;

    m_relative = shot.relative_to_parent && radial_angles(shot.target, scene, &m_base_yaw, &m_base_pitch);
    if (!m_relative) {
        m_base_yaw = 0.0;
        m_base_pitch = 0.0;
    }
    m_yaw = wrap_pi(abs_yaw - m_base_yaw);
    m_pitch = abs_pitch - m_base_pitch;

    if (duration <= 0.0 || m_distance <= 0.0) {
        m_target = m_prev_target = shot.target;
        m_yaw = wrap_pi(shot.yaw);
        m_pitch = shot.pitch;
        m_distance = shot.distance;
        m_transition = 1.0;
        m_animate_view = false;
        return;
    }

    m_prev_target = m_target;
    m_target = shot.target;
    m_transition = 0.0;
    m_transition_duration = duration;
    m_animate_view = true;

    m_yaw_from = m_yaw;
    m_yaw_to = m_yaw + wrap_pi(shot.yaw - m_yaw); // shortest way round
    m_pitch_from = m_pitch;
    m_pitch_to = shot.pitch;
    m_log_distance_from = std::log(m_distance);
    m_log_distance_to = std::log(shot.distance);
}

void OrbitCamera::angles_for_direction(const glm::dvec3& dir_world, double* yaw, double* pitch) const
{
    const glm::dvec3 local = glm::transpose(m_basis) * glm::normalize(dir_world);
    *yaw = std::atan2(local.y, local.x);
    *pitch = std::clamp(std::asin(std::clamp(local.z, -1.0, 1.0)), -kPitchLimit, kPitchLimit);
}

void OrbitCamera::look_from(const glm::dvec3& dir_world)
{
    angles_for_direction(dir_world, &m_yaw, &m_pitch);
    m_relative = false;
    m_base_yaw = m_base_pitch = 0.0;
    m_animate_view = false;
}

void OrbitCamera::rotate(double d_yaw, double d_pitch)
{
    m_animate_view = false;
    m_yaw = wrap_pi(m_yaw + d_yaw);
    m_pitch = std::clamp(m_pitch + m_base_pitch + d_pitch, -kPitchLimit, kPitchLimit) - m_base_pitch;
}

void OrbitCamera::zoom(double steps)
{
    m_animate_view = false;
    m_distance *= std::pow(0.85, steps);
}

void OrbitCamera::update(double real_dt, const Scene& scene)
{
    const glm::dvec3 target_pos = scene.bodies[static_cast<size_t>(m_target)].world_position;
    if (m_relative) {
        radial_angles(m_target, scene, &m_base_yaw, &m_base_pitch);
    }

    if (m_transition < 1.0) {
        m_transition = std::min(1.0, m_transition + real_dt / m_transition_duration);
        const double s = smoothstep(m_transition);
        const glm::dvec3 from = scene.bodies[static_cast<size_t>(m_prev_target)].world_position;
        m_focus_point = from + (target_pos - from) * s;

        if (m_animate_view) {
            // Drift keeps acting on the animated yaw so the motion stays continuous.
            m_yaw_to += m_drift * real_dt;
            m_yaw = wrap_pi(m_yaw_from + (m_yaw_to - m_yaw_from) * s);
            m_pitch = m_pitch_from + (m_pitch_to - m_pitch_from) * s;
            m_distance = std::exp(m_log_distance_from + (m_log_distance_to - m_log_distance_from) * s);
        } else {
            m_yaw = wrap_pi(m_yaw + m_drift * real_dt);
        }
    } else {
        m_focus_point = target_pos;
        m_animate_view = false;
        m_yaw = wrap_pi(m_yaw + m_drift * real_dt);
    }

    m_distance = std::clamp(m_distance, min_distance(m_target, scene), kMaxDistanceKm);

    // Collision avoidance: push the camera out of any body it would be inside.
    m_position = m_focus_point + orbit_offset();
    for (const Body& body : scene.bodies) {
        const double clearance = clearance_radius(body);
        const glm::dvec3 rel = m_position - body.world_position;
        const double d = glm::length(rel);
        if (d < clearance) {
            m_position = body.world_position + (d > 0.0 ? rel / d : m_basis[2]) * clearance;
        }
    }
}

glm::dvec3 OrbitCamera::orbit_offset() const
{
    const double yaw = m_yaw + m_base_yaw;
    const double pitch = std::clamp(m_pitch + m_base_pitch, -kPitchLimit, kPitchLimit);
    const glm::dvec3 dir_local(std::cos(pitch) * std::cos(yaw), std::cos(pitch) * std::sin(yaw), std::sin(pitch));
    return m_basis * dir_local * m_distance;
}

CameraView OrbitCamera::view(int width_px, int height_px) const
{
    // Look at the focus from the (collision-adjusted) camera position.
    glm::dvec3 offset = m_position - m_focus_point;
    if (glm::length(offset) <= 0.0) {
        offset = orbit_offset();
    }
    const glm::dvec3 up = m_basis[2];

    CameraView v;
    v.position = m_focus_point + offset;
    v.viewport = glm::vec2(static_cast<float>(std::max(width_px, 1)), static_cast<float>(std::max(height_px, 1)));

    // Rotation-only view: the camera is at the render-space origin.
    v.view = glm::mat4(glm::lookAt(glm::dvec3(0.0), -offset, up));

    // Reversed-Z infinite projection (D3D clip z in [0, 1]): z_ndc = near / distance.
    const float aspect = v.viewport.x / v.viewport.y;
    const float f = static_cast<float>(1.0 / std::tan(m_fov_y * 0.5));
    v.proj = glm::mat4(0.0f);
    v.proj[0][0] = f / aspect;
    v.proj[1][1] = f;
    v.proj[2][3] = -1.0f;
    v.proj[3][2] = static_cast<float>(m_near_km);
    v.focal_y = f;

    v.view_proj = v.proj * v.view;
    v.sky_view_proj = v.view_proj;
    return v;
}

} // namespace astraxis
