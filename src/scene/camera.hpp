#pragma once

#include "scene/scene.hpp"

#include <glm/mat3x3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

namespace astraxis {

// Per-frame camera data for rendering. The camera sits at the origin of render
// space: world positions must be made camera-relative (in double) before use.
struct CameraView {
    glm::mat4 view{1.0f};      // rotation only
    glm::mat4 proj{1.0f};      // reversed-Z, infinite far plane
    glm::mat4 view_proj{1.0f};
    glm::mat4 sky_view_proj{1.0f}; // for ICRF directions (stars): view_proj * display-from-ICRF rotation
    glm::dvec3 position{0.0};  // world (km)
    glm::vec2 viewport{1.0f};  // pixels
    float focal_y = 1.0f;      // proj[1][1] = 1 / tan(fov_y / 2)
};

// Where the camera should be: orbiting `target` at the given angles (radians,
// in the up-axis basis) and distance (km). With `relative_to_parent`, the
// angles are offsets from the direction parent -> target, so the parent stays
// in the same place behind a moving moon.
struct CameraShot {
    int target = 0;
    double yaw = 0.0;
    double pitch = 0.0;
    double distance = 1.0;
    bool relative_to_parent = false;
};

// Orbit camera around a focus body, with smooth transitions between shots.
class OrbitCamera {
public:
    // Axis that points "up" on screen (e.g. the primary's spin pole).
    void set_up_axis(const glm::dvec3& up);

    // Changes the focus body, keeping the viewing angles; animated.
    void focus(int body, const Scene& scene);
    // Animated move to a full shot over `duration` seconds (0 = immediate).
    void fly_to(const CameraShot& shot, double duration, const Scene& scene);
    // Yaw/pitch that place the camera on the side of `dir_world` from the focus.
    void angles_for_direction(const glm::dvec3& dir_world, double* yaw, double* pitch) const;
    // Places the camera on the side of `dir_world` (unit vector from the focus).
    void look_from(const glm::dvec3& dir_world);

    // User input: these cancel any running angle/distance animation.
    void rotate(double d_yaw, double d_pitch);
    void zoom(double steps);

    // Slow automatic orbiting (rad/s), used by the auto tour.
    void set_drift(double yaw_rate) { m_drift = yaw_rate; }

    int target() const { return m_target; }
    bool transitioning() const { return m_transition < 1.0; }
    double default_distance(int body, const Scene& scene) const;

    void update(double real_dt, const Scene& scene);
    CameraView view(int width_px, int height_px) const;

private:
    double min_distance(int body, const Scene& scene) const;
    // Offset from the focus point to the camera for the current angles/distance.
    glm::dvec3 orbit_offset() const;
    // Angles of the parent -> body direction; false if the body has no parent.
    bool radial_angles(int body, const Scene& scene, double* yaw, double* pitch) const;

    glm::dmat3 m_basis{1.0};
    int m_target = 0;
    int m_prev_target = 0;

    // Current view parameters. In relative mode, yaw/pitch are offsets from
    // m_base_yaw/m_base_pitch (the parent -> target direction, updated per frame).
    double m_yaw = -1.2;
    double m_pitch = 0.35;
    double m_distance = 0.0; // 0 = not yet initialized
    bool m_relative = false;
    double m_base_yaw = 0.0;
    double m_base_pitch = 0.0;

    // Transition state (focus always; angles/distance unless cancelled by input).
    double m_transition = 1.0; // 0..1
    double m_transition_duration = 1.0;
    bool m_animate_view = false;
    double m_yaw_from = 0.0, m_yaw_to = 0.0;
    double m_pitch_from = 0.0, m_pitch_to = 0.0;
    double m_log_distance_from = 0.0, m_log_distance_to = 0.0;

    double m_drift = 0.0;

    double m_fov_y = 45.0 * 3.14159265358979323846 / 180.0;
    double m_near_km = 1e-3;

    glm::dvec3 m_focus_point{0.0};
    glm::dvec3 m_position{0.0}; // camera position after collision avoidance
};

} // namespace astraxis
