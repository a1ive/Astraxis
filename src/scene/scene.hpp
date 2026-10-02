#pragma once

#include "ephem/motion.hpp"

#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace astraxis {

// Values match the STYLE_* constants of shaders/body.frag.hlsl (2 and 3 are
// used there for stars and black holes).
enum class SurfaceStyle : int {
    Solid = 0,
    GasGiantBands = 1,
    DeathStar = 4, // procedural battle station: panels, equatorial trench, superlaser dish
};

enum class BodyKind {
    Planet,     // drawn as a lit ellipsoid (planets and moons)
    Star,       // self-luminous; the first star is the scene's light source
    Spacecraft, // drawn as a marker only
    BlackHole,  // drawn as its (black) event horizon
    Barycenter, // invisible reference point
};

enum class TrailMode {
    None,
    Orbit,   // osculating orbit around the parent (closed ellipse)
    History, // path actually travelled over the last trail_history_days
};

struct Body {
    std::string name;
    int parent = -1; // index into Scene::bodies; must precede this body
    BodyKind kind = BodyKind::Planet;

    double equatorial_radius_km = 1.0;
    double polar_radius_km = 1.0;
    double gm_km3_s2 = 0.0; // enables osculating-orbit trails of its children
    double spin = 0.0;      // black holes: a / M

    // Black holes: illustrative thin accretion disk in the equatorial plane
    // (normal = pole), from the ISCO to disk_outer_m (units of M). 0 = no disk.
    double disk_outer_m = 0.0;
    double disk_temperature_k = 8000.0; // peak of the temperature profile
    double disk_brightness = 3.0;       // radiance scale (linear HDR)

    // Stars: effective temperature (color) and luminosity (apparent brightness).
    double temperature_k = 5772.0;
    double luminosity_solar = 1.0;

    // IAU rotation model: spin pole in ICRF and prime meridian angle
    // W = pm_w0_deg + pm_rate_deg_per_day * days since J2000.
    glm::dvec3 pole{0.0, 0.0, 1.0};
    double pole_ra = 0.0; // radians; with pole_dec defines the IAU node direction
    double pole_dec = 1.5707963267948966;
    double pm_w0_deg = 0.0;
    double pm_rate_deg_per_day = 0.0;

    glm::vec3 color{1.0f}; // sRGB
    SurfaceStyle style = SurfaceStyle::Solid;
    glm::vec3 orbit_color{0.5f}; // sRGB
    TrailMode trail = TrailMode::Orbit;
    double trail_history_days = 0.0; // History trails: 0 = everything since the motion became valid
    // Draw the line of apsides at each periapsis within the trail, so the
    // apsidal advance (GR precession) shows as a fan of rotating lines.
    bool mark_periapsides = false;

    // Optional equirectangular surface map (path relative to the asset
    // directory). Longitudes increase eastward across the map, or westward if
    // texture_west_positive; texture_left_lon_deg is the left edge's longitude
    // in that same convention.
    std::string texture;
    double texture_left_lon_deg = -180.0;
    bool texture_west_positive = false;

    // Motion relative to the parent; null keeps the body at the scene origin.
    std::unique_ptr<MotionSource> motion;

    // Updated by Scene::update.
    bool visible = true;
    glm::dvec3 icrf_position{0.0};  // relative to the scene origin, ICRF axes
    glm::dvec3 world_position{0.0}; // in the active display frame (what is rendered)
    glm::dmat3 orientation{1.0};    // body-fixed axes in the display frame (columns)
};

// How positions are presented. Inertial frames keep ICRF axes around an origin
// body; rotating frames turn with the primary -> secondary line (x axis) and
// the secondary's orbit normal (z axis), e.g. the Sun-Earth frame in which
// JWST's halo orbit around L2 is a closed loop.
struct DisplayFrame {
    enum class Type { Inertial, Rotating };
    std::string name;
    Type type = Type::Inertial;
    int origin = 0;
    int primary = -1; // -1 = the sun
    int secondary = -1;
};

// display = axes^T * (icrf - origin)
struct FrameTransform {
    glm::dvec3 origin{0.0};
    glm::dmat3 axes{1.0};

    glm::dvec3 to_display(const glm::dvec3& icrf) const { return glm::transpose(axes) * (icrf - origin); }
    glm::dvec3 direction_to_display(const glm::dvec3& icrf_dir) const { return glm::transpose(axes) * icrf_dir; }
};

// Collinear Lagrange point of a primary/secondary pair.
struct Marker {
    std::string name;
    int point = 2;      // 1 or 2
    int primary = -1;   // -1 = the sun
    int secondary = -1;
    double mass_ratio = 0.0; // m_secondary / (m_primary + m_secondary)

    glm::dvec3 world_position{0.0}; // updated by Scene::update (display frame)
};

// Background sky. The catalog stars are always drawn; the optional Milky Way
// map (equirectangular in ICRF right ascension / declination, path relative to
// the asset directory) is drawn behind them.
struct SceneSky {
    std::string milky_way;
    double milky_way_brightness = 1.0; // linear scale of the map's values
};

// A moment worth jumping to (e.g. a flyby).
struct SceneEvent {
    std::string name;
    double t_tdb = 0.0;
    int focus = -1;
    double warp = 0.0;          // 0 = keep
    double distance_km = 0.0;   // 0 = keep
    int frame = -1;             // -1 = keep
    int from_orbit_normal = -1; // >= 0: view along this body's orbit normal (face-on orbit)
    int from_body = -1;         // >= 0: view from this body's direction (e.g. the Sun: as seen from Earth)
};

struct SceneView {
    int focus = 0;
    double distance_km = 0.0; // 0 = automatic
    int frame = 0;
    bool start_now = true;
    double start_tdb = 0.0;
    double warp = 3600.0;
    glm::dvec3 up_icrf{0.0}; // zero = the root body's pole
};

class Scene {
public:
    std::string name;
    std::vector<Body> bodies;
    std::vector<DisplayFrame> frames; // frames[0] is inertial around the root
    std::vector<Marker> markers;
    std::vector<SceneEvent> events;
    SceneView view;
    SceneSky sky;

    // Heliocentric motion of the scene origin, used for the sun position when
    // the scene has no Star body.
    std::unique_ptr<MotionSource> origin_heliocentric;

    void update(double t_tdb);
    // Time of the last update().
    double time() const { return m_time; }

    int find(std::string_view body_name) const;
    int star_index() const;

    void set_active_frame(int index);
    int active_frame() const { return m_active_frame; }
    const DisplayFrame& frame() const { return frames[static_cast<size_t>(m_active_frame)]; }
    bool frame_is_rotating() const { return frame().type == DisplayFrame::Type::Rotating; }
    FrameTransform frame_transform(double t_tdb) const;
    const FrameTransform& current_transform() const { return m_transform; }

    // State of a body relative to the scene origin (ICRF) at any time.
    State icrf_state_at(int body, double t_tdb) const;
    // State of the sun relative to the scene origin (ICRF).
    State sun_icrf_state_at(double t_tdb) const;

    // Sun position in the display frame.
    const glm::dvec3& sun_position() const { return m_sun_position; }

    // "Up" for the camera, in the display frame.
    glm::dvec3 up_axis() const;

    // Distance from the root at which the whole (visible) system fits.
    double system_extent_km() const { return m_system_extent_km; }

    // Trail of a body in the display frame at t: head (current position) first,
    // going back in time, with fade 0 at the head and 1 at the tail.
    void trail(int body, double t_tdb, int max_points, std::vector<glm::dvec3>& points,
               std::vector<float>& fades) const;

    // Times in [t0, t1] at which a body is closest to its parent (radial
    // velocity crossing zero upward), oldest first.
    std::vector<double> periapsis_times(int body, double t0, double t1) const;

    // Unit normal of a body's orbit around its parent (angular momentum
    // direction) at the current time, in the display frame.
    glm::dvec3 orbit_normal(int body) const;

    // Osculating Keplerian elements of a body around its parent at the current
    // time; false if it has no motion, the parent no GM, or the orbit is unbound.
    bool osculating_apsides(int body, double* periapsis_km, double* apoapsis_km, double* period_s) const;

private:
    double m_time = 0.0;
    int m_active_frame = 0;
    FrameTransform m_transform;
    glm::dvec3 m_sun_position{1.0, 0.0, 0.0};
    double m_system_extent_km = 1.0;
    mutable std::vector<double> m_scratch_times;
    mutable std::vector<glm::dvec3> m_scratch_points;
};

// Distance of a collinear Lagrange point from the secondary, in units of the
// primary-secondary separation (circular restricted three-body problem):
// the root of Szebehely's quintic for L1 (point 1) or L2 (point 2).
double lagrange_gamma(double mass_ratio, int point);

// Ecliptic north pole (J2000, IAU 2006 obliquity) in ICRF.
glm::dvec3 ecliptic_pole_icrf();

} // namespace astraxis
