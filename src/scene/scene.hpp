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
    Saturn = 5,    // procedural pale bands and the north polar hexagon
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

// One component of a faint planetary ring, in the planet's equatorial plane.
struct RingBand {
    std::string name;
    double inner_km = 0.0; // from the planet's center
    double outer_km = 0.0;
    double optical_depth = 0.0; // normal optical depth
    double thickness_km = 0.0;  // vertical extent: limits the brightening seen edge-on
};

// Planetary rings, drawn by single scattering for any optical depth (thin
// dusty rings like Jupiter's, or Saturn's opaque B ring). `gain` scales the
// physical brightness (I/F), which is far too faint to see for Jupiter's rings.
struct RingSystem {
    std::vector<RingBand> bands; // as listed in the scene file (if any)
    glm::vec3 color{1.0f};       // sRGB tint (single-scattering albedo)
    double gain = 1.0;
    double phase_g = 0.0; // Henyey-Greenstein asymmetry (> 0: forward scattering)

    // Radial profile used for rendering, sampled evenly over [inner_km,
    // outer_km]: x = normal optical depth, y = lower bound for the |cos| of the
    // viewing angle (band thickness / width; 0 for a thin sheet). Loaded from a
    // measured profile or rasterized from `bands`. Empty: no rings.
    std::vector<glm::vec2> profile;
    double inner_km = 0.0;
    double outer_km = 0.0;

    // Normal optical depth at radius r (linear between samples, 0 outside).
    double optical_depth(double r_km) const;
};

// Fills rings.profile from rings.bands: overlapping bands add up, and each
// band's edges are softened over 2% of its width.
void rasterize_ring_bands(RingSystem& rings, int samples);

struct Body {
    std::string name;
    int parent = -1; // index into Scene::bodies; must precede this body
    BodyKind kind = BodyKind::Planet;

    double equatorial_radius_km = 1.0;   // largest (a, along the prime meridian for triaxial bodies)
    double equatorial_radius_b_km = 1.0; // second equatorial semi-axis (b); = a for spheroids
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
    // W = pm_w0_deg + pm_rate_deg_per_day * days since J2000, plus (optional)
    // linear pole drift and the periodic terms of the IAU/PCK model:
    //   RA  += ra_rate * T + sum nut_prec_ra[i]  * sin(theta_i)
    //   Dec += dec_rate * T + sum nut_prec_dec[i] * cos(theta_i)
    //   W   += sum nut_prec_pm[i] * sin(theta_i)
    // with T in Julian centuries since J2000 and theta_i = a + b T from
    // `nut_prec_angles` of the nearest ancestor that has them (the planet).
    glm::dvec3 pole{0.0, 0.0, 1.0}; // at J2000, without periodic terms
    double pole_ra = 0.0; // radians; with pole_dec defines the IAU node direction
    double pole_dec = 1.5707963267948966;
    double pm_w0_deg = 0.0;
    double pm_rate_deg_per_day = 0.0;
    glm::dvec2 pole_rate_deg_per_century{0.0}; // RA, Dec
    std::vector<glm::dvec2> nut_prec_angles;   // planets: (a deg, b deg per century) for each theta_i
    std::vector<double> nut_prec_ra;           // degrees, one per theta_i (trailing ones may be left out)
    std::vector<double> nut_prec_dec;
    std::vector<double> nut_prec_pm;
    int nut_prec_source = -1; // body whose nut_prec_angles apply (-1: none)

    glm::vec3 color{1.0f}; // sRGB
    SurfaceStyle style = SurfaceStyle::Solid;
    glm::vec3 orbit_color{0.5f}; // sRGB
    TrailMode trail = TrailMode::Orbit;
    double trail_history_days = 0.0; // History trails: 0 = everything since the motion became valid
    double trail_linger_days = 0.0;  // History trails: how long the trail stays after the motion ends
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

    RingSystem rings; // no bands: no rings

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
    double pitch_rad = 0.35;  // camera elevation above the frame's reference plane
    double yaw_rad = -1.2;    // camera azimuth from the frame's x axis
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

    // The body that `body` is shown orbiting: its parent, except around a
    // barycenter, where the first child is the primary (shown orbiting the
    // barycenter's own host) and later children orbit that primary. -1 if none.
    int satellite_host(int body) const;

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
    // The star that lights `body`: its nearest star ancestor (e.g. Proxima for
    // Proxima b in the alpha Cen scene), or -1 for the scene's star / the sun.
    int lighting_star(int body) const;
    // Position (display frame) of the star that lights `body`.
    glm::dvec3 light_position(int body) const;

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

// How much of each body's label, marker and orbit line to show (0..1) for a
// camera at `camera` (display frame): a body whose distance from its
// satellite_host spans fewer than hide_px pixels on screen is hidden, and more
// than show_px fully shown, so that moons merge into their planet when the
// camera pulls back. Satellites of a hidden body are hidden too.
// `px_per_radian` is the screen scale (pixels per radian at the view center).
void satellite_fades(const Scene& scene, const glm::dvec3& camera, double px_per_radian, double hide_px,
                     double show_px, std::vector<float>& fades);

// Ecliptic north pole (J2000, IAU 2006 obliquity) in ICRF.
glm::dvec3 ecliptic_pole_icrf();

} // namespace astraxis
