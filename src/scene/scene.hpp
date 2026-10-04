#pragma once

#include "ephem/motion.hpp"
#include "scene/shape_model.hpp"
#include "scene/star_catalog.hpp"

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

// Fills rings.profile from rings.bands: each sample is the optical depth
// averaged over its cell (overlapping bands add up), so a ring narrower than a
// sample (e.g. Uranus' 2 km rings) keeps its equivalent width, tau x width.
void rasterize_ring_bands(RingSystem& rings, int samples);

// Effective wavelengths (nm) of the display channels R, G, B, at which
// atmospheric optical depths are evaluated.
inline constexpr double kChannelWavelengthsNm[3] = {680.0, 550.0, 440.0};

// An atmosphere drawn by single scattering in a shell around the body: gas
// (Rayleigh scattering, optical depth ~ lambda^-4) and haze or dust (Henyey-
// Greenstein phase function), each with an exponential density profile. An
// opaque cloud or haze deck (Venus, Titan) is drawn in the body's color instead
// of its surface, and the shell starts on top of it. Optical depths are
// vertical, from the base of the shell, per display channel.
struct Atmosphere {
    bool enabled = false;
    double deck_altitude_km = 0.0; // > 0: the opaque deck, above the surface
    double height_km = 0.0;        // top of the shell above its base
    glm::dvec3 rayleigh_depth{0.0};
    double rayleigh_scale_height_km = 1.0;
    glm::dvec3 haze_depth{0.0}; // extinction
    double haze_scale_height_km = 1.0;
    glm::dvec3 haze_albedo{1.0}; // single-scattering albedo
    double haze_g = 0.0;         // Henyey-Greenstein asymmetry (> 0: forward scattering)
    double gain = 1.0;           // scales the scattered light (1 = physical)
};

struct Body {
    std::string name;
    int parent = -1; // index into Scene::bodies; must precede this body
    BodyKind kind = BodyKind::Planet;

    double equatorial_radius_km = 1.0;   // largest (a, along the prime meridian for triaxial bodies)
    double equatorial_radius_b_km = 1.0; // second equatorial semi-axis (b); = a for spheroids
    double polar_radius_km = 1.0;
    double gm_km3_s2 = 0.0; // enables osculating-orbit trails of its children
    // GM of the body alone (its mass, shown in the info panel): body_gm_km3_s2
    // in the scene file, else gm_km3_s2 (which may include the satellites) or
    // the body's N-body GM. 0 = unknown.
    double body_gm_km3_s2 = 0.0;
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

    // Optional shape model (irregular bodies), drawn instead of the ellipsoid;
    // the radii still serve as its overall size (shadows, camera distances).
    std::shared_ptr<const ShapeModel> shape;

    RingSystem rings; // no bands: no rings
    Atmosphere atmosphere;

    // Pulsars: two radio beams from the magnetic poles, inclined to the spin
    // axis (the orbit normal of `spin_axis_orbit_of`) and sweeping around it.
    // The true spin (milliseconds) cannot be shown; the beams turn once every
    // display_period_s of real time while the clock runs.
    struct Pulsar {
        bool enabled = false;
        int spin_axis_orbit_of = -1;
        double display_period_s = 2.0;
        double magnetic_inclination_deg = 45.0;
        double beam_half_angle_deg = 8.0;
        glm::vec3 beam_color{0.7f, 0.8f, 1.0f}; // sRGB
    } pulsar;

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
    // Caps every history trail while this frame is active (0: no cap), e.g. a
    // few days of orbits in a planet-centered view of a years-long tour.
    double trail_history_days = 0.0;
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
    // Stars that replace the catalog (e.g. a globular cluster seen from inside);
    // empty: the catalog sky as seen from Earth.
    std::vector<CatalogStar> stars;
};

// A cloud of small bodies drawn as points (the main asteroid belt, the Kuiper
// belt): two-body heliocentric orbits from one epoch, J2000 ecliptic elements
// (tools/belts/make_belts.py). The renderer propagates them on the GPU;
// belt_position() is the same computation on the CPU.
struct SceneBelt {
    std::string name;
    glm::vec3 color{1.0f};   // sRGB
    double brightness = 1.0; // artistic radiance scale of one point
    double point_size_px = 1.5;
    double epoch_tdb = 0.0;  // TDB seconds since J2000
    // kStride floats per object: a (au), e, i, node, arg_peri, mean anomaly
    // at the epoch (rad), absolute magnitude H.
    static constexpr size_t kStride = 7;
    std::vector<float> elements;

    size_t size() const { return elements.size() / kStride; }
};

// Heliocentric ICRF position (km) of object `index` of a belt at t_tdb.
glm::dvec3 belt_position(const SceneBelt& belt, size_t index, double t_tdb);

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
    double phase_deg = 0.0;     // with from_body: turned this far around the up axis (180 = from behind)
};

// Osculating two-body orbit.
struct OrbitElements {
    double eccentricity = 0.0;
    double periapsis_km = 0.0;
    double apoapsis_km = 0.0; // 0 if unbound
    double period_s = 0.0;    // 0 if unbound
};

// Elements of the orbit with state `relative` (to the central body) and `gm`
// (of the central body plus the orbiting one); false if gm or r is not positive.
bool osculating_elements(const State& relative, double gm, OrbitElements* out);

// A star lighting a body, with its irradiance relative to the brightest one.
struct StarLight {
    int star = -1;
    double relative_flux = 1.0;
};
inline constexpr int kMaxStarLights = 2;
inline constexpr double kMinRelativeStarFlux = 1e-3;

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
    std::vector<SceneBelt> belts;

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

    // What a body is described as orbiting (info panel): its satellite_host,
    // except that a body farther from its barycenter than a massive sibling of
    // the host (a circumbinary planet) orbits the barycenter. `gm` is that of
    // the two-body problem (central mass plus the body's), 0 if unknown.
    struct OrbitCenter {
        int body = -1;
        double gm = 0.0;
    };
    OrbitCenter orbit_center(int body) const;

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
    // The stars that light `body`, brightest first (by L / d^2): up to
    // kMaxStarLights of the visible stars, leaving out those fainter than
    // kMinRelativeStarFlux of the brightest (e.g. both suns of a circumbinary
    // planet, but not the Sun far away). None without star bodies: the sun.
    int lighting_stars(int body, StarLight* out, int max_lights) const;
    // The brightest of them, or -1 for the scene's star / the sun.
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
