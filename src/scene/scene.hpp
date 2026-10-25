#pragma once

#include "core/time.hpp"
#include "ephem/motion.hpp"
#include "scene/aurora.hpp"
#include "scene/corona.hpp"
#include "scene/pinwheel.hpp"
#include "scene/roche.hpp"
#include "scene/shape_model.hpp"
#include "scene/star_catalog.hpp"

#include <glm/mat3x3.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

#include <limits>
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
    Spacecraft, // drawn as a marker, or as its shape model if it has one
    BlackHole,  // drawn as its (black) event horizon
    Barycenter, // invisible reference point
    // A comparison point (e.g. where a star would be without relativity):
    // drawn as a hollow marker with its orbit line; no mass, light or surface.
    // Ghosts pair up among themselves (satellite_host), not with real bodies.
    Ghost,
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

    // Narrow eccentric or inclined rings (Uranus'): inner_km and outer_km are
    // then the semi-major axes of the edges, ellipses with the planet at a
    // focus. Both edges share the ring's apse and node, which precess at
    // constant rates (self-gravity keeps a ring's apses aligned). Angles count
    // from the ascending node of the planet's equator on the ICRF equator, in
    // the sense of the planet's rotation, at RingSystem::epoch_tdb. The optical
    // depth is the mean; where the ring narrows it rises (tau x width is kept).
    double inner_ae_km = 0.0; // a * e of each edge
    double outer_ae_km = 0.0;
    double peri_deg = 0.0; // longitude of periapsis
    double peri_rate_deg_per_day = 0.0;
    double a_sin_i_km = 0.0; // height above the equator 90 deg past the node
    double node_deg = 0.0;   // longitude of the ascending node
    double node_rate_deg_per_day = 0.0;

    // Arcs (Neptune's Adams arcs, Quaoar's dense Q1R): the band only spans
    // arc_length_deg around a longitude that moves at arc_rate_deg_per_day
    // (same convention as above), with soft ends. 0: a complete ring.
    double arc_length_deg = 0.0;
    double arc_center_deg = 0.0;
    double arc_rate_deg_per_day = 0.0;
    // Optional fading: the optical depth falls linearly to 0 over this span
    // (TDB seconds since J2000; equal: never fades).
    double fade_from_tdb = 0.0;
    double fade_to_tdb = 0.0;

    bool shaped() const
    {
        return inner_ae_km != 0.0 || outer_ae_km != 0.0 || a_sin_i_km != 0.0 || arc_length_deg > 0.0 ||
               fade_to_tdb > fade_from_tdb;
    }
    // Fraction of the optical depth left at t_tdb.
    double fade(double t_tdb) const;
};
inline constexpr int kMaxShapedRingBands = 16;

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

    // Epoch of the band elements (TDB seconds since J2000).
    double epoch_tdb = 0.0;

    // Normal optical depth at radius r (linear between samples, 0 outside).
    double optical_depth(double r_km) const;
    // Any band eccentric or inclined: the bands are then drawn one by one in
    // their true shapes (the profile, circular, only shades the planet).
    bool shaped() const;
};

// Fills rings.profile from rings.bands: each sample is the optical depth
// averaged over its cell (overlapping bands add up), so a ring narrower than a
// sample (e.g. Uranus' 2 km rings) keeps its equivalent width, tau x width.
// Shaped bands are taken as circles of their semi-major axes.
void rasterize_ring_bands(RingSystem& rings, int samples);

// Orientation of a shaped band at some time, in the planet's body-fixed frame
// (x, y in the equator, z along the IAU north pole).
struct RingBandShape {
    glm::dvec2 periapsis{1.0, 0.0}; // unit vector toward the periapsis
    // Toward the ascending node, of length a sin i: a point of the ring in
    // direction u (unit) lies at height z = cross(node, u) = node.x u.y - node.y u.x.
    glm::dvec2 node{0.0};
    double arc_center = 0.0; // angle from the body x axis toward y (radians, -pi..pi)
};

// Radius of an edge (semi-major axis a_km, a e = ae_km) at an angle
// `from_periapsis` past its periapsis.
double ring_edge_radius(double a_km, double ae_km, double from_periapsis);

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

// A plume of particles above a vent (or a row of vents), always erupting,
// drawn by single scattering. Positions are planetocentric (east longitude).
//   umbrella: Io's volcanic plumes. Particles on ballistic paths fill the
//     canopy, altitude < H (1 - (d / R)^2) at surface distance d from the
//     vent, concentrated in a shell under it;
//   jets: `count` cones evenly spaced from `lat_lon_deg` to `end_lat_lon_deg`
//     (Enceladus' tiger stripes), their density falling off as exp(-z / H);
//   geyser: a narrow column up to height H, whose material drifts downwind
//     in a tail at that height (Triton).
struct Plume {
    enum class Type { Umbrella, Jets, Geyser };
    std::string name;
    Type type = Type::Umbrella;
    glm::dvec2 lat_lon_deg{0.0};
    glm::dvec2 end_lat_lon_deg{0.0}; // jets: the other end of the row
    int count = 1;                   // jets
    double height_km = 1.0;          // umbrella: canopy top; jets: scale height; geyser: column top
    double radius_km = 1.0;          // umbrella: canopy radius at the ground; geyser: column radius
    double spread_deg = 10.0;        // jets: half-width of each cone
    double tail_km = 0.0;            // geyser
    double tail_azimuth_deg = 0.0;   // geyser: direction of the tail, from north through east
    // Umbrella: vertical, through the canopy above the vent; jets: across a
    // jet at its base; geyser: across the column.
    double optical_depth = 0.0;
    glm::dvec3 albedo{1.0}; // single-scattering albedo
    double g = 0.0;         // Henyey-Greenstein asymmetry
};
inline constexpr int kMaxPlumeJets = 12;

// A rotation fitted to attitude data, for a spin that changes (comet 67P,
// which turned faster after perihelion): knots of the pole and of the prime
// meridian angle W with its rate. Between knots the pole is interpolated
// linearly and W by cubic Hermite; before the first knot and after the last,
// the pole stays and W turns at the end rate (tools/philae/c67p_spin.py).
struct SpinTable {
    struct Knot {
        double t_tdb = 0.0;
        double pole_ra = 0.0; // radians
        double pole_dec = 0.0;
        double w_deg = 0.0; // counted on without wrapping
        double rate_deg_per_day = 0.0;
    };
    std::vector<Knot> knots; // ascending, at least one

    // Pole, W and dW/dt at t_tdb.
    Knot at(double t_tdb) const;
};

// A pulsating star (a Cepheid): its radius, effective temperature and
// luminosity over the cycle, tabulated by phase, and an ephemeris that puts
// table phase `epoch_phase` at epoch_tdb with the period then and its linear
// change: the cycle count since the epoch is the integral of dt / P(t),
// P(t) = period_s + period_change (t - epoch_tdb).
struct Pulsation {
    struct Sample {
        double phase = 0.0; // 0..1, ascending
        double radius_km = 1.0;
        double temperature_k = 5772.0;
        double luminosity_solar = 1.0;
    };
    std::vector<Sample> samples;
    double epoch_tdb = 0.0;
    double epoch_phase = 0.0;
    double period_s = 1.0;
    double period_change = 0.0; // dP/dt (seconds per second)

    // Cycles since the epoch (not wrapped).
    double cycles(double t_tdb) const;
    // The star at t_tdb (linear between samples, wrapping around the cycle).
    Sample at(double t_tdb) const;
    double max_temperature_k() const;
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
    // Oblateness of the gravity field (unnormalized zonal harmonics, about the
    // pole): zonal_harmonics[n] is J_n at reference radius gravity_radius_km.
    // Only an event's ghost coasts in it (CoastGravity); empty: a point mass.
    std::vector<double> zonal_harmonics;
    double gravity_radius_km = 0.0;
    double spin = 0.0;      // black holes: a / M

    // Black holes: illustrative thin accretion disk in the equatorial plane
    // (normal = pole), from disk_inner_m (0: the ISCO) to disk_outer_m (units
    // of M). 0 = no disk. Its temperature follows the thin-disk profile, or is
    // the same everywhere (disk_uniform, e.g. Interstellar's cooled disk).
    double disk_inner_m = 0.0;
    double disk_outer_m = 0.0;
    bool disk_uniform = false;
    double disk_temperature_k = 8000.0; // peak of the temperature profile, or the uniform one
    double disk_brightness = 3.0;       // radiance scale (linear HDR)

    // Stars: effective temperature (color) and luminosity (apparent brightness).
    double temperature_k = 5772.0;
    double luminosity_solar = 1.0;
    // A pulsating star: Scene::update sets the radii, temperature and
    // luminosity from it (the radii given in the file are the mean).
    std::shared_ptr<const Pulsation> pulsation;

    // IAU rotation model: spin pole in ICRF and prime meridian angle
    // W = pm_w0_deg + pm_rate_deg_per_day * d + pm_accel_deg_per_day2 * d^2
    // (d days since J2000), plus (optional) linear pole drift and the periodic
    // terms of the IAU/PCK model:
    //   RA  += ra_rate * T + sum nut_prec_ra[i]  * sin(theta_i)
    //   Dec += dec_rate * T + sum nut_prec_dec[i] * cos(theta_i)
    //   W   += sum nut_prec_pm[i] * sin(theta_i)
    // with T in Julian centuries since J2000 and theta_i = a + b T + c T^2 from
    // `nut_prec_angles` of the body itself or else the nearest ancestor that
    // has them (the planet).
    glm::dvec3 pole{0.0, 0.0, 1.0}; // at J2000, without periodic terms
    double pole_ra = 0.0; // radians; with pole_dec defines the IAU node direction
    double pole_dec = 1.5707963267948966;
    double pm_w0_deg = 0.0;
    double pm_rate_deg_per_day = 0.0;
    double pm_accel_deg_per_day2 = 0.0;        // Phobos' tidal acceleration
    glm::dvec2 pole_rate_deg_per_century{0.0}; // RA, Dec
    std::vector<glm::dvec3> nut_prec_angles;   // (a deg, b deg per century, c deg per century^2) for each theta_i
    std::vector<double> nut_prec_ra;           // degrees, one per theta_i (trailing ones may be left out)
    std::vector<double> nut_prec_dec;
    std::vector<double> nut_prec_pm;
    int nut_prec_source = -1; // body whose nut_prec_angles apply (-1: none)
    // Instead of the IAU model: the pole and W from a table. pole_ra, pole_dec
    // and pm_rate_deg_per_day then hold the first knot (the nominal rotation).
    std::shared_ptr<const SpinTable> spin_table;

    // Instead of the IAU model: torque-free precession of an axially symmetric
    // body (comet Halley's long-axis mode; 'Oumuamua's short-axis mode, body z
    // its disc's normal, spinning backward). The angular momentum M is fixed
    // along `pole`; the body z axis is inclined to it by `nutation` and turns
    // around it at precession_rate while the body rotates about z at spin_rate
    // (z-x-z Euler angles phi, theta, psi):
    //   body -> ICRF = iau_pole_frame(M) * Rz(phi) * Rx(theta) * Rz(psi),
    //   phi = phi0 + precession_rate (t - epoch), psi = psi0 + spin_rate (t - epoch).
    struct FreePrecession {
        bool enabled = false;
        double nutation_deg = 0.0;
        double precession_deg_per_day = 0.0;
        double spin_deg_per_day = 0.0;
        double epoch_tdb = 0.0;
        double phi0_deg = 0.0;
        double psi0_deg = 0.0;
    } free_precession;

    glm::vec3 color{1.0f}; // sRGB
    SurfaceStyle style = SurfaceStyle::Solid;
    glm::vec3 orbit_color{0.5f}; // sRGB
    bool label = true; // false: no name on screen unless focused (e.g. a crowd of stars)
    int event = -1;    // shown only while this event is active: its ghost (SceneEvent::Ghost) or a body
                       // it shows (SceneEvent's `show` in the file, e.g. a planet claim since refuted)
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
    // A second shape model turning about the body z axis (a spacecraft's rotating
    // ring), once every period_s of simulated time; on screen at most
    // Simulation::kRotorMaxTurnsPerSecond, so that it does not strobe at high warp.
    struct Rotor {
        std::shared_ptr<const ShapeModel> shape; // null: none
        double period_s = 0.0;
    } rotor;
    // A star filling its Roche lobe around `synchronous_with` (Algol B): its
    // surface, also drawn as `shape`; the radii are then the distances from its
    // center to L1, the sides and the poles.
    std::shared_ptr<const RocheSurface> roche;
    // Instead of a rotation model: turns with its orbit around this body,
    // x toward it and z along the orbit normal (-1: no).
    int synchronous_with = -1;

    RingSystem rings; // no bands: no rings
    Atmosphere atmosphere;
    // City lights on the night side: a gray map of emitted light in the same
    // projection as `texture` (decoded as sRGB), times `color` (sRGB) and
    // `gain` (I/F of a white pixel). They come on as the sun sets. No map: none.
    struct NightLights {
        std::string texture;
        glm::vec3 color{1.0f};
        double gain = 0.0;
    } night_lights;
    // Earth's aurora (scene/aurora.hpp), driven by this Kp record; null: none.
    std::shared_ptr<const KpRecord> aurora;
    std::vector<Plume> plumes;

    // Comets: a gas coma (C2, Haser model), an ion tail and dust (a coma and a
    // Finson-Probstein tail), all glowing in proportion to the comet's total
    // brightness law. The ion tail's width and length are near the nucleus and
    // along the tail; the dust parameters are illustrative (see comet.hpp).
    struct Comet {
        bool enabled = false;
        double m1 = 0.0; // total magnitude law (comet_total_magnitude)
        double k1 = 0.0;
        double ion_tail_diameter_km = 0.0;  // near the nucleus
        double ion_tail_length_km = 0.0;    // brightness e-folding length along the tail
        double solar_wind_km_s = 400.0;     // aberration of the tail
        double dust_tail_days = 40.0;  // age of the oldest grains drawn
        double dust_beta_min = 0.003;  // range of radiation pressure / gravity
        double dust_beta_max = 0.9;
        double dust_speed_km_s = 0.3;  // ejection speed at beta = 1, 1 au; ~ sqrt(beta / r_h)
        glm::vec3 coma_color{0.49f, 1.0f, 0.69f};     // sRGB (C2 Swan bands)
        glm::vec3 ion_tail_color{0.31f, 0.48f, 1.0f}; // sRGB (CO+)
        glm::vec3 dust_color{1.0f, 0.93f, 0.82f};     // sRGB (sunlight, slightly reddened)
    } comet;

    // Pulsars: two radio beams from the magnetic poles, inclined to the spin
    // axis, optionally precessing about the orbit normal of `spin_axis_orbit_of`.
    // The true spin (milliseconds) cannot be shown; the beams turn once every
    // display_period_s of real time while the clock runs.
    struct Pulsar {
        bool enabled = false;
        int spin_axis_orbit_of = -1;
        double display_period_s = 2.0;
        double magnetic_inclination_deg = 45.0;
        double beam_half_angle_deg = 8.0;
        bool beam_cutoff = false; // taper the outer 20% to zero at the cone edge
        struct Precession {
            bool enabled = false;
            double misalignment_deg = 0.0;
            double epoch_tdb = 0.0; // WT02 T0: spin closest to the direction away from Earth
            double rate_rad_s = 0.0; // Barker & O'Connell, from the timing orbit and masses
            glm::dvec3 axis_icrf{0.0, 0.0, 1.0};
            glm::dvec3 x_icrf{1.0, 0.0, 0.0}; // toward the projection of observer-to-pulsar
        } precession;
        glm::vec3 beam_color{0.7f, 0.8f, 1.0f}; // sRGB
    } pulsar;

    // Motion relative to the parent; null keeps the body at the scene origin.
    std::unique_ptr<MotionSource> motion;

    // Updated by Scene::update.
    bool visible = true;
    glm::dvec3 icrf_position{0.0};  // relative to the scene origin, ICRF axes
    glm::dvec3 icrf_velocity{0.0};
    glm::dvec3 world_position{0.0}; // in the active display frame (what is rendered)
    glm::dmat3 orientation{1.0};    // body-fixed axes in the display frame (columns)
    // Size of the body's orbit around its satellite_host, for satellite_fades:
    // the osculating semi-major axis around its orbit_center, or the current
    // distance from the host if the orbit is unbound or unknown.
    double fade_radius_km = 0.0;
};

// Body-fixed axes of a body in ICRF (columns) at t_tdb from its rotation
// model: IAU pole and prime meridian (periodic terms with `nut_prec_angles`,
// those of its nut_prec_source; empty: none), free precession or a spin table.
glm::dmat3 body_rotation(const Body& body, const std::vector<glm::dvec3>& nut_prec_angles, double t_tdb);

// Rate of the prime meridian (deg/day, negative if retrograde) at t_tdb: the
// spin table's, else the IAU model's secular rate (without periodic terms).
double prime_meridian_rate(const Body& body, double t_tdb);

// Shape of one of the planet's shaped ring bands at t_tdb; body_to_icrf holds
// the planet's body-fixed axes in ICRF at that time (columns).
RingBandShape ring_band_shape(const Body& planet, const RingBand& band, const glm::dmat3& body_to_icrf,
                              double t_tdb);

// How positions are presented. Inertial frames keep ICRF axes around an origin
// body; rotating frames turn with the primary -> secondary line (x axis) and
// the secondary's orbit normal (z axis), e.g. the Sun-Earth frame in which
// JWST's halo orbit around L2 is a closed loop; body-fixed frames turn with
// the origin body's own axes (its IAU rotation model: z along the north pole,
// x through the prime meridian), so that its surface stands still, e.g. the
// Pluto-fixed frame in which Charon hangs over one hemisphere.
struct DisplayFrame {
    enum class Type { Inertial, Rotating, BodyFixed };
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

// Background sky: stars, and behind them a Milky Way map (equirectangular in
// ICRF right ascension / declination, path relative to the asset directory;
// empty for none). By default the sky seen from Earth.
struct SceneSky {
    std::string milky_way = "textures/milky_way.jpg";
    double milky_way_brightness = 0.25; // linear scale of the map's values
    // Where the sky is seen from (ICRF, pc from the Sun).
    glm::dvec3 viewer_pc{0.0};
    // Stars that replace the catalog seen from Earth: the catalog from
    // viewer_pc, or a model sky made for it (tools/sky/make_galaxy_sky.py),
    // plus the members of a cluster seen from inside. Empty: from Earth.
    // Magnitudes are as displayed: a model sky's stars are shifted fainter by
    // display_mag_offset where the eye adapts to a very bright sky.
    std::vector<CatalogStar> stars;
    double display_mag_offset = 0.0;
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
    bool from_sun = false;      // view from the sun's direction in a scene without a Sun body
                                // (centered on a planet; from_body = "Sun" in the file)
    bool from_earth = false;   // from -sky.viewer_pc, with a fixed line of sight
    double phase_deg = 0.0;     // with from_body: turned this far around the up axis (180 = from behind)
    double elevation_deg = std::numeric_limits<double>::quiet_NaN(); // set: camera this far above the
                                                                     // plane normal to the up axis
    int observer = -1; // camera anchored to this body's center plus an ICRF offset
    glm::dvec3 observer_offset_km{0.0};
    double fov_deg = 45.0;
    bool corona_enhanced = false; // white-light enhancement near Parker
    // Black holes' disks: Doppler and gravitational frequency shifts (color and
    // brightness) while the event is active. Off: the disk as filmed for
    // Interstellar, which left them out.
    bool disk_shifts = true;
    // A pinwheel's dust in visible light, as the eye would see it, while the
    // event is active; otherwise in false-color near infrared.
    bool dust_visible = false;
    int flyby = -1;             // >= 0: the body whose pass of the focus the event shows (-1: looked for)
    std::string caption;        // a museum-style note on why the view is worth seeing; may be empty

    // Trail settings that replace a body's own while the event is active
    // (Scene::set_active_event), e.g. a short trail for a flyby close-up in a
    // planet-centered frame, where the earlier flybys would cross the planet,
    // or a long one where the frame caps trails short. The frame's cap does
    // not apply to them.
    struct Trail {
        int body = -1;
        TrailMode mode = TrailMode::History;
        double history_days = 0.0;
    };
    std::vector<Trail> trails;

    // The path a spacecraft would have taken without a burn (e.g. the
    // hyperbola it was on before an orbit insertion), to compare with the
    // real one: `ghost` (a Ghost body, shown only while this event is active)
    // coasts from `body`'s state at t_tdb, just before the burn, in the
    // center's gravity (its zonal harmonics, if given) and the tides of the
    // bodies the event names (CoastGravity). Its path is drawn whole, from
    // start_tdb to end_tdb.
    struct Ghost {
        int body = -1; // -1: no ghost
        int ghost = -1;
        double burn_tdb = 0.0; // the burn's start
        double t_tdb = 0.0;
        double start_tdb = 0.0;
        double end_tdb = 0.0;
        // The burn's delta-v so far (km/s, unexplained_delta_v) at even steps
        // from t_tdb over delta_v_span_s, which covers the burn; empty if the
        // burn's length is not given.
        double delta_v_span_s = 0.0;
        std::vector<double> delta_v_km_s;
    } ghost;
};

// The scene's stars as a distant observer sees them, summed into one point of
// light (scene/light_curve.hpp computes it; the hosts plot it). Each star is a
// disk of its radius, radiating as a black body at its effective temperature
// at one wavelength, with linear limb darkening; nearer stars hide those
// behind them.
struct LightCurve {
    int observer = -1;          // the body it is seen from (the Sun: from Earth); -1: no light curve
    std::vector<int> stars;     // every Star body but the observer
    std::string band = "V";     // shown with the plot
    double wavelength_nm = 550.0;
    double limb_darkening = 0.6; // u in I(mu) = 1 - u (1 - mu)
    double span_days = 1.0;      // shown at once
    double depth_mag = 1.0;      // the plot's range, from the deepest dip (set on load)
    // Magnitudes count from this light: every star as bright as it gets
    // (a pulsating one at its brightest), nothing hidden (set on load).
    double reference_flux = 0.0;

    bool enabled() const { return observer >= 0 && !stars.empty(); }
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
    int initial_event = -1; // named event selected when the scene opens
    int focus = 0;
    double distance_km = 0.0; // 0 = automatic
    double pitch_rad = 0.35;  // camera elevation above the frame's reference plane
    double yaw_rad = -1.2;    // camera azimuth from the frame's x axis
    int frame = 0;
    bool start_now = true;
    double start_tdb = 0.0;
    double warp = 3600.0;
    glm::dvec3 up_icrf{0.0}; // zero = the root body's pole
    // Star surfaces: all equally bright (only their color tells their
    // temperature), or as bright relative to each other as they are at 550 nm
    // (black bodies), the hottest at the usual brightness. Relative surfaces
    // let an eclipse show which star is in front (a cool star before a hot
    // one); a lone cool star would look almost black with them.
    bool relative_star_surfaces = false;
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
    // How times are shown and entered (top-level `calendar`: "gregorian" or
    // "kerbin"); also the calendar of the scene's own start and event times.
    Calendar calendar = Calendar::Gregorian;
    std::vector<SceneBelt> belts;
    LightCurve light_curve;
    // A colliding-wind binary's dust spiral (WR 104); disabled in most scenes.
    Pinwheel pinwheel;
    Corona corona;

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

    // Switching to another frame ends the active event (its trails were set
    // for the event's own frame).
    void set_active_frame(int index);
    int active_frame() const { return m_active_frame; }
    // The event whose trail settings apply (-1: none).
    void set_active_event(int index);
    int active_event() const { return m_active_event; }
    const DisplayFrame& frame() const { return frames[static_cast<size_t>(m_active_frame)]; }
    // Whether the axes of the active frame turn (rotating or body-fixed): its
    // z axis is then "up", and orbits are drawn as history trails.
    bool frame_is_rotating() const { return frame().type != DisplayFrame::Type::Inertial; }
    FrameTransform frame_transform(double t_tdb) const;
    // How fast the active frame's axes turn (rad per simulated second; 0 when inertial).
    double frame_turn_rate(double t_tdb) const;
    // Body-fixed axes of a body in ICRF (columns) at t_tdb: its IAU rotation
    // model, its free precession or spin table, or its orbit if it turns with
    // it (Body::synchronous_with). `body` must be one of `bodies`.
    glm::dmat3 body_to_icrf(const Body& body, double t_tdb) const;
    const FrameTransform& current_transform() const { return m_transform; }

    // State of a body relative to the scene origin (ICRF) at any time.
    State icrf_state_at(int body, double t_tdb) const;
    // State of the sun relative to the scene origin (ICRF).
    State sun_icrf_state_at(double t_tdb) const;
    // GM of the sun: the Star body's, or the Sun's in a scene centered on a planet.
    double sun_gm() const;

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

    // A body's trail settings: its own, or the active event's for it. The
    // trail drawn may still differ (a rotating frame turns orbits into history
    // trails, and a frame may cap the length of those not set by the event).
    TrailMode trail_mode(int body) const;
    // The active event's burn delta-v so far for `body` (its ghost's
    // SceneEvent::Ghost::delta_v_km_s at the current time), from a minute
    // before the burn on; false otherwise. `total_km_s`, if given, gets the
    // whole burn's.
    bool burn_delta_v(int body, double* km_s, double* total_km_s = nullptr) const;
    double trail_history_days(int body) const;

    // Trail of a body in the display frame at t: head (current position) first,
    // going back in time, with fade 0 at the head and 1 at the tail. A
    // body-fixed frame shows only the trails of its origin's system.
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
    // Whether a body belongs to the system of the active frame's origin: the
    // origin and its satellites, or everything around its barycenter.
    bool in_frame_system(int body) const;
    // The active event's trail settings for a body, or null.
    const SceneEvent::Trail* event_trail(int body) const;

    double m_time = 0.0;
    int m_active_frame = 0;
    int m_active_event = -1;
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
// camera at `camera` (display frame): a body whose orbit (fade_radius_km)
// around its satellite_host spans fewer than hide_px pixels on screen is hidden, and more
// than show_px fully shown, so that moons merge into their planet when the
// camera pulls back. Satellites of a hidden body are hidden too.
// `px_per_radian` is the screen scale (pixels per radian at the view center).
void satellite_fades(const Scene& scene, const glm::dvec3& camera, double px_per_radian, double hide_px,
                     double show_px, std::vector<float>& fades);

// Ecliptic north pole (J2000, IAU 2006 obliquity) in ICRF.
glm::dvec3 ecliptic_pole_icrf();

// Point of a body's ellipsoid (body-fixed frame, km) in the planetocentric
// direction (latitude, east longitude), and the outward unit normal there.
glm::dvec3 body_surface_point(const Body& body, double lat_deg, double lon_deg);
glm::dvec3 body_surface_normal(const Body& body, const glm::dvec3& point);

} // namespace astraxis
