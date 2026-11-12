#pragma once

#include <glm/mat3x3.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <cstdint>
#include <filesystem>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace astraxis {

// Azimuthal terms of a tabulated surface in solar wind axes.
struct BoundaryAzimuth {
    glm::dvec2 rotation{1.0, 0.0}; // cosine and sine of the turn about +x
    double north_south = 1.0;      // +1 for the table's season, -1 for its reflection across z = 0
    double strength = 1.0;        // 0 for the azimuthal mean at each angle from +x
};

// The solar wind upstream of a planet.
struct SolarWind {
    double pressure_npa = 2.0;  // dynamic pressure
    double bz_nt = 0.0;         // the interplanetary field's north-south part (GSM); Earth's model only
    double speed_km_s = 430.0;  // sets the aberration of the tail by the planet's own motion
    // Where the boundaries were met, in place of what the model gives (0: the model's). The
    // heliosphere has no formula for them here, only these.
    double standoff_radii = 0.0;
    double shock_radii = 0.0;
    // Tabulated: where the planet's turn has brought its dipole, 0 the table's first state
    // (Neptune "Earth-like") .. 1 its second ("pole-on"); dipole_phase. Half: the middle of a turn.
    double phase = 0.5;
    BoundaryAzimuth azimuth;
};

// Surfaces known only as curves in two planes through the Sun's line, in two states of the
// planet's turn (Neptune: Mejnertsen et al. 2016, JGR Space Physics 121, 7497, Figure 4, baked
// by tools/geomag/neptune_boundaries_fig4.py). Each has four half profiles per state, the
// distance against the angle from the Sun's line toward +y, +z, -y and -z; between them the
// distance goes round as a + b cos(phi) + c sin(phi) + d cos(2 phi), the least that passes
// through the four. Beyond the last angle a profile goes on as r (1 + cos theta)^flaring = const.
class BoundaryTable {
public:
    enum Surface { Shock = 0, Magnetopause = 1 };

    bool load(const std::filesystem::path& file, std::string* error);
    // Distance to a surface along a unit direction (x toward the Sun), the two states blended.
    double radius(Surface surface, double blend, const glm::dvec3& direction,
                  const BoundaryAzimuth& azimuth = {}) const;

private:
    struct Profile {
        double flaring = 0.6;
        std::vector<double> radii; // at 0, step, 2 step ..
    };
    double at(const Profile& profile, double theta) const;

    double m_step = 0.0;       // radians
    Profile m_profiles[2][2][4]; // surface, state, azimuth
};

// The solar wind at the Earth hour by hour (NASA's OMNI2: King & Papitashvili 2005, JGR 110,
// A02104), shifted to the bow shock's nose, as baked by tools/geomag/solar_wind_bake.py into
// assets/geomag/solar_wind.bin. Many hours have no measurement, most of those before 1995.
class SolarWindRecord {
public:
    bool load(const std::filesystem::path& file, std::string* error);

    // Simulation times of the first hour's start and the last one's end.
    double start_tdb() const;
    double end_tdb() const;
    // The pressure, Bz and speed averaged over [t0, t1], each hour weighted by its overlap;
    // an hour without a measurement, or outside the record, counts as `fallback`'s. The
    // rest of the result is `fallback`'s.
    SolarWind mean(double t0_tdb, double t1_tdb, const SolarWind& fallback) const;
    // The share of [t0, t1] with a measured pressure.
    double coverage(double t0_tdb, double t1_tdb) const;

private:
    struct Hour {
        uint16_t pressure; // 0.01 nPa
        int8_t bz;         // 0.5 nT
        uint8_t speed;     // 150 + 8 n km/s
    };
    // Hours since the record's start (fractional) at a simulation time.
    double hour(double t_tdb) const;

    double m_start_jd_utc = 0.0;
    std::vector<Hour> m_hours;
};

// Empirical models of where the solar wind stops a planet's field (the magnetopause) and of
// the shock standing upstream of it.
enum class MagnetosphereModel {
    None,
    // Standoff from pressure balance, R = (B0^2 / (2 mu0 P))^(1 / compressibility), with the
    // shape r = R (2 / (1 + cos theta))^flaring (Jasinski et al. 2024 for Uranus).
    PressureBalance,
    // Shue et al. 1998, JGR 103, 17691: the Earth, from the pressure and Bz.
    Shue1998,
    // Joy et al. 2002, JGR 107, 1309: Jupiter, quadric surfaces whose coefficients depend on
    // the pressure; polar flattening and a dawn-dusk asymmetry included. For use within
    // 250 radii down the tail (the bow shock: 100).
    Joy2002,
    // Fits of spacecraft crossings to Shue's shape: standoff = standoff_at_1npa
    // P^(-1 / compressibility) and flaring + flaring_slope P (Saturn: Kanani et al. 2010;
    // Mercury: Winslow et al. 2013), or a fixed standoff_radii where the pressure's part is
    // not known. The shock may have its own power law (Saturn: Went et al. 2011).
    PowerLaw,
    // The Sun's own: the heliopause, where the solar wind meets the interstellar gas, in
    // the shape of Parker 1961 (ApJ 134, 20) for incompressible flow from a point source in
    // a uniform stream, r = R / cos(theta / 2), with the nose toward the interstellar wind;
    // inside it the termination shock, a sphere. Sizes are given, not computed. Lengths are
    // usually in au.
    Heliosphere,
    // Surfaces from a table (BoundaryTable), neither axisymmetric nor steady: they change with
    // the planet's turn (Neptune, whose dipole swings toward the Sun and away every 16 hours).
    // No dependence on the wind.
    Tabulated,
};

// A planet's magnetosphere ([bodies.magnetosphere] in a scene file). Lengths in `radius_km`.
struct BodyMagnetosphere {
    MagnetosphereModel model = MagnetosphereModel::None;
    double radius_km = 1.0;
    SolarWind typical;              // what is drawn without an event's own
    // The wind hour by hour (the Earth), for a state that follows it; null: none.
    std::shared_ptr<const SolarWindRecord> record;
    bool aberration = false;        // the axis leans into the planet's motion (Earth: about 4 deg)
    double tail_radii = 60.0;       // drawn this far down the tail
    // PressureBalance.
    double surface_field_nt = 0.0;
    double compressibility = 6.0;
    double flaring = 0.6;
    // Where the tail stops widening and goes on as a cylinder (0: it keeps flaring).
    double flare_end_radii = 0.0;
    // PowerLaw (0: unused). The shock's own law replaces shock_standoff_ratio.
    double standoff_at_1npa = 0.0;
    double flaring_slope = 0.0;
    double shock_at_1npa = 0.0;
    double shock_compressibility = 0.0;
    // The shock's focus, this far up the axis from the center (Mercury).
    double shock_focus_radii = 0.0;
    // Tabulated: the table, and what its distances are multiplied by (a simulation's surfaces
    // brought to where a spacecraft met them). Needs the body's magnetic_field, for the phase.
    std::shared_ptr<const BoundaryTable> boundaries;
    double standoff_scale = 1.0;
    double shock_scale = 1.0;
    // The center of both surfaces from the body's, in the axes below (Mercury's dipole is
    // north of the planet's center).
    glm::dvec3 center_offset_radii{0.0};
    // Heliosphere, and PowerLaw without a pressure law: the nose's distance. Heliosphere:
    // its direction in ICRF (a unit vector).
    double standoff_radii = 0.0;
    glm::dvec3 nose_icrf{1.0, 0.0, 0.0};
    // The bow shock of the first two models: a conic of revolution with its focus at the
    // planet, its nose this many magnetopause standoffs out. Heliosphere: the termination
    // shock, inside the heliopause (a ratio below 1).
    double shock_standoff_ratio = 1.275;
    double shock_eccentricity = 1.0;
    glm::vec3 color{0.55f, 0.82f, 1.0f};        // sRGB
    glm::vec3 shock_color{1.0f, 0.72f, 0.45f};  // sRGB

    bool enabled() const { return model != MagnetosphereModel::None; }
    // PressureBalance: the pressure that puts the nose at `radii`.
    double pressure_for_standoff(double radii) const;
};

// One drawing of a body's boundaries (an event's `magnetosphere` list).
struct MagnetosphereState {
    int body = -1;
    SolarWind wind;
    bool recorded = false;    // the body's record of the wind at the time; `wind` where it has none
    bool faint = false;       // for comparison with the others
    double tail_radii = 0.0;  // 0: the body's
};

// A boundary in the solar wind's axes: x upstream, z north of the orbit, y toward dusk.
// Distances in planetary radii.
struct BoundarySurface {
    enum class Shape { Flared, Conic, Quadric, Tabulated };
    Shape shape = Shape::Flared;
    // Flared: r = standoff (2 / (1 + cos theta))^flaring, a cylinder of tail_radius beyond
    // x = -flare_end (if set). Conic: r' = (standoff - focus) (1 + e) / (1 + e cos theta')
    // about a focus at x = focus; standoff is still the nose's distance from the center.
    double standoff = 1.0;
    double flaring = 0.6;
    double flare_end = 0.0;
    double tail_radius = 0.0;
    double eccentricity = 1.0;
    double focus = 0.0;
    // Quadric: z^2 = A + B x + C x^2 + D y + E y^2 + F x y, in units of `scale` radii.
    double coefficients[6]{};
    double scale = 1.0;
    // Tabulated: the two states and azimuthal terms of the table, times `scale`.
    std::shared_ptr<const BoundaryTable> table;
    BoundaryTable::Surface which = BoundaryTable::Magnetopause;
    double blend = 0.5;
    BoundaryAzimuth azimuth;

    // Distance to the surface along a unit direction; infinite where it is open.
    double radius(const glm::dvec3& direction) const;
    bool contains(const glm::dvec3& point) const;
    double nose() const { return radius({1.0, 0.0, 0.0}); }
    // Half the width across the tail (along y) at x = -downtail.
    double tail_half_width(double downtail) const;
};

BoundarySurface magnetopause(const BodyMagnetosphere& body, const SolarWind& wind);
BoundarySurface bow_shock(const BodyMagnetosphere& body, const SolarWind& wind);

// Where a planet's turn has brought its dipole, between the two phases at which the dipole
// is in the plane of the spin axis and the wind: 1 where it is nearest the wind's line
// (pole-on), 0 half a turn later. All unit vectors; the cosine of the angle between dipole
// and wind goes evenly from one to the other. An aligned dipole does not turn: 1 if it is
// nearer the wind's line than across it. The mean over a turn is one half.
double dipole_phase(const glm::dvec3& dipole, const glm::dvec3& spin, const glm::dvec3& upstream);

// The axes above (columns) in the frame of `from_sun` (the planet from the Sun) and
// `velocity` (its heliocentric velocity, km/s). wind_speed_km_s > 0: x is where the wind
// comes from as the planet sees it, turned toward the planet's motion.
glm::dmat3 solar_wind_axes(const glm::dvec3& from_sun, const glm::dvec3& velocity, double wind_speed_km_s);

// Lines over the surface down to x = -tail: `meridians` from the nose, and rings around
// the axis, the last one at the cut. xyz in radii, w = 0 on the day side .. 1 at the cut.
void boundary_wireframe(const BoundarySurface& surface, double tail, int meridians,
                        std::vector<std::vector<glm::dvec4>>& lines);

} // namespace astraxis
