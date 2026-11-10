#pragma once

#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <limits>
#include <vector>

namespace astraxis {

// The solar wind upstream of a planet.
struct SolarWind {
    double pressure_npa = 2.0;  // dynamic pressure
    double bz_nt = 0.0;         // the interplanetary field's north-south part (GSM); Earth's model only
    double speed_km_s = 430.0;  // sets the aberration of the tail by the planet's own motion
    // Where the boundaries were met, in place of what the model gives (0: the model's). The
    // heliosphere has no formula for them here, only these.
    double standoff_radii = 0.0;
    double shock_radii = 0.0;
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
    // The Sun's own: the heliopause, where the solar wind meets the interstellar gas, in
    // the shape of Parker 1961 (ApJ 134, 20) for incompressible flow from a point source in
    // a uniform stream, r = R / cos(theta / 2), with the nose toward the interstellar wind;
    // inside it the termination shock, a sphere. Sizes are given, not computed. Lengths are
    // usually in au.
    Heliosphere,
};

// A planet's magnetosphere ([bodies.magnetosphere] in a scene file). Lengths in `radius_km`.
struct BodyMagnetosphere {
    MagnetosphereModel model = MagnetosphereModel::None;
    double radius_km = 1.0;
    SolarWind typical;              // what is drawn without an event's own
    bool aberration = false;        // the axis leans into the planet's motion (Earth: about 4 deg)
    double tail_radii = 60.0;       // drawn this far down the tail
    // PressureBalance.
    double surface_field_nt = 0.0;
    double compressibility = 6.0;
    double flaring = 0.6;
    // Where the tail stops widening and goes on as a cylinder (0: it keeps flaring).
    double flare_end_radii = 0.0;
    // Heliosphere: the nose's distance, and its direction in ICRF (a unit vector).
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
    bool faint = false;       // for comparison with the others
    double tail_radii = 0.0;  // 0: the body's
};

// A boundary in the solar wind's axes: x upstream, z north of the orbit, y toward dusk.
// Distances in planetary radii.
struct BoundarySurface {
    enum class Shape { Flared, Conic, Quadric };
    Shape shape = Shape::Flared;
    // Flared: r = standoff (2 / (1 + cos theta))^flaring, a cylinder of tail_radius beyond
    // x = -flare_end (if set). Conic: r = standoff (1 + e) / (1 + e cos theta).
    double standoff = 1.0;
    double flaring = 0.6;
    double flare_end = 0.0;
    double tail_radius = 0.0;
    double eccentricity = 1.0;
    // Quadric: z^2 = A + B x + C x^2 + D y + E y^2 + F x y, in units of `scale` radii.
    double coefficients[6]{};
    double scale = 1.0;

    // Distance to the surface along a unit direction; infinite where it is open.
    double radius(const glm::dvec3& direction) const;
    bool contains(const glm::dvec3& point) const;
    double nose() const { return radius({1.0, 0.0, 0.0}); }
    // Half the width across the tail (along y) at x = -downtail.
    double tail_half_width(double downtail) const;
};

BoundarySurface magnetopause(const BodyMagnetosphere& body, const SolarWind& wind);
BoundarySurface bow_shock(const BodyMagnetosphere& body, const SolarWind& wind);

// The axes above (columns) in the frame of `from_sun` (the planet from the Sun) and
// `velocity` (its heliocentric velocity, km/s). wind_speed_km_s > 0: x is where the wind
// comes from as the planet sees it, turned toward the planet's motion.
glm::dmat3 solar_wind_axes(const glm::dvec3& from_sun, const glm::dvec3& velocity, double wind_speed_km_s);

// Lines over the surface down to x = -tail: `meridians` from the nose, and rings around
// the axis, the last one at the cut. xyz in radii, w = 0 on the day side .. 1 at the cut.
void boundary_wireframe(const BoundarySurface& surface, double tail, int meridians,
                        std::vector<std::vector<glm::dvec4>>& lines);

} // namespace astraxis
