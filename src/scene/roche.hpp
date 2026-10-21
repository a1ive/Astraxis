#pragma once

#include "scene/shape_model.hpp"

#include <glm/vec3.hpp>

#include <vector>

namespace astraxis {

// A star that fills its Roche lobe in a circular binary and turns with the
// orbit (Algol B): its surface is the equipotential through the inner
// Lagrangian point L1 of the potential, in units of the separation a and
// G M_star / a, of the star (origin) and its companion (at x = 1), turning
// about the z axis with the orbit:
//   Omega = 1 / r1 + q (1 / r2 - x) + (1 + q) / 2 (x^2 + y^2),  q = M_companion / M_star
// (Kopal's normalization). Local gravity is |grad Omega|; the surface
// temperature follows gravity darkening, T^4 ~ g^gravity_darkening, scaled so
// that the star gives out the luminosity of its effective temperature, and is
// raised where the companion shines on it: T^4 += albedo T_c^4 (R_c / s)^2 cos(gamma),
// the companion taken as a point at distance s, gamma its angle from the normal.
struct RocheStar {
    double mass_ratio = 1.0;    // q = M_companion / M_star
    double separation_km = 1.0; // a
    double temperature_k = 5772.0;  // effective, before irradiation
    double gravity_darkening = 0.32; // exponent g (0.32 for a convective envelope, Lucy 1967)
    double albedo = 0.5;             // bolometric, for the companion's light
    double companion_radius_km = 0.0;
    double companion_temperature_k = 0.0;
};

// Wavelength (nm) of RocheSurface's vertex brightness and facet radiance:
// that at which star surfaces are compared (view/scene_renderer.cpp).
inline constexpr double kRocheRadianceWavelengthNm = 550.0;

// The surface as a mesh in the star's frame: x toward the companion, z along
// the orbital angular momentum (km).
struct RocheSurface {
    struct Facet {
        glm::dvec3 center{0.0}; // km
        glm::dvec3 normal{0.0}; // unit, outward
        double area_km2 = 0.0;
        double temperature_k = 0.0;
        double radiance = 0.0; // planck_shape at kRocheRadianceWavelengthNm
    };

    // Vertex albedo: the surface radiance at kRocheRadianceWavelengthNm over
    // that of the effective temperature (gravity darkening and irradiation).
    ShapeModel mesh;
    std::vector<Facet> facets;
    std::vector<double> vertex_temperature_k;
    // From the center (km): to L1 (the farthest point), the far side (-x),
    // the sides (y) and the poles (z); the radius of a sphere of equal volume.
    double l1_km = 0.0;
    double back_km = 0.0;
    double side_km = 0.0;
    double pole_km = 0.0;
    double volume_radius_km = 0.0;
};

// Distance of L1 from the star, in units of the separation.
double roche_l1(double mass_ratio);
// Distance from the star's center to its critical surface along `direction`
// (unit, in the star's frame), in units of the separation.
double roche_radius(double mass_ratio, const glm::dvec3& direction);

// The critical surface on `rings` circles of latitude about the x axis (L1 is
// the vertex at one end) with `segments` vertices each.
RocheSurface make_roche_surface(const RocheStar& star, int rings = 64, int segments = 128);

} // namespace astraxis
