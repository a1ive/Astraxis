#pragma once

#include <glm/vec3.hpp>

namespace astraxis {

// Haser scale lengths of C2 at 1 au, parent and daughter (km), both growing as
// r_h^2: A'Hearn et al. 1995 (Icarus 118, 223), as tabulated in Gilbert et al.
// 2010 (arXiv:0910.0416, Table 3), with an outflow speed of 1 km/s.
inline constexpr double kC2ParentScaleKm = 2.2e4;
inline constexpr double kC2DaughterScaleKm = 6.6e4;

// Total visual magnitude of a comet seen from 1 au, at `r_au` from the sun:
// m = M1 + 5 log10(delta) + K1 log10(r) with delta = 1 (JPL SBDB's comet
// magnitude parameters).
double comet_total_magnitude(double m1, double k1, double r_au);

// Direction of an ion (type I) tail: along the solar wind as the comet sees it,
// the radial wind at `solar_wind_km_s` minus the comet's own heliocentric
// velocity (the aberration of the tail). Heliocentric position and velocity in
// any one set of axes; the result is a unit vector in the same axes.
glm::dvec3 ion_tail_direction(const glm::dvec3& helio_position, const glm::dvec3& helio_velocity,
                              double solar_wind_km_s);

} // namespace astraxis
