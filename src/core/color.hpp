#pragma once

#include <glm/vec3.hpp>

namespace astraxis {

// Effective temperature (K) from Johnson B-V color index.
// Ballesteros (2012), EPL 97, 34008: T = 4600 (1/(0.92 BV + 1.7) + 1/(0.92 BV + 0.62)).
double temperature_from_bv(double bv);

// Linear sRGB color of a blackbody, normalized so the largest component is 1.
// Planckian locus in CIE 1931 xy from the cubic spline fit of Kim et al. (2002),
// valid 1667-25000 K (input is clamped), then XYZ -> linear sRGB (IEC 61966-2-1).
glm::dvec3 blackbody_linear_srgb(double kelvin);

// Share of a black body's power (sigma T^4) emitted between lo_nm and hi_nm:
// (15 / pi^4) times the integral of x^3 / (e^x - 1) over x = c2 / (lambda T)
// (c2 = hc / k, CODATA 2018).
double blackbody_band_share(double kelvin, double lo_nm, double hi_nm);

} // namespace astraxis
