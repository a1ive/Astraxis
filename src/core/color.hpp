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

// Linear sRGB color (largest component 1) of a spectrum given as F_lambda by `spectrum`
// (any units), sampled 380-780 nm against the CIE 1931 2-degree color matching
// functions in the multi-lobe Gaussian fit of Wyman, Sloan & Shirley 2013 (JCGT 2, 1).
glm::dvec3 spectrum_linear_srgb(double (*spectrum)(double lambda_nm, double param), double param);

// Color of a power-law spectrum F_nu ~ nu^-alpha (optically thin synchrotron).
glm::dvec3 power_law_linear_srgb(double alpha);

// Share of a black body's power (sigma T^4) emitted between lo_nm and hi_nm:
// (15 / pi^4) times the integral of x^3 / (e^x - 1) over x = c2 / (lambda T)
// (c2 = hc / k, CODATA 2018).
double blackbody_band_share(double kelvin, double lo_nm, double hi_nm);

// Linear sRGB of a color written in sRGB (IEC 61966-2-1).
glm::dvec3 srgb_to_linear(const glm::dvec3& c);

} // namespace astraxis
