#pragma once

#include <glm/vec2.hpp>

namespace astraxis {

// White-light K corona around the scene's Sun. Density structures are
// illustrative; the radial profile and finite-disk scattering have physical units.
struct Corona {
    bool enabled = false;
    double display_gamma = 0.65; // compresses the radial contrast for display
    double display_gain = 8.0e5; // mean solar brightness -> linear HDR, enhanced view
};

inline constexpr double kCoronaOuterRadius = 32.0; // integration cutoff, solar radii
inline constexpr int kCoronaRadialCount = 128;

// Leblanc, Dulk & Bougeret (1998), electron density in cm^-3, r in solar radii.
double corona_electron_density(double r);

// Coefficients A+B*cos(scattering_angle)^2, per solar radius of path, in mean
// solar brightness. Integrates a finite solar disk with linear limb darkening.
glm::dvec2 corona_scattering(double r, double limb_darkening = 0.6);

// Illustrative interpolation through independently published cycle dates;
// 0 = minimum, 1 = maximum. Outside the anchor span the endpoints are held.
double corona_activity(double t_tdb);

// Fraction of a limb-darkened disk visible behind one circular occulter.
// Radii and separation use the same units (angles for a distant eclipse).
double corona_visible_fraction(double sun_radius, double occulter_radius, double separation);

} // namespace astraxis
