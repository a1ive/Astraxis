#pragma once

#include <glm/geometric.hpp>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <cmath>

namespace astraxis {

inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kTwoPi = 2.0 * kPi;
inline constexpr double kDegToRad = kPi / 180.0;
inline constexpr double kRadToDeg = 180.0 / kPi;

// Astronomical unit in km (IAU 2012 Resolution B2, exact).
inline constexpr double kAuKm = 149597870.7;

// Nominal solar radius in km (IAU 2015 Resolution B3).
inline constexpr double kSunRadiusKm = 695700.0;

// Speed of light in km/s (exact, SI definition).
inline constexpr double kSpeedOfLightKmS = 299792.458;

// Parsec in km: 648000 / pi au (IAU 2015 Resolution B2).
inline constexpr double kParsecKm = 648000.0 / 3.14159265358979323846 * kAuKm;

// Heliocentric gravitational constant in km^3/s^2 (JPL DE440, https://ssd.jpl.nasa.gov/astro_par.html).
inline constexpr double kSunGmKm3S2 = 1.32712440041279419e11;

// Wraps an angle to [0, 2*pi).
inline double wrap_two_pi(double a)
{
    a = std::fmod(a, kTwoPi);
    return a < 0.0 ? a + kTwoPi : a;
}

// Wraps an angle to [-pi, pi).
inline double wrap_pi(double a)
{
    return wrap_two_pi(a + kPi) - kPi;
}

// Active rotation of vectors by angle `a` (radians) about the x axis.
inline glm::dmat3 rotation_x(double a)
{
    const double c = std::cos(a);
    const double s = std::sin(a);
    return glm::dmat3(1.0, 0.0, 0.0, // column 0
                      0.0, c, s,     // column 1
                      0.0, -s, c);   // column 2
}

// Active rotation of vectors by angle `a` (radians) about the z axis.
inline glm::dmat3 rotation_z(double a)
{
    const double c = std::cos(a);
    const double s = std::sin(a);
    return glm::dmat3(c, s, 0.0,  // column 0
                      -s, c, 0.0, // column 1
                      0.0, 0.0, 1.0);
}

// Unit vector from right ascension / declination (radians).
inline glm::dvec3 unit_from_ra_dec(double ra, double dec)
{
    return {std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec)};
}

// Rotation from Galactic (x toward the Galactic centre, z toward the north
// Galactic pole) to ICRF axes, built from the Hipparcos definition (ESA SP-1200,
// vol. 1, sec. 1.5.3): north Galactic pole at RA 192.85948 deg, Dec 27.12825 deg,
// and Galactic longitude of the north celestial pole 122.93192 deg.
inline glm::dmat3 galactic_to_icrf()
{
    constexpr double kDeg = 3.14159265358979323846 / 180.0;
    const double ra = 192.85948 * kDeg;
    const double dec = 27.12825 * kDeg;
    const double l_ncp = 122.93192 * kDeg;
    const glm::dvec3 z(std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec));
    // The celestial pole's projection on the Galactic plane lies at longitude l_ncp.
    const glm::dvec3 p = glm::normalize(glm::dvec3(0.0, 0.0, 1.0) - z * z.z);
    const glm::dvec3 q = glm::cross(z, p);
    const glm::dvec3 x = std::cos(l_ncp) * p - std::sin(l_ncp) * q;
    return glm::dmat3(x, glm::cross(z, x), z);
}

// IAU body-fixed frame before rotation by W: z = pole at (ra, dec), x = the
// ascending node of the body equator on the ICRF equator, at RA ra + 90 deg.
// Unlike frame_from_pole, this stays defined for a pole at Dec = +-90 deg
// (e.g. Earth, alpha0 = 0, delta0 = 90: the node is at RA 90 deg).
inline glm::dmat3 iau_pole_frame(double ra, double dec)
{
    const glm::dvec3 z(std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec));
    const glm::dvec3 x(-std::sin(ra), std::cos(ra), 0.0);
    return glm::dmat3(x, glm::cross(z, x), z);
}

// Orthonormal frame whose z axis is `pole` and whose x axis is the ascending
// node of the pole's equator on the ICRF equator (the IAU convention used for
// body rotation and for Laplace-plane orbital elements). Columns are the axes.
inline glm::dmat3 frame_from_pole(const glm::dvec3& pole)
{
    glm::dvec3 x = glm::cross(glm::dvec3(0.0, 0.0, 1.0), pole);
    const double len = glm::length(x);
    x = len > 1e-12 ? x / len : glm::dvec3(1.0, 0.0, 0.0);
    const glm::dvec3 y = glm::cross(pole, x);
    return glm::dmat3(x, y, pole);
}

} // namespace astraxis
