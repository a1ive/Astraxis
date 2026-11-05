#pragma once

#include "core/math.hpp"

#include <glm/geometric.hpp>
#include <glm/vec3.hpp>

#include <algorithm>
#include <array>
#include <cmath>

namespace astraxis {

// Where a moon is in Jupiter's magnetosphere. Positions are in System III axes, the
// scene's body frame of Jupiter: z along the spin, east longitude = atan2(y, x).

// The dipole of VIP4 (Connerney et al. 1998), 9.5 deg toward 159.2 E (200.8 W). Marzok et
// al. 2022's magnetic latitudes are this dipole's (test_ganymede_aurora, their table 1).
inline constexpr double kJovianDipoleTiltDeg = 9.5;
inline constexpr double kJovianDipoleEastDeg = 159.2;

inline double jovian_magnetic_latitude_deg(const glm::dvec3& system_iii)
{
    const double tilt = kJovianDipoleTiltDeg * kDegToRad;
    const double east = kJovianDipoleEastDeg * kDegToRad;
    const glm::dvec3 axis(std::sin(tilt) * std::cos(east), std::sin(tilt) * std::sin(east), std::cos(tilt));
    return std::asin(std::clamp(glm::dot(glm::normalize(system_iii), axis), -1.0, 1.0)) / kDegToRad;
}

// Height above the plasma sheet's centrifugal equator (units of Jupiter's 71,492 km radius),
// positive to the north: the tangent approximation of Milby 2026 (Caltech thesis) Eqs. 1.6
// and 1.11, after Phipps & Bagenal 2021, good to 2%. The sign follows his tables 1 and 3.1.
inline double centrifugal_equator_height_rj(const glm::dvec3& system_iii)
{
    constexpr double kJupiterRadiusKm = 71492.0;
    const double distance = glm::length(system_iii);
    const double r = distance / kJupiterRadiusKm;
    const double east = std::atan2(system_iii.y, system_iii.x);
    const double offset_deg = (1.66 * std::tanh(0.131 * r - 1.62) + 7.76) * std::cos(east - kJovianDipoleEastDeg * kDegToRad);
    return r * std::tan(std::asin(std::clamp(system_iii.z / distance, -1.0, 1.0)) + offset_deg * kDegToRad);
}

// The electrons around the equator, relative to on it: Gaussian (Hill & Michel 1976).
inline double plasma_sheet_density(double height_rj, double scale_rj)
{
    return std::exp(-(height_rj * height_rj) / (scale_rj * scale_rj));
}

// Its mean over System III longitude, for a moon on a circular equatorial orbit.
inline double plasma_sheet_mean_density(double radius_km, double scale_rj)
{
    constexpr int kSteps = 72;
    double sum = 0.0;
    for (int k = 0; k < kSteps; ++k) {
        const double east = kTwoPi * (k + 0.5) / kSteps;
        sum += plasma_sheet_density(centrifugal_equator_height_rj(radius_km * glm::dvec3(std::cos(east), std::sin(east), 0.0)), scale_rj);
    }
    return sum / kSteps;
}

// How much of the moon's plasma-sheet state to draw, 0..1. Jupiter takes a moon through the
// sheet twice per its roughly 10.5 hours relative to the moon; past a cycle per second of
// real time that would be flicker, so the views fall back on the long-term mean.
inline double plasma_sheet_steadiness(double warp)
{
    constexpr double kCycleSeconds = 5.25 * 3600.0;
    const double x = std::clamp((std::abs(warp) / kCycleSeconds - 0.25) / 0.75, 0.0, 1.0);
    return 1.0 - x * x * (3.0 - 2.0 * x);
}

// Shares of Marzok et al. 2022's maps of Ganymede above, inside and below the current sheet.
// Their boundary is 6 deg of magnetic latitude, section 3.2. It is not a sharp edge in the
// aurora, so the views blend over two degrees (half an hour of the moon's time).
inline std::array<double, 3> ganymede_sheet_shares(double magnetic_latitude_deg)
{
    const auto beyond = [](double latitude) {
        const double x = std::clamp((latitude - 5.0) / 2.0, 0.0, 1.0);
        return x * x * (3.0 - 2.0 * x);
    };
    const double above = beyond(magnetic_latitude_deg);
    const double below = beyond(-magnetic_latitude_deg);
    return {above, 1.0 - above - below, below};
}

} // namespace astraxis
