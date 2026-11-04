#pragma once

#include <cmath>

namespace astraxis {

// Assumed exponential oxygen layer: H=50 km, truncated at 300 km.
// Integrate the visible shell (including off-limb emission) divided by the
// projected solid disk area. This converts a zenith column to aperture R.
inline double europa_aurora_disk_factor(double radius_km)
{
    constexpr int steps = 4096;
    constexpr double dh = 300.0 / steps;
    double sum = 0.0;
    for (int i = 0; i < steps; ++i) {
        const double h = (i + 0.5) * dh;
        const double r = 1.0 + h / radius_km;
        sum += 2.0 * r * r * (1.0 + std::sqrt(1.0 - 1.0 / (r * r))) *
               std::exp(-h / 50.0) / (50.0 * (1.0 - std::exp(-6.0))) * dh;
    }
    return sum;
}

} // namespace astraxis
