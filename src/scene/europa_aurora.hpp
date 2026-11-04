#pragma once

#include <cmath>

namespace astraxis {

// Exponential emitting layer; defaults are Europa's assumed H=50 km and top=300 km.
// Integrate the visible shell (including off-limb emission) divided by the
// projected solid disk area. This converts a zenith column to aperture R.
inline double europa_aurora_disk_factor(double radius_km, double scale_km = 50.0, double top_km = 300.0)
{
    constexpr int steps = 4096;
    const double dh = top_km / steps;
    double sum = 0.0;
    for (int i = 0; i < steps; ++i) {
        const double h = (i + 0.5) * dh;
        const double r = 1.0 + h / radius_km;
        sum += 2.0 * r * r * (1.0 + std::sqrt(1.0 - 1.0 / (r * r))) *
               std::exp(-h / scale_km) / (scale_km * (1.0 - std::exp(-top_km / scale_km))) * dh;
    }
    return sum;
}

} // namespace astraxis
