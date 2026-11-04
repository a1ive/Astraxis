#pragma once

#include <cmath>

namespace astraxis {
// Knutsen et al. 2025: 101 +/- 14 R at 20 deg elevation; their dust-aware
// transfer model predicts 2.1 kR at the limb from orbit. Not a zenith column.
struct MarsAurora {
    static constexpr double observed_rayleigh = 101.0;
    static constexpr double limb_rayleigh = 2100.0;
    // Assumed diffuse SEP profile, not a measured 2024 height distribution.
    static constexpr double bottom_km = 20.0, top_km = 140.0;
    static constexpr double peak_km = 60.0, sigma_km = 15.0;
    static constexpr double display_gain = 3.0e7;
    static double profile(double h) {
        if (h < bottom_km || h > top_km) return 0.0;
        const double z = (h - peak_km) / sigma_km;
        return std::exp(-0.5 * z * z);
    }
    // Peak tangent radiance of the assumed spherical reference profile.
    static double limb_integral() {
        static const double value = [] {
            constexpr double R = 3396.19;
            double maximum = 0.0;
            for (int i = 0; i <= 240; ++i) {
                const double b = R + bottom_km + (top_km - bottom_km) * i / 240.0;
                const double end = std::sqrt((R + top_km) * (R + top_km) - b * b);
                double sum = 0.0;
                for (int j = 0; j < 2048; ++j) {
                    const double z = end * (j + 0.5) / 2048.0;
                    sum += profile(std::sqrt(b*b + z*z) - R);
                }
                const double integral = 2.0 * end * sum / 2048.0;
                if (integral > maximum) maximum = integral;
            }
            return maximum;
        }();
        return value;
    }
};
} // namespace astraxis
