#pragma once

#include "core/time.hpp"

#include <array>
#include <cmath>
#include <utility>

namespace astraxis {
// Knutsen et al. 2025: 101 +/- 14 R at 20 deg elevation; their dust-aware
// transfer model predicts 2.1 kR at the limb from orbit. Not a zenith column.
struct MarsAurora {
    static constexpr double observed_rayleigh = 101.0;
    static constexpr double limb_rayleigh = 2100.0;
    // The emission behind the CO2+ UVD limb curve of Schneider et al. 2018 Fig. 2, a Gaussian
    // with a sigma per side fitted through the line-of-sight integral
    // (tools/geomag/mars_aurora_profile.py). Green line assumed to have the same shape, as in
    // Knutsen's own model (they use the 297.2 nm line, same event).
    static constexpr double bottom_km = 20.0, top_km = 140.0;
    static constexpr double peak_km = 63.0, sigma_below_km = 12.0, sigma_above_km = 17.0;
    static constexpr double display_gain = 3.0e7;
    static double profile(double h) {
        if (h < bottom_km || h > top_km) return 0.0;
        const double z = (h - peak_km) / (h < peak_km ? sigma_below_km : sigma_above_km);
        return std::exp(-0.5 * z * z);
    }
    // Peak tangent radiance of the reference profile.
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

    // September 2017, Schneider 2018 Fig. 3: (day of the month, CO Cameron limb kR), the brightest
    // of each 6 hours, noise floor 0.04 kR (tools/geomag/mars_aurora_profile.py).
    static constexpr std::array<std::pair<double, double>, 42> kSeptember2017{{
        {7.01, 0.040}, {8.43, 0.040}, {10.40, 0.100}, {11.21, 0.099}, {11.27, 0.047}, {11.58, 0.083},
        {11.88, 0.287}, {12.01, 0.069}, {12.50, 0.103}, {12.99, 0.786}, {13.12, 3.906}, {13.43, 0.885},
        {13.61, 1.503}, {13.86, 0.868}, {14.17, 0.220}, {14.47, 0.283}, {14.60, 0.139}, {14.91, 0.049},
        {15.21, 0.083}, {15.34, 0.142}, {16.69, 0.040}, {17.18, 0.051}, {17.43, 0.070}, {17.87, 0.040},
        {18.17, 0.091}, {18.54, 0.064}, {18.92, 0.073}, {19.28, 0.053}, {20.39, 0.050}, {20.89, 0.099},
        {21.14, 0.145}, {21.26, 0.299}, {21.50, 0.312}, {22.00, 0.172}, {22.38, 0.043}, {22.74, 0.078},
        {23.48, 0.040}, {23.73, 0.040}, {24.22, 0.071}, {24.47, 0.040}, {25.21, 0.057}, {25.58, 0.040},
    }};
    static constexpr double kSeptember2017PeakKr = 3.906; // orbit 5730, the limb scan of Fig. 2
    // Knutsen's 2024 event is 10 times fainter than this peak at 297.2 nm: the green line is
    // inferred from the ultraviolet, nobody saw this one in visible light.
    static constexpr double kSeptember2017PeakFactor = 10.0;

    // Green line at the limb at this time (R). The reference observation's level at any other
    // date, reached over a day on each side of the record.
    static double limb_rayleigh_at(double t_tdb) {
        static const double month_start =
            (jd_from_calendar({2017, 9, 1, 0, 0, 0}) - 1.0 - kJ2000Jd) * kSecondsPerDay; // day 0
        const double day = (t_tdb - month_start) / kSecondsPerDay;
        const auto& record = kSeptember2017;
        const double peak_rayleigh = limb_rayleigh * kSeptember2017PeakFactor;
        const double reference_kr = kSeptember2017PeakKr / kSeptember2017PeakFactor;
        const auto log_mix = [](double day_now, double d0, double a, double d1, double b) {
            return a * std::pow(b / a, (day_now - d0) / (d1 - d0));
        };
        const double first = record.front().first - 1.0;
        const double last = record.back().first + 1.0;
        if (day <= first || day >= last) return limb_rayleigh;
        double kr = 0.0;
        if (day < record.front().first) {
            kr = log_mix(day, first, reference_kr, record.front().first, record.front().second);
        } else if (day >= record.back().first) {
            kr = log_mix(day, record.back().first, record.back().second, last, reference_kr);
        } else {
            for (size_t k = 0; k + 1 < record.size(); ++k) {
                if (day < record[k + 1].first) {
                    kr = log_mix(day, record[k].first, record[k].second, record[k + 1].first, record[k + 1].second);
                    break;
                }
            }
        }
        return peak_rayleigh * kr / kSeptember2017PeakKr;
    }
};
} // namespace astraxis
