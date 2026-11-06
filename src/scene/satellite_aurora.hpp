#pragma once

#include "scene/aurora.hpp"
#include "scene/aurora_model.hpp"

#include <glm/geometric.hpp>
#include <glm/vec3.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace astraxis {

struct SatelliteAuroraLine {
    double wavelength_nm;
    double rayleigh;
    double uncertainty;
    glm::dvec3 rgb_per_kr;
};

// Milby et al. 2025, Table 2, 2023-08-09. RGB: aurora_colors.py.
inline const std::array<SatelliteAuroraLine, 6> kIoOxygenSulfurLines{{
    {406.9, 310.0, 18.0, {1.2049e-11, 0.0, 5.1308e-11}},
    {407.6, 87.0, 11.0, {1.3207e-11, 0.0, 5.6305e-11}},
    {458.9, 100.0, 8.0, {3.1611e-10, 0.0, 4.5313e-9}},
    {557.7, 461.0, 12.0, kGreenPerKr},
    {630.0, 6090.0, 110.0, kRedPerKr},
    {636.4, 2200.0, 20.0, {3.3823e-9, 0.0, 2.6238e-10}}
}};
inline const std::array<SatelliteAuroraLine, 2> kIoSodiumLines{{
    {589.0, 12100.0, 300.0, {8.0494e-9, 2.3081e-9, 0.0}},
    {589.6, 5070.0, 140.0, {8.1386e-9, 2.2241e-9, 0.0}}
}};
// Cassini's BL1 filter (390-500 nm) on 2001-01-05, 38 min after ingress: Milby et al. 2025
// Table 3 (from Geissler et al. 2004) gives 25 +- 10 kR against 0.5 kR from the lines above,
// which Geissler ascribes to SO2's electron-impact quasi-continuum. RGB: aurora_colors.py,
// a flat band.
inline const SatelliteAuroraLine kIoSo2Band{445.0, 24500.0, 10000.0, {0.0, 1.1788e-10, 4.6749e-9}};
// On the plasma sheet's equator: Milby 2026 (Caltech thesis), the fit of his Figure 5.5 to
// the 630.0 nm exposures of 2025-11-29, 19 +- 2 R. 636.4 nm has 1/3.09 of it (Wiese et al.
// 1996, as his Table 5.3 allows, 8.2 and 3.2 R). The other lines are upper limits.
inline constexpr double kCallistoPlasmaScaleRj = 3.6; // +- 0.7, Bagenal & Delamere 2011 predict 3.7
inline const std::array<SatelliteAuroraLine, 2> kCallistoOxygenLines{{
    {630.0, 19.0, 2.0, kRedPerKr},
    {636.4, 19.0 / 3.09, 2.0 / 3.09, {3.3823e-9, 0.0, 2.6238e-10}}
}};
// May 20 2021 disk-average oxygen lines (geomag/SOURCES.md). 636.4 nm uses the
// nearby 630 nm RGB chromaticity and photon energy ratio.
inline const std::array<SatelliteAuroraLine, 3> kEuropaOxygenLines{{
    {630.0, 573.0, 14.0, kRedPerKr},
    {636.4, 190.0, 7.0, kRedPerKr * (630.0 / 636.4)},
    {557.7, 37.0, 4.0, kGreenPerKr}
}};

template<std::size_t N>
inline glm::dvec3 satellite_line_rgb(const std::array<SatelliteAuroraLine, N>& lines)
{
    glm::dvec3 rgb(0.0);
    for (const auto& line : lines) rgb += line.rgb_per_kr * (line.rayleigh / 1000.0);
    return rgb;
}

// Exponential emitting layer of scale_km below top_km. Integrates the visible
// shell (including off-limb emission) divided by the projected solid disk
// area: this converts a zenith column to aperture R.
inline double aurora_disk_factor(double radius_km, double scale_km, double top_km)
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

struct SatelliteVisibleAurora {
    double scale_km;
    double top_km;
    double sodium_scale_km;
    glm::dvec3 column_rgb;
    glm::dvec3 sodium_column_rgb;
    glm::dvec3 so2_column_rgb; // follows the same pattern and height as column_rgb
    double display_gain;
};

// Europa's brightness follows the electrons of the plasma sheet: de Kleer et al. 2023 fit a
// scale of 1.67 RJ over single eclipses, Roth et al. 2016 1.6 RJ in the UV, Milby 2026 uses
// Bagenal & Delamere's 1.7 RJ (his fits of single nights give 1.24-1.70). The lines above
// are the mean over the eclipse of 2021-05-20, when Europa went from the equator to
// 0.683 RJ below it (Milby 2026, table 3.1).
inline constexpr double kEuropaPlasmaScaleRj = 1.7;

// Europa's brightness over the surface. Roth et al. 2016 from 71 HST/STIS images, fitted by
// tools/geomag/europa_pattern_fit.py to their ratios (Figures 7 and 8), in the 135.6 nm line:
// - brighter on the afternoon side, 1.59 between the hemispheres, always (de Kleer & Brown
//   2018 and de Kleer et al. 2023 see the same in the optical lines). The side is fixed to
//   the sun, not to the body: in eclipse the sun is behind Jupiter and this is the trailing side.
// - brighter on the side that faces the plasma sheet, up to about 4:1 seen from the trailing
//   (upstream) side and 2:1 from the leading (downstream) one.
// - gathered along the field line through the moon, which tilts with the magnetic latitude.
// In exp: sheet n.z (kSheet + kUpstream n.y) + kFieldTilt ((n.b)^2 - n.z^2); sheet = +1 when
// farthest below the magnetic equator, -1 when farthest above.
inline constexpr double kEuropaAfternoon = 0.425;
inline constexpr double kEuropaSheet = 0.575;
inline constexpr double kEuropaUpstream = 0.838;
inline constexpr double kEuropaFieldTilt = 1.925;

// shaders/aurora.frag.hlsl does the same (PATTERN_EUROPA). All vectors are in the moon's axes;
// field: jovian_field_line_at_moon, to_sun: a unit vector. Before europa_pattern_mean.
inline double europa_pattern(const glm::dvec3& n, const glm::dvec3& field, double sheet, const glm::dvec3& to_sun)
{
    const glm::dvec3 afternoon(-to_sun.y, to_sun.x, 0.0); // z x sun: east of the sub-solar point
    const double across = glm::length(afternoon);
    const double dusk = across > 1e-6 ? glm::dot(n, afternoon) / across : 0.0;
    const double along = glm::dot(n, field);
    return (1.0 + kEuropaAfternoon * dusk) *
           std::exp(sheet * n.z * (kEuropaSheet + kEuropaUpstream * n.y) + kEuropaFieldTilt * (along * along - n.z * n.z));
}

// Its mean over the sphere, so that the pattern moves the light without changing the total,
// which the plasma's density sets. Without the afternoon term, which averages out but for a
// few percent where the sheet's terms correlate with it.
inline double europa_pattern_mean(const glm::dvec3& field, double sheet)
{
    constexpr int kSamples = 256;
    constexpr double kGoldenAngle = 2.39996322972865332;
    double sum = 0.0;
    for (int k = 0; k < kSamples; ++k) {
        const double z = 1.0 - 2.0 * (k + 0.5) / kSamples;
        const double rho = std::sqrt(1.0 - z * z);
        const glm::dvec3 n(rho * std::cos(k * kGoldenAngle), rho * std::sin(k * kGoldenAngle), z);
        const double along = glm::dot(n, field);
        sum += std::exp(sheet * n.z * (kEuropaSheet + kEuropaUpstream * n.y) + kEuropaFieldTilt * (along * along - n.z * n.z));
    }
    return sum / kSamples;
}
inline double europa_reference_density()
{
    constexpr int kSteps = 64;
    constexpr double kFarthestRj = 0.683;
    double sum = 0.0;
    for (int k = 0; k < kSteps; ++k) {
        const double z = kFarthestRj * (k + 0.5) / kSteps;
        sum += std::exp(-(z * z) / (kEuropaPlasmaScaleRj * kEuropaPlasmaScaleRj));
    }
    return sum / kSteps;
}

// Io's sodium aurora, relative to before the eclipse. Schmidt et al. 2023: it falls with an
// e-folding time of 10 min after ingress, to a steady 27% of its mean over the first 10 min
// within about 30 min, and recovers with one of 2 hours after egress.
inline constexpr double kIoSodiumDeclineSeconds = 600.0;
inline constexpr double kIoSodiumRecoverySeconds = 7200.0;
inline double io_sodium_level(const MoonEclipse& eclipse, double t)
{
    // s = 0.27 (s + (1 - s) (1 - 1/e)), the mean of the decay over one e-folding time.
    const double first = 1.0 - std::exp(-1.0);
    const double steady = 0.27 * first / (1.0 - 0.27 + 0.27 * first);
    const auto in_shadow = [&](double time) {
        return steady + (1.0 - steady) * std::exp(-(time - eclipse.ingress) / kIoSodiumDeclineSeconds);
    };
    if (t < eclipse.egress) return in_shadow(t);
    return 1.0 - (1.0 - in_shadow(eclipse.egress)) * std::exp(-(t - eclipse.egress) / kIoSodiumRecoverySeconds);
}

// Io's SO2, relative to before the eclipse. de Pater et al. 2020 (ALMA): the disk-integrated
// flux falls exponentially at ingress, about 10 min (Schmidt et al. 2023), to 2-3 times less
// than in sunlight, and returns linearly within about 10 min of egress. Whether the auroral
// SO2 emission follows the bulk gas is open (geomag/SOURCES.md).
inline constexpr double kIoSo2Steady = 1.0 / 2.5;
inline constexpr double kIoSo2Seconds = 600.0;
inline double io_so2_level(const MoonEclipse& eclipse, double t)
{
    const auto in_shadow = [&](double time) {
        return kIoSo2Steady + (1.0 - kIoSo2Steady) * std::exp(-(time - eclipse.ingress) / kIoSo2Seconds);
    };
    if (t < eclipse.egress) return in_shadow(t);
    const double end = in_shadow(eclipse.egress);
    return std::min(1.0, end + (1.0 - end) * (t - eclipse.egress) / kIoSo2Seconds);
}

// The levels of the reference observations: the mean of io_sodium_level over the seven
// exposures of 2023-08-09 (Milby et al. 2025, test_io_eclipse_levels recomputes it), and
// io_so2_level when Cassini's BL1 image was taken, 38.25 min after ingress.
inline constexpr double kIoSodiumReferenceLevel = 0.2073;
inline double io_so2_reference_level()
{
    return kIoSo2Steady + (1.0 - kIoSo2Steady) * std::exp(-38.25 * 60.0 / kIoSo2Seconds);
}

inline const SatelliteVisibleAurora& europa_visible_aurora()
{
    // Roth et al. 2016 section 4.2.4: 150 km is the favored height of the 135.6 nm line
    // above the limb, about O2's. The top is assumed.
    static const SatelliteVisibleAurora model{
        150.0, 600.0, 150.0,
        satellite_line_rgb(kEuropaOxygenLines) / aurora_disk_factor(1560.8, 150.0, 600.0),
        glm::dvec3(0.0), glm::dvec3(0.0), 3.0e7
    };
    return model;
}

inline const SatelliteVisibleAurora& io_visible_aurora()
{
    // Assumed heights. The equatorial pattern averages to one over a hemisphere
    // because its x^2 and z^2 terms cancel.
    static const SatelliteVisibleAurora model{
        50.0, 600.0, 100.0,
        satellite_line_rgb(kIoOxygenSulfurLines) / aurora_disk_factor(1821.49, 50.0, 600.0),
        satellite_line_rgb(kIoSodiumLines) / aurora_disk_factor(1821.49, 100.0, 600.0),
        kIoSo2Band.rgb_per_kr * (kIoSo2Band.rayleigh / 1000.0) / aurora_disk_factor(1821.49, 50.0, 600.0),
        1.0e6
    };
    return model;
}

inline const SatelliteVisibleAurora& callisto_visible_aurora()
{
    static const SatelliteVisibleAurora model{
        50.0, 300.0, 50.0,
        satellite_line_rgb(kCallistoOxygenLines) / aurora_disk_factor(2410.3, 50.0, 300.0),
        glm::dvec3(0.0), glm::dvec3(0.0), 1.5e9
    };
    return model;
}

} // namespace astraxis
