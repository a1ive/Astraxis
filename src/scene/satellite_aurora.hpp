#pragma once

#include "scene/europa_aurora.hpp"

#include <glm/vec3.hpp>
#include <array>
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
    {557.7, 461.0, 12.0, {1.7361e-9, 6.3813e-9, 0.0}},
    {630.0, 6090.0, 110.0, {4.4807e-9, 0.0, 2.9771e-10}},
    {636.4, 2200.0, 20.0, {3.3823e-9, 0.0, 2.6238e-10}}
}};
inline const std::array<SatelliteAuroraLine, 2> kIoSodiumLines{{
    {589.0, 12100.0, 300.0, {8.0494e-9, 2.3081e-9, 0.0}},
    {589.6, 5070.0, 140.0, {8.1386e-9, 2.2241e-9, 0.0}}
}};
// de Kleer et al. 2023, Table 2, 2021-07-04. Other lines are upper limits.
inline const std::array<SatelliteAuroraLine, 2> kCallistoOxygenLines{{
    {630.0, 6.1, 1.7, {4.4807e-9, 0.0, 2.9771e-10}},
    {636.4, 3.3, 1.2, {3.3823e-9, 0.0, 2.6238e-10}}
}};

template<std::size_t N>
inline glm::dvec3 satellite_line_rgb(const std::array<SatelliteAuroraLine, N>& lines)
{
    glm::dvec3 rgb(0.0);
    for (const auto& line : lines) rgb += line.rgb_per_kr * (line.rayleigh / 1000.0);
    return rgb;
}

struct SatelliteVisibleAurora {
    double scale_km;
    double top_km;
    double sodium_scale_km;
    glm::dvec3 column_rgb;
    glm::dvec3 sodium_column_rgb;
    double display_gain;
};

inline const SatelliteVisibleAurora& io_visible_aurora()
{
    // Assumed heights. The equatorial pattern averages to one in a
    // sub-Jovian aperture because its y^2 and z^2 terms cancel.
    static const SatelliteVisibleAurora model{
        50.0, 600.0, 100.0,
        satellite_line_rgb(kIoOxygenSulfurLines) / europa_aurora_disk_factor(1821.49, 50.0, 600.0),
        satellite_line_rgb(kIoSodiumLines) / europa_aurora_disk_factor(1821.49, 100.0, 600.0),
        1.0e6
    };
    return model;
}

inline const SatelliteVisibleAurora& callisto_visible_aurora()
{
    static const SatelliteVisibleAurora model{
        50.0, 300.0, 50.0,
        satellite_line_rgb(kCallistoOxygenLines) / europa_aurora_disk_factor(2410.3),
        glm::dvec3(0.0), 3.0e9
    };
    return model;
}

} // namespace astraxis
