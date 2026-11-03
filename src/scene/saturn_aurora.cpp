#include "scene/saturn_aurora.hpp"

#include "core/math.hpp"

#include <glm/geometric.hpp>
#include <cmath>

namespace astraxis {

glm::dvec4 SaturnAurora::profile(double local_hours)
{
    const double hour = std::remainder(local_hours, 24.0);
    const double amplitude[2][4] = {{15.1, 1.4, 0.2, 0.1}, {15.9, 2.4, 0.8, 0.3}};
    const double phase[2][3] = {{4.0, 20.7, 1.9}, {4.3, 21.9, 2.8}};
    glm::dvec4 p{};
    for (int h = 0; h < 2; ++h) {
        double colat = amplitude[h][0];
        for (int n = 1; n <= 3; ++n) {
            colat += amplitude[h][n] * std::cos(n * (hour - phase[h][n - 1]) * kPi / 12.0);
        }
        p[h] = colat * kDegToRad;
    }
    // Smooth interpolants through the reported dawn/dusk widths (~4/8 deg
    // FWHM) and 05 h brightness maximum. These are not additional fits.
    p.z = (6.0 - 2.0 * std::sin(hour * kPi / 12.0)) * kDegToRad / 2.354820045;
    p.w = std::exp(3.0 * (std::cos((hour - 5.0) * kPi / 12.0) - 1.0));
    return p;
}

glm::dvec3 SaturnAurora::noon(const glm::dvec3& sun_body)
{
    const glm::dvec3 equator(sun_body.x, sun_body.y, 0.0);
    return glm::length(equator) > 1e-12 ? glm::normalize(equator) : glm::dvec3(1, 0, 0);
}

} // namespace astraxis
