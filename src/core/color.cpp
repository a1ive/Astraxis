#include "core/color.hpp"

#include <glm/common.hpp>

#include <algorithm>

namespace astraxis {

double temperature_from_bv(double bv)
{
    return 4600.0 * (1.0 / (0.92 * bv + 1.7) + 1.0 / (0.92 * bv + 0.62));
}

glm::dvec3 blackbody_linear_srgb(double kelvin)
{
    const double t = std::clamp(kelvin, 1667.0, 25000.0);
    const double t2 = t * t;
    const double t3 = t2 * t;

    double x = 0.0;
    if (t <= 4000.0) {
        x = -0.2661239e9 / t3 - 0.2343589e6 / t2 + 0.8776956e3 / t + 0.179910;
    } else {
        x = -3.0258469e9 / t3 + 2.1070379e6 / t2 + 0.2226347e3 / t + 0.240390;
    }

    const double x2 = x * x;
    const double x3 = x2 * x;
    double y = 0.0;
    if (t <= 2222.0) {
        y = -1.1063814 * x3 - 1.34811020 * x2 + 2.18555832 * x - 0.20219683;
    } else if (t <= 4000.0) {
        y = -0.9549476 * x3 - 1.37418593 * x2 + 2.09137015 * x - 0.16748867;
    } else {
        y = 3.0817580 * x3 - 5.87338670 * x2 + 3.75112997 * x - 0.37001483;
    }

    // xyY (Y = 1) -> XYZ -> linear sRGB.
    const double X = x / y;
    const double Y = 1.0;
    const double Z = (1.0 - x - y) / y;
    glm::dvec3 rgb(3.2406 * X - 1.5372 * Y - 0.4986 * Z,
                   -0.9689 * X + 1.8758 * Y + 0.0415 * Z,
                   0.0557 * X - 0.2040 * Y + 1.0570 * Z);
    rgb = glm::max(rgb, glm::dvec3(0.0));
    const double m = std::max({rgb.r, rgb.g, rgb.b});
    return m > 0.0 ? rgb / m : glm::dvec3(1.0);
}

} // namespace astraxis
