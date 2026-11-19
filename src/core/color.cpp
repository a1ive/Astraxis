#include "core/color.hpp"

#include "core/math.hpp"

#include <glm/common.hpp>

#include <algorithm>
#include <cmath>

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

namespace {

// Piecewise Gaussian lobe of Wyman et al. 2013: sigma1 below the peak, sigma2 above.
double cie_lobe(double lambda, double mu, double sigma1, double sigma2)
{
    const double t = (lambda - mu) / (lambda < mu ? sigma1 : sigma2);
    return std::exp(-0.5 * t * t);
}

} // namespace

glm::dvec3 spectrum_linear_srgb(double (*spectrum)(double lambda_nm, double param), double param)
{
    glm::dvec3 xyz(0.0);
    for (double lambda = 380.0; lambda <= 780.0; lambda += 1.0) {
        const double f = spectrum(lambda, param);
        xyz.x += f * (1.056 * cie_lobe(lambda, 599.8, 37.9, 31.0) + 0.362 * cie_lobe(lambda, 442.0, 16.0, 26.7) -
                      0.065 * cie_lobe(lambda, 501.1, 20.4, 26.2));
        xyz.y += f * (0.821 * cie_lobe(lambda, 568.8, 46.9, 40.5) + 0.286 * cie_lobe(lambda, 530.9, 16.3, 31.1));
        xyz.z += f * (1.217 * cie_lobe(lambda, 437.0, 11.8, 36.0) + 0.681 * cie_lobe(lambda, 459.0, 26.0, 13.8));
    }
    glm::dvec3 rgb(3.2406 * xyz.x - 1.5372 * xyz.y - 0.4986 * xyz.z,
                   -0.9689 * xyz.x + 1.8758 * xyz.y + 0.0415 * xyz.z,
                   0.0557 * xyz.x - 0.2040 * xyz.y + 1.0570 * xyz.z);
    rgb = glm::max(rgb, glm::dvec3(0.0));
    const double m = std::max({rgb.r, rgb.g, rgb.b});
    return m > 0.0 ? rgb / m : glm::dvec3(1.0);
}

glm::dvec3 power_law_linear_srgb(double alpha)
{
    // F_nu ~ nu^-alpha is F_lambda ~ lambda^(alpha - 2).
    return spectrum_linear_srgb([](double lambda, double a) { return std::pow(lambda / 550.0, a - 2.0); }, alpha);
}

double blackbody_band_share(double kelvin, double lo_nm, double hi_nm)
{
    constexpr double kC2NmK = 1.438776877e7; // second radiation constant hc / k (nm K), CODATA 2018
    constexpr double kNorm = 15.0 / (kPi * kPi * kPi * kPi);
    const double x_lo = kC2NmK / (hi_nm * kelvin);
    const double x_hi = std::min(kC2NmK / (lo_nm * kelvin), x_lo + 200.0); // e^-200: nothing left
    auto f = [](double x) { return x > 1e-8 ? x * x * x / std::expm1(x) : x * x; };
    // Simpson's rule on an even number of steps.
    constexpr int kSteps = 2000;
    const double h = (x_hi - x_lo) / kSteps;
    double sum = f(x_lo) + f(x_hi);
    for (int i = 1; i < kSteps; ++i) {
        sum += (i % 2 ? 4.0 : 2.0) * f(x_lo + i * h);
    }
    return std::max(0.0, kNorm * sum * h / 3.0);
}

glm::dvec3 srgb_to_linear(const glm::dvec3& c)
{
    auto channel = [](double v) { return v <= 0.04045 ? v / 12.92 : std::pow((v + 0.055) / 1.055, 2.4); };
    return {channel(c.r), channel(c.g), channel(c.b)};
}

} // namespace astraxis
