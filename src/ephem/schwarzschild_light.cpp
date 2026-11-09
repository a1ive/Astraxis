#include "ephem/schwarzschild_light.hpp"

#include "core/math.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace astraxis {

namespace {

double tangent_impact(double mass, double r)
{
    return r / std::sqrt(1.0 - 2.0 * mass / r);
}

// Poutanen (2020), Eq. 6, over one monotonic radial leg starting at r.
// x = r / r' = 1 - t² gives a denominator sqrt(a + b*t² + O(t^4)).
// Integrating in q = sqrt(a + b*t²) also smooths out the nearly grazing
// case a -> 0. Eight Gauss-Legendre samples on [0, 1] suffice at the
// compactness of the stars in these scenes.
double leg_angle(double mass, double r, double impact, double span, double a)
{
    if (span <= 0.0 || impact == 0.0) {
        return 0.0;
    }
    const double u = 2.0 * mass / r;
    const double w = impact / r;
    const double w2 = w * w;
    const double q0 = std::sqrt(a);
    const double q1 = std::sqrt(a + w2 * (2.0 - 3.0 * u) * span);
    constexpr std::array nodes{0.0198550718, 0.1016667613, 0.2372337950, 0.4082826788,
                               0.5917173212, 0.7627662050, 0.8983332387, 0.9801449282};
    constexpr std::array weights{0.0506142681, 0.1111905172, 0.1568533229, 0.1813418917,
                                 0.1813418917, 0.1568533229, 0.1111905172, 0.0506142681};
    double sum = 0.0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        const double q = std::lerp(q0, q1, nodes[i]);
        const double t2 = span * nodes[i] * (q + q0) / (q1 + q0);
        const double correction = w2 * (t2 / q) * (t2 / q) * (3.0 * u - 1.0 - u * t2);
        sum += weights[i] / std::sqrt(1.0 + correction);
    }
    return 2.0 * w * span * sum / (q1 + q0);
}

} // namespace

double SchwarzschildStarLight::limb_angle() const
{
    return std::asin(std::min(1.0, tangent_impact(mass, radius) / tangent_impact(mass, observer)));
}

double SchwarzschildStarLight::surface_angle(double emission_cosine) const
{
    const double mu2 = emission_cosine * emission_cosine;
    const double impact = tangent_impact(mass, radius) * std::sqrt(1.0 - mu2);
    return leg_angle(mass, radius, impact, 1.0 - radius / observer, mu2);
}

double SchwarzschildStarLight::sky_deflection(double angle) const
{
    const double impact = tangent_impact(mass, observer) * std::sin(angle);
    double swept;
    if (angle >= 0.5 * kPi) {
        const double mu = std::cos(angle);
        swept = leg_angle(mass, observer, impact, 1.0, mu * mu);
    } else {
        // b = r / sqrt(1 - 2M/r) is monotonic outside the photon sphere.
        double lo = radius;
        double hi = std::min(observer, impact);
        for (int i = 0; i < 48; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (tangent_impact(mass, mid) < impact) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        const double closest = 0.5 * (lo + hi);
        swept = leg_angle(mass, closest, impact, 1.0, 0.0) +
                leg_angle(mass, closest, impact, 1.0 - closest / observer, 0.0);
    }
    return swept - (kPi - angle);
}

} // namespace astraxis
