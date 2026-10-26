#include "scene/refraction.hpp"

#include "scene/scene.hpp"

#include <algorithm>
#include <cmath>

namespace astraxis {

RefractedRay atmosphere_refracted_ray(const Atmosphere& a, double radius, double height)
{
    RefractedRay ray;
    if (height < 0.0 || radius <= 0.0) return ray;
    const double h = a.rayleigh_scale_height_km;
    const double tangent = radius + height;
    const double nu_t = a.refractivity * std::exp(-height / h);
    ray.impact_km = tangent * (1.0 + nu_t);
    // r = r_t + x^2 removes the tangent-point square-root singularity.
    // 32 scale heights make the omitted exponential tail negligible.
    constexpr int steps = 512;
    const double dx = std::sqrt(32.0 * std::max(h, a.haze_scale_height_km)) / steps;
    glm::dvec3 tau(0.0);
    for (int i = 0; i < steps; ++i) {
        const double x = (i + 0.5) * dx;
        const double dr = x * x;
        const double r = tangent + dr;
        const double decay = std::exp(-dr / h);
        const double nu = nu_t * decay;
        const double nr = r * (1.0 + nu);
        // Factored difference retains accuracy near the turning point.
        const double difference = dr + nu_t * (tangent * std::expm1(-dr / h) + dr * decay);
        const double weight = 4.0 * x * dx / std::sqrt(difference * (nr + ray.impact_km));
        ray.bending_rad += weight * ray.impact_km * nu / (h * (1.0 + nu));
        tau += weight * nr *
               (a.rayleigh_depth * (std::exp(-(height + dr) / h) / h) +
                a.haze_depth * (std::exp(-(height + dr) / a.haze_scale_height_km) / a.haze_scale_height_km));
    }
    ray.transmission = glm::exp(-tau);
    return ray;
}

RefractionProfile atmosphere_refraction_profile(const Atmosphere& a, double radius)
{
    RefractionProfile profile;
    profile.min_impact = 1.0 + a.refractivity;
    profile.max_impact = 1.0 + a.height_km / radius;
    for (int i = 0; i < kRefractionSamples; ++i) {
        const double b = radius * (profile.min_impact +
            (profile.max_impact - profile.min_impact) * i / (kRefractionSamples - 1));
        double lo = 0.0, hi = a.height_km;
        for (int j = 0; j < 48; ++j) {
            const double mid = 0.5 * (lo + hi);
            const double nr = (radius + mid) * (1.0 + a.refractivity * std::exp(-mid / a.rayleigh_scale_height_km));
            if (nr < b) lo = mid;
            else hi = mid;
        }
        const auto ray = atmosphere_refracted_ray(a, radius, 0.5 * (lo + hi));
        profile.samples[static_cast<size_t>(i)] = glm::vec4(ray.transmission, ray.bending_rad);
    }
    return profile;
}

} // namespace astraxis
