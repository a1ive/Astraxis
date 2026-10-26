#pragma once

#include <array>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

namespace astraxis {

struct Atmosphere;
inline constexpr int kRefractionSamples = 128;

struct RefractedRay {
    double impact_km = 0.0; // asymptotic impact parameter n(r_t) r_t
    double bending_rad = 0.0;
    glm::dvec3 transmission{0.0};
};

// Spherical exponential atmosphere, with n-1 proportional to gas density.
// Integrates bending and extinction along the curved ray. See assets/corona/SOURCES.md.
RefractedRay atmosphere_refracted_ray(const Atmosphere& atmosphere, double radius_km, double tangent_km);

struct RefractionProfile {
    double min_impact = 1.0; // planet radii
    double max_impact = 1.0;
    // Uniform in asymptotic impact parameter; RGB transmission, W bending.
    std::array<glm::vec4, kRefractionSamples> samples{};
};

RefractionProfile atmosphere_refraction_profile(const Atmosphere& atmosphere, double radius_km);

} // namespace astraxis
