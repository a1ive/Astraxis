#pragma once

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <string>

namespace astraxis {

// A body's relief map (`[bodies.relief]`): the slopes of its surface and its open water, in
// the projection of its albedo map, baked from a height map by tools/textures/relief_map.py.
// Red and green hold the slopes toward the east and the north (rise over run) as
// 0.5 + 0.5 sign(s) sqrt(|s| / kReliefMaxSlope), blue the share of the pixel that is water.
// body.frag.hlsl tilts the normal by the slopes: the relief shows near the terminator, where
// a few degrees of slope are a large part of the sun's elevation. There are no cast shadows.
struct SurfaceRelief {
    std::string texture;       // empty: none
    double exaggeration = 1.0; // of the slopes (1: as measured)
};

// Sunlight reflected off the water of the relief map (`[bodies.sea]`): the facets of a
// wind-roughened surface (Cox & Munk 1954) that happen to mirror the sun toward the viewer.
struct SeaSurface {
    double wind_m_s = -1.0; // negative: no reflection
    double refractive_index = 1.0;
};

// Clouds over the surface (`[bodies.clouds]`): a gray map of their cover in the projection
// of the albedo map, drawn with the atmosphere (not at all on a body without one). Where
// they are, the surface gives way to a smooth white layer of this albedo, which also hides
// the sea's glint and the city lights.
// No shadows, no height. The picture is a still one, set flowing along the parallels for
// the look of it: westward in the tropics and eastward in the middle latitudes, at the two
// rates below (degrees of longitude per day, east positive; both 0: the picture stays
// still on the ground).
struct CloudLayer {
    std::string texture; // empty: none
    double albedo = 0.0;
    double tropical_drift_deg_per_day = 0.0;
    double midlatitude_drift_deg_per_day = 0.0;
};

// The flow would wind a still picture up without end, so it runs in cycles of this many
// days of drift: two copies of the picture half a cycle apart, each sheared from half a
// cycle back to half a cycle ahead and faded in and out so that it is unseen when it starts
// over. Simulation::cloud_phase() is the place in the cycle; it follows the clock, but no
// faster on screen than kCloudMaxDaysPerSecond of drift per second.
inline constexpr double kCloudCycleDays = 3.0;
inline constexpr double kCloudMaxDaysPerSecond = 1.0 / 12.0;

// The drift's rate at a latitude, degrees of longitude per day: the tropical rate at the
// equator, the mid-latitude one at 45 degrees, none at the poles. As in body.frag.hlsl.
inline double cloud_drift_deg_per_day(const CloudLayer& clouds, double lat_rad)
{
    const double c = std::cos(lat_rad);
    const double s2 = std::sin(2.0 * lat_rad);
    const double c2 = c * c;
    return clouds.midlatitude_drift_deg_per_day * s2 * s2 + clouds.tropical_drift_deg_per_day * c2 * c2 * c2;
}

// One of the two copies (0 or 1) at a place in the cycle (0..1): how far east its picture
// has moved at this latitude (degrees) and its weight. The weights of the two add up to 1.
struct CloudCopy {
    double shift_deg = 0.0;
    double weight = 0.0;
};
inline CloudCopy cloud_copy(const CloudLayer& clouds, double lat_rad, double phase, int copy)
{
    double p = phase + 0.5 * copy;
    p -= std::floor(p);
    return {cloud_drift_deg_per_day(clouds, lat_rad) * kCloudCycleDays * (p - 0.5), 1.0 - std::abs(2.0 * p - 1.0)};
}

// What a cloud cover (0..1, the map's value) leaves of the sunlight mirrored by the sea
// below it, which crosses the layer twice, and of the light of the cities, which crosses
// it once. As in body.frag.hlsl.
inline double cloud_glint_share(double cover)
{
    return (1.0 - cover) * (1.0 - cover);
}
inline double cloud_lights_share(double cover)
{
    return 1.0 - cover;
}

inline constexpr double kReliefMaxSlope = 1.0;

// The slope from a channel of the relief map (0..1).
inline double relief_decode(double encoded)
{
    const double s = 2.0 * encoded - 1.0;
    return kReliefMaxSlope * s * std::abs(s);
}

// The surface normal tilted by the slopes toward the east and the north: n is the smooth
// surface's unit normal and pole the body's z axis, in the same axes. As in body.frag.hlsl.
inline glm::dvec3 relief_normal(const glm::dvec3& n, const glm::dvec3& pole, double east_slope, double north_slope)
{
    const glm::dvec3 across = glm::cross(pole, n);
    const double len = glm::length(across);
    if (len < 1e-6) {
        return n; // at a pole east and north mean nothing
    }
    const glm::dvec3 east = across / len;
    const glm::dvec3 north = glm::cross(n, east);
    return glm::normalize(n - east_slope * east - north_slope * north);
}

// Mean square slope of the sea (both components together) under a wind of `wind_m_s` at
// 12.5 m: Cox & Munk 1954, J. Opt. Soc. Am. 44, 838, as given by Mobley's Ocean Optics Web
// Book ("Cox-Munk sea surface slope statistics", Eq. 4): 0.003 + 5.12e-3 U, +-0.004.
inline double sea_slope_variance(double wind_m_s)
{
    return 0.003 + 5.12e-3 * std::max(wind_m_s, 0.0);
}

// Reflectance of a smooth surface of refractive index n at normal incidence (Fresnel).
inline double fresnel_normal_reflectance(double n)
{
    return (n - 1.0) * (n - 1.0) / ((n + 1.0) * (n + 1.0));
}

// Smith's shadowing term of a surface with Gaussian slopes of variance `variance`, seen or
// lit at cosine mu from its mean normal: the rational fit of Walter et al. 2007 (EGSR,
// "Microfacet Models for Refraction through Rough Surfaces", Eq. 27).
inline double sea_shadowing(double mu, double variance)
{
    const double sine = std::sqrt(std::max(1.0 - mu * mu, 0.0));
    if (mu <= 0.0) {
        return 0.0;
    }
    const double a = mu / std::max(std::sqrt(variance) * sine, 1e-9);
    return a < 1.6 ? (3.535 * a + 2.181 * a * a) / (1.0 + 2.276 * a + 2.577 * a * a) : 1.0;
}

// The glint of the sun off the sea, as the albedo of a Lambert surface that would be as
// bright (times the cosine of the sun's zenith angle, like the shader's diffuse term):
// pi f mu_sun, with f the reflectance of a surface of mirror facets whose slopes are
// Gaussian (Cox & Munk), f = F p(slope) G / (4 mu_sun mu_view cos^4 beta). n, to_sun and
// to_viewer are unit vectors. As in body.frag.hlsl.
inline double sea_glint(const glm::dvec3& n, const glm::dvec3& to_sun, const glm::dvec3& to_viewer, double variance,
                        double normal_reflectance)
{
    const double mu_sun = glm::dot(n, to_sun);
    const double mu_view = glm::dot(n, to_viewer);
    if (mu_sun <= 0.0 || mu_view <= 0.0 || variance <= 0.0) {
        return 0.0;
    }
    const glm::dvec3 half = glm::normalize(to_sun + to_viewer);
    const double c = std::max(glm::dot(n, half), 1e-4); // cos beta, the facet's tilt
    const double c2 = c * c;
    const double facets = std::exp(-(1.0 - c2) / (c2 * variance)) / (variance * c2 * c2); // pi p / cos^4 beta
    // Schlick's approximation of the Fresnel reflectance at the facet.
    const double fresnel =
        normal_reflectance + (1.0 - normal_reflectance) * std::pow(1.0 - std::min(glm::dot(to_viewer, half), 1.0), 5.0);
    return fresnel * facets * sea_shadowing(mu_sun, variance) * sea_shadowing(mu_view, variance) / (4.0 * mu_view);
}

} // namespace astraxis
