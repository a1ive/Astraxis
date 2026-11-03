#pragma once

#include <glm/mat3x3.hpp>
#include <glm/vec2.hpp>

#include <array>
#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// Head et al. 2024, A&A 688 A205: Juno-UVS PJ1-54 reference contours.
// Each radius is a planetocentric angle around its hemisphere's reference
// point, sampled uniformly in tangent-plane azimuth. See geomag/SOURCES.md.
struct JovianAurora {
    static constexpr int kSamples = 256;
    std::array<glm::dvec2, kSamples> radii{}; // north, south (rad)

    bool load(const std::filesystem::path& file, std::string* error);
    double radius(bool south, double azimuth) const;
    // Columns: local east, north, center, in right-handed body coordinates.
    static glm::dmat3 frame(bool south);

    // JRM33 + CON2020, at 900 km above the one-bar spheroid. A periodic
    // longitude grid retains the actual moon's radius and latitude.
    static constexpr int kLongitudes = 360;
    static constexpr std::array<const char*, 3> kMoons{"Io", "Europa", "Ganymede"};
    struct Footpath {
        double min_radius = 0.0, max_radius = 0.0; // in 71492 km radii
        double min_latitude = 0.0, max_latitude = 0.0; // radians
        std::vector<std::array<glm::dvec3, 2>> points;
    };
    std::array<Footpath, 3> footpaths;
    bool load_footprints(const std::filesystem::path& file, std::string* error);
    // Positive eastward equatorial lead (rad), Hue 2023 table 1. The few
    // noninteger fitted frequencies are blended across the longitude seam.
    static double lead_angle(int moon, bool south, double west_longitude);
    // False outside the baked domain; never extrapolate a magnetic trace.
    bool footprint(int moon, bool south, const glm::dvec3& moon_body_km, glm::dvec3& direction) const;
};

} // namespace astraxis
