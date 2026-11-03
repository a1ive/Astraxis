#pragma once

#include <glm/mat3x3.hpp>
#include <glm/vec2.hpp>

#include <array>
#include <filesystem>
#include <string>

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
};

} // namespace astraxis
