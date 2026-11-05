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
    // Where the field line through the moon lands now, with no lead: not the bright spot.
    bool instant_footprint(int moon, bool south, const glm::dvec3& moon_body_km, glm::dvec3& direction) const;

    // The lines from the moons to the 900 km spheroid, same field as the footprints, at the
    // moons' nominal radii: by moon, then longitude, hemisphere (north, south) and point. xyz
    // in 71492 km, turned by minus the line's longitude; w = the part of the Alfven wave's
    // travel time spent by that point. See geomag/SOURCES.md.
    static constexpr int kTubeLongitudes = 72;
    static constexpr int kTubePoints = 40;
    std::array<std::vector<glm::vec4>, 3> flux_tubes;
    bool load_flux_tubes(const std::filesystem::path& file, std::string* error);
    // A line from the moon to one hemisphere, in km in the body's axes. alfven_wing false: the
    // field line through the moon now, which ends at instant_footprint. True: where the wave
    // the moon launched is now. Each part of it left the moon when the moon was farther east
    // along its way through the field, so it runs across the field lines and ends at the
    // bright spot, footprint(). False where footprint() is.
    bool flux_tube(int moon, bool south, const glm::dvec3& moon_body_km, bool alfven_wing,
                   std::vector<glm::dvec3>& line) const;

private:
    bool footprint_at(int moon, bool south, const glm::dvec3& moon_body_km, double east, glm::dvec3& direction) const;
    glm::dvec4 tube_point(int moon, bool south, double east, int point) const;
};

} // namespace astraxis
