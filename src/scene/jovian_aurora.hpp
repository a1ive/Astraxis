#pragma once

#include <glm/mat3x3.hpp>
#include <glm/vec2.hpp>
#include <glm/vec4.hpp>

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
    // The spots of a footprint (Bonfond et al. 2008, GRL 35, L05107; 2009, JGR 114, A07224).
    // Main: where the Alfven wing the moon launched toward this hemisphere arrives (the MAW
    // spot). Beam: electrons the other hemisphere's wing sent back along the field cross the
    // plasma sheet in seconds, so they land on the field line of that wing's foot (the TEB
    // spot): ahead of the main spot in the hemisphere the moon is farther from, behind it in
    // the other. Reflected: the wing that left toward the other hemisphere, turned back at
    // the sheet's far edge (the RAW spot; seen at Io).
    enum class Spot { Main, Beam, Reflected };
    // The eastward equatorial lead (rad) of a spot: lead_angle of this hemisphere, of the
    // other one, or of this one plus twice the other's (the reflected wave crosses the other
    // half of the sheet twice, and nearly all of the delay is in the sheet).
    static double spot_lead_angle(int moon, bool south, double west_longitude, Spot spot);
    // False outside the baked domain; never extrapolate a magnetic trace.
    bool footprint(int moon, bool south, const glm::dvec3& moon_body_km, glm::dvec3& direction, Spot spot = Spot::Main) const;
    // Where the field line through the moon lands now, with no lead: not the bright spot.
    bool instant_footprint(int moon, bool south, const glm::dvec3& moon_body_km, glm::dvec3& direction) const;

    // The main spot and the downstream wake along its mapped track. Each point contains
    // a unit direction and brightness relative to the main spot. False outside the table.
    static constexpr int kTailPoints = 20;
    using Tail = std::array<glm::dvec4, kTailPoints>;
    bool footprint_tail(int moon, bool south, const glm::dvec3& moon_body_km, Tail& tail) const;

    // The lines from the moons to the 900 km spheroid, same field as the footprints, at the
    // moons' nominal radii: by moon, then longitude, hemisphere (north, south) and point. xyz
    // in 71492 km, turned by minus the line's longitude; w = the part of the Alfven wave's
    // travel time spent by that point. See geomag/SOURCES.md.
    static constexpr int kTubeLongitudes = 72;
    static constexpr int kTubePoints = 40;
    std::array<std::vector<glm::vec4>, 3> flux_tubes;
    // The same field (nT) at each moon's nominal radius in the equatorial plane, by longitude,
    // turned like the lines.
    std::array<std::vector<glm::vec3>, 3> moon_fields;
    bool load_flux_tubes(const std::filesystem::path& file, std::string* error);
    // Jupiter's field around a moon (nT, in Jupiter's axes), taken at the moon's longitude on
    // its nominal orbit: what the moon's own field stands in (Ganymede). False without the table.
    bool field_at_moon(int moon, const glm::dvec3& moon_body_km, glm::dvec3& field_nt) const;
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
