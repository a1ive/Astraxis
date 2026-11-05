#pragma once

#include "scene/aurora.hpp"
#include "scene/aurora_map.hpp"
#include "scene/jovian_aurora.hpp"

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <array>
#include <memory>

namespace astraxis {

// A body's aurora: at most one model, chosen by its [aurora] table. Only the
// models with data of their own use the pointers.
enum class AuroraModel { None, Earth, Jupiter, Saturn, Ganymede, Europa, Io, Callisto, Mars, Uranus, Neptune };

struct BodyAurora {
    AuroraModel model = AuroraModel::None;
    std::shared_ptr<const KpRecord> kp;             // Earth
    std::shared_ptr<const JovianAurora> jovian;     // Jupiter
    std::shared_ptr<const AuroraMap> map;           // Ganymede, Uranus, Neptune
};

// Whether a model draws under the current view: Jupiter's and Ganymede's data
// are ultraviolet, Uranus's and Neptune's infrared, the others also have visible
// light. Each view's flag stands alone: pass visible with either false color.
bool aurora_shown(AuroraModel model, bool visible, bool ultraviolet, bool infrared);

// What shaders/aurora.frag.hlsl does for a shell, as three independent choices
// (the values are its PATTERN_, VERTICAL_ and LIGHTING_ constants).
enum class AuroraPattern { // brightness across the surface
    EarthOval = 0,         // oval, diffuse aurora and arcs at magnetic local time
    JovianOval = 1,        // the statistical contour, per hemisphere
    JovianFootprints = 2,  // spots at the moons' footprints
    SaturnOval = 3,        // the statistical oval at solar local time
    BodyMapShares = 4,     // a brightness map fixed to the body (AuroraMap), its four layers
                           // in the shares map_sun
    Uniform = 5,
    EuropaTrailing = 6,    // assumed smooth trailing enhancement
    IoFlanks = 7,          // broad equatorial flanks
    BodyMapSeasons = 8,    // AuroraMap's two layers mixed by the sun at the latitude's noon:
                           // map_sun = zenith angles (rad) where the sunlit layer's share is
                           // 1, the terminator's and 0, the terminator's share
    BodyMapDusk = 9,       // an AuroraMap plus a gradient fixed to the sun: map_sun.x, in
                           // map units, from -1 at the dawn limb to +1 at the dusk limb
};
enum class AuroraVertical {  // brightness with height; AuroraProfile::layers reads, for each:
    EarthChapman = 0,        // low peak, low scale, high peak, high scale (km); the high
                             // layer's color is per kR, by the local red-to-blue ratio
    ChapmanPair = 1,         // low peak, low scale, high peak, high scale (km)
    Exponential = 2,         // low scale, layer top, high scale (km); the high layer
                             // ignores the pattern (Io's sodium corona)
    Gaussian = 3,            // peak, sigma, bottom, top (km), one layer
    UniformColumn = 4,       // thickness (km), one layer
    Table = 5,               // first height, step (km), one layer: AuroraProfile::table
};
enum class AuroraLighting { // when it shows
    Always = 0,
    NightSide = 1,          // fades out as the sun rises at the point below
    NightOrEclipsed = 2,    // also in the shadow of the body's parent
};

inline constexpr int kAuroraTableSize = 16;

// What a moon's aurora depends on right now (Simulation::moon_aurora). The defaults
// reproduce the reference observation.
struct MoonAuroraState {
    // Ganymede: the shares of the long-term map, and of it above, inside and below Jupiter's
    // current sheet (AuroraMap's layers).
    glm::dvec4 map_shares{1.0, 0.0, 0.0, 0.0};
};

// A shell of emission, in linear RGB with its display exposure.
struct AuroraProfile {
    AuroraPattern pattern = AuroraPattern::Uniform;
    AuroraVertical vertical = AuroraVertical::ChapmanPair;
    AuroraLighting lighting = AuroraLighting::Always;
    double bottom_km = 0.0; // above the equatorial radius
    double top_km = 0.0;    // above the polar radius
    glm::dvec3 low_rgb{0.0};
    glm::dvec3 high_rgb{0.0};
    glm::dvec4 layers{0.0}; // see AuroraVertical
    double blue_kr = 0.0;       // EarthChapman: 427.8 nm of a column through the oval
    double diffuse_flux = 0.0;  // EarthOval: the diffuse aurora's flux relative to the oval's
    double map_gain = 0.0;      // BodyMap: low_rgb's fraction per unit of the map
    glm::dvec4 map_sun{0.0};    // see AuroraPattern
    double sigma_km = 0.0;      // JovianOval, JovianFootprints: Gaussian width across the contour
    // At 24 k / kOvalSamples hours. EarthOval: colatitudes (rad) of the diffuse
    // aurora, the oval's equatorward and poleward edges. SaturnOval: SaturnAurora::profile.
    std::array<glm::vec4, AuroraEmission::kOvalSamples> oval{};
    // Table: the column's share per km at layers.x + k * layers.y, linear between
    // them and zero beyond.
    std::array<float, kAuroraTableSize> table{};
};

// The model's shell. ultraviolet: Saturn's false-color view instead of its
// visible one. kp: Earth's activity. moon: the moons' state.
AuroraProfile aurora_profile(AuroraModel model, bool ultraviolet, double kp, const MoonAuroraState& moon = {});

// Jupiter's second shell: its moons' footprints at the height they are mapped
// to, a display gain per moon (Io, Europa, Ganymede).
AuroraProfile jovian_footprint_profile();

// The pattern of the map-based models at a point of the body, as map_pattern and
// its callers in shaders/aurora.frag.hlsl: to_sun is a unit vector in the body's frame.
double aurora_map_pattern(const AuroraProfile& profile, const AuroraMap& map, double east_longitude,
                          double latitude, const glm::dvec3& to_sun);

// How much of a body's reflected light its false-color view keeps. Uranus and
// Neptune are dark at the 3-4 um of the H3+ lines (methane absorbs), and their
// low-contrast glow covers the whole disk: on the visible planet it would vanish.
double aurora_surface_gain(AuroraModel model, bool infrared);
inline constexpr std::array<double, 3> kJovianFootprintGains{6.0, 2.0, 4.0};

} // namespace astraxis
