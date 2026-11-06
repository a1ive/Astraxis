#include "scene/aurora_model.hpp"

#include "core/color.hpp"
#include "core/math.hpp"
#include "scene/mars_aurora.hpp"
#include "scene/satellite_aurora.hpp"
#include "scene/saturn_aurora.hpp"

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr double kGanymedeTopKm = 70.0; // Marzok 2022 section 2.3, uniform column

// The aurora is ~1e-8 of the sunlit ground (I/F; 1 kR of the green line is
// 5e-9): boosted like a long exposure, with the night lights' gain (scene
// files) chosen to match. Its shell (km above the equatorial radius) holds
// both emission layers.
constexpr double kEarthExposure = 2.5e6;
constexpr double kEarthBottomKm = 80.0;
constexpr double kEarthTopKm = 450.0;

// The planets' ultraviolet false color; the exposures are arbitrary.
glm::dvec3 uv_false_color(double exposure)
{
    return srgb_to_linear(glm::dvec3(0.3, 0.55, 1.0)) * exposure;
}

// The infrared one: warm, to stand apart from the UV blue and from the ice
// giants' own color.
glm::dvec3 ir_false_color(double exposure)
{
    return srgb_to_linear(glm::dvec3(1.0, 0.6, 0.15)) * exposure;
}

// Melin et al. 2025's maps are of H3+ total emission, uW m-2 sr-1 (Uranus); the
// display exposure of 1 uW m-2 sr-1.
constexpr double kUranusExposure = 0.28;
// Neptune's map is relative to the disk's background; the background's exposure.
constexpr double kNeptuneExposure = 0.095;
// Melin et al. 2025 (Neptune) sees the H3+ rise from dawn to dusk across the disk, likely
// the build-up of photoionization over the day; they give no number. Fitted to their
// spaxels, relative to the background (tools/geomag/ice_giant_h3p_bake.py prints it).
constexpr double kNeptuneDuskGradient = 0.228;

// How much of Uranus's sunlit emission a latitude keeps, by the sun's zenith angle at its
// noon. Non-auroral H3+ comes from solar photoionization and lives about 1 h (Melin et al.
// 2025, after Moore et al. 2019), far less than the 17.2 h rotation, so it follows the sun.
// Fully lit at the zenith angle of the map's last trusted latitude when it was observed
// (8 N, sub-solar 58.83 N, the scene at Melin's dates, test_ice_giant_aurora).
constexpr double kUranusFullSunZenithDeg = 58.83 - 8.0;
// At the terminator, the column near it (Tiranti et al. 2026, 6.6e14 m-2) over the sunlit
// disk's median (Melin et al. 2025, 2.27e15 m-2).
constexpr double kUranusTerminatorShare = 6.6e14 / 2.27e15;
// None once the sun is below the horizon of the altitude of the primary photoelectron
// heating near the terminator (4000 km, Tiranti et al. 2026), above the 1 bar radius.
constexpr double kUranusPhotoionizationKm = 4000.0;
constexpr double kUranusEquatorialKm = 25559.0;
// Tiranti et al. 2026: the mean of their 36 longitudes of H3+ total emission at
// 475 .. 5025 km, every 350 km (tools/geomag/ice_giant_h3p_bake.py prints it).
// One step below it, zero: nothing observed lower, 300 km is the homopause. One
// step above it, zero: the tail is unreliable higher up.
constexpr double kUranusTableBottomKm = 125.0;
constexpr double kUranusTableStepKm = 350.0;
constexpr std::array<float, kAuroraTableSize> kUranusTable{
    0.0f, 2.9405e-04f, 3.2810e-04f, 3.7366e-04f, 3.7617e-04f, 3.5086e-04f, 3.0463e-04f, 2.4104e-04f,
    1.6730e-04f, 1.1350e-04f, 8.6082e-05f, 7.1455e-05f, 5.9887e-05f, 4.9504e-05f, 4.0907e-05f, 0.0f};

AuroraProfile moon_profile(const SatelliteVisibleAurora& model, AuroraPattern pattern, const MoonAuroraState& moon)
{
    AuroraProfile p;
    p.pattern = pattern;
    p.vertical = AuroraVertical::Exponential;
    p.lighting = AuroraLighting::NightOrEclipsed;
    p.top_km = model.top_km;
    p.layers = glm::dvec4(model.scale_km, model.top_km, model.sodium_scale_km, 0.0);
    p.low_rgb = (model.column_rgb * moon.plasma + model.so2_column_rgb * moon.so2) * model.display_gain;
    p.high_rgb = model.sodium_column_rgb * moon.sodium * model.display_gain;
    return p;
}

} // namespace

bool aurora_shown(AuroraModel model, bool visible, bool ultraviolet, bool infrared)
{
    switch (model) {
    case AuroraModel::None:
        return false;
    case AuroraModel::Jupiter:
    case AuroraModel::Ganymede:
        return ultraviolet;
    case AuroraModel::Uranus:
    case AuroraModel::Neptune:
        return infrared;
    default:
        return visible;
    }
}

AuroraProfile aurora_profile(AuroraModel model, bool ultraviolet, double kp, const MoonAuroraState& moon)
{
    AuroraProfile p;
    switch (model) {
    case AuroraModel::None:
        break;
    case AuroraModel::Earth: {
        const AuroraEmission emission = aurora_emission(kp);
        p.pattern = AuroraPattern::EarthOval;
        p.vertical = AuroraVertical::EarthChapman;
        p.lighting = AuroraLighting::NightSide;
        p.bottom_km = kEarthBottomKm;
        p.top_km = kEarthTopKm;
        p.low_rgb = emission.low_rgb * kEarthExposure;
        p.high_rgb = kRedPerKr * kEarthExposure;
        p.layers = glm::dvec4(AuroraEmission::kLowLayerKm, AuroraEmission::kLowScaleKm,
                              AuroraEmission::kHighLayerKm, AuroraEmission::kHighScaleKm);
        p.blue_kr = emission.blue_kr;
        p.diffuse_flux = AuroraEmission::kDiffuseFlux;
        for (size_t k = 0; k < p.oval.size(); ++k) {
            p.oval[k] = glm::vec4(emission.oval[k], 0.0f);
        }
        break;
    }
    case AuroraModel::Jupiter:
        p.pattern = AuroraPattern::JovianOval;
        p.bottom_km = 150.0;
        p.top_km = 1000.0;
        // False color and display gain. The column brightness is illustrative.
        p.low_rgb = uv_false_color(2.0);
        // Head 2024 uses a 400 km projection altitude. Scale height is assumed.
        p.layers = glm::dvec4(400.0, 80.0, 400.0, 80.0);
        p.sigma_km = 200.0; // illustrative
        break;
    case AuroraModel::Saturn:
        p.pattern = AuroraPattern::SaturnOval;
        if (ultraviolet) {
            p.bottom_km = SaturnAurora::kBottomKm;
            p.top_km = SaturnAurora::kTopKm;
            // A 40 kR column gets this arbitrary false color / exposure.
            p.low_rgb = uv_false_color(2.0);
            p.layers = glm::dvec4(SaturnAurora::kPeakKm, SaturnAurora::kScaleKm,
                                  SaturnAurora::kPeakKm, SaturnAurora::kScaleKm);
        } else {
            // Dyudina 2015 Fig. 8: illustrative local colors and Chapman
            // layers for the observed pink lower / purple upper curtains.
            p.lighting = AuroraLighting::NightSide;
            p.bottom_km = 200.0;
            p.top_km = 2000.0;
            p.low_rgb = srgb_to_linear(glm::dvec3(1.0, 0.35, 0.55)) * 0.35;
            p.high_rgb = srgb_to_linear(glm::dvec3(0.65, 0.22, 1.0)) * 0.35;
            p.layers = glm::dvec4(450.0, 140.0, 1100.0, 180.0);
        }
        for (size_t k = 0; k < p.oval.size(); ++k) {
            p.oval[k] = glm::vec4(SaturnAurora::profile(24.0 * static_cast<double>(k) / static_cast<double>(p.oval.size())));
        }
        break;
    case AuroraModel::Ganymede:
        p.pattern = AuroraPattern::BodyMapShares;
        p.map_sun = moon.map_shares;
        p.vertical = AuroraVertical::UniformColumn;
        p.top_km = kGanymedeTopKm;
        // Display exposure per 100 R, independent of Jupiter's much brighter UV.
        p.low_rgb = uv_false_color(0.2);
        p.layers.x = kGanymedeTopKm;
        p.map_gain = 1.0 / 100.0;
        break;
    case AuroraModel::Uranus:
        p.pattern = AuroraPattern::BodyMapSeasons;
        p.map_sun = glm::dvec4(kUranusFullSunZenithDeg * kDegToRad, kPi / 2.0,
                               kPi / 2.0 + std::acos(kUranusEquatorialKm / (kUranusEquatorialKm + kUranusPhotoionizationKm)),
                               kUranusTerminatorShare);
        p.vertical = AuroraVertical::Table;
        p.bottom_km = kUranusTableBottomKm;
        p.top_km = kUranusTableBottomKm + kUranusTableStepKm * (kAuroraTableSize - 1);
        p.low_rgb = ir_false_color(kUranusExposure);
        p.layers = glm::dvec4(kUranusTableBottomKm, kUranusTableStepKm, 0.0, 0.0);
        p.map_gain = 1.0;
        p.table = kUranusTable;
        break;
    case AuroraModel::Neptune:
        // No vertical profile is observed. Assumed: Uranus's shape, as a Chapman
        // layer; models put the ion above 1000 km (Melin et al. 2025).
        p.pattern = AuroraPattern::BodyMapDusk;
        p.map_sun = glm::dvec4(kNeptuneDuskGradient, 0.0, 0.0, 0.0);
        p.bottom_km = 300.0;
        p.top_km = 5000.0;
        p.low_rgb = ir_false_color(kNeptuneExposure);
        p.layers = glm::dvec4(1500.0, 700.0, 1500.0, 700.0);
        p.map_gain = 1.0;
        break;
    case AuroraModel::Europa:
        p = moon_profile(europa_visible_aurora(), AuroraPattern::Europa, moon);
        p.low_rgb /= europa_pattern_mean(moon.field, moon.sheet);
        p.map_sun = glm::dvec4(moon.sheet * kEuropaSheet, moon.sheet * kEuropaUpstream, kEuropaFieldTilt, kEuropaAfternoon);
        break;
    case AuroraModel::Io:
        p = moon_profile(io_visible_aurora(), AuroraPattern::IoFlanks, moon);
        break;
    case AuroraModel::Callisto:
        p = moon_profile(callisto_visible_aurora(), AuroraPattern::Uniform, moon); // unresolved glow
        break;
    case AuroraModel::Mars:
        p.vertical = AuroraVertical::Gaussian;
        p.lighting = AuroraLighting::NightSide;
        p.bottom_km = MarsAurora::bottom_km;
        p.top_km = MarsAurora::top_km;
        p.layers = glm::dvec4(MarsAurora::peak_km, MarsAurora::sigma_km, MarsAurora::bottom_km, MarsAurora::top_km);
        p.low_rgb = kGreenPerKr * (MarsAurora::limb_rayleigh / 1000.0 * MarsAurora::display_gain / MarsAurora::limb_integral());
        break;
    }
    return p;
}

double aurora_map_pattern(const AuroraProfile& profile, const AuroraMap& map, double east_longitude,
                          double latitude, const glm::dvec3& to_sun)
{
    const double value = map.brightness(east_longitude, latitude) * profile.map_gain;
    switch (profile.pattern) {
    case AuroraPattern::BodyMapSeasons: {
        const glm::dvec4& s = profile.map_sun;
        const double zenith = std::abs(latitude - std::asin(std::clamp(to_sun.z, -1.0, 1.0)));
        const double sunlit = zenith <= s.y ? std::lerp(1.0, s.w, std::clamp((zenith - s.x) / (s.y - s.x), 0.0, 1.0))
                                            : s.w * std::clamp((s.z - zenith) / (s.z - s.y), 0.0, 1.0);
        return std::lerp(map.dark_brightness(east_longitude, latitude) * profile.map_gain, value, sunlit);
    }
    case AuroraPattern::BodyMapShares: {
        double sum = 0.0;
        for (int layer = 0; layer < 4; ++layer) {
            sum += profile.map_sun[layer] * map.layer_brightness(static_cast<size_t>(layer), east_longitude, latitude);
        }
        return sum * profile.map_gain;
    }
    case AuroraPattern::BodyMapDusk: {
        const double across = std::hypot(to_sun.x, to_sun.y);
        if (across < 1e-9) return value;
        const double dusk = std::cos(latitude) * (std::sin(east_longitude) * to_sun.x - std::cos(east_longitude) * to_sun.y) / across;
        return std::max(value + profile.map_sun.x * dusk, 0.0);
    }
    default:
        return value;
    }
}

double aurora_surface_gain(AuroraModel model, bool infrared)
{
    return infrared && (model == AuroraModel::Uranus || model == AuroraModel::Neptune) ? 0.03 : 1.0;
}

AuroraProfile jovian_footprint_profile()
{
    AuroraProfile p = aurora_profile(AuroraModel::Jupiter, true, 0.0);
    p.pattern = AuroraPattern::JovianFootprints;
    p.top_km = 3500.0;
    p.layers = glm::dvec4(900.0, 366.0, 900.0, 366.0);
    // Illustrative spot size and relative brightness, not photometry.
    p.sigma_km = 250.0;
    return p;
}

} // namespace astraxis
