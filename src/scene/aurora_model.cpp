#include "scene/aurora_model.hpp"

#include "core/color.hpp"
#include "scene/mars_aurora.hpp"
#include "scene/satellite_aurora.hpp"
#include "scene/saturn_aurora.hpp"

namespace astraxis {

namespace {

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

AuroraProfile moon_profile(const SatelliteVisibleAurora& model, AuroraPattern pattern)
{
    AuroraProfile p;
    p.pattern = pattern;
    p.vertical = AuroraVertical::Exponential;
    p.lighting = AuroraLighting::NightOrEclipsed;
    p.top_km = model.top_km;
    p.layers = glm::dvec4(model.scale_km, model.top_km, model.sodium_scale_km, 0.0);
    p.low_rgb = model.column_rgb * model.display_gain;
    p.high_rgb = model.sodium_column_rgb * model.display_gain;
    return p;
}

} // namespace

bool aurora_shown(AuroraModel model, bool visible, bool ultraviolet)
{
    switch (model) {
    case AuroraModel::None:
        return false;
    case AuroraModel::Jupiter:
    case AuroraModel::Ganymede:
        return ultraviolet;
    default:
        return visible;
    }
}

AuroraProfile aurora_profile(AuroraModel model, bool ultraviolet, double kp)
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
        p.pattern = AuroraPattern::GanymedeMap;
        p.vertical = AuroraVertical::UniformColumn;
        p.top_km = GanymedeAurora::kTopKm;
        // Display exposure per 100 R, independent of Jupiter's much brighter UV.
        p.low_rgb = uv_false_color(0.2);
        p.layers.x = GanymedeAurora::kTopKm;
        p.map_gain = GanymedeAurora::kTextureRange / 100.0;
        break;
    case AuroraModel::Europa:
        p = moon_profile(europa_visible_aurora(), AuroraPattern::EuropaTrailing);
        break;
    case AuroraModel::Io:
        p = moon_profile(io_visible_aurora(), AuroraPattern::IoFlanks);
        break;
    case AuroraModel::Callisto:
        p = moon_profile(callisto_visible_aurora(), AuroraPattern::Uniform); // unresolved glow
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
