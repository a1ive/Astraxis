#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// A brightness map fixed to a body, in the body's frame: south-to-north rows of
// west longitude cell centers, in the units of its model (scene/aurora_model.hpp).
// Baked by tools/geomag/ganymede_aurora_bake.py and ice_giant_h3p_bake.py. The header
// AXAURORAMAPn has n layers one after the other: 1, or 2 (the second is what is left without
// sunlight, Uranus: AuroraPattern::BodyMapSeasons), or 4 (the long-term map, then Ganymede
// above, inside and below Jupiter's current sheet: AuroraPattern::BodyMapShares).
struct AuroraMap {
    int width = 0;
    int height = 0;
    double range = 0.0; // the largest value the texture holds
    std::vector<double> values;
    std::vector<std::vector<double>> more_layers; // after the first

    bool load(const std::filesystem::path& file, std::string* error = nullptr);
    double brightness(double east_longitude, double latitude) const; // radians
    double layer_brightness(size_t layer, double east_longitude, double latitude) const; // layer 0 is values
    double dark_brightness(double east_longitude, double latitude) const { return layer_brightness(1, east_longitude, latitude); }
    // The layers in red, green, blue and alpha. With one layer it is in the color channels,
    // with two the blue channel repeats the first, and alpha is opaque below four.
    std::vector<uint8_t> texture_rgba() const;
};

} // namespace astraxis
