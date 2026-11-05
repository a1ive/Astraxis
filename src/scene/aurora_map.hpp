#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// A brightness map fixed to a body, in the body's frame: south-to-north rows of
// west longitude cell centers, in the units of its model (scene/aurora_model.hpp).
// Baked by tools/geomag/ganymede_aurora_bake.py and ice_giant_h3p_bake.py. The header
// AXAURORAMAP1 has one layer; AXAURORAMAP2 has a second after it, what is left without
// sunlight (Uranus: AuroraPattern::BodyMapSeasons).
struct AuroraMap {
    int width = 0;
    int height = 0;
    double range = 0.0; // the largest value the texture holds
    std::vector<double> values;
    std::vector<double> dark; // empty with one layer

    bool load(const std::filesystem::path& file, std::string* error = nullptr);
    double brightness(double east_longitude, double latitude) const; // radians
    double dark_brightness(double east_longitude, double latitude) const;
    // The layers in red and green; with one layer it is in all the channels.
    std::vector<uint8_t> texture_rgba() const;
};

} // namespace astraxis
