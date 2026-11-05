#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// A brightness map fixed to a body, in the body's frame: south-to-north rows of
// west longitude cell centers, in the units of its model (scene/aurora_model.hpp).
// Baked by tools/geomag/ganymede_aurora_bake.py and ice_giant_h3p_bake.py.
struct AuroraMap {
    int width = 0;
    int height = 0;
    double range = 0.0; // the largest value the texture holds
    std::vector<double> values;

    bool load(const std::filesystem::path& file, std::string* error = nullptr);
    double brightness(double east_longitude, double latitude) const; // radians
    std::vector<uint8_t> texture_rgba() const;
};

} // namespace astraxis
