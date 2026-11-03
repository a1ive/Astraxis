#include "scene/jovian_aurora.hpp"

#include "core/math.hpp"

#include <cmath>
#include <fstream>
#include <locale>

namespace astraxis {

bool JovianAurora::load(const std::filesystem::path& file, std::string* error)
{
    std::ifstream input(file);
    input.imbue(std::locale::classic());
    std::string magic;
    int count = 0;
    std::array<glm::dvec2, kSamples> values{};
    bool valid = bool(input >> magic >> count) && magic == "AXJOVIAN1" && count == kSamples;
    for (auto& value : values) {
        if (!valid) break;
        valid = bool(input >> value.x >> value.y);
        for (int h = 0; h < 2; ++h) {
            valid = valid && std::isfinite(value[h]) && value[h] > 0.01 && value[h] < 0.7;
        }
    }
    std::string extra;
    if (!valid || (input >> extra)) {
        if (error) *error = "invalid Jovian reference oval: " + file.string();
        return false;
    }
    radii = values;
    return true;
}

double JovianAurora::radius(bool south, double azimuth) const
{
    const double turns = azimuth / kTwoPi;
    const double x = (turns - std::floor(turns)) * kSamples;
    const int i = static_cast<int>(x) % kSamples;
    const int h = south ? 1 : 0;
    return radii[static_cast<size_t>(i)][h] + (x - std::floor(x)) *
        (radii[static_cast<size_t>((i + 1) % kSamples)][h] - radii[static_cast<size_t>(i)][h]);
}

glm::dmat3 JovianAurora::frame(bool south)
{
    // Head 2024 section 4.1 gives WEST longitude; body y points east.
    const double lon = -(south ? 32.0 : 185.0) * kDegToRad;
    const double lat = (south ? -82.0 : 74.0) * kDegToRad;
    return {{-std::sin(lon), std::cos(lon), 0.0},
            {-std::sin(lat) * std::cos(lon), -std::sin(lat) * std::sin(lon), std::cos(lat)},
            {std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat)}};
}

} // namespace astraxis
