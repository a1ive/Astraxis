#include "scene/aurora_map.hpp"

#include "core/math.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <locale>

namespace astraxis {

bool AuroraMap::load(const std::filesystem::path& file, std::string* error)
{
    std::ifstream input(file);
    input.imbue(std::locale::classic());
    std::string magic, extra;
    int w = 0, h = 0;
    double max_value = 0.0;
    bool valid = bool(input >> magic >> w >> h >> max_value) && magic == "AXAURORAMAP1" && w > 0 && h > 0 &&
                 w <= 4096 && h <= 2048 && std::isfinite(max_value) && max_value > 0.0;
    std::vector<double> cells(valid ? static_cast<size_t>(w) * static_cast<size_t>(h) : 0u);
    for (double& v : cells) {
        if (!valid) break;
        valid = bool(input >> v) && std::isfinite(v) && v >= 0.0 && v <= max_value;
    }
    if (!valid || (input >> extra)) {
        if (error) *error = "invalid aurora map: " + file.string();
        return false;
    }
    width = w;
    height = h;
    range = max_value;
    values = std::move(cells);
    return true;
}

double AuroraMap::brightness(double east_longitude, double latitude) const
{
    if (!std::isfinite(east_longitude) || !std::isfinite(latitude) || values.empty()) return 0.0;
    const double west_turns = -east_longitude / kTwoPi;
    const double x = (west_turns - std::floor(west_turns)) * width - 0.5;
    const double y = std::clamp((latitude / kPi + 0.5) * height - 0.5, 0.0, double(height - 1));
    const int i = (static_cast<int>(std::floor(x)) + width) % width;
    const int j = static_cast<int>(y);
    auto at = [&](int row, int col) { return values[static_cast<size_t>(row * width + col)]; };
    const double a = std::lerp(at(j, i), at(j, (i + 1) % width), x - std::floor(x));
    const int next = std::min(j + 1, height - 1);
    const double b = std::lerp(at(next, i), at(next, (i + 1) % width), x - std::floor(x));
    return std::lerp(a, b, y - j);
}

std::vector<uint8_t> AuroraMap::texture_rgba() const
{
    std::vector<uint8_t> pixels(values.size() * 4);
    for (size_t i = 0; i < values.size(); ++i) {
        const auto v = static_cast<uint8_t>(std::lround(std::clamp(values[i] / range, 0.0, 1.0) * 255.0));
        pixels[4 * i] = pixels[4 * i + 1] = pixels[4 * i + 2] = v;
        pixels[4 * i + 3] = 255;
    }
    return pixels;
}

} // namespace astraxis
