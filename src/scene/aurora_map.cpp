#include "scene/aurora_map.hpp"

#include "core/math.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <locale>

namespace astraxis {

namespace {

double sample_layer(const std::vector<double>& layer, int width, int height, double east_longitude, double latitude)
{
    if (!std::isfinite(east_longitude) || !std::isfinite(latitude) || layer.empty()) return 0.0;
    const double west_turns = -east_longitude / kTwoPi;
    const double x = (west_turns - std::floor(west_turns)) * width - 0.5;
    const double y = std::clamp((latitude / kPi + 0.5) * height - 0.5, 0.0, double(height - 1));
    const int i = (static_cast<int>(std::floor(x)) + width) % width;
    const int j = static_cast<int>(y);
    auto at = [&](int row, int col) { return layer[static_cast<size_t>(row * width + col)]; };
    const double a = std::lerp(at(j, i), at(j, (i + 1) % width), x - std::floor(x));
    const int next = std::min(j + 1, height - 1);
    const double b = std::lerp(at(next, i), at(next, (i + 1) % width), x - std::floor(x));
    return std::lerp(a, b, y - j);
}

} // namespace

bool AuroraMap::load(const std::filesystem::path& file, std::string* error)
{
    std::ifstream input(file);
    input.imbue(std::locale::classic());
    std::string magic, extra;
    int w = 0, h = 0;
    double max_value = 0.0;
    bool valid = bool(input >> magic >> w >> h >> max_value) &&
                 (magic == "AXAURORAMAP1" || magic == "AXAURORAMAP2" || magic == "AXAURORAMAP4") &&
                 w > 0 && h > 0 && w <= 4096 && h <= 2048 && std::isfinite(max_value) && max_value > 0.0;
    const size_t count = valid ? static_cast<size_t>(w) * static_cast<size_t>(h) : 0u;
    const auto read_layer = [&](std::vector<double>& layer) {
        layer.resize(count);
        for (double& v : layer) {
            if (!valid) break;
            valid = bool(input >> v) && std::isfinite(v) && v >= 0.0 && v <= max_value;
        }
    };
    std::vector<double> cells;
    read_layer(cells);
    std::vector<std::vector<double>> more(valid ? static_cast<size_t>(magic.back() - '1') : 0u);
    for (auto& layer : more) read_layer(layer);
    if (!valid || (input >> extra)) {
        if (error) *error = "invalid aurora map: " + file.string();
        return false;
    }
    width = w;
    height = h;
    range = max_value;
    values = std::move(cells);
    more_layers = std::move(more);
    return true;
}

double AuroraMap::brightness(double east_longitude, double latitude) const
{
    return sample_layer(values, width, height, east_longitude, latitude);
}

double AuroraMap::layer_brightness(size_t layer, double east_longitude, double latitude) const
{
    if (layer == 0) return brightness(east_longitude, latitude);
    if (layer > more_layers.size()) return 0.0;
    return sample_layer(more_layers[layer - 1], width, height, east_longitude, latitude);
}

std::vector<uint8_t> AuroraMap::texture_rgba() const
{
    const auto byte = [&](double v) { return static_cast<uint8_t>(std::lround(std::clamp(v / range, 0.0, 1.0) * 255.0)); };
    std::vector<uint8_t> pixels(values.size() * 4);
    for (size_t i = 0; i < values.size(); ++i) {
        pixels[4 * i] = byte(values[i]);
        for (size_t channel = 1; channel < 3; ++channel) {
            pixels[4 * i + channel] = channel <= more_layers.size() ? byte(more_layers[channel - 1][i]) : pixels[4 * i];
        }
        pixels[4 * i + 3] = more_layers.size() >= 3 ? byte(more_layers[2][i]) : uint8_t{255};
    }
    return pixels;
}

} // namespace astraxis
