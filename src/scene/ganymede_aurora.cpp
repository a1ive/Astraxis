#include "scene/ganymede_aurora.hpp"

#include "core/math.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <locale>

namespace astraxis {

bool GanymedeAurora::load(const std::filesystem::path& file, std::string* error)
{
    std::ifstream input(file);
    input.imbue(std::locale::classic());
    std::string magic, extra;
    int width = 0, height = 0;
    std::array<double, kWidth * kHeight> values{};
    bool valid = bool(input >> magic >> width >> height) && magic == "AXGANYMEDE1" &&
                 width == kWidth && height == kHeight;
    for (double& v : values) {
        if (!valid) break;
        valid = bool(input >> v) && std::isfinite(v) && v >= 0.0 && v <= kTextureRange;
    }
    if (!valid || (input >> extra)) {
        if (error) *error = "invalid Ganymede aurora map: " + file.string();
        return false;
    }
    rayleigh = values;
    return true;
}

double GanymedeAurora::brightness(double east_longitude, double latitude) const
{
    if (!std::isfinite(east_longitude) || !std::isfinite(latitude)) return 0.0;
    const double west_turns = -east_longitude / kTwoPi;
    const double x = (west_turns - std::floor(west_turns)) * kWidth - 0.5;
    const double y = std::clamp((latitude / kPi + 0.5) * kHeight - 0.5, 0.0, double(kHeight - 1));
    const int i = (static_cast<int>(std::floor(x)) + kWidth) % kWidth;
    const int j = static_cast<int>(y);
    auto at = [&](int row, int col) { return rayleigh[static_cast<size_t>(row * kWidth + col)]; };
    const double a = std::lerp(at(j, i), at(j, (i + 1) % kWidth), x - std::floor(x));
    const int next = std::min(j + 1, kHeight - 1);
    const double b = std::lerp(at(next, i), at(next, (i + 1) % kWidth), x - std::floor(x));
    return std::lerp(a, b, y - j);
}

std::array<uint8_t, GanymedeAurora::kWidth * GanymedeAurora::kHeight * 4> GanymedeAurora::texture_rgba() const
{
    std::array<uint8_t, kWidth * kHeight * 4> pixels{};
    for (size_t i = 0; i < rayleigh.size(); ++i) {
        const auto v = static_cast<uint8_t>(std::lround(std::clamp(rayleigh[i] / kTextureRange, 0.0, 1.0) * 255.0));
        pixels[4 * i] = pixels[4 * i + 1] = pixels[4 * i + 2] = v;
        pixels[4 * i + 3] = 255;
    }
    return pixels;
}

} // namespace astraxis
