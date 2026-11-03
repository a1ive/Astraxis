#include "scene/jovian_aurora.hpp"

#include "core/math.hpp"

#include <glm/geometric.hpp>

#include <algorithm>
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

bool JovianAurora::load_footprints(const std::filesystem::path& file, std::string* error)
{
    std::ifstream input(file);
    input.imbue(std::locale::classic());
    std::string magic;
    int lon = 0, radial = 0, latitude = 0;
    bool valid = bool(input >> magic >> lon >> radial >> latitude) && magic == "AXJOVIANFP1" &&
                 lon == kLongitudes && radial == 3 && latitude == 3;
    std::array<Footpath, 3> paths;
    for (size_t m = 0; m < paths.size() && valid; ++m) {
        auto& path = paths[m];
        std::string name;
        valid = bool(input >> name >> path.min_radius >> path.max_radius >>
                     path.min_latitude >> path.max_latitude) && name == kMoons[m] &&
                path.min_radius > 1.0 && path.max_radius < 20.0 && path.min_radius < path.max_radius &&
                path.min_latitude >= -5.0 && path.max_latitude <= 5.0 &&
                path.min_latitude < path.max_latitude;
        path.min_latitude *= kDegToRad;
        path.max_latitude *= kDegToRad;
        if (!valid) break;
        path.points.resize(9 * kLongitudes);
        for (auto& pair : path.points) {
            for (int h = 0; h < 2; ++h) {
                auto& p = pair[static_cast<size_t>(h)];
                valid = valid && bool(input >> p.x >> p.y >> p.z) &&
                        std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z) &&
                        std::abs(glm::length(p) - 1.0) < 1e-6 && (h == 0 ? p.z > 0.5 : p.z < -0.5);
            }
            if (!valid) break;
        }
    }
    std::string extra;
    if (!valid || (input >> extra)) {
        if (error) *error = "invalid Jovian footprint table: " + file.string();
        return false;
    }
    footpaths = std::move(paths);
    return true;
}

double JovianAurora::lead_angle(int moon, bool south, double west_longitude)
{
    if (moon < 0 || moon >= 3 || !std::isfinite(west_longitude)) return 0.0;
    const auto fit = [moon, south](double x) {
        if (moon == 0) {
            return south ? 4.14 - 1.89*std::cos(x) - .70*std::sin(x) + .22*std::cos(2*x) - .12*std::sin(2*x)
                         : 4.26 + 2.64*std::cos(x) + .50*std::sin(x) + .20*std::cos(2*x) + .126*std::sin(2*x);
        }
        if (moon == 1) return south ? 7.06 - 3.37*std::cos(x) + .15*std::sin(x)
                                   : 5.78 + 3.29*std::cos(.99*x) + 1.41*std::sin(.99*x);
        return south ? 12.97 - 5.39*std::cos(.99*x) - 4.50*std::sin(.99*x)
                     : 12.45 + 6.16*std::cos(.98*x) + 3.64*std::sin(.98*x);
    };
    double x = west_longitude - std::floor(west_longitude / kTwoPi) * kTwoPi;
    constexpr double width = 15.0 * kDegToRad;
    if (x < width || x > kTwoPi - width) {
        if (x > kTwoPi - width) x -= kTwoPi;
        const double u = (x + width) / (2.0 * width);
        const double w = u*u*(3.0 - 2.0*u);
        return ((1.0 - w) * fit(x + kTwoPi) + w * fit(x)) * kDegToRad;
    }
    return fit(x) * kDegToRad;
}

bool JovianAurora::footprint(int moon, bool south, const glm::dvec3& p, glm::dvec3& direction) const
{
    direction = glm::dvec3(0.0);
    if (moon < 0 || moon >= 3) return false;
    const auto& path = footpaths[static_cast<size_t>(moon)];
    const double r = glm::length(p) / 71492.0;
    if (path.points.empty() || !std::isfinite(r) || r < path.min_radius || r > path.max_radius) return false;
    const double lat = std::atan2(p.z, std::hypot(p.x, p.y));
    if (lat < path.min_latitude || lat > path.max_latitude) return false;
    const double east = std::atan2(p.y, p.x);
    // The co-rotating plasma overtakes the moon: the emitting tube leads
    // eastward, i.e. toward smaller System III WEST longitude.
    const double turns = (east + lead_angle(moon, south, -east)) / kTwoPi;
    const double x = (turns - std::floor(turns)) * kLongitudes;
    const double y = 2.0 * (lat - path.min_latitude) / (path.max_latitude - path.min_latitude);
    const double z = 2.0 * (r - path.min_radius) / (path.max_radius - path.min_radius);
    const int ix = static_cast<int>(std::floor(x)) % kLongitudes;
    const int iy = std::min(static_cast<int>(std::floor(y)), 1);
    const int iz = std::min(static_cast<int>(std::floor(z)), 1);
    const glm::dvec3 f(x - std::floor(x), y - iy, z - iz);
    for (int dz = 0; dz < 2; ++dz) {
        for (int dy = 0; dy < 2; ++dy) {
            for (int dx = 0; dx < 2; ++dx) {
                const size_t index = static_cast<size_t>(((iz+dz)*3+iy+dy)*kLongitudes + (ix+dx)%kLongitudes);
                const double w = (dx ? f.x : 1-f.x) * (dy ? f.y : 1-f.y) * (dz ? f.z : 1-f.z);
                direction += w * path.points[index][south ? 1 : 0];
            }
        }
    }
    direction = glm::normalize(direction);
    return true;
}

} // namespace astraxis
