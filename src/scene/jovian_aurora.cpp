#include "scene/jovian_aurora.hpp"

#include "core/math.hpp"
#include "scene/jovian_plasma_sheet.hpp"

#include <glm/common.hpp>
#include <glm/geometric.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>

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

double JovianAurora::spot_lead_angle(int moon, bool south, double west_longitude, Spot spot)
{
    const double here = lead_angle(moon, south, west_longitude);
    const double there = lead_angle(moon, !south, west_longitude);
    return spot == Spot::Main ? here : spot == Spot::Beam ? there : here + 2.0 * there;
}

bool JovianAurora::footprint(int moon, bool south, const glm::dvec3& p, glm::dvec3& direction, Spot spot) const
{
    // The co-rotating plasma overtakes the moon: the emitting tube leads
    // eastward, i.e. toward smaller System III WEST longitude.
    const double east = std::atan2(p.y, p.x);
    return footprint_at(moon, south, p, east + spot_lead_angle(moon, south, -east, spot), direction);
}

bool JovianAurora::field_at_moon(int moon, const glm::dvec3& p, glm::dvec3& field_nt) const
{
    field_nt = glm::dvec3(0.0);
    if (moon < 0 || moon >= 3 || moon_fields[static_cast<size_t>(moon)].empty()) return false;
    const auto& table = moon_fields[static_cast<size_t>(moon)];
    const double east = std::atan2(p.y, p.x);
    const double turns = east / kTwoPi;
    const double x = (turns - std::floor(turns)) * kTubeLongitudes;
    const int i = static_cast<int>(std::floor(x)) % kTubeLongitudes;
    const glm::dvec3 a(table[static_cast<size_t>(i)]);
    const glm::dvec3 b(table[static_cast<size_t>((i + 1) % kTubeLongitudes)]);
    // Interpolated in the moon's own frame, as the lines are.
    const glm::dvec3 f = a + (x - std::floor(x)) * (b - a);
    const double c = std::cos(east), s = std::sin(east);
    field_nt = {c * f.x - s * f.y, s * f.x + c * f.y, f.z};
    return true;
}

bool JovianAurora::instant_footprint(int moon, bool south, const glm::dvec3& p, glm::dvec3& direction) const
{
    return footprint_at(moon, south, p, std::atan2(p.y, p.x), direction);
}

bool JovianAurora::footprint_tail(int moon, bool south, const glm::dvec3& p, Tail& tail) const
{
    tail = {};
    glm::dvec3 head;
    if (!footprint(moon, south, p, head)) return false;

    // Bonfond et al. 2017, sections 2 and 4.2. The Ganymede observation gives a lower
    // bound on the visible extent. Use it as a nominal decay scale for the display.
    constexpr double kDecayDeg[] = {40.0, 14.0, 24.0};
    constexpr double kHeadGain[] = {0.32, 0.20, 0.22};
    double gain = kHeadGain[moon];
    if (moon != 0) {
        const double x = std::log10(glm::length(p) / 71492.0);
        // Bagenal & Delamere 2011, equation 6 and Table 1.
        const double scale = std::pow(10.0, -0.116 + x * (2.14 + x * (-2.05 + x * (0.491 + 0.126 * x))));
        gain *= plasma_sheet_density(centrifugal_equator_height_rj(p), scale);
    }
    const double east = std::atan2(p.y, p.x);
    const double start = east + lead_angle(moon, south, -east);
    tail.front() = glm::dvec4(head, gain);
    for (int i = 1; i < kTailPoints; ++i) {
        const double distance = 4.0 * i / (kTailPoints - 1);
        glm::dvec3 direction;
        footprint_at(moon, south, p, start + distance * kDecayDeg[moon] * kDegToRad, direction);
        const double fade = 1.0 - glm::smoothstep(3.0, 4.0, distance);
        tail[static_cast<size_t>(i)] = glm::dvec4(direction, gain * std::exp(-distance) * fade);
    }
    return true;
}

bool JovianAurora::footprint_at(int moon, bool south, const glm::dvec3& p, double east, glm::dvec3& direction) const
{
    direction = glm::dvec3(0.0);
    if (moon < 0 || moon >= 3) return false;
    const auto& path = footpaths[static_cast<size_t>(moon)];
    const double r = glm::length(p) / 71492.0;
    if (path.points.empty() || !std::isfinite(r) || r < path.min_radius || r > path.max_radius) return false;
    const double lat = std::atan2(p.z, std::hypot(p.x, p.y));
    if (lat < path.min_latitude || lat > path.max_latitude) return false;
    const double turns = east / kTwoPi;
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

bool JovianAurora::load_flux_tubes(const std::filesystem::path& file, std::string* error)
{
    const auto failed = [&](const char* what) {
        if (error) *error = "invalid Jovian flux tube table: " + file.string() + ": " + what;
        return false;
    };
    std::ifstream input(file, std::ios::binary);
    if (!input) return failed("cannot open");
    std::string header;
    std::getline(input, header);
    std::istringstream row(header);
    row.imbue(std::locale::classic());
    std::string magic;
    int longitudes = 0, points = 0;
    if (!(row >> magic >> longitudes >> points) || magic != "AXFLUXTUBE2" || longitudes != kTubeLongitudes ||
        points != kTubePoints) {
        return failed("bad header");
    }
    std::array<std::vector<glm::vec4>, 3> tubes;
    std::array<std::vector<glm::vec3>, 3> fields;
    const size_t count = static_cast<size_t>(kTubeLongitudes) * 2 * kTubePoints;
    for (size_t moon = 0; moon < tubes.size(); ++moon) {
        auto& tube = tubes[moon];
        tube.resize(count);
        const auto bytes = static_cast<std::streamsize>(count * sizeof(glm::vec4));
        input.read(reinterpret_cast<char*>(tube.data()), bytes);
        if (input.gcount() != bytes) return failed("too short");
        for (size_t i = 0; i < count; ++i) {
            const glm::vec4& p = tube[i];
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || !(p.w >= 0.0f && p.w <= 1.0f)) {
                return failed("bad point");
            }
        }
        auto& field = fields[moon];
        field.resize(kTubeLongitudes);
        const auto field_bytes = static_cast<std::streamsize>(field.size() * sizeof(glm::vec3));
        input.read(reinterpret_cast<char*>(field.data()), field_bytes);
        if (input.gcount() != field_bytes) return failed("too short");
        for (const glm::vec3& b : field) {
            if (!std::isfinite(b.x) || !std::isfinite(b.y) || !std::isfinite(b.z) || !(glm::length(b) > 1.0f)) {
                return failed("bad field");
            }
        }
    }
    if (input.peek() != std::char_traits<char>::eof()) return failed("too long");
    flux_tubes = std::move(tubes);
    moon_fields = std::move(fields);
    return true;
}

glm::dvec4 JovianAurora::tube_point(int moon, bool south, double east, int point) const
{
    const double turns = east / kTwoPi;
    const double x = (turns - std::floor(turns)) * kTubeLongitudes;
    const int i = static_cast<int>(std::floor(x)) % kTubeLongitudes;
    const int next = (i + 1) % kTubeLongitudes;
    const auto& tube = flux_tubes[static_cast<size_t>(moon)];
    const auto at = [&](int longitude) {
        const size_t index = (static_cast<size_t>(longitude) * 2 + (south ? 1 : 0)) * kTubePoints + static_cast<size_t>(point);
        return glm::dvec4(tube[index]);
    };
    // Interpolated in the line's own frame, so that the corner between two lines is not cut.
    const glm::dvec4 p = at(i) + (x - std::floor(x)) * (at(next) - at(i));
    const double c = std::cos(east), s = std::sin(east);
    return {c * p.x - s * p.y, s * p.x + c * p.y, p.z, p.w};
}

bool JovianAurora::flux_tube(int moon, bool south, const glm::dvec3& p, bool alfven_wing,
                             std::vector<glm::dvec3>& line) const
{
    line.clear();
    if (moon < 0 || moon >= 3 || flux_tubes[static_cast<size_t>(moon)].empty()) return false;
    const double east = std::atan2(p.y, p.x);
    const double lead = alfven_wing ? lead_angle(moon, south, -east) : 0.0;
    glm::dvec3 foot_direction;
    if (!footprint_at(moon, south, p, east + lead, foot_direction)) return false;
    // On the same spheroid as the table.
    constexpr double kShellA = 71492.0 + 900.0, kShellB = 66854.0 + 900.0, kRadiusKm = 71492.0;
    const glm::dvec3 foot = foot_direction / std::sqrt((foot_direction.x * foot_direction.x + foot_direction.y * foot_direction.y) / (kShellA * kShellA) +
                                                       foot_direction.z * foot_direction.z / (kShellB * kShellB));

    // The table's lines start at the moon's nominal place: move them to the moon and to the
    // foot, the moon's offset fading as the line comes in, as a dipole's field line would.
    const glm::dvec3 start = glm::dvec3(tube_point(moon, south, east, 0)) * kRadiusKm;
    const glm::dvec3 end = glm::dvec3(tube_point(moon, south, east + lead, kTubePoints - 1)) * kRadiusKm;
    const double start_radius = glm::length(start);
    const double end_radius = glm::length(glm::dvec3(tube_point(moon, south, east, kTubePoints - 1))) * kRadiusKm;
    const glm::dvec3 moon_offset = p - start;
    const glm::dvec3 foot_offset = foot - end;
    std::array<glm::dvec3, kTubePoints> knots;
    for (int k = 0; k < kTubePoints; ++k) {
        const glm::dvec4 nominal = tube_point(moon, south, east, k);
        // A part of the wave k along its way left when the moon was on a line lead * w to the east.
        const glm::dvec3 here = glm::dvec3(alfven_wing ? tube_point(moon, south, east + lead * nominal.w, k) : nominal) * kRadiusKm;
        const double near = std::clamp((glm::length(glm::dvec3(nominal)) * kRadiusKm - end_radius) / (start_radius - end_radius), 0.0, 1.0);
        const double far = 1.0 - near;
        const double far2 = far * far;
        const double far4 = far2 * far2;
        knots[static_cast<size_t>(k)] = here + near * moon_offset + far4 * far4 * foot_offset;
    }
    knots.front() = p;
    knots.back() = foot;
    // Far out the table's points are 2 radii apart: a Catmull-Rom curve through them, not corners.
    constexpr int kSubdivisions = 3;
    line.reserve((kTubePoints - 1) * kSubdivisions + 1);
    // The curve leaves the moon along the field there (southward on the equator, so the
    // northern line goes against it), not along the chord to the next point, which in the
    // stretched field of the current sheet is up to 15 deg off: a moon's own open lines
    // (Ganymede's) run on into it.
    glm::dvec3 leaving(0.0);
    if (glm::dvec3 field; field_at_moon(moon, p, field)) {
        leaving = (south ? 1.0 : -1.0) * glm::normalize(field);
    }
    for (size_t k = 0; k + 1 < knots.size(); ++k) {
        const glm::dvec3& a = knots[k];
        const glm::dvec3& b = knots[k + 1];
        const glm::dvec3 before = k > 0                       ? knots[k - 1]
                                  : leaving != glm::dvec3(0.0) ? b - 2.0 * glm::length(b - a) * leaving
                                                               : 2.0 * a - b;
        const glm::dvec3 after = k + 2 < knots.size() ? knots[k + 2] : 2.0 * b - a;
        for (int i = 0; i < kSubdivisions; ++i) {
            const double t = static_cast<double>(i) / kSubdivisions;
            line.push_back(a + 0.5 * t * ((b - before) + t * ((2.0 * before - 5.0 * a + 4.0 * b - after) +
                                                              t * (3.0 * (a - b) + after - before))));
        }
    }
    line.push_back(foot);
    return true;
}

} // namespace astraxis
