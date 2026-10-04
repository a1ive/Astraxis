#include "scene/star_catalog.hpp"

#include "core/math.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>

namespace astraxis {

namespace {

// Splits on commas without allocating per field beyond the output strings.
int split_csv(const std::string& line, std::string* fields, int max_fields)
{
    int count = 0;
    size_t start = 0;
    while (count < max_fields) {
        const size_t comma = line.find(',', start);
        fields[count++] = line.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        if (comma == std::string::npos) {
            break;
        }
        start = comma + 1;
    }
    return count;
}

} // namespace

bool load_star_catalog(const std::filesystem::path& path, std::vector<CatalogStar>& out, std::string* error)
{
    std::ifstream file(path);
    if (!file) {
        if (error) {
            *error = "cannot open " + path.string();
        }
        return false;
    }

    out.clear();
    std::string line;
    bool header_seen = false;
    int line_number = 0;
    std::string fields[6];
    while (std::getline(file, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        if (!header_seen) {
            header_seen = true; // "hr,ra_deg,dec_deg,vmag,bv,name"
            continue;
        }
        if (split_csv(line, fields, 6) < 5) {
            if (error) {
                *error = path.string() + ":" + std::to_string(line_number) + ": too few fields";
            }
            return false;
        }
        CatalogStar star;
        star.ra_deg = std::strtod(fields[1].c_str(), nullptr);
        star.dec_deg = std::strtod(fields[2].c_str(), nullptr);
        star.vmag = std::strtod(fields[3].c_str(), nullptr);
        star.has_bv = !fields[4].empty();
        star.bv = star.has_bv ? std::strtod(fields[4].c_str(), nullptr) : 0.0;
        out.push_back(star);
    }
    return true;
}

bool load_cluster_stars(const std::filesystem::path& path, const ClusterView& view, std::vector<CatalogStar>& out,
                        std::string* error)
{
    std::ifstream file(path);
    if (!file) {
        if (error) {
            *error = "cannot open " + path.string();
        }
        return false;
    }

    // Tangent-plane axes at the centre, and the viewer on them (gnomonic projection).
    const double ra0 = view.center_ra_deg * kDegToRad;
    const double dec0 = view.center_dec_deg * kDegToRad;
    const glm::dvec3 east(-std::sin(ra0), std::cos(ra0), 0.0);
    const glm::dvec3 north(-std::sin(dec0) * std::cos(ra0), -std::sin(dec0) * std::sin(ra0), std::cos(dec0));
    const glm::dvec3 away = unit_from_ra_dec(ra0, dec0);
    const glm::dvec3 viewer_dir = unit_from_ra_dec(view.viewer_ra_deg * kDegToRad, view.viewer_dec_deg * kDegToRad);
    const double along = glm::dot(viewer_dir, away);
    const glm::dvec3 viewer(glm::dot(viewer_dir, east) / along * view.distance_pc,
                            glm::dot(viewer_dir, north) / along * view.distance_pc, view.viewer_depth_pc);

    out.clear();
    std::string line;
    bool header_seen = false;
    int line_number = 0;
    std::string fields[5];
    while (std::getline(file, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        if (!header_seen) {
            header_seen = true; // "east_pc,north_pc,away_pc,abs_vmag,bv"
            continue;
        }
        if (split_csv(line, fields, 5) < 5) {
            if (error) {
                *error = path.string() + ":" + std::to_string(line_number) + ": too few fields";
            }
            return false;
        }
        const glm::dvec3 star(std::strtod(fields[0].c_str(), nullptr), std::strtod(fields[1].c_str(), nullptr),
                              std::strtod(fields[2].c_str(), nullptr));
        const glm::dvec3 d = star - viewer;
        const double distance = std::max(glm::length(d), 1e-3); // pc
        const glm::dvec3 dir = (east * d.x + north * d.y + away * d.z) / distance;
        CatalogStar s;
        s.ra_deg = wrap_two_pi(std::atan2(dir.y, dir.x)) * kRadToDeg;
        s.dec_deg = std::asin(std::clamp(dir.z, -1.0, 1.0)) * kRadToDeg;
        s.vmag = std::strtod(fields[3].c_str(), nullptr) + 5.0 * std::log10(distance / 10.0);
        s.bv = std::strtod(fields[4].c_str(), nullptr);
        s.has_bv = true;
        out.push_back(s);
    }
    return true;
}

} // namespace astraxis
