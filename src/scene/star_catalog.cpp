#include "scene/star_catalog.hpp"

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

} // namespace astraxis
