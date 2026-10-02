#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

struct CatalogStar {
    double ra_deg = 0.0; // J2000
    double dec_deg = 0.0;
    double vmag = 0.0;
    double bv = 0.0;
    bool has_bv = false;
};

// Loads the CSV written by tools/stars/convert_bsc5.py
// (header "hr,ra_deg,dec_deg,vmag,bv,name"; lines starting with '#' are comments).
bool load_star_catalog(const std::filesystem::path& path, std::vector<CatalogStar>& out, std::string* error);

} // namespace astraxis
