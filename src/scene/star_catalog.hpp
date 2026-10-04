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

// A star cluster seen from inside: member positions in pc relative to the
// cluster centre on the tangent-plane axes at the centre (east, north, away
// from us) with absolute V magnitudes and B - V, as written by
// tools/stars/make_m4.py. The viewer sits at sky position viewer_ra/dec_deg,
// viewer_depth_pc behind the centre's distance.
struct ClusterView {
    double center_ra_deg = 0.0;
    double center_dec_deg = 0.0;
    double distance_pc = 0.0;
    double viewer_ra_deg = 0.0;
    double viewer_dec_deg = 0.0;
    double viewer_depth_pc = 0.0;
};

// The members as the viewer sees them: ICRF directions and apparent V.
bool load_cluster_stars(const std::filesystem::path& path, const ClusterView& view, std::vector<CatalogStar>& out,
                        std::string* error);

} // namespace astraxis
