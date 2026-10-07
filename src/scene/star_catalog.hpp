#pragma once

#include <glm/vec3.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// The faintest stars drawn from the catalog (the naked-eye limit).
constexpr double kNakedEyeMag = 6.5;

// Catalog stars closer than this to the viewer are left out of its sky: they
// belong to the viewer's own system, which the scene draws as bodies.
constexpr double kOwnSystemPc = 0.5;

struct CatalogStar {
    double ra_deg = 0.0; // J2000
    double dec_deg = 0.0;
    double vmag = 0.0;
    double bv = 0.0;
    bool has_bv = false;
    double distance_pc = 0.0; // from the observer vmag refers to; 0: unknown (infinitely far)
};

// Stars as seen from the Sun, with distances.
struct StarCatalog {
    std::vector<CatalogStar> stars;
    // Complete to kNakedEyeMag for viewers within this distance of the Sun.
    double max_viewer_pc = 0.0;
};

// Loads the CSV written by tools/stars/convert_hyg.py
// (header "hip,hr,ra_deg,dec_deg,dist_pc,vmag,bv,name"; lines starting with '#'
// are comments, one of them "# max_viewer_pc = <pc>").
bool load_star_catalog(const std::filesystem::path& path, StarCatalog& out, std::string* error);

// The catalog's sky from viewer_pc (ICRF, pc from the Sun): directions and
// apparent V of the stars brighter than kNakedEyeMag there, without those
// within kOwnSystemPc of the viewer. Stars of unknown distance keep their
// direction and magnitude.
std::vector<CatalogStar> catalog_sky(const StarCatalog& catalog, const glm::dvec3& viewer_pc);

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
