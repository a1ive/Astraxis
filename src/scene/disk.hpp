#pragma once

#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <array>
#include <limits>
#include <string>
#include <vector>

namespace astraxis {

struct Nebula;

inline constexpr int kDiskMaxComponents = 4;
inline constexpr int kDiskRadialCount = 128;
inline constexpr int kDiskViewCount = 4;
inline constexpr int kDiskMaxClumps = 2;
// The table of optical depth from the star (see disk_build_star_depth).
inline constexpr int kDiskDepthWidth = 256;
inline constexpr int kDiskDepthHeight = 128;

// What a disk is seen in: starlight scattered by small grains, their thermal
// glow in the mid-infrared (false color), the millimeter glow of the large
// grains that stay near their parent bodies (false color), or the rotational
// line of its CO (colored by the gas's velocity along the line of sight).
enum class DiskView { Scattered = 0, Infrared = 1, Millimeter = 2, Gas = 3 };

// One part of a disk: dust in its own plane, axisymmetric, with a column
// density given by radius and a vertical shape of two Gaussians (the second is
// wider). Their sigmas are constant inside flare_radius_au and grow as
// r^flare_index outside it. Lengths are in au here: only ratios of them reach
// the light, and the shader works in single precision.
struct DiskComponent {
    std::string name;
    glm::dmat3 to_icrf{1.0};    // columns: in-plane axes, then the normal
    double sigma_au = 1.0;
    double sigma2_au = 1.0;
    double share2 = 0.0;        // the second Gaussian's share of the column
    double flare_radius_au = 1.0;
    double flare_index = 1.0;
    // In an optically thick disk: the optical depth straight through the disk,
    // per unit of the column below (0: it blocks nothing and is not lit).
    double tau_per_column = 0.0;
    // Column at disk_radial_au(i), relative: only the disk's own ratios matter.
    std::array<double, kDiskRadialCount> column{};
    std::array<double, kDiskViewCount> gain{}; // 0: not seen in that view
    glm::dvec3 color{1.0};      // scattered light (linear, largest component 1)
    double temperature_k = 100.0; // grains' temperature temperature_radius_au from the star
    double temperature_radius_au = 100.0;
    double temperature_index = 0.5; // T ~ d^-index
};

// Something small that shines in a disk (the glow of a planet accreting, a
// circumplanetary disk): a Gaussian spot, in a view's own brightness relative
// to the disk's reference. Not attenuated by the disk. Never drawn narrower than
// a pixel or so, keeping its peak: it is a marker when too small to resolve.
struct DiskSource {
    int body = -1;
    std::array<double, kDiskViewCount> sigma_au{}; // its size in each view
    std::array<double, kDiskViewCount> gain{};
    glm::dvec3 color{1.0}; // linear, in scattered light; the other views have their false color
};

// A concentration in a disk's gas (the Gas view's components): its column is
// multiplied by 1 + gain * radial Gaussian * azimuthal shape. The clump turns at
// its own rate about the plane's normal, in the sense of the plane's orbits; the
// gas it releases lags behind and fades exponentially with the angle, a Gaussian
// bounds it ahead.
struct DiskClump {
    double radius_au = 85.0;
    double sigma_au = 15.0;
    double azimuth_deg = 0.0;   // from the plane's first axis (the node), at epoch_tdb
    double epoch_tdb = 0.0;
    double period_s = 1.0;
    double sigma_deg = 15.0;    // ahead of the clump
    double tail_deg = 60.0;     // behind it
    double gain = 1.0;
};

// A disk of dust around a star. In a debris disk (optically thin, seen through)
// each sample along a line of sight adds its light, with no extinction and no
// shadow:
//   Scattered:  column * vertical shape * phase function / d^2
//   Infrared:   column * vertical shape * a black body's light in the band
//   Millimeter: column * vertical shape * T (Rayleigh-Jeans)
//   Gas:        column * vertical shape * clumps
// A protoplanetary disk is `thick` in scattered light: its components with an
// opacity lie in one plane, and each sample scatters
//   extinction * (gain * phase function * exp(-optical depth from the star) + D) / d^2
// and dims what lies behind it. With an albedo this is the "truncated multiple
// scattering" of Nakajima & Tanaka 1988 (JQSRT 40, 51): the forward share g of the
// scattering counts as no scattering at all (the similarity relations, so the
// optical depths are (1 - albedo g) times the true ones and the rest scatters
// isotropically with a smaller albedo); the first scattering keeps the true phase
// function with gain = albedo / (1 - albedo g); D is the isotropic source of all the
// later ones, on a table (disk_build_diffuse). The other views stay optically thin.
// `reference` is what a ray sees when it passes display_reference_au from the
// star along the first component's midplane (in a thick disk, straight through
// the disk there), which sets the display: the range of radiance relative to it
// is raised to display_gamma.
struct Disk {
    int center = -1;            // the star
    double inner_au = 0.0;
    double outer_au = 0.0;
    double extent_au = 0.0;     // how far the disk shows well (the director's system extent)
    std::array<double, 3> phase_g{}; // Henyey-Greenstein mixture
    std::array<double, 3> phase_w{};
    double display_reference_au = 100.0;
    double display_gamma = 0.5;
    // A camera closer to the star than this sees the disk dimmed in proportion
    // (the eye adapting), so the dust it sits in does not wash out the view.
    double display_adapt_au = 1000.0;
    std::array<double, kDiskViewCount> reference{};
    std::vector<DiskComponent> components;
    std::vector<DiskSource> sources;
    std::vector<DiskClump> clumps;
    double gas_star_gm_km3_s2 = 0.0; // the rotation the Gas view shows
    double gas_velocity_km_s = 5.0;  // the line-of-sight speed at full color
    double gas_display_gain = 1.0;   // on the Gas view's display brightness
    bool thick = false;
    // A thick disk's grains, from [disk]'s albedo and asymmetry. tau_per_column is
    // scaled already by the similarity relation; 0 albedo: single scattering, gain 1.
    double albedo = 0.0;
    double asymmetry = 0.0;
    double single_scatter_gain = 1.0;
    // ln(1 + optical depth) from the star, kDiskDepthWidth x kDiskDepthHeight
    // (see disk_build_star_depth). Empty unless thick.
    std::vector<float> star_depth;
    // ln(D + kDiskDiffuseFloor) on the same grid (see disk_build_diffuse): empty
    // without an albedo.
    std::vector<float> diffuse;

    bool enabled() const { return center >= 0 && !components.empty(); }
};

inline constexpr double kDiskDiffuseFloor = 1e-9;

// Radius of table entry i: even steps in ln r from inner to outer.
double disk_radial_au(const Disk& disk, int i);

// A component's column at radius r (linear between the table's entries).
double disk_column(const Disk& disk, const DiskComponent& component, double r_au);

// The two Gaussians' common growth with radius.
double disk_flare(const DiskComponent& component, double r_au);

// The vertical shape at (r, z), per au: its integral over z is 1.
double disk_vertical(const DiskComponent& component, double r_au, double z_au);

// Grains' temperature d from the star.
double disk_temperature(const DiskComponent& component, double d_au);

// The mixture of Henyey-Greenstein functions relative to isotropic scattering.
double disk_phase(const Disk& disk, double cos_theta);

// A component's column times its gain in a view at radius r (the shader's table).
double disk_table_light(const Disk& disk, const DiskComponent& component, DiskView view, double r_au);

// The rest of a sample's light in the thermal views, d from the star: the
// grains' black-body light in the infrared band, or their temperature at
// millimeter wavelengths (Rayleigh-Jeans). Scattered light takes the phase
// function and 1 / d^2 instead.
double disk_thermal_light(const DiskComponent& component, DiskView view, double d_au);

// What the clumps multiply the Gas view's column by at radius r and azimuth phi
// (radians, from the plane's first axis) at time t.
double disk_clump_factor(const Disk& disk, double r_au, double phi, double t_tdb);

// A clump's azimuth at time t (radians, not wrapped).
double disk_clump_azimuth(const DiskClump& clump, double t_tdb);

// The speed along `direction` (positive: moving away along it) of gas on a circular
// orbit at `position` (relative to the star, any length unit with position_au in au).
double disk_line_of_sight_km_s(const DiskComponent& component, double gm_km3_s2, const glm::dvec3& position_au,
                               const glm::dvec3& direction);

// Light along a ray (au, ICRF axes, relative to the star; `direction` a unit
// vector, the way the camera looks), integrated by the shader's steps
// (shaders/disk.hlsli) as far as t_max, with the clumps as at t_tdb.
double disk_radiance(const Disk& disk, DiskView view, const glm::dvec3& origin, const glm::dvec3& direction,
                     double t_max = std::numeric_limits<double>::infinity(), double t_tdb = 0.0);

// A thick disk's extinction at (r, z), per au.
double disk_extinction(const Disk& disk, double r_au, double z_au);

// Fills disk.star_depth: the optical depth along the straight line from the
// star, on even steps of ln r (inner to outer) and of the square root of the
// elevation angle over 90 degrees, in the plane of the components with opacity.
void disk_build_star_depth(Disk& disk);

// Optical depth from the star to (r, z), from disk.star_depth (bilinear in
// ln(1 + depth), as the shader samples it).
double disk_star_depth(const Disk& disk, double r_au, double z_au);

// Sets the similarity relations from disk.albedo and disk.asymmetry: the scaled
// opacity of the components and single_scatter_gain. Call once, before the star
// depth.
void disk_apply_similarity(Disk& disk);

// The albedo after the similarity relation: what scatters isotropically.
double disk_scaled_albedo(const Disk& disk);

// Fills disk.diffuse: the source D of the scatterings after the first, relative to
// the unattenuated starlight (units of 1 / d^2), by two more orders of scattering on
// a polar grid and a geometric series for the rest. Each order sums the source along
// a set of rays from each grid point. Needs disk.star_depth.
void disk_build_diffuse(Disk& disk);

// D at (r, z), as the shader samples it (0 with no table).
double disk_diffuse(const Disk& disk, double r_au, double z_au);

// Optical depth of a thick disk along a ray (as disk_radiance), which stops
// counting once little light is left (~5). 0 for a thin disk.
double disk_optical_depth(const Disk& disk, const glm::dvec3& origin, const glm::dvec3& direction,
                          double t_max = std::numeric_limits<double>::infinity());

// Sets disk.reference (once the components, and a thick disk's star_depth, are
// complete).
void disk_set_references(Disk& disk);

// A burst of dust from one collision (see build_dust_stream): grains leave a
// parent body on a circular orbit with its velocity plus a small kick, and
// radiation pressure takes the share `beta` of the star's gravity away.
struct DustStream {
    double gm_km3_s2 = 0.0;
    double radius_km = 0.0;      // of the parent's orbit
    double argument_deg = 0.0;   // where on it, from the ascending node
    glm::dmat3 orbit_to_icrf{1.0}; // columns: to the node, 90 deg ahead, the normal
    double beta_lo = 0.5;
    double beta_hi = 10.0;
    double age_s = 0.0;
    double kick_share = 0.0;     // RMS of the kick, per axis, relative to the orbital speed
    double weight_index = 0.0;   // a grain's light ~ beta^-weight_index
    // Drag by the interstellar medium: acceleration = beta_hi * (beta / beta_hi)^index
    // * ism_drag_per_km * |w| * w with w the wind's velocity relative to the grain.
    // Index 1 if the two forces both go with a grain's area over its mass, as for large
    // grains. Grains well under the wavelength absorb less light as they shrink, so beta
    // levels off while the drag per mass goes on rising: the index exceeds 1 at the tip.
    glm::dvec3 ism_velocity_km_s{0.0};
    double ism_drag_per_km = 0.0;
    double ism_beta_index = 1.0;
};

// What ism_drag_per_km is: drag rho v^2 pi s^2, radiation pressure's Q L s^2 /
// (4 c d^2), so 4 pi G M c rho / (Q L), with rho of hydrogen atoms at ism_h_per_cm3.
double dust_stream_drag_per_km(double gm_km3_s2, double luminosity_solar, double radiation_q, double ism_h_per_cm3);

// Fills nebula.grains with the stream's dust where it has reached after age_s
// (RK4; radiation pressure and gravity only), evenly in beta and in that order,
// their weights averaging 1. Also sets the nebula's extent and mean distance, and
// its patch size to patch_sigma_km.
void build_dust_stream(const DustStream& stream, int count, double patch_sigma_km, Nebula& nebula);

} // namespace astraxis
