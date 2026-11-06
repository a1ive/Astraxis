#pragma once

#include <glm/vec3.hpp>

#include <string>
#include <vector>

namespace astraxis {

// What a comet fragment did to a gas giant: Shoemaker-Levy 9 at Jupiter, July 1994.
//
// [HD01]  Harrington & Deming 2001, ApJ 561, 455, Table 2: the ballistic plume.
// [HAM95] Hammel et al. 1995, Science 267, 1288: sites, plumes, rings; site classes (Table 2).
// [HAM96] Hammel 1996, in "The Collision of Comet Shoemaker-Levy 9 and Jupiter", pp. 111-120.
// [BEE96] Beebe 1996, same book, pp. 307-328: how the sites spread out.
// [SL98]  Sanchez-Lavega et al. 1998, Icarus 131, 341: two years of the debris (Tables III, IV).
//
// All of it depends on the time alone, so jumps and reversed time need no history.
// shaders/impact_plume.frag.hlsl and impact_scar.frag.hlsl hold the same formulas: change
// them together.

inline constexpr int kImpactSiteLimit = 24; // a uniform array in the scar shader

// [HAM95] Table 2: 1 = dark region over 10,000 km, large ejecta, probably several rings; 2a = 4,000
// to 8,000 km, medium ejecta; 2b = slightly smaller, one ring; 2c = under 6,000 km, on earlier sites;
// 3 = under 3,000 km, no ejecta, no ring; 4 = not detected.
enum class ImpactClass { C1, C2a, C2b, C2c, C3, C4 };

struct ImpactClassParams {
    double ejecta_km;       // outer edge of the crescent; 0: no ejecta [HAM95]
    double core_km;         // sigma of the dark core
    double weight;          // the debris' optical depth, relative to the largest (display choice)
    double plume_weight;    // the plume's mass, relative to the largest (display choice)
    int rings;              // 0, 1 main ring, 2 and the inner one
    double fade_days;       // an e-folding time for the small sites; 0: they follow [SL98]
    double expansion_km_s;  // [SL98] Table III: outward over the first two days (K, L, G only)
};
const ImpactClassParams& impact_class_params(ImpactClass klass);

struct ImpactSite {
    std::string name;
    ImpactClass klass = ImpactClass::C4;
    double t_tdb = 0.0;       // at Jupiter, not at the Earth
    double lat_deg = 0.0;     // planetocentric
    double west_lon_deg = 0.0; // System III
};

// --- The plume [HD01] ---------------------------------------------------------
// A swarm thrown out at t = 0 from one point into a cone, then falling freely: M(>v) = k v^-alpha
// between v_min and v_max, the mass above v_max gathered in a thin shell at v_max (the "vanguard").
inline constexpr double kPlumeVMinKmS = 4.5;
inline constexpr double kPlumeVMaxKmS = 11.81;
inline constexpr double kPlumeAlpha = 1.55;
inline constexpr double kPlumeGravityKmS2 = 0.02325;
inline constexpr double kPlumeConeDeg = 75.0;      // opening angle
inline constexpr double kPlumeTiltDeg = 30.0;      // of the axis from the vertical
inline constexpr double kPlumeAzimuthDeg = 145.0;  // of the axis, north through east
inline constexpr double kPlumeVanguardWidth = 0.08; // the shell's sigma, of v_max (it was a delta function)

// Unit vector of the cone's axis: x east, y north, z up.
glm::dvec3 impact_plume_axis();
double impact_plume_vanguard_share();
// The top of the swarm above the 100-mbar level, and when it has all come down.
double impact_plume_top_km(double age_s);
double impact_plume_end_s();
// Mass per km^3, as a share of the plume's mass, at (x east, y north, altitude) km.
double impact_plume_density(const glm::dvec3& p, double age_s);

// --- The flash [HAM95] --------------------------------------------------------
// Light of the bolide and the fireball, as a share of Jupiter's integrated flux: up in about 2 s,
// over in 25 to 50 s; G reached 3% outside the methane bands.
inline constexpr double kImpactFlashShare = 0.03;
double impact_flash(double age_s);

// --- Rings [HAM95] ------------------------------------------------------------
// Waves in the atmosphere: 454 m/s from 586 km, the inner one about 189 m/s. Seen 1 to 2.5 h after an
// impact, and not after 10 h.
double impact_ring_radius_km(double age_s, bool inner);

// --- The debris [HAM95, BEE96, SL98] -------------------------------------------
// Optical depth at 890 nm in the core of a large site, from [SL98] Table IV.
double impact_core_depth(double age_s);
// The band that the debris spreads into, in visible light, with `days` since the first impact.
double impact_band_depth(double days);
// How much of a band of that depth lies at this latitude: 1 inside the main part.
double impact_band_profile(double lat_deg, double days);
// Share of a site's own pattern at this age (it goes into the band in the end).
double impact_site_fade(const ImpactClassParams& params, double age_s);
// Eastward wind relative to the site at `north_km` from it [SL98] Table III, km/s.
double impact_wind_km_s(double north_km);
// The site's pattern, 1 at its darkest: at (east, north) km from the impact after the wind has
// been undone. `seed` varies the rays and the waviness.
double impact_scar_shape(const ImpactClassParams& params, double east_km, double north_km, double age_s, double seed);
// A seed for a site; the same on the CPU and the GPU.
double impact_seed(const ImpactSite& site);
// Optical depth at 890 nm over all the sites at a place on the planet.
double impact_scar_depth(const std::vector<ImpactSite>& sites, double lat_deg, double west_lon_deg, double t_tdb);

// The distance over which latitude and longitude turn into km near the sites.
inline constexpr double kImpactRadiusKm = 69300.0;
// Time after which nothing is left to draw (about 2000: [SL98] still saw some in 1996).
inline constexpr double kImpactBandDays = 2000.0;

} // namespace astraxis
