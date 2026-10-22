#pragma once

#include <glm/vec3.hpp>

#include <vector>

namespace astraxis {

class Scene;

// One dust grain of a pinwheel (see Pinwheel): the direction it coasts in
// (unit, ICRF) and when it leaves, phase + turn periods after the reference
// time (0 <= phase < 1, turn = 0 .. turns - 1), and again every `turns`
// periods; `tangent` (unit) points along its ring toward the next grain, and
// `weight` (0..1) is the share of dust it carries relative to the heaviest
// grain (the production rate over the orbit and around the cone).
struct PinwheelGrain {
    glm::vec3 direction{1.0f, 0.0f, 0.0f};
    float phase = 0.0f;
    glm::vec3 tangent{0.0f, 1.0f, 0.0f};
    float turn = 0.0f;
    float weight = 1.0f;
};

// The dust plume of a colliding-wind binary (WR 104's "pinwheel", WR 140's
// shells), as the geometric model of White et al. 2025 (arXiv:2507.14610,
// Appendix B, "xenomorph"; Han et al. 2022 for WR 140): dust made where the
// winds collide leaves as a ring on a hollow cone whose axis points from the
// star with the stronger wind (`primary`) to the other one (`secondary`), and
// moves outward from the pair's barycenter (`center`). As the pair orbits, the
// cone turns and the rings wind into a spiral (an Archimedean one for a
// circular orbit at constant speed). An eccentric binary may make dust only
// around periastron (`windowed`: from true anomaly window_on to window_off),
// at a rate that varies along the orbit and around the cone.
//
// The plume is a fixed set of grains. The orbit repeats every period, so a
// grain that leaves `phase` periods (plus whole periods) after the reference
// time always leaves in the same direction. Each grain flies for `turns`
// periods and then starts again at the center (faded out there, unborn
// here): at time t it is
//   age = turns * period * frac(((t - t_ref) / period - phase - turn) / turns)
// old, at distance pinwheel_radius(age). Positions are a function of time
// alone, so jumps and running backward cost nothing; the renderer computes
// them the same way on the GPU.
struct Pinwheel {
    int center = -1;
    int primary = -1;
    int secondary = -1;
    double period_s = 0.0;
    double t_ref_tdb = 0.0;
    // Speed at launch; from the nucleation radius on, a constant acceleration
    // up to acceleration_end_km, falling off as 1/r^2 beyond (radiation
    // pressure on the dust, Han et al. 2022 eq. 1). No acceleration: constant
    // speed.
    double speed_km_s = 0.0;
    double acceleration_km_s2 = 0.0;
    double acceleration_end_km = 0.0;
    double opening_angle_deg = 90.0; // full angle of the cone
    // Dust condenses at this distance from the center (none closer), at this
    // temperature, and cools as T = T_nuc (r / r_nuc)^-temperature_index.
    double nucleation_radius_km = 0.0;
    double nucleation_temperature_k = 2000.0;
    double temperature_index = 0.45;
    // Gaussian thickness of the cone's walls, across the expanding shell (km).
    double wall_sigma_km = 0.0;
    // Share of the visible starlight a grain scatters (its albedo there).
    double visible_albedo = 0.0;
    // Dust production: all around the orbit (rings evenly spread in time), or
    // only between two true anomalies of the secondary around the primary
    // (rings evenly spread in true anomaly). Within the window the rate may
    // peak at both ends as half-Gaussians in true anomaly (modulation_sigma_deg;
    // 0: constant), and around each ring peak at the trailing edge of the cone
    // as a Gaussian in the angle around the axis (trailing_sigma_deg; 0: even).
    bool windowed = false;
    double window_on_deg = -180.0;
    double window_off_deg = 180.0;
    double modulation_sigma_deg = 0.0;
    double trailing_sigma_deg = 0.0;
    // Display: a face-on wall is as bright at display_reference_km as the
    // renderer's reference brightness (0: at the nucleation radius), and its
    // fall-off with distance (r^-4) is raised to display_gamma, which keeps
    // the outer turns in view.
    double display_reference_km = 0.0;
    double display_gamma = 0.6;
    int turns = 0;
    int rings_per_turn = 0;
    int grains_per_ring = 0;
    std::vector<PinwheelGrain> grains;
    // When each ring leaves (fraction of the period after t_ref) and the dust
    // it carries relative to the heaviest ring (its grains' weights also vary
    // around the cone).
    struct Ring {
        double phase = 0.0;
        double weight = 1.0;
    };
    std::vector<Ring> rings;

    bool enabled() const { return center >= 0 && !grains.empty(); }
    // Distance between the outermost two turns (km): speed * period without
    // acceleration.
    double step_km() const;
    // Angle the cone's axis sweeps while dust forms (radians; 2 pi all around).
    double window_rad() const;
};

// Fills pinwheel.grains (and sets rings_per_turn and grains_per_ring):
// `rings_per_turn` rings over one period (evenly spread in time, or in true
// anomaly over the window), `grains_per_ring` grains around each, repeated for
// every turn. The cone's axis at each ring's time comes from the scene's
// motions (icrf_state_at), so the scene need not be updated; a windowed
// pinwheel needs the masses of the pair (body_gm_km3_s2) for its orbit.
void build_pinwheel_grains(const Scene& scene, Pinwheel& pinwheel, int rings_per_turn, int grains_per_ring);

// Distance from the center (km) of dust `age_s` seconds after it left.
double pinwheel_radius(const Pinwheel& pinwheel, double age_s);

// Age of a grain at t_tdb (seconds; 0 .. turns * period).
double pinwheel_grain_age(const Pinwheel& pinwheel, const PinwheelGrain& grain, double t_tdb);

// Position of a grain relative to the center (ICRF, km) at t_tdb.
glm::dvec3 pinwheel_grain_offset(const Pinwheel& pinwheel, const PinwheelGrain& grain, double t_tdb);

// Dust production rate per unit time at true anomaly nu_deg (relative to the
// peak; 0 outside the window), and its share at angle psi_deg around the
// cone's axis from the leading edge (1 at the trailing edge, psi = 180).
double pinwheel_orbit_rate(const Pinwheel& pinwheel, double nu_deg);
double pinwheel_azimuth_rate(const Pinwheel& pinwheel, double psi_deg);

// The brightest a face-on wall of the youngest dust is at t_tdb, relative to
// one at display_reference_km (no closer than the nucleation radius, where
// the dust appears): the largest ring weight times
// ((reference / r)^4)^display_gamma over each ring's last launch. A pinwheel
// making dust all the time (WR 104, reference at the nucleation radius) is
// always at 1; WR 140's newest shell, near the stars, at thousands.
double pinwheel_peak_brightness(const Pinwheel& pinwheel, double t_tdb);

// The exposure the eye would adapt to at t_tdb (the dust's brightness factor,
// 0..1) from camera_km off the center: walls brighter than twice the
// reference bring everything down so that the brightest stays at that, but
// by no more than 20 times, so the older shells stay in sight while new dust
// forms (WR 140). Inside the innermost turn (step_km), where the old shells
// are out of the picture or far off, the cap opens up as the distance squared.
double pinwheel_exposure(const Pinwheel& pinwheel, double t_tdb, double camera_km);

// Dust temperature (K) at distance r_km from the center (that of the
// nucleation radius closer in).
double pinwheel_dust_temperature(const Pinwheel& pinwheel, double r_km);

// Light of one grain at distance r_km from the center, relative to a grain at
// the nucleation radius: the starlight it takes up, (r_nuc / r)^2, so that the
// face-on surface brightness falls as r^-3. It fades in over the first 25%
// beyond the nucleation radius (0 closer) and out over the last turn, where
// grains wrap around to the center. A grain's own light is this times its
// weight.
double pinwheel_grain_brightness(const Pinwheel& pinwheel, double r_km);

// Band of "visible light" for the eye's view of the dust (nm).
inline constexpr double kVisibleLoNm = 400.0;
inline constexpr double kVisibleHiNm = 700.0;

// The visible light of a grain at distance r_km from the center, as shares of
// the starlight it takes up: what it absorbs and re-emits as a black body at
// its temperature, and what it scatters (visible_albedo times the stars'
// visible share).
struct PinwheelVisibleLight {
    double thermal = 0.0;
    double scattered = 0.0;
};
PinwheelVisibleLight pinwheel_visible_light(const Pinwheel& pinwheel, double r_km, double star_visible_share);

// Share of the pair's light in the visible band (black bodies at the two
// stars' temperatures, weighted by their luminosities).
double pinwheel_star_visible_share(const Scene& scene, const Pinwheel& pinwheel);

} // namespace astraxis
