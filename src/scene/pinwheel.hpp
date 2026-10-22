#pragma once

#include <glm/vec3.hpp>

#include <vector>

namespace astraxis {

class Scene;

// One dust grain of a pinwheel (see Pinwheel): the direction it coasts in
// (unit, ICRF) and when it leaves, phase + turn periods after the reference
// time (0 <= phase < 1, turn = 0 .. turns - 1), and again every `turns`
// periods; `tangent` (unit) points along its ring toward the next grain.
// 32 bytes; the renderer uploads these as they are.
struct PinwheelGrain {
    glm::vec3 direction{1.0f, 0.0f, 0.0f};
    float phase = 0.0f;
    glm::vec3 tangent{0.0f, 1.0f, 0.0f};
    float turn = 0.0f;
};
static_assert(sizeof(PinwheelGrain) == 32, "PinwheelGrain must match the shader's raw buffer loads");

// The dust plume of a colliding-wind binary (WR 104's "pinwheel"), as the
// geometric model of White et al. 2025 (arXiv:2507.14610, Appendix B,
// "xenomorph"): dust made where the winds collide leaves as a ring on a
// hollow cone whose axis points from the star with the stronger wind
// (`primary`) to the other one (`secondary`), and coasts outward from the
// pair's barycenter (`center`) at constant speed. As the pair orbits, the cone
// turns and the rings wind into a spiral (an Archimedean one for a circular
// orbit: speed * period between turns).
//
// The plume is a fixed set of grains. The orbit repeats every period, so a
// grain that leaves `phase` periods (plus whole periods) after the reference
// time always leaves in the same direction. Each grain coasts for `turns`
// periods and then starts again at the center (faded out there, unborn
// here): at time t it is
//   age = turns * period * frac(((t - t_ref) / period - phase - turn) / turns)
// old. Positions are a function of time alone, so jumps and running backward
// cost nothing; the renderer computes them the same way on the GPU.
struct Pinwheel {
    int center = -1;
    int primary = -1;
    int secondary = -1;
    double period_s = 0.0;
    double t_ref_tdb = 0.0;
    double speed_km_s = 0.0;
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
    int turns = 0;
    int rings_per_turn = 0;
    int grains_per_ring = 0;
    std::vector<PinwheelGrain> grains;

    bool enabled() const { return center >= 0 && !grains.empty(); }
    // Distance between successive turns of the spiral (km).
    double step_km() const { return speed_km_s * period_s; }
};

// Fills pinwheel.grains (and sets rings_per_turn and grains_per_ring):
// `rings_per_turn` rings evenly spread over one period, `grains_per_ring`
// grains around each, repeated for every turn. The cone's axis at each ring's
// time comes from the scene's motions (icrf_state_at), so the scene need not
// be updated.
void build_pinwheel_grains(const Scene& scene, Pinwheel& pinwheel, int rings_per_turn, int grains_per_ring);

// Age of a grain at t_tdb (seconds; 0 .. turns * period).
double pinwheel_grain_age(const Pinwheel& pinwheel, const PinwheelGrain& grain, double t_tdb);

// Position of a grain relative to the center (ICRF, km) at t_tdb.
glm::dvec3 pinwheel_grain_offset(const Pinwheel& pinwheel, const PinwheelGrain& grain, double t_tdb);

// Dust temperature (K) at distance r_km from the center (that of the
// nucleation radius closer in).
double pinwheel_dust_temperature(const Pinwheel& pinwheel, double r_km);

// Light of one grain at distance r_km from the center, relative to a grain at
// the nucleation radius: the starlight it takes up, (r_nuc / r)^2, so that the
// face-on surface brightness falls as r^-3. It fades in over the first 25%
// beyond the nucleation radius (0 closer) and out over the last turn, where
// grains wrap around to the center.
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
