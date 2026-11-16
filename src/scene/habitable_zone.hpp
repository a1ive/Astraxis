#pragma once

#include <glm/vec3.hpp>

#include <array>
#include <vector>

namespace astraxis {

class Scene;
struct Body;

// Habitable zones of main-sequence stars and of binaries of them, where an
// Earth-like planet could keep liquid water on its surface.
//
// [KO14] Kopparapu et al. 2014, ApJL 787, L29 (arXiv:1404.5292): Eq. (4)-(5) and
//        Table 1 (1 Earth mass), for 2600 K <= T_eff <= 7200 K. The conservative
//        zone runs from the runaway greenhouse to the maximum greenhouse limit; the
//        optimistic one from "recent Venus" to "early Mars" (the insolation those
//        planets got when they last had surface water).
// [HK13] Haghighipour & Kaltenegger 2013, ApJ 777, 166: in a binary, a limit
//        lies where the stars' insolations, each divided by its own star's S_eff,
//        add up to 1 (S_eff depends on the spectrum, so on T_eff).

// Insolation at each limit, relative to the Earth's (S_eff of [KO14]).
struct HabitableZoneFlux {
    double recent_venus = 0.0;
    double runaway_greenhouse = 0.0;
    double maximum_greenhouse = 0.0;
    double early_mars = 0.0;
};

// [KO14] Eq. (4) at `teff_k`.
HabitableZoneFlux habitable_zone_flux(double teff_k);

// Range of T_eff with a zone: [KO14]'s, extended 100 K down so that TRAPPIST-1
// (2566 K) is in, as the papers on its planets apply these limits.
inline constexpr double kHabitableZoneMinTeffK = 2500.0;
inline constexpr double kHabitableZoneMaxTeffK = 7200.0;
// Main-sequence stars in that range are smaller; evolved ones (Algol B,
// Polaris, delta Cephei) larger.
inline constexpr double kHabitableZoneMaxRadiusSolar = 2.0;

// Whether a star gets a zone: a main-sequence star (by its radius) with its
// T_eff in range and with sourced values (hr_diagram), not pulsating. Its
// visibility is not checked.
bool has_habitable_zone(const Body& body);

inline constexpr int kMaxHabitableZoneStars = 4;

// One zone, drawn on a plane: the stars whose light shapes it, and the annulus
// around `center` that holds its optimistic part.
struct HabitableZone {
    int count = 0;
    std::array<int, kMaxHabitableZoneStars> stars{};
    // Per star: L / S_eff for each limit (au^2): the square of the distance of
    // that limit from the star, alone.
    std::array<HabitableZoneFlux, kMaxHabitableZoneStars> weights_au2{};
    glm::dvec3 center{0.0}; // display frame (km)
    glm::dvec3 normal{0.0, 0.0, 1.0};
    double inner_km = 0.0;
    double outer_km = 0.0;
    double width_km = 0.0; // of the conservative zone (along the radius, roughly)
};

// The zones of the visible stars at the scene's current time (after
// Scene::update). Stars whose zones would overlap share one (a planet around
// both, Kepler-47's); a star bright enough nearby to shift a zone, but without
// one of its own (a hot or evolved companion), leaves it out. A zone lies in
// the orbital plane of its planet nearest the zone's middle, else of its stars
// around each other, else across the scene's up axis.
void find_habitable_zones(const Scene& scene, std::vector<HabitableZone>& out);

// Each limit's sum of insolation over S_eff at a point (display frame, km):
// above 1, closer in than the limit (too hot for the inner two); below 1,
// farther out. Within the conservative zone runaway_greenhouse < 1 <
// maximum_greenhouse.
HabitableZoneFlux habitable_zone_insolation(const Scene& scene, const HabitableZone& zone,
                                            const glm::dvec3& point);

} // namespace astraxis
