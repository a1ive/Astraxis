#include "scene/habitable_zone.hpp"

#include "core/math.hpp"
#include "scene/scene.hpp"

#include <glm/geometric.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace astraxis {

namespace {

// [KO14] Eq. (4): S_eff = S_eff_sun + a T + b T^2 + c T^3 + d T^4, T = T_eff - 5780 K.
struct Limit {
    double s_sun, a, b, c, d;
};
// [KO14] Table 1, 1 Earth mass.
constexpr Limit kRecentVenus{1.776, 2.136e-4, 2.533e-8, -1.332e-11, -3.097e-15};
constexpr Limit kRunawayGreenhouse{1.107, 1.332e-4, 1.58e-8, -8.308e-12, -1.931e-15};
constexpr Limit kMaximumGreenhouse{0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16};
constexpr Limit kEarlyMars{0.32, 5.547e-5, 1.526e-9, -2.874e-12, -5.011e-16};

// [RA18] Table 1, with the corrected Early Mars cubic coefficient.
constexpr Limit kHotRecentVenus{1.768, 1.3151e-4, 5.8695e-10, -2.8895e-12, 3.2174e-16};
constexpr Limit kHotRunawayGreenhouse{1.1105, 1.1921e-4, 9.5932e-9, -2.6189e-12, 1.3710e-16};
constexpr Limit kHotMaximumGreenhouse{0.3587, 5.8087e-5, 1.5393e-9, -8.3547e-13, 1.0319e-16};
constexpr Limit kHotEarlyMars{0.3246, 5.213e-5, 4.5245e-10, -1.0223e-12, 9.6376e-17};

constexpr double kKopparapuMaxTeffK = 7200.0;
constexpr double kTemperatureOriginK = 5780.0;

constexpr double s_eff(const Limit& l, double t)
{
    return l.s_sun + t * (l.a + t * (l.b + t * (l.c + t * l.d)));
}

double hot_s_eff(const Limit& cool, const Limit& hot, double t)
{
    // Take the fractional change along [RA18]'s curve from the [KO14] endpoint.
    constexpr double join = kKopparapuMaxTeffK - kTemperatureOriginK;
    return s_eff(cool, join) * (s_eff(hot, t) / s_eff(hot, join));
}

// Another star's insolation on a zone (L_sun / au^2) below which it is left
// out: 1% of the smallest S_eff (early Mars, ~0.2 at 2600 K).
constexpr double kNegligibleInsolation = 0.002;
// Two stars share a zone where one would add more than ~5% to the insolation
// at the other's outer edge: within this many of its own outer radii of it.
constexpr double kShareOuterRadii = 4.5;
// The annulus beyond the limits found from the stars' spread around its center.
constexpr double kInnerMargin = 0.9;
constexpr double kOuterMargin = 1.15;

struct Candidate {
    int body = -1;
    glm::dvec3 position{0.0};
    HabitableZoneFlux weights;
    double outer_au = 0.0; // early Mars, alone
};

HabitableZoneFlux weights_au2(const Body& body)
{
    const HabitableZoneFlux s = habitable_zone_flux(body.temperature_k);
    const double l = body.luminosity_solar;
    return {l / s.recent_venus, l / s.runaway_greenhouse, l / s.maximum_greenhouse, l / s.early_mars};
}

int find_root(std::vector<int>& parent, int i)
{
    while (parent[static_cast<size_t>(i)] != i) {
        parent[static_cast<size_t>(i)] = parent[static_cast<size_t>(parent[static_cast<size_t>(i)])];
        i = parent[static_cast<size_t>(i)];
    }
    return i;
}

} // namespace

HabitableZoneFlux habitable_zone_flux(double teff_k)
{
    const double t = teff_k - kTemperatureOriginK;
    if (teff_k > kKopparapuMaxTeffK) {
        return {hot_s_eff(kRecentVenus, kHotRecentVenus, t), hot_s_eff(kRunawayGreenhouse, kHotRunawayGreenhouse, t),
                hot_s_eff(kMaximumGreenhouse, kHotMaximumGreenhouse, t), hot_s_eff(kEarlyMars, kHotEarlyMars, t)};
    }
    return {s_eff(kRecentVenus, t), s_eff(kRunawayGreenhouse, t), s_eff(kMaximumGreenhouse, t),
            s_eff(kEarlyMars, t)};
}

bool has_habitable_zone(const Body& body)
{
    const bool inferred_main_sequence = body.temperature_k <= kKopparapuMaxTeffK &&
        body.equatorial_radius_km >= kHabitableZoneMinRadiusSolar * kSunRadiusKm &&
        body.equatorial_radius_km <= kHabitableZoneMaxRadiusSolar * kSunRadiusKm;
    return body.kind == BodyKind::Star && !body.neutron_star && !body.pulsation && body.hr_diagram &&
           (body.white_dwarf || body.main_sequence.value_or(inferred_main_sequence)) && body.luminosity_solar > 0.0 &&
           body.temperature_k >= kHabitableZoneMinTeffK &&
           body.temperature_k <= kHabitableZoneMaxTeffK;
}

void find_habitable_zones(const Scene& scene, std::vector<HabitableZone>& out)
{
    out.clear();
    std::vector<Candidate> candidates;
    for (size_t i = 0; i < scene.bodies.size(); ++i) {
        const Body& body = scene.bodies[i];
        if (body.visible && has_habitable_zone(body)) {
            Candidate c;
            c.body = static_cast<int>(i);
            c.position = body.world_position;
            c.weights = weights_au2(body);
            c.outer_au = std::sqrt(c.weights.early_mars);
            candidates.push_back(c);
        }
    }
    const int n = static_cast<int>(candidates.size());
    std::vector<int> parent(static_cast<size_t>(n));
    std::iota(parent.begin(), parent.end(), 0);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const Candidate& a = candidates[static_cast<size_t>(i)];
            const Candidate& b = candidates[static_cast<size_t>(j)];
            const double d = glm::length(a.position - b.position) / kAuKm;
            if (d < std::max(a.outer_au + kShareOuterRadii * b.outer_au, b.outer_au + kShareOuterRadii * a.outer_au)) {
                parent[static_cast<size_t>(find_root(parent, i))] = find_root(parent, j);
            }
        }
    }

    for (int root = 0; root < n; ++root) {
        if (find_root(parent, root) != root) {
            continue;
        }
        HabitableZone zone;
        bool too_many = false;
        double luminosity = 0.0;
        HabitableZoneFlux sum;
        for (int i = 0; i < n; ++i) {
            if (find_root(parent, i) != root) {
                continue;
            }
            if (zone.count == kMaxHabitableZoneStars) {
                too_many = true;
                break;
            }
            const Candidate& c = candidates[static_cast<size_t>(i)];
            zone.stars[static_cast<size_t>(zone.count)] = c.body;
            zone.weights_au2[static_cast<size_t>(zone.count)] = c.weights;
            ++zone.count;
            const double l = scene.bodies[static_cast<size_t>(c.body)].luminosity_solar;
            zone.center += l * c.position;
            luminosity += l;
            sum.recent_venus += c.weights.recent_venus;
            sum.runaway_greenhouse += c.weights.runaway_greenhouse;
            sum.maximum_greenhouse += c.weights.maximum_greenhouse;
            sum.early_mars += c.weights.early_mars;
        }
        if (too_many) {
            continue;
        }
        zone.center /= luminosity;
        const int members = zone.count;
        double spread_au = 0.0;
        for (int k = 0; k < members; ++k) {
            const glm::dvec3& p = scene.bodies[static_cast<size_t>(zone.stars[static_cast<size_t>(k)])].world_position;
            spread_au = std::max(spread_au, glm::length(p - zone.center) / kAuKm);
        }
        const double outer_au = std::sqrt(sum.early_mars) + spread_au;
        zone.inner_km = kInnerMargin * std::max(0.0, std::sqrt(sum.recent_venus) - spread_au) * kAuKm;
        zone.outer_km = kOuterMargin * outer_au * kAuKm;
        zone.width_km = (std::sqrt(sum.maximum_greenhouse) - std::sqrt(sum.runaway_greenhouse)) * kAuKm;

        // The light of other stars: those with zones of their own add to it,
        // those without leave the limits unknown.
        bool blocked = false;
        for (size_t k = 0; k < scene.bodies.size() && !blocked; ++k) {
            const Body& star = scene.bodies[k];
            if (star.kind != BodyKind::Star || !star.visible ||
                std::find(zone.stars.begin(), zone.stars.begin() + members, static_cast<int>(k)) !=
                    zone.stars.begin() + members) {
                continue;
            }
            const double gap_au = std::max(glm::length(star.world_position - zone.center) / kAuKm - outer_au, 1e-3);
            if (star.luminosity_solar / (gap_au * gap_au) < kNegligibleInsolation) {
                continue;
            }
            if (!has_habitable_zone(star)) {
                blocked = true;
            } else if (zone.count < kMaxHabitableZoneStars) {
                zone.stars[static_cast<size_t>(zone.count)] = static_cast<int>(k);
                zone.weights_au2[static_cast<size_t>(zone.count)] = weights_au2(star);
                ++zone.count;
            }
        }
        if (blocked) {
            continue;
        }

        // The plane: of the planet around these stars nearest the zone's middle
        // (always shown, not one an event shows), else of the stars' own orbits.
        // A planet far beyond the zone (PSR B1620-26 b, 300 times as far as its
        // white dwarf's) says nothing of where a planet in it would orbit.
        const auto is_member = [&](int body) {
            return std::find(zone.stars.begin(), zone.stars.begin() + members, body) != zone.stars.begin() + members;
        };
        const double middle_km = std::sqrt(std::sqrt(sum.runaway_greenhouse * sum.maximum_greenhouse)) * kAuKm;
        double nearest = std::numeric_limits<double>::infinity();
        bool found = false;
        for (size_t p = 0; p < scene.bodies.size(); ++p) {
            const Body& planet = scene.bodies[p];
            if (planet.kind != BodyKind::Planet || !planet.visible || planet.event >= 0 || !planet.motion) {
                continue;
            }
            const int center = scene.orbit_center(static_cast<int>(p)).body;
            if (center < 0) {
                continue;
            }
            const bool around = is_member(center) ||
                (scene.bodies[static_cast<size_t>(center)].kind == BodyKind::Barycenter &&
                 std::any_of(zone.stars.begin(), zone.stars.begin() + members,
                             [&](int s) { return scene.bodies[static_cast<size_t>(s)].parent == center; }));
            const double r = glm::length(planet.world_position - zone.center);
            if (r > kHabitableZonePlaneReach * outer_au * kAuKm) {
                continue;
            }
            const double off = r > 0.0 ? std::abs(std::log(r / middle_km)) : nearest;
            if (around && off < nearest) {
                nearest = off;
                zone.normal = scene.orbit_normal(static_cast<int>(p));
                found = true;
            }
        }
        if (!found) {
            // A pair's orbit, or a lone star's around its barycenter (the up
            // axis without a motion).
            zone.normal = scene.orbit_normal(zone.stars[static_cast<size_t>(members - 1)]);
        }
        out.push_back(zone);
    }
}

HabitableZoneFlux habitable_zone_insolation(const Scene& scene, const HabitableZone& zone, const glm::dvec3& point)
{
    HabitableZoneFlux sum;
    for (int k = 0; k < zone.count; ++k) {
        const glm::dvec3& star = scene.bodies[static_cast<size_t>(zone.stars[static_cast<size_t>(k)])].world_position;
        const double r = glm::length(point - star) / kAuKm;
        const double inv_r2 = 1.0 / std::max(r * r, 1e-30);
        const HabitableZoneFlux& w = zone.weights_au2[static_cast<size_t>(k)];
        sum.recent_venus += w.recent_venus * inv_r2;
        sum.runaway_greenhouse += w.runaway_greenhouse * inv_r2;
        sum.maximum_greenhouse += w.maximum_greenhouse * inv_r2;
        sum.early_mars += w.early_mars * inv_r2;
    }
    return sum;
}

} // namespace astraxis
