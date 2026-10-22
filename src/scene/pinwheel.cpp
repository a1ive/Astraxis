#include "scene/pinwheel.hpp"

#include "core/color.hpp"
#include "core/math.hpp"
#include "scene/scene.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

void build_pinwheel_grains(const Scene& scene, Pinwheel& pinwheel, int rings_per_turn, int grains_per_ring)
{
    pinwheel.rings_per_turn = rings_per_turn;
    pinwheel.grains_per_ring = grains_per_ring;
    pinwheel.grains.clear();
    pinwheel.grains.reserve(static_cast<size_t>(pinwheel.turns) * static_cast<size_t>(rings_per_turn) *
                            static_cast<size_t>(grains_per_ring));
    const double half = 0.5 * pinwheel.opening_angle_deg * kDegToRad;
    // Golden-ratio offsets of the grains around each ring, so that grains of
    // neighboring rings and turns do not line up into streaks.
    constexpr double kGolden = 0.6180339887498949;
    for (int k = 0; k < rings_per_turn; ++k) {
        const double phase = (k + 0.5) / rings_per_turn;
        const double t = pinwheel.t_ref_tdb + phase * pinwheel.period_s;
        const State p = scene.icrf_state_at(pinwheel.primary, t);
        const State s = scene.icrf_state_at(pinwheel.secondary, t);
        // The cone (White et al. 2025, eqs. B9-B12): its axis from the primary
        // to the secondary, the ring around it in the plane of the axis' normal
        // in the orbit plane and the orbit normal.
        const glm::dvec3 rel = s.position - p.position;
        const glm::dvec3 axis = glm::normalize(rel);
        const glm::dvec3 normal = glm::normalize(glm::cross(rel, s.velocity - p.velocity));
        const glm::dvec3 side = glm::cross(normal, axis);
        for (int n = 0; n < pinwheel.turns; ++n) {
            const double offset = std::fmod((k + n * rings_per_turn) * kGolden, 1.0);
            for (int j = 0; j < grains_per_ring; ++j) {
                const double psi = kTwoPi * (j + offset) / grains_per_ring;
                const glm::dvec3 d = axis * std::cos(half) +
                                     (side * std::sin(psi) + normal * std::cos(psi)) * std::sin(half);
                PinwheelGrain g;
                g.direction = glm::vec3(d);
                g.tangent = glm::vec3(side * std::cos(psi) - normal * std::sin(psi));
                g.phase = static_cast<float>(phase);
                g.turn = static_cast<float>(n);
                pinwheel.grains.push_back(g);
            }
        }
    }
}

double pinwheel_grain_age(const Pinwheel& pinwheel, const PinwheelGrain& grain, double t_tdb)
{
    const double n = static_cast<double>(pinwheel.turns);
    const double y = ((t_tdb - pinwheel.t_ref_tdb) / pinwheel.period_s - static_cast<double>(grain.phase) -
                      static_cast<double>(grain.turn)) /
                     n;
    return n * pinwheel.period_s * (y - std::floor(y));
}

glm::dvec3 pinwheel_grain_offset(const Pinwheel& pinwheel, const PinwheelGrain& grain, double t_tdb)
{
    return glm::dvec3(grain.direction) * (pinwheel.speed_km_s * pinwheel_grain_age(pinwheel, grain, t_tdb));
}

double pinwheel_dust_temperature(const Pinwheel& pinwheel, double r_km)
{
    const double r = std::max(r_km, pinwheel.nucleation_radius_km);
    return pinwheel.nucleation_temperature_k * std::pow(r / pinwheel.nucleation_radius_km, -pinwheel.temperature_index);
}

double pinwheel_grain_brightness(const Pinwheel& pinwheel, double r_km)
{
    const double r_nuc = pinwheel.nucleation_radius_km;
    if (r_km <= r_nuc) {
        return 0.0;
    }
    const double fade_in = std::min((r_km - r_nuc) / (0.25 * r_nuc), 1.0);
    const double end = pinwheel.turns * pinwheel.step_km();
    const double fade_out = std::clamp((end - r_km) / pinwheel.step_km(), 0.0, 1.0);
    const double q = r_nuc / r_km;
    return q * q * fade_in * fade_out;
}

PinwheelVisibleLight pinwheel_visible_light(const Pinwheel& pinwheel, double r_km, double star_visible_share)
{
    const double t = pinwheel_dust_temperature(pinwheel, r_km);
    PinwheelVisibleLight light;
    light.thermal = (1.0 - pinwheel.visible_albedo) * blackbody_band_share(t, kVisibleLoNm, kVisibleHiNm);
    light.scattered = pinwheel.visible_albedo * star_visible_share;
    return light;
}

double pinwheel_star_visible_share(const Scene& scene, const Pinwheel& pinwheel)
{
    double visible = 0.0;
    double total = 0.0;
    for (int star : {pinwheel.primary, pinwheel.secondary}) {
        const Body& b = scene.bodies[static_cast<size_t>(star)];
        visible += b.luminosity_solar * blackbody_band_share(b.temperature_k, kVisibleLoNm, kVisibleHiNm);
        total += b.luminosity_solar;
    }
    return total > 0.0 ? visible / total : 0.0;
}

} // namespace astraxis
