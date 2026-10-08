#include "scene/comet.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/kepler.hpp"
#include "scene/scene.hpp"

#include <glm/geometric.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

// Deterministic random numbers in [0, 1) for a grain (SplitMix64).
double grain_random(int64_t slot, int k, int which)
{
    uint64_t z = static_cast<uint64_t>(slot) * 0x9E3779B97F4A7C15ull + static_cast<uint64_t>(k) * 0xBF58476D1CE4E5B9ull +
                 static_cast<uint64_t>(which) * 0x94D049BB133111EBull + 0x2545F4914F6CDD1Dull;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z ^= z >> 31;
    return static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0);
}

double smoothstep(double edge0, double edge1, double x)
{
    const double t = std::clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

} // namespace

double comet_total_magnitude(double m1, double k1, double r_au)
{
    return m1 + k1 * std::log10(r_au);
}

glm::dvec3 ion_tail_direction(const glm::dvec3& helio_position, const glm::dvec3& helio_velocity,
                              double solar_wind_km_s)
{
    return glm::normalize(solar_wind_km_s * glm::normalize(helio_position) - helio_velocity);
}

State dust_grain_state(const State& comet_at_ejection, const glm::dvec3& ejection_velocity, double beta,
                       double gm_sun, double dt, double* universal_anomaly)
{
    const State start{comet_at_ejection.position, comet_at_ejection.velocity + ejection_velocity};
    return propagate_kepler(start, (1.0 - beta) * gm_sun, dt, universal_anomaly);
}

DustTail::Slot DustTail::make_slot(const Scene& scene, int body, int64_t index, double slot_s) const
{
    const Body& comet_body = scene.bodies[static_cast<size_t>(body)];
    const Body::Comet& comet = comet_body.comet;
    const double log_beta_min = std::log(comet.dust_beta_min);
    const double log_beta_step = (std::log(comet.dust_beta_max) - log_beta_min) / kBetas;
    const double production_exponent = std::max(0.4 * comet.k1 - 2.0, 0.0);

    Slot slot;
    slot.index = index;
    slot.grains.resize(kBetas * kDirections);
    for (int k = 0; k < kBetas * kDirections; ++k) {
        Grain& g = slot.grains[static_cast<size_t>(k)];
        g.emitted_tdb = (static_cast<double>(index) + grain_random(index, k, 0)) * slot_s;
        if (!comet_body.motion || !comet_body.motion->valid_at(g.emitted_tdb)) {
            continue;
        }
        g.beta = std::exp(log_beta_min + log_beta_step * (k / kDirections + grain_random(index, k, 1)));
        const State sun = scene.sun_icrf_state_at(g.emitted_tdb);
        const State c = scene.icrf_state_at(body, g.emitted_tdb);
        const glm::dvec3 r = c.position - sun.position;
        const double r_au = glm::length(r) / kAuKm;
        // Toward the sunlit hemisphere: the sun direction plus a random unit vector.
        const double cos_t = 2.0 * grain_random(index, k, 2) - 1.0;
        const double phi = kTwoPi * grain_random(index, k, 3);
        const double sin_t = std::sqrt(std::max(0.0, 1.0 - cos_t * cos_t));
        const glm::dvec3 random(sin_t * std::cos(phi), sin_t * std::sin(phi), cos_t);
        const glm::dvec3 dir = glm::normalize(-glm::normalize(r) + random * 0.999);
        const double speed = comet.dust_speed_km_s * std::sqrt(g.beta / r_au);
        g.start = {r, c.velocity - sun.velocity + speed * dir};
        g.ejection_km_s = speed;
        g.random = grain_random(index, k, 4);
        g.production = std::pow(r_au, -production_exponent) * std::sqrt(g.beta);
        g.valid = true;
    }
    return slot;
}

void DustTail::update(const Scene& scene, int body, double t_tdb)
{
    if (body == m_body && t_tdb == m_last_tdb) {
        return;
    }
    m_grains.clear();
    m_last_tdb = t_tdb;
    const Body& comet_body = scene.bodies[static_cast<size_t>(body)];
    const Body::Comet& comet = comet_body.comet;
    if (body != m_body) {
        m_slots.clear();
        m_body = body;
    }
    const double span_s = comet.dust_tail_days * kSecondsPerDay;
    const double slot_s = span_s / kSlots;
    const int64_t first = static_cast<int64_t>(std::floor((t_tdb - span_s) / slot_s));
    const int64_t last = static_cast<int64_t>(std::floor(t_tdb / slot_s));

    // Keep the cached slots that are still in the span; add the missing ones.
    if (!m_slots.empty() && (m_slots.back().index < first || m_slots.front().index > last)) {
        m_slots.clear();
    }
    while (!m_slots.empty() && m_slots.front().index < first) {
        m_slots.pop_front();
    }
    while (!m_slots.empty() && m_slots.back().index > last) {
        m_slots.pop_back();
    }
    if (m_slots.empty()) {
        m_slots.push_back(make_slot(scene, body, last, slot_s));
    }
    while (m_slots.front().index > first) {
        m_slots.push_front(make_slot(scene, body, m_slots.front().index - 1, slot_s));
    }
    while (m_slots.back().index < last) {
        m_slots.push_back(make_slot(scene, body, m_slots.back().index + 1, slot_s));
    }

    const double gm_sun = scene.sun_gm(); // the parent may be the planet a scene is centered on
    const State sun = scene.sun_icrf_state_at(t_tdb);
    const double r_comet = glm::length(scene.icrf_state_at(body, t_tdb).position - sun.position);
    double total = 0.0;
    for (Slot& slot : m_slots) {
        for (Grain& g : slot.grains) {
            const double age = t_tdb - g.emitted_tdb;
            if (!g.valid || age <= 0.0 || age > span_s) {
                continue;
            }
            const State s = dust_grain_state(g.start, glm::dvec3(0.0), g.beta, gm_sun, age, &g.anomaly);
            const double fade = 1.0 - smoothstep(0.9 * span_s, span_s, age);
            const double weight = g.production * fade;
            total += weight;
            const double r_grain = glm::length(s.position);
            DustGrain out;
            out.position = sun.position + s.position;
            out.weight = weight * (r_comet * r_comet) / (r_grain * r_grain);
            out.age_s = age;
            out.ejection_km_s = g.ejection_km_s;
            out.random = g.random;
            m_grains.push_back(out);
        }
    }
    if (total > 0.0) {
        for (DustGrain& g : m_grains) {
            g.weight /= total;
        }
    }
}

} // namespace astraxis
