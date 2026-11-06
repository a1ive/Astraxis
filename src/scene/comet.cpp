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

namespace {

glm::dmat3 gravity_gradient(const glm::dvec3& r, double gm)
{
    const double d = glm::length(r);
    const glm::dvec3 u = r / d;
    const double k = gm / (d * d * d);
    glm::dmat3 g(0.0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            g[i][j] = k * (3.0 * u[i] * u[j] - (i == j ? 1.0 : 0.0));
        }
    }
    return g;
}

// One grain's deviation for one source term: d'' = G d + f.
struct Deviation {
    glm::dvec3 d{0.0};
    glm::dvec3 v{0.0};
};

Deviation rk4_step(const Deviation& y, const glm::dmat3& g0, const glm::dmat3& g1, const glm::dmat3& g2, const glm::dvec3& f0,
                   const glm::dvec3& f1, const glm::dvec3& f2, double h)
{
    const glm::dvec3 k1d = y.v;
    const glm::dvec3 k1v = g0 * y.d + f0;
    const glm::dvec3 k2d = y.v + 0.5 * h * k1v;
    const glm::dvec3 k2v = g1 * (y.d + 0.5 * h * k1d) + f1;
    const glm::dvec3 k3d = y.v + 0.5 * h * k2v;
    const glm::dvec3 k3v = g1 * (y.d + 0.5 * h * k2d) + f1;
    const glm::dvec3 k4d = y.v + h * k3v;
    const glm::dvec3 k4v = g2 * (y.d + h * k3d) + f2;
    return {y.d + h / 6.0 * (k1d + 2.0 * k2d + 2.0 * k3d + k4d), y.v + h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)};
}

} // namespace

void BoundDustTail::build(const Scene& scene, int body)
{
    if (m_body == body && !m_times.empty()) {
        return;
    }
    m_body = body;
    m_last_tdb = std::numeric_limits<double>::quiet_NaN();
    m_times.clear();
    m_release.clear();
    m_first_row.clear();
    m_rows.clear();
    m_seeds.clear();
    const Body& comet_body = scene.bodies[static_cast<size_t>(body)];
    const MotionSource& motion = *comet_body.motion;
    const int parent = comet_body.parent;
    const double gm_parent = scene.bodies[static_cast<size_t>(parent)].gm_km3_s2;
    const double gm_sun = scene.sun_gm();

    // The span of the comet's table (a fragment that split off later starts later).
    const double from = comet_body.comet.bound_tail_from_tdb;
    const double to = comet_body.comet.bound_tail_to_tdb;
    constexpr double kSearchS = 4.0 * kDaysPerJulianYear * kSecondsPerDay;
    double start = from;
    if (!motion.valid_at(start)) {
        double lo = from;
        double hi = from + kSearchS;
        while (hi - lo > 1.0) {
            const double mid = 0.5 * (lo + hi);
            (motion.valid_at(mid) ? hi : lo) = mid;
        }
        start = hi;
    }
    double lo = start;
    double hi = start + kSearchS;
    while (hi - lo > 1.0) {
        const double mid = 0.5 * (lo + hi);
        (motion.valid_at(mid) ? lo : hi) = mid;
    }
    const double end = lo;
    if (end <= start) {
        return;
    }

    // Steps by the path: a fiftieth of its time scale, at most 3 days. The sources at each step's ends and middle.
    const auto relative = [&](double t) {
        return scene.icrf_state_at(body, t).position - scene.icrf_state_at(parent, t).position;
    };
    const auto sources = [&](double t, glm::dmat3* g, glm::dvec3* f) {
        const glm::dvec3 r = relative(t);
        const glm::dvec3 from_sun = scene.icrf_state_at(body, t).position - scene.sun_icrf_state_at(t).position;
        *g = gravity_gradient(r, gm_parent) + gravity_gradient(from_sun, gm_sun);
        const double s = glm::length(from_sun);
        *f = gm_sun * from_sun / (s * s * s); // per unit beta: away from the Sun
    };
    std::vector<glm::dmat3> g_ends;
    std::vector<glm::dmat3> g_mids;
    std::vector<glm::dvec3> f_ends;
    std::vector<glm::dvec3> f_mids;
    double t = start;
    for (;;) {
        glm::dmat3 g;
        glm::dvec3 f;
        sources(t, &g, &f);
        m_times.push_back(t);
        g_ends.push_back(g);
        f_ends.push_back(f);
        if (t >= end) {
            break;
        }
        const State s = motion.eval(t);
        const double speed = std::max(glm::length(s.velocity), 1e-3);
        const double h = std::clamp(0.02 * glm::length(s.position) / speed, 120.0, 3.0 * kSecondsPerDay);
        const double next = std::min(t + h, end);
        sources(0.5 * (t + next), &g, &f);
        g_mids.push_back(g);
        f_mids.push_back(f);
        t = next;
    }

    // One set of rows per time of release, each from the step at or before it.
    const double first_release = std::max(from, start);
    const double last_release = std::max(to, start + 30.0 * kSecondsPerDay);
    const size_t n = m_times.size();
    m_rows.assign(static_cast<size_t>(kReleases) * n, Row{});
    for (int k = 0; k < kReleases; ++k) {
        const double when = first_release + (last_release - first_release) * k / (kReleases - 1);
        const size_t first = static_cast<size_t>(std::upper_bound(m_times.begin(), m_times.end(), when) - m_times.begin()) - 1;
        m_release.push_back(m_times[first]);
        m_first_row.push_back(first);
        Deviation columns[4];
        columns[0].v = glm::dvec3(1.0, 0.0, 0.0);
        columns[1].v = glm::dvec3(0.0, 1.0, 0.0);
        columns[2].v = glm::dvec3(0.0, 0.0, 1.0);
        for (size_t i = first; i + 1 < n; ++i) {
            const double h = m_times[i + 1] - m_times[i];
            for (int c = 0; c < 3; ++c) {
                columns[c] = rk4_step(columns[c], g_ends[i], g_mids[i], g_ends[i + 1], glm::dvec3(0.0), glm::dvec3(0.0), glm::dvec3(0.0), h);
            }
            columns[3] = rk4_step(columns[3], g_ends[i], g_mids[i], g_ends[i + 1], f_ends[i], f_mids[i], f_ends[i + 1], h);
            Row& out = m_rows[static_cast<size_t>(k) * n + i + 1];
            for (int c = 0; c < 3; ++c) {
                for (int e = 0; e < 3; ++e) {
                    out.p[c * 3 + e] = static_cast<float>(columns[c].d[e]);
                }
                out.b[c] = static_cast<float>(columns[3].d[c]);
            }
        }
    }

    // The grains, the same for all time.
    const double log_min = std::log(kBetaMin);
    const double log_step = (std::log(kBetaMax) - log_min) / kGrains;
    double total = 0.0;
    m_seeds.resize(kGrains);
    for (int k = 0; k < kGrains; ++k) {
        Seed& seed = m_seeds[static_cast<size_t>(k)];
        seed.release = k % kReleases;
        seed.beta = std::exp(log_min + log_step * (k + grain_random(body, k, 0)));
        const double cos_t = 2.0 * grain_random(body, k, 1) - 1.0;
        const double phi = kTwoPi * grain_random(body, k, 2);
        const double sin_t = std::sqrt(std::max(0.0, 1.0 - cos_t * cos_t));
        const double speed = kSpeedMinKmS + (kSpeedMaxKmS - kSpeedMinKmS) * grain_random(body, k, 3);
        seed.velocity = speed * glm::dvec3(sin_t * std::cos(phi), sin_t * std::sin(phi), cos_t);
        seed.weight = seed.beta;
        seed.random = grain_random(body, k, 4);
        total += seed.weight;
    }
    for (Seed& seed : m_seeds) {
        seed.weight /= total;
    }
}

glm::dvec3 BoundDustTail::offset(int release, double beta, const glm::dvec3& velocity, double t_tdb) const
{
    if (m_times.empty() || t_tdb <= m_release[static_cast<size_t>(release)]) {
        return glm::dvec3(0.0);
    }
    const size_t hi = std::clamp<size_t>(static_cast<size_t>(std::upper_bound(m_times.begin(), m_times.end(), t_tdb) - m_times.begin()), 1,
                                         m_times.size() - 1);
    const double s = std::clamp((t_tdb - m_times[hi - 1]) / (m_times[hi] - m_times[hi - 1]), 0.0, 1.0);
    const Row& a = row(release, hi - 1);
    const Row& b = row(release, hi);
    glm::dvec3 out(0.0);
    for (int e = 0; e < 3; ++e) {
        for (int c = 0; c < 3; ++c) {
            out[e] += velocity[c] * (a.p[c * 3 + e] + s * (b.p[c * 3 + e] - a.p[c * 3 + e]));
        }
        out[e] += beta * (a.b[e] + s * (b.b[e] - a.b[e]));
    }
    return out;
}

void BoundDustTail::update(const Scene& scene, int body, double t_tdb)
{
    if (body == m_body && t_tdb == m_last_tdb) {
        return;
    }
    build(scene, body);
    m_last_tdb = t_tdb;
    m_grains.clear();
    if (m_times.empty() || t_tdb <= m_times.front() || t_tdb >= m_times.back()) {
        return;
    }
    const glm::dvec3 here = scene.icrf_state_at(body, t_tdb).position;
    const double from_parent = glm::length(here - scene.icrf_state_at(scene.bodies[static_cast<size_t>(body)].parent, t_tdb).position);
    for (const Seed& seed : m_seeds) {
        const double released = m_release[static_cast<size_t>(seed.release)];
        if (t_tdb <= released) {
            continue;
        }
        const glm::dvec3 d = offset(seed.release, seed.beta, seed.velocity, t_tdb);
        const double reach = glm::length(d);
        if (reach > 0.25 * from_parent) {
            continue;
        }
        DustGrain out;
        out.position = here + d;
        out.weight = seed.weight;
        out.age_s = t_tdb - released;
        out.ejection_km_s = glm::length(seed.velocity);
        out.random = seed.random;
        m_grains.push_back(out);
    }
}

} // namespace astraxis
