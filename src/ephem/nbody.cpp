#include "ephem/nbody.hpp"

#include <glm/geometric.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

void accelerations(const std::vector<NBodyParticle>& p, std::vector<glm::dvec3>& acc)
{
    acc.assign(p.size(), glm::dvec3(0.0));
    for (size_t i = 0; i < p.size(); ++i) {
        for (size_t j = i + 1; j < p.size(); ++j) {
            const glm::dvec3 d = p[j].position - p[i].position;
            const double r2 = glm::dot(d, d);
            const glm::dvec3 f = d / (r2 * std::sqrt(r2));
            acc[i] += f * p[j].gm;
            acc[j] -= f * p[i].gm;
        }
    }
}

void yoshida4_step(std::vector<NBodyParticle>& p, double dt, std::vector<glm::dvec3>& acc)
{
    const double cbrt2 = std::cbrt(2.0);
    const double w1 = 1.0 / (2.0 - cbrt2);
    const double w0 = -cbrt2 * w1;
    const double c[4] = {w1 / 2.0, (w0 + w1) / 2.0, (w0 + w1) / 2.0, w1 / 2.0};
    const double d[3] = {w1, w0, w1};

    for (int stage = 0; stage < 4; ++stage) {
        for (NBodyParticle& q : p) {
            q.position += q.velocity * (c[stage] * dt);
        }
        if (stage < 3) {
            accelerations(p, acc);
            for (size_t i = 0; i < p.size(); ++i) {
                p[i].velocity += acc[i] * (d[stage] * dt);
            }
        }
    }
}

void yoshida_step(std::vector<NBodyParticle>& p, double dt, int order, std::vector<glm::dvec3>& acc)
{
    if (order == 6) {
        const double z1 = 1.0 / (2.0 - std::pow(2.0, 0.2));
        const double z0 = -std::pow(2.0, 0.2) * z1;
        yoshida4_step(p, z1 * dt, acc);
        yoshida4_step(p, z0 * dt, acc);
        yoshida4_step(p, z1 * dt, acc);
    } else {
        yoshida4_step(p, dt, acc);
    }
}

EphemerisTable::Knot knot(const std::vector<NBodyParticle>& p, size_t i, double t, const NBodyOptions& options)
{
    glm::dvec3 cp(0.0);
    glm::dvec3 cv(0.0);
    if (!options.origin.empty()) {
        double gm = 0.0;
        for (size_t k : options.origin) {
            cp += p[k].position * p[k].gm;
            cv += p[k].velocity * p[k].gm;
            gm += p[k].gm;
        }
        cp /= gm;
        cv /= gm;
    }
    return {t, p[i].position - cp, p[i].velocity - cv};
}

void run(std::vector<NBodyParticle> p, double t0, double t_end, double step_s, const NBodyOptions& options,
         std::vector<std::vector<EphemerisTable::Knot>>& out)
{
    const double sign = t_end >= t0 ? 1.0 : -1.0;
    std::vector<glm::dvec3> acc;
    double t = t0;
    for (long long steps = 1; sign * (t_end - t) > 1e-6; ++steps) {
        const double dt = std::min(step_s, sign * (t_end - t));
        yoshida_step(p, sign * dt, options.order, acc);
        t += sign * dt;
        const bool last = sign * (t_end - t) <= 1e-6;
        for (size_t i = 0; i < p.size(); ++i) {
            const long long stride = options.stride.empty() ? 1 : std::max(1, options.stride[i]);
            if (last || steps % stride == 0) {
                out[i].push_back(knot(p, i, t, options));
            }
        }
    }
}

} // namespace

std::vector<std::vector<EphemerisTable::Knot>> integrate_nbody(const std::vector<NBodyParticle>& particles,
                                                               double t_epoch, double t_begin, double t_end,
                                                               double step_s, const NBodyOptions& options)
{
    std::vector<std::vector<EphemerisTable::Knot>> backward(particles.size());
    std::vector<std::vector<EphemerisTable::Knot>> forward(particles.size());
    run(particles, t_epoch, std::min(t_begin, t_epoch), step_s, options, backward);
    run(particles, t_epoch, std::max(t_end, t_epoch), step_s, options, forward);

    std::vector<std::vector<EphemerisTable::Knot>> out(particles.size());
    for (size_t i = 0; i < particles.size(); ++i) {
        out[i].reserve(backward[i].size() + 1 + forward[i].size());
        out[i].assign(backward[i].rbegin(), backward[i].rend());
        out[i].push_back(knot(particles, i, t_epoch, options));
        out[i].insert(out[i].end(), forward[i].begin(), forward[i].end());
    }
    return out;
}

double nbody_energy(const std::vector<NBodyParticle>& p)
{
    // Energy in units where particle "mass" is gm (consistent for conservation checks).
    double e = 0.0;
    for (size_t i = 0; i < p.size(); ++i) {
        e += 0.5 * p[i].gm * glm::dot(p[i].velocity, p[i].velocity);
        for (size_t j = i + 1; j < p.size(); ++j) {
            e -= p[i].gm * p[j].gm / glm::length(p[j].position - p[i].position);
        }
    }
    return e;
}

} // namespace astraxis
