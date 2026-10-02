#include "ephem/kepler.hpp"

#include "core/math.hpp"

#include <cmath>

namespace astraxis {

double solve_kepler(double mean_anomaly, double e)
{
    const double m = wrap_pi(mean_anomaly);
    double E = e < 0.8 ? m : (m < 0.0 ? -kPi : kPi);
    for (int iter = 0; iter < 50; ++iter) {
        const double f = E - e * std::sin(E) - m;
        const double dE = f / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-14) {
            break;
        }
    }
    return E;
}

glm::dmat3 perifocal_to_reference(const KeplerElements& el)
{
    return rotation_z(el.node) * rotation_x(el.i) * rotation_z(el.arg_peri);
}

State kepler_state(const KeplerElements& el, double mean_motion)
{
    const double E = solve_kepler(el.mean_anomaly, el.e);
    const double cos_E = std::cos(E);
    const double sin_E = std::sin(E);
    const double b_over_a = std::sqrt(1.0 - el.e * el.e);

    const glm::dvec3 pos(el.a * (cos_E - el.e), el.a * b_over_a * sin_E, 0.0);
    const double k = mean_motion * el.a / (1.0 - el.e * cos_E);
    const glm::dvec3 vel(-k * sin_E, k * b_over_a * cos_E, 0.0);

    const glm::dmat3 r = perifocal_to_reference(el);
    return {r * pos, r * vel};
}

bool sample_osculating_ellipse(const State& s, double mu, int count, std::vector<glm::dvec3>& out)
{
    out.clear();
    const glm::dvec3& r = s.position;
    const glm::dvec3& v = s.velocity;
    const double rn = glm::length(r);
    const glm::dvec3 h = glm::cross(r, v);
    const double hn = glm::length(h);
    if (mu <= 0.0 || rn <= 0.0 || hn <= 0.0 || count < 2) {
        return false;
    }
    const double energy = 0.5 * glm::dot(v, v) - mu / rn;
    if (energy >= 0.0) {
        return false; // parabolic or hyperbolic
    }
    const double a = -mu / (2.0 * energy);
    const glm::dvec3 e_vec = ((glm::dot(v, v) - mu / rn) * r - glm::dot(r, v) * v) / mu;
    const double e = glm::length(e_vec);
    if (e >= 1.0) {
        return false;
    }

    // Perifocal basis: P toward periapsis (any in-plane axis if circular), Q 90 deg ahead.
    const glm::dvec3 w = h / hn;
    const glm::dvec3 p = e > 1e-10 ? e_vec / e : r / rn;
    const glm::dvec3 q = glm::cross(w, p);
    const double b = a * std::sqrt(1.0 - e * e);
    const double start_E = std::atan2(glm::dot(r, q) / b, glm::dot(r, p) / a + e);

    out.reserve(static_cast<size_t>(count));
    for (int k = 0; k < count; ++k) {
        const double E = start_E - kTwoPi * static_cast<double>(k) / static_cast<double>(count - 1);
        out.push_back(p * (a * (std::cos(E) - e)) + q * (b * std::sin(E)));
    }
    return true;
}

void sample_ellipse(const KeplerElements& el, double start_E, int count, std::vector<glm::dvec3>& out)
{
    out.clear();
    if (count < 2) {
        return;
    }
    out.reserve(static_cast<size_t>(count));

    const glm::dmat3 r = perifocal_to_reference(el);
    const double b = el.a * std::sqrt(1.0 - el.e * el.e);
    for (int k = 0; k < count; ++k) {
        const double E = start_E - kTwoPi * static_cast<double>(k) / static_cast<double>(count - 1);
        out.push_back(r * glm::dvec3(el.a * (std::cos(E) - el.e), b * std::sin(E), 0.0));
    }
}

} // namespace astraxis
