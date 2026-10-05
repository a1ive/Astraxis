#include "ephem/kepler.hpp"

#include "core/math.hpp"

#include <algorithm>
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

namespace {

// Stumpff functions C(z) and S(z), with series near z = 0.
double stumpff_c(double z)
{
    if (z > 1e-6) {
        return (1.0 - std::cos(std::sqrt(z))) / z;
    }
    if (z < -1e-6) {
        return (std::cosh(std::sqrt(-z)) - 1.0) / -z;
    }
    return 0.5 - z / 24.0 + z * z / 720.0;
}

double stumpff_s(double z)
{
    if (z > 1e-6) {
        const double s = std::sqrt(z);
        return (s - std::sin(s)) / (s * s * s);
    }
    if (z < -1e-6) {
        const double s = std::sqrt(-z);
        return (std::sinh(s) - s) / (s * s * s);
    }
    return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0;
}

} // namespace

State propagate_kepler(const State& s, double mu, double dt, double* universal_anomaly)
{
    const double r0 = glm::length(s.position);
    if (dt == 0.0 || r0 <= 0.0) {
        return s;
    }
    if (mu <= 0.0) {
        return {s.position + s.velocity * dt, s.velocity};
    }
    const double sqrt_mu = std::sqrt(mu);
    const double vr0 = glm::dot(s.position, s.velocity) / r0;
    const double alpha = 2.0 / r0 - glm::dot(s.velocity, s.velocity) / mu; // 1 / a

    // Newton iteration for the universal anomaly x (Curtis Algorithm 3.3).
    double x = sqrt_mu * std::abs(alpha) * dt;
    if (std::abs(alpha) < 1e-12) {
        x = sqrt_mu * dt / r0; // near-parabolic start
    }
    if (universal_anomaly && *universal_anomaly != 0.0) {
        x = *universal_anomaly;
    }
    for (int iter = 0; iter < 50; ++iter) {
        const double x2 = x * x;
        const double z = alpha * x2;
        const double c = stumpff_c(z);
        const double sz = stumpff_s(z);
        const double f = r0 * vr0 / sqrt_mu * x2 * c + (1.0 - alpha * r0) * x2 * x * sz + r0 * x - sqrt_mu * dt;
        const double df = r0 * vr0 / sqrt_mu * x * (1.0 - z * sz) + (1.0 - alpha * r0) * x2 * c + r0;
        const double dx = f / df;
        x -= dx;
        if (std::abs(dx) <= 1e-12 * std::max(1.0, std::abs(x))) {
            break;
        }
    }

    if (universal_anomaly) {
        *universal_anomaly = x;
    }

    // Lagrange coefficients (Curtis Algorithm 3.4).
    const double x2 = x * x;
    const double z = alpha * x2;
    const double c = stumpff_c(z);
    const double sz = stumpff_s(z);
    const double f = 1.0 - x2 / r0 * c;
    const double g = dt - x2 * x / sqrt_mu * sz;
    State out;
    out.position = f * s.position + g * s.velocity;
    const double r = glm::length(out.position);
    const double fdot = sqrt_mu / (r * r0) * (z * x * sz - x);
    const double gdot = 1.0 - x2 / r * c;
    out.velocity = fdot * s.position + gdot * s.velocity;
    return out;
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
