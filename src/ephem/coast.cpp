#include "ephem/coast.hpp"

#include "ephem/kepler.hpp"

#include <glm/geometric.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

bool CoastGravity::point_mass() const
{
    const bool zonal_terms = std::any_of(zonal.begin(), zonal.end(), [](double j) { return j != 0.0; });
    return !zonal_terms && tides.empty();
}

glm::dvec3 CoastGravity::acceleration(const glm::dvec3& r, double t_tdb) const
{
    const double rn = glm::length(r);
    glm::dvec3 a = -mu / (rn * rn * rn) * r;
    // Zonal harmonics: the potential mu / r (1 - sum J_n (R / r)^n P_n(s)),
    // s = u . pole, has the gradient
    //   -mu J_n R^n / r^(n+2) (P_n'(s) (pole - s u) - (n + 1) P_n(s) u).
    if (zonal.size() > 2 && radius_km > 0.0) {
        const glm::dvec3 u = r / rn;
        const double s = glm::dot(u, pole);
        // Legendre polynomials and their derivatives by recurrence.
        double p_prev = 1.0, p = s; // P_0, P_1
        double dp = 1.0;            // P_1'
        double scale = mu / (rn * rn) * (radius_km / rn);
        for (size_t n = 2; n < zonal.size(); ++n) {
            const double nd = static_cast<double>(n);
            const double p_next = ((2.0 * nd - 1.0) * s * p - (nd - 1.0) * p_prev) / nd;
            const double dp_next = nd * p + s * dp; // P_n' = n P_(n-1) + s P_(n-1)'
            p_prev = p;
            p = p_next;
            dp = dp_next;
            scale *= radius_km / rn;
            if (zonal[n] != 0.0) {
                a -= zonal[n] * scale * (dp * (pole - s * u) - (nd + 1.0) * p * u);
            }
        }
    }
    for (const Tide& tide : tides) {
        const glm::dvec3 d = tide.position(t_tdb);
        const glm::dvec3 to = d - r;
        const double to_n = glm::length(to);
        const double d_n = glm::length(d);
        a += tide.gm * (to / (to_n * to_n * to_n) - d / (d_n * d_n * d_n));
    }
    return a;
}

namespace {

// Step for the local dynamical time (see coast).
double coast_step(const CoastGravity& gravity, const State& s)
{
    const double r = glm::length(s.position);
    const double v = glm::length(s.velocity);
    double scale = std::sqrt(r * r * r / gravity.mu);
    if (v > 0.0) {
        scale = std::min(scale, r / v);
    }
    return scale / 500.0;
}

State rk4_step(const CoastGravity& gravity, const State& s, double t, double h)
{
    const glm::dvec3 a1 = gravity.acceleration(s.position, t);
    const glm::dvec3 r2 = s.position + 0.5 * h * s.velocity;
    const glm::dvec3 v2 = s.velocity + 0.5 * h * a1;
    const glm::dvec3 a2 = gravity.acceleration(r2, t + 0.5 * h);
    const glm::dvec3 r3 = s.position + 0.5 * h * v2;
    const glm::dvec3 v3 = s.velocity + 0.5 * h * a2;
    const glm::dvec3 a3 = gravity.acceleration(r3, t + 0.5 * h);
    const glm::dvec3 r4 = s.position + h * v3;
    const glm::dvec3 v4 = s.velocity + h * a3;
    const glm::dvec3 a4 = gravity.acceleration(r4, t + h);
    State out;
    out.position = s.position + h / 6.0 * (s.velocity + 2.0 * v2 + 2.0 * v3 + v4);
    out.velocity = s.velocity + h / 6.0 * (a1 + 2.0 * a2 + 2.0 * a3 + a4);
    return out;
}

// Integrates from t over dt (either sign), calling visit(t, state) after each step.
template <typename Visit>
State integrate(const CoastGravity& gravity, State s, double t, double dt, Visit&& visit)
{
    const double end = t + dt;
    const double dir = dt < 0.0 ? -1.0 : 1.0;
    while (dir * (end - t) > 0.0) {
        double h = coast_step(gravity, s);
        // Land on the end without a sliver of a step.
        if (h >= dir * (end - t) * 0.999) {
            h = dir * (end - t);
        }
        s = rk4_step(gravity, s, t, dir * h);
        t = h == dir * (end - t) ? end : t + dir * h;
        visit(t, s);
    }
    return s;
}

} // namespace

State coast(const CoastGravity& gravity, const State& s, double t_tdb, double dt)
{
    if (gravity.point_mass()) {
        return propagate_kepler(s, gravity.mu, dt);
    }
    return integrate(gravity, s, t_tdb, dt, [](double, const State&) {});
}

std::vector<EphemerisTable::Knot> integrate_coast(const CoastGravity& gravity, const State& state, double epoch_tdb,
                                                  double t0, double t1)
{
    std::vector<EphemerisTable::Knot> back;
    integrate(gravity, state, epoch_tdb, t0 - epoch_tdb,
              [&](double t, const State& s) { back.push_back({t, s.position, s.velocity}); });
    std::vector<EphemerisTable::Knot> knots(back.rbegin(), back.rend());
    knots.push_back({epoch_tdb, state.position, state.velocity});
    integrate(gravity, state, epoch_tdb, t1 - epoch_tdb,
              [&](double t, const State& s) { knots.push_back({t, s.position, s.velocity}); });
    return knots;
}

CoastMotion::CoastMotion(const State& state, const CoastGravity& gravity, double epoch_tdb, double t0, double t1)
    : m_state(state)
    , m_mu(gravity.mu)
    , m_epoch(epoch_tdb)
    , m_t0(t0)
    , m_t1(t1)
{
    if (!gravity.point_mass()) {
        m_table = std::make_shared<EphemerisTable>();
        m_table->set_knots(
            integrate_coast(gravity, state, epoch_tdb, std::min(t0, epoch_tdb), std::max(t1, epoch_tdb)));
    }
}

State CoastMotion::eval(double t_tdb) const
{
    return m_table ? m_table->eval(t_tdb) : propagate_kepler(m_state, m_mu, t_tdb - m_epoch);
}

bool CoastMotion::sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const
{
    return sample_osculating_ellipse(eval(t_tdb), m_mu, count, out);
}

void CoastMotion::history_times(double t0, double t1, int max_points, std::vector<double>& out) const
{
    if (m_table) {
        m_table->history_times(t0, t1, max_points, out);
    } else {
        MotionSource::history_times(t0, t1, max_points, out);
    }
}

std::vector<double> unexplained_delta_v(const MotionSource& motion, const CoastGravity& gravity, double t0,
                                        double duration, int steps)
{
    steps = std::max(steps, 1);
    std::vector<double> times;
    motion.history_times(t0, t0 + duration, steps + 1, times);
    std::vector<double> sums{0.0};
    State s = motion.eval(times.front());
    for (size_t i = 1; i < times.size(); ++i) {
        const State next = motion.eval(times[i]);
        const State coasted = coast(gravity, s, times[i - 1], times[i] - times[i - 1]);
        sums.push_back(sums.back() + glm::length(next.velocity - coasted.velocity));
        s = next;
    }

    std::vector<double> total;
    size_t i = 0;
    for (int k = 0; k <= steps; ++k) {
        const double t = t0 + duration * k / steps;
        while (i + 2 < times.size() && times[i + 1] < t) {
            ++i;
        }
        const double span = times[i + 1] - times[i];
        const double x = span > 0.0 ? std::clamp((t - times[i]) / span, 0.0, 1.0) : 1.0;
        total.push_back(sums[i] + (sums[i + 1] - sums[i]) * x);
    }
    return total;
}

} // namespace astraxis
