#include "ephem/kerr_orbit.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr double kC = kSpeedOfLightKmS;

} // namespace

KerrOrbitMotion::KerrOrbitMotion(const KerrOrbitSetup& setup)
    : m_setup(setup)
    , m_kerr(setup.gm_km3_s2 / (kC * kC), setup.spin * setup.gm_km3_s2 / (kC * kC))
    , m_plane(visual_orbit_frame(setup.elements))
    , m_mass_km(setup.gm_km3_s2 / (kC * kC))
{
}

KerrOrbitMotion::KerrOrbitMotion(const KerrOrbitSetup& setup, double t_begin, double t_end)
    : KerrOrbitMotion(setup)
{
    const PlaneState s0 = initial_state(setup.elements_epoch_tdb);
    std::vector<EphemerisTable::Knot> backward;
    std::vector<EphemerisTable::Knot> forward;
    integrate(s0, std::min(t_begin, s0.t), backward);
    integrate(s0, std::max(t_end, s0.t), forward);

    std::vector<EphemerisTable::Knot> knots(backward.rbegin(), backward.rend());
    knots.push_back(knot(s0));
    knots.insert(knots.end(), forward.begin(), forward.end());
    m_table.set_knots(std::move(knots));
}

std::unique_ptr<KerrOrbitMotion> KerrOrbitMotion::rolling(const KerrOrbitSetup& setup, double window_days)
{
    std::unique_ptr<KerrOrbitMotion> m(new KerrOrbitMotion(setup));
    m->m_rolling = true;
    m->m_window_s = window_days * kSecondsPerDay;
    return m;
}

KerrOrbitMotion::PlaneState KerrOrbitMotion::initial_state(double epoch_tdb) const
{
    const State s = visual_orbit_state(m_setup.elements, m_setup.gm_km3_s2, epoch_tdb);
    const glm::dvec3 p = m_plane[0];
    const glm::dvec3 q = m_plane[1];
    PlaneState out;
    out.t = epoch_tdb;
    out.x = glm::dvec2(glm::dot(s.position, p), glm::dot(s.position, q));
    const glm::dvec2 v(glm::dot(s.velocity, p) / kC, glm::dot(s.velocity, q) / kC);
    out.u = m_kerr.u_from_velocity(out.x, v);
    return out;
}

EphemerisTable::Knot KerrOrbitMotion::knot(const PlaneState& s) const
{
    glm::dvec2 dxdt;
    glm::dvec2 dudt;
    m_kerr.derivatives(s.x, s.u, dxdt, dudt);
    const glm::dvec3 p = m_plane[0];
    const glm::dvec3 q = m_plane[1];
    return {s.t, p * s.x.x + q * s.x.y, (p * dxdt.x + q * dxdt.y) * kC};
}

KerrOrbitMotion::PlaneState KerrOrbitMotion::integrate(PlaneState s, double t_end,
                                                       std::vector<EphemerisTable::Knot>& out) const
{
    const double sign = t_end >= s.t ? 1.0 : -1.0;
    const double plunge = 2.0 * m_kerr.horizon_radius();
    while (sign * (t_end - s.t) > 1e-6) {
        const double r = glm::length(s.x);
        if (r < plunge) {
            break; // captured: the trajectory ends here
        }
        // Adaptive step ~ local dynamical time (geometric units: km of time).
        double dt_s = m_setup.step_eta * std::sqrt(r * r * r / m_mass_km) / kC;
        dt_s = std::min(dt_s, sign * (t_end - s.t));
        m_kerr.step(s.x, s.u, sign * dt_s * kC);
        s.t += sign * dt_s;
        out.push_back(knot(s));
    }
    return s;
}

void KerrOrbitMotion::restart(double t_tdb) const
{
    // Hypothetical orbits: re-seed from the elements, keeping the periastron
    // phase (an integer number of Keplerian periods away).
    const double a = m_setup.elements.a_km;
    const double period = kTwoPi * std::sqrt(a * a * a / m_setup.gm_km3_s2);
    KerrOrbitSetup shifted = m_setup;
    shifted.elements.t_peri_tdb += std::round((t_tdb - m_setup.elements.t_peri_tdb) / period) * period;
    const KerrOrbitMotion seed(shifted);

    const PlaneState s0 = seed.initial_state(t_tdb);
    const double half = 0.5 * m_window_s;
    std::vector<EphemerisTable::Knot> backward;
    std::vector<EphemerisTable::Knot> forward;
    m_front = integrate(s0, t_tdb - half, backward);
    m_back = integrate(s0, t_tdb + half, forward);
    m_front_valid = m_back_valid = true;

    std::vector<EphemerisTable::Knot> knots(backward.rbegin(), backward.rend());
    knots.push_back(knot(s0));
    knots.insert(knots.end(), forward.begin(), forward.end());
    m_table.set_knots(std::move(knots));
}

void KerrOrbitMotion::ensure_covered(double t_tdb) const
{
    if (!m_rolling || m_table.covers(t_tdb)) {
        return;
    }
    auto& knots = m_table.mutable_knots();
    const double chunk = 0.25 * m_window_s;

    if (!knots.empty() && m_back_valid && t_tdb > m_table.end() && t_tdb - m_table.end() < m_window_s) {
        std::vector<EphemerisTable::Knot> forward;
        m_back = integrate(m_back, t_tdb + chunk, forward);
        knots.insert(knots.end(), forward.begin(), forward.end());
        // Drop knots far behind, keeping the window (the front state moves with them).
        const auto keep = std::lower_bound(knots.begin(), knots.end(), t_tdb - m_window_s,
                                           [](const EphemerisTable::Knot& k, double v) { return k.t < v; });
        if (keep != knots.begin()) {
            knots.erase(knots.begin(), keep);
            m_front_valid = false;
        }
    } else if (!knots.empty() && m_front_valid && t_tdb < m_table.start() && m_table.start() - t_tdb < m_window_s) {
        std::vector<EphemerisTable::Knot> backward;
        m_front = integrate(m_front, t_tdb - chunk, backward);
        knots.insert(knots.begin(), backward.rbegin(), backward.rend());
        const auto keep = std::upper_bound(knots.begin(), knots.end(), t_tdb + m_window_s,
                                           [](double v, const EphemerisTable::Knot& k) { return v < k.t; });
        if (keep != knots.end()) {
            knots.erase(keep, knots.end());
            m_back_valid = false;
        }
    } else {
        // Far jump, or the needed end was trimmed: start over around t.
        restart(t_tdb);
    }
}

State KerrOrbitMotion::eval(double t_tdb) const
{
    ensure_covered(t_tdb);
    return m_table.eval(t_tdb);
}

bool KerrOrbitMotion::valid_at(double t_tdb) const
{
    return m_rolling || m_table.covers(t_tdb);
}

void KerrOrbitMotion::history_times(double t0, double t1, int max_points, std::vector<double>& out) const
{
    ensure_covered(t1);
    m_table.history_times(std::max(t0, m_table.start()), t1, max_points, out);
}

} // namespace astraxis
