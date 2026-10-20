#include "ephem/kerr_orbit.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr double kC = kSpeedOfLightKmS;
constexpr int kMaxStepHalvings = 4;

// Hole frame -> ICRF: z along the spin axis, x toward the ascending node of
// the hole's equator on the ICRF equator.
glm::dmat3 hole_frame_from_pole(const glm::dvec3& pole)
{
    const glm::dvec3 z = glm::normalize(pole);
    const glm::dvec3 node = glm::cross(glm::dvec3(0.0, 0.0, 1.0), z);
    const glm::dvec3 x = glm::length(node) > 1e-12 ? glm::normalize(node) : glm::dvec3(1.0, 0.0, 0.0);
    return glm::dmat3(x, glm::cross(z, x), z);
}

// One RK4 step of proper time h (units of M) for the particle and its coordinate time t.
void particle_step(Photon& s, double& t, double a, double h)
{
    glm::dvec3 kx[4], kp[4];
    double kt[4];
    Photon y = s;
    for (int i = 0; i < 4; ++i) {
        if (i > 0) {
            const double c = i == 3 ? h : 0.5 * h;
            y.x = s.x + c * kx[i - 1];
            y.p = s.p + c * kp[i - 1];
        }
        geodesic_derivatives(y, a, kx[i], kp[i], kt[i]);
    }
    s.x += h / 6.0 * (kx[0] + 2.0 * kx[1] + 2.0 * kx[2] + kx[3]);
    s.p += h / 6.0 * (kp[0] + 2.0 * kp[1] + 2.0 * kp[2] + kp[3]);
    t += h / 6.0 * (kt[0] + 2.0 * kt[1] + 2.0 * kt[2] + kt[3]);
}

double particle_step_length(const Photon& s, double a, double eta)
{
    const double r = kerr_schild_r(s.x, a);
    return eta * r * std::sqrt(r);
}


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
        // Retry a failed Newton solve with shorter steps; if that fails too, end
        // the trajectory here rather than store a diverged state.
        glm::dvec2 x = s.x;
        glm::dvec2 u = s.u;
        bool ok = false;
        for (int halving = 0; !ok && halving <= kMaxStepHalvings; ++halving) {
            if (halving > 0) {
                dt_s *= 0.5;
                x = s.x;
                u = s.u;
            }
            ok =m_kerr.step(x, u, sign * dt_s * kC);
        }
        if (!ok) {
            break;
        }
        s.x = x;
        s.u = u;
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

KerrCircularMotion::KerrCircularMotion(double gm_km3_s2, double spin, double r_m, const glm::dvec3& pole,
                                       double phase_deg)
{
    const double m = gm_km3_s2 / (kC * kC);
    m_radius_km = m * std::sqrt(r_m * r_m + spin * spin);
    m_omega = kC / (m * (r_m * std::sqrt(r_m) + spin));
    m_phase = phase_deg * kDegToRad;
    const glm::dmat3 frame = hole_frame_from_pole(pole);
    m_e1 = frame[0];
    m_e2 = frame[1];
}

State KerrCircularMotion::eval(double t_tdb) const
{
    const double phi = m_phase + m_omega * t_tdb;
    const glm::dvec3 radial = std::cos(phi) * m_e1 + std::sin(phi) * m_e2;
    const glm::dvec3 along = -std::sin(phi) * m_e1 + std::cos(phi) * m_e2;
    return {m_radius_km * radial, m_radius_km * m_omega * along};
}

bool KerrCircularMotion::sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const
{
    out.clear();
    const int n = std::max(count, 2);
    for (int i = 0; i < n; ++i) {
        out.push_back(eval(t_tdb - period_s() * static_cast<double>(i) / static_cast<double>(n - 1)).position);
    }
    return true;
}

double KerrCircularMotion::period_s() const
{
    return kTwoPi / m_omega;
}


KerrApoapsisMotion::KerrApoapsisMotion(const KerrApoapsisSetup& setup, double t_begin, double t_end)
    : m_setup(setup)
    , m_frame(hole_frame_from_pole(setup.pole))
    , m_mass_km(setup.gm_km3_s2 / (kC * kC))
{
    const double i = setup.inclination_deg * kDegToRad;
    const double node = setup.node_deg * kDegToRad;
    const double w = setup.arg_apoapsis_deg * kDegToRad;
    const glm::dvec3 n_hat(std::cos(node), std::sin(node), 0.0);
    const glm::dvec3 m_hat(-std::cos(i) * std::sin(node), std::cos(i) * std::cos(node), std::sin(i));
    m_p = std::cos(w) * n_hat + std::sin(w) * m_hat;
    m_q = -std::sin(w) * n_hat + std::cos(w) * m_hat;

    // Periapsis grows with the speed at the apoapsis: bisect between a plunge
    // and the circular speed.
    double lo = 0.0;
    double hi = 1.0 / std::sqrt(setup.apoapsis_m);
    for (int k = 0; k < 60; ++k) {
        const double mid = 0.5 * (lo + hi);
        (first_periapsis(mid) < setup.periapsis_m ? lo : hi) = mid;
    }
    m_speed = 0.5 * (lo + hi);
    m_start = particle(m_speed);

    const double t0_m = 0.0;
    const double to_m = kC / m_mass_km; // seconds -> units of M
    std::vector<EphemerisTable::Knot> forward;
    std::vector<EphemerisTable::Knot> backward;
    integrate(m_start, t0_m, (t_end - setup.apoapsis_tdb) * to_m, forward);
    integrate(m_start, t0_m, (t_begin - setup.apoapsis_tdb) * to_m, backward);
    std::vector<EphemerisTable::Knot> knots(backward.rbegin(), backward.rend());
    knots.insert(knots.end(), forward.begin() + 1, forward.end());
    m_table.set_knots(std::move(knots));
}

Photon KerrApoapsisMotion::particle(double speed) const
{
    return make_particle(m_setup.apoapsis_m * m_p, speed * m_q, m_setup.spin);
}

double KerrApoapsisMotion::first_periapsis(double speed) const
{
    const double a = m_setup.spin;
    const double capture = horizon_radius(a) * 1.05;
    Photon s = particle(speed);
    double t = 0.0;
    double r_prev = kerr_schild_r(s.x, a);
    bool inward = false;
    for (int n = 0; n < 2000000; ++n) {
        particle_step(s, t, a, particle_step_length(s, a, m_setup.step_eta));
        const double r = kerr_schild_r(s.x, a);
        if (r < capture) {
            return 0.0;
        }
        if (r < r_prev) {
            inward = true;
        } else if (inward) {
            return r_prev;
        }
        r_prev = r;
    }
    return r_prev;
}

void KerrApoapsisMotion::integrate(Photon s, double t_m, double t_end_m, std::vector<EphemerisTable::Knot>& out) const
{
    const double a = m_setup.spin;
    const double sign = t_end_m >= t_m ? 1.0 : -1.0;
    const double capture = horizon_radius(a) * 1.05;
    auto knot = [&](const Photon& p, double t) {
        glm::dvec3 dx;
        glm::dvec3 dp;
        double dt = 1.0;
        geodesic_derivatives(p, a, dx, dp, dt);
        return EphemerisTable::Knot{m_setup.apoapsis_tdb + t * m_mass_km / kC, m_frame * p.x * m_mass_km,
                                    m_frame * (dx / dt) * kC};
    };
    out.push_back(knot(s, t_m));
    while (sign * (t_end_m - t_m) > 0.0) {
        particle_step(s, t_m, a, sign * particle_step_length(s, a, m_setup.step_eta));
        if (kerr_schild_r(s.x, a) < capture) {
            break; // plunged: the trajectory ends here
        }
        out.push_back(knot(s, t_m));
    }
}

State KerrApoapsisMotion::eval(double t_tdb) const
{
    return m_table.eval(std::clamp(t_tdb, m_table.start(), m_table.end()));
}

bool KerrApoapsisMotion::valid_at(double t_tdb) const
{
    return m_table.covers(t_tdb);
}

void KerrApoapsisMotion::history_times(double t0, double t1, int max_points, std::vector<double>& out) const
{
    m_table.history_times(std::max(t0, m_table.start()), std::min(t1, m_table.end()), max_points, out);
}

} // namespace astraxis
