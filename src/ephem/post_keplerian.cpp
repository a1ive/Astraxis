#include "ephem/post_keplerian.hpp"

#include "core/math.hpp"
#include "ephem/kepler.hpp"

#include <cmath>

namespace astraxis {

namespace {

double semi_major_axis(double gm, double period_s)
{
    const double n = kTwoPi / period_s;
    return std::cbrt(gm / (n * n));
}

double current_period(const PostKeplerianOrbit& orbit, double t_tdb)
{
    return orbit.period_s + orbit.pb_dot * (t_tdb - orbit.elements.t_peri_tdb);
}

// Perifocal frame (P, Q, W) of the elements with the periastron at `arg_peri_rad`.
glm::dmat3 frame_at(const VisualOrbit& elements, double arg_peri_rad)
{
    VisualOrbit o = elements;
    o.arg_peri_deg = arg_peri_rad * kRadToDeg;
    return visual_orbit_frame(o);
}

} // namespace

PostKeplerianOrbit make_post_keplerian(const VisualOrbit& elements, double period_s, double omega_dot_rad_s,
                                       double pb_dot, double gm)
{
    PostKeplerianOrbit orbit;
    orbit.elements = elements;
    orbit.elements.a_km = semi_major_axis(gm, period_s);
    orbit.period_s = period_s;
    orbit.omega_dot_rad_s = omega_dot_rad_s;
    orbit.pb_dot = pb_dot;
    return orbit;
}

PostKeplerianPhase post_keplerian_phase(const PostKeplerianOrbit& orbit, double gm, double t_tdb)
{
    const double e = orbit.elements.e;
    const double orbits = (t_tdb - orbit.elements.t_peri_tdb) / orbit.period_s;
    const double turns = orbits - 0.5 * orbit.pb_dot * orbits * orbits; // M / (2 pi)
    const double whole = std::floor(turns);

    PostKeplerianPhase p;
    p.period_s = current_period(orbit, t_tdb);
    p.a_km = semi_major_axis(gm, p.period_s);
    p.mean_anomaly = kTwoPi * (turns - whole);
    p.ecc_anomaly = solve_kepler(p.mean_anomaly, e); // in [-pi, pi]
    const double nu = 2.0 * std::atan2(std::sqrt(1.0 + e) * std::sin(0.5 * p.ecc_anomaly),
                                       std::sqrt(1.0 - e) * std::cos(0.5 * p.ecc_anomaly));
    p.true_anomaly = wrap_two_pi(nu);
    // A_e: the true anomaly counted continuously from T0.
    const double k = orbit.omega_dot_rad_s * orbit.period_s / kTwoPi;
    const double a_e = p.true_anomaly + kTwoPi * whole;
    p.arg_peri_rad = orbit.elements.arg_peri_deg * kDegToRad + k * a_e;
    return p;
}

State post_keplerian_state(const PostKeplerianOrbit& orbit, double gm, double t_tdb)
{
    const PostKeplerianPhase p = post_keplerian_phase(orbit, gm, t_tdb);
    const double e = orbit.elements.e;
    const double a = p.a_km;
    const double cos_e = std::cos(p.ecc_anomaly);
    const double sin_e = std::sin(p.ecc_anomaly);
    const double b_over_a = std::sqrt(1.0 - e * e);

    // In the perifocal frame of the moment: Keplerian motion (mean motion dM/dt),
    // plus the turning of the frame itself, d omega / dt = k d nu / dt.
    const double n = kTwoPi / orbit.period_s *
                     (1.0 - orbit.pb_dot * (t_tdb - orbit.elements.t_peri_tdb) / orbit.period_s);
    const glm::dvec3 pos(a * (cos_e - e), a * b_over_a * sin_e, 0.0);
    const double speed = n * a / (1.0 - e * cos_e);
    glm::dvec3 vel(-speed * sin_e, speed * b_over_a * cos_e, 0.0);
    const double k = orbit.omega_dot_rad_s * orbit.period_s / kTwoPi;
    const double nu_dot = n * b_over_a / ((1.0 - e * cos_e) * (1.0 - e * cos_e));
    vel += k * nu_dot * glm::dvec3(-pos.y, pos.x, 0.0);

    const glm::dmat3 frame = frame_at(orbit.elements, p.arg_peri_rad);
    return {frame * pos, frame * vel};
}

VisualOrbit newtonian_continuation(const PostKeplerianOrbit& orbit, double gm, double t_ref_tdb)
{
    const PostKeplerianPhase p = post_keplerian_phase(orbit, gm, t_ref_tdb);
    VisualOrbit o = orbit.elements;
    o.a_km = p.a_km;
    o.arg_peri_deg = p.arg_peri_rad * kRadToDeg;
    // visual_orbit_state takes the mean motion from a through Kepler's third
    // law, i.e. the period of the moment.
    o.t_peri_tdb = t_ref_tdb - p.mean_anomaly / (kTwoPi / p.period_s);
    return o;
}

bool PostKeplerianMotion::valid_at(double t_tdb) const
{
    return current_period(m_orbit, t_tdb) > 0.0;
}

bool PostKeplerianMotion::sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const
{
    const PostKeplerianPhase p = post_keplerian_phase(m_orbit, m_gm, t_tdb);
    KeplerElements k;
    k.a = p.a_km;
    k.e = m_orbit.elements.e;
    sample_ellipse(k, p.ecc_anomaly, count, out);
    const glm::dmat3 frame = frame_at(m_orbit.elements, p.arg_peri_rad);
    for (glm::dvec3& point : out) {
        point = frame * point;
    }
    return true;
}

} // namespace astraxis
