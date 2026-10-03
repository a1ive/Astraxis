#include "ephem/mean_element_orbit.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

namespace astraxis {

namespace {

double rate_from_period(double period_seconds)
{
    return period_seconds != 0.0 ? kTwoPi / period_seconds : 0.0;
}

} // namespace

MeanElementOrbit::MeanElementOrbit(const MeanElements& el)
    : m_el(el)
{
    const glm::dvec3 pole = unit_from_ra_dec(el.pole_ra_deg * kDegToRad, el.pole_dec_deg * kDegToRad);
    m_plane_to_icrf = frame_from_pole(pole);

    constexpr double kSecondsPerYear = kDaysPerJulianYear * kSecondsPerDay;
    m_mean_motion = rate_from_period(el.period_days * kSecondsPerDay);
    m_apsis_rate = rate_from_period(el.apsis_period_years * kSecondsPerYear);
    m_node_rate = rate_from_period(el.node_period_years * kSecondsPerYear);
}

KeplerElements MeanElementOrbit::elements_at(double t_tdb) const
{
    const double dt = t_tdb - m_el.epoch_tdb;
    KeplerElements k;
    k.a = m_el.a_km;
    k.e = m_el.e;
    k.i = m_el.i_deg * kDegToRad;
    k.node = wrap_two_pi(m_el.node_deg * kDegToRad + m_node_rate * dt);
    k.arg_peri = wrap_two_pi(m_el.arg_peri_deg * kDegToRad + m_apsis_rate * dt);
    const double days = dt / kSecondsPerDay;
    k.mean_anomaly = wrap_two_pi((m_el.mean_anomaly_deg + m_el.mean_anomaly_accel_deg_per_day2 * days * days) * kDegToRad +
                                 m_mean_motion * dt);
    return k;
}

State MeanElementOrbit::eval(double t_tdb) const
{
    const double accel = m_el.mean_anomaly_accel_deg_per_day2 * kDegToRad / (kSecondsPerDay * kSecondsPerDay);
    const State s = kepler_state(elements_at(t_tdb), m_mean_motion + 2.0 * accel * (t_tdb - m_el.epoch_tdb));
    return {m_plane_to_icrf * s.position, m_plane_to_icrf * s.velocity};
}

bool MeanElementOrbit::sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const
{
    const KeplerElements k = elements_at(t_tdb);
    sample_ellipse(k, solve_kepler(k.mean_anomaly, k.e), count, out);
    for (glm::dvec3& p : out) {
        p = m_plane_to_icrf * p;
    }
    return true;
}

} // namespace astraxis
