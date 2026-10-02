#include "ephem/jpl_approx_planet.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <cmath>

namespace astraxis {

namespace {

// Obliquity of the ecliptic at J2000 used by Standish's tables: 23.43928 deg.
constexpr double kObliquityJ2000 = 23.43928 * kDegToRad;

glm::dmat3 ecliptic_to_icrf()
{
    return rotation_x(kObliquityJ2000);
}

} // namespace

KeplerElements JplApproxPlanetOrbit::ecliptic_elements_at(double t_tdb, double* mean_motion) const
{
    const double T = julian_centuries_from_tdb(t_tdb);
    const ApproxPlanetElements& p = m_el;

    const double a = p.a_au + p.a_rate * T;
    const double e = p.e + p.e_rate * T;
    const double i = p.i_deg + p.i_rate * T;
    const double L = p.mean_longitude_deg + p.mean_longitude_rate * T;
    const double lon_peri = p.lon_peri_deg + p.lon_peri_rate * T;
    const double node = p.node_deg + p.node_rate * T;
    const double fT = p.f * T * kDegToRad;
    const double M = L - lon_peri + p.b * T * T + p.c * std::cos(fT) + p.s * std::sin(fT);

    KeplerElements k;
    k.a = a * kAuKm;
    k.e = e;
    k.i = i * kDegToRad;
    k.node = node * kDegToRad;
    k.arg_peri = (lon_peri - node) * kDegToRad;
    k.mean_anomaly = M * kDegToRad;

    if (mean_motion) {
        constexpr double kSecondsPerCentury = kDaysPerJulianCentury * kSecondsPerDay;
        *mean_motion = p.mean_longitude_rate * kDegToRad / kSecondsPerCentury;
    }
    return k;
}

State JplApproxPlanetOrbit::eval(double t_tdb) const
{
    double n = 0.0;
    const KeplerElements k = ecliptic_elements_at(t_tdb, &n);
    const State s = kepler_state(k, n);
    const glm::dmat3 r = ecliptic_to_icrf();
    return {r * s.position, r * s.velocity};
}

bool JplApproxPlanetOrbit::sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const
{
    const KeplerElements k = ecliptic_elements_at(t_tdb, nullptr);
    sample_ellipse(k, solve_kepler(k.mean_anomaly, k.e), count, out);
    const glm::dmat3 r = ecliptic_to_icrf();
    for (glm::dvec3& p : out) {
        p = r * p;
    }
    return true;
}

} // namespace astraxis
