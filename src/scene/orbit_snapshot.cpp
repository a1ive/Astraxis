#include "scene/orbit_snapshot.hpp"

#include "core/math.hpp"

#include <cmath>

namespace astraxis {

glm::dmat3 ecliptic_j2000_to_icrf()
{
    return rotation_x(84381.448 / 3600.0 * kDegToRad);
}

double relativistic_apsidal_rate(double gm_km3_s2, double a_km, double e)
{
    const double c2 = kSpeedOfLightKmS * kSpeedOfLightKmS;
    return 3.0 * std::pow(gm_km3_s2, 1.5) / (c2 * std::pow(a_km, 2.5) * (1.0 - e * e));
}

KeplerElements newtonian_elements(const OrbitSnapshot& snapshot, double gm_km3_s2, double t_from_tdb)
{
    KeplerElements el = snapshot.elements;
    el.arg_peri -= relativistic_apsidal_rate(gm_km3_s2, el.a, el.e) * (snapshot.t_tdb - t_from_tdb);
    return el;
}

glm::dvec3 orbit_point_icrf(const KeplerElements& el, double nu)
{
    const double r = el.a * (1.0 - el.e * el.e) / (1.0 + el.e * std::cos(nu));
    return ecliptic_j2000_to_icrf() * (perifocal_to_reference(el) * glm::dvec3(r * std::cos(nu), r * std::sin(nu), 0.0));
}

} // namespace astraxis
