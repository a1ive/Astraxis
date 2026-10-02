#include "ephem/visual_orbit.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/kepler.hpp"

#include <cmath>

namespace astraxis {

namespace {

// Direction (ICRF) of the tangent-plane point at true-anomaly-plus-omega angle u.
glm::dvec3 sky_direction(const VisualOrbit& o, double u)
{
    const double ra = o.ra_deg * kDegToRad;
    const double dec = o.dec_deg * kDegToRad;
    const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
    const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
    const glm::dvec3 away(std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec));

    const double node = o.node_deg * kDegToRad;
    const double i = o.i_deg * kDegToRad;
    const double n = std::cos(u) * std::cos(node) - std::sin(u) * std::sin(node) * std::cos(i);
    const double e = std::cos(u) * std::sin(node) + std::sin(u) * std::cos(node) * std::cos(i);
    const double d = std::sin(u) * std::sin(i);
    return north * n + east * e + away * d;
}

} // namespace

glm::dmat3 visual_orbit_frame(const VisualOrbit& orbit)
{
    const double w = orbit.arg_peri_deg * kDegToRad;
    const glm::dvec3 p = sky_direction(orbit, w);
    const glm::dvec3 q = sky_direction(orbit, w + kPi / 2.0);
    return glm::dmat3(p, q, glm::cross(p, q));
}

State visual_orbit_state(const VisualOrbit& orbit, double gm, double t_tdb)
{
    const double n = std::sqrt(gm / (orbit.a_km * orbit.a_km * orbit.a_km));
    KeplerElements k;
    k.a = orbit.a_km;
    k.e = orbit.e;
    k.mean_anomaly = n * (t_tdb - orbit.t_peri_tdb);
    // In the perifocal frame (node = i = arg_peri = 0), then rotate to ICRF.
    const State s = kepler_state(k, n);
    const glm::dmat3 frame = visual_orbit_frame(orbit);
    return {frame * s.position, frame * s.velocity};
}

VisualOrbit visual_orbit_from_transit(const TransitOrbit& transit, double gm)
{
    VisualOrbit o;
    const double n = kTwoPi / (transit.period_days * kSecondsPerDay);
    o.a_km = std::cbrt(gm / (n * n));
    o.e = transit.e;
    o.i_deg = transit.i_deg;
    o.node_deg = transit.node_deg;
    o.arg_peri_deg = transit.arg_peri_deg + 270.0 - transit.transit_u_deg;
    o.ra_deg = transit.ra_deg;
    o.dec_deg = transit.dec_deg;

    // Mid-transit at u = 270 deg: true, eccentric and mean anomaly there.
    const double f = 1.5 * kPi - o.arg_peri_deg * kDegToRad;
    const double ecc_anomaly =
        2.0 * std::atan2(std::sqrt(1.0 - o.e) * std::sin(0.5 * f), std::sqrt(1.0 + o.e) * std::cos(0.5 * f));
    const double mean_anomaly = ecc_anomaly - o.e * std::sin(ecc_anomaly);
    o.t_peri_tdb = transit.t_transit_tdb - mean_anomaly / n;
    return o;
}

} // namespace astraxis
