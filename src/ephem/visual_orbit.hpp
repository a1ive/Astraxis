#pragma once

#include "ephem/motion.hpp"

#include <glm/mat3x3.hpp>

namespace astraxis {

// Campbell elements of an astrometric (visual binary) orbit: the companion's
// orbit relative to the primary, oriented on the sky at the primary's position.
//
// Conventions (e.g. Heintz, "Double Stars", 1978): position angles and the node
// Omega are measured from north through east; the tangent-plane offsets are
//   north = r [cos(nu + omega) cos Omega - sin(nu + omega) sin Omega cos i]
//   east  = r [cos(nu + omega) sin Omega + sin(nu + omega) cos Omega cos i]
//   depth = r  sin(nu + omega) sin i
// with depth measured away from the observer (the ascending node is where the
// companion recedes). Astrometry alone cannot tell the depth sign; it comes
// from radial velocities, and only mirrors the orbit through the sky plane.
struct VisualOrbit {
    double a_km = 0.0;
    double e = 0.0;
    double i_deg = 0.0;
    double node_deg = 0.0;     // Omega
    double arg_peri_deg = 0.0; // omega (of the companion)
    double t_peri_tdb = 0.0;   // time of periastron (TDB seconds since J2000)
    double ra_deg = 0.0;       // sky position of the system (ICRF)
    double dec_deg = 0.0;
};

// Perifocal frame in ICRF: columns P (toward periastron), Q (90 deg ahead in
// the direction of motion) and W (orbit normal).
glm::dmat3 visual_orbit_frame(const VisualOrbit& orbit);

// Keplerian state of the companion relative to the primary (ICRF, km, km/s);
// `gm` is the total G(m1 + m2) in km^3/s^2.
State visual_orbit_state(const VisualOrbit& orbit, double gm, double t_tdb);

// Osculating elements of a transiting planet as transit-timing fits publish
// them: the period and a time of mid-transit take the place of a and the time
// of periastron. Fitting codes differ in where they put the observer, which
// shows up as the argument of latitude u = omega + f at mid-transit:
//   90 deg  - the usual exoplanet convention (Winn 2010), e.g. Phodymm (Mills et al. 2016)
//   270 deg - NbodyGradient (Agol et al. 2021), whose z axis points away from the observer
// Either way the published angles describe the same orbit once omega is
// shifted by 270 deg - transit_u_deg, which puts the planet in front of the
// star (depth -r sin i in VisualOrbit's convention) at mid-transit.
struct TransitOrbit {
    double period_days = 0.0;
    double t_transit_tdb = 0.0; // a time of mid-transit (TDB seconds since J2000)
    double e = 0.0;
    double arg_peri_deg = 0.0;  // omega, in the fit's own convention
    double i_deg = 90.0;
    double node_deg = 0.0;      // Omega on the sky; transits do not constrain it
    double transit_u_deg = 90.0;
    double ra_deg = 0.0;        // sky position of the system (ICRF)
    double dec_deg = 0.0;
};

// The equivalent visual orbit, for a total G(m1 + m2) of `gm` (sets a through
// Kepler's third law, so the orbit keeps the published period).
VisualOrbit visual_orbit_from_transit(const TransitOrbit& transit, double gm);

// A fixed two-body (Keplerian) orbit on the sky: the companion relative to the
// primary. With ScaledMotion it places both around their barycenter, which is
// how hierarchical multiples are nested (pair of pairs of stars).
class VisualOrbitMotion final : public MotionSource {
public:
    VisualOrbitMotion(const VisualOrbit& orbit, double gm)
        : m_orbit(orbit)
        , m_gm(gm)
        , m_frame(visual_orbit_frame(orbit))
    {
    }
    State eval(double t_tdb) const override { return visual_orbit_state(m_orbit, m_gm, t_tdb); }
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override;
    const VisualOrbit& orbit() const { return m_orbit; }

private:
    VisualOrbit m_orbit;
    double m_gm;
    glm::dmat3 m_frame;
};

} // namespace astraxis
