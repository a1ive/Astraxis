#pragma once

#include "ephem/motion.hpp"
#include "ephem/visual_orbit.hpp"

namespace astraxis {

// A binary pulsar's orbit as relativistic timing models describe it (the
// orbital part of the Damour-Deruelle "DD" model; Damour & Deruelle 1986,
// Ann. Inst. Henri Poincare 44, 263): a Keplerian ellipse whose periastron
// advances and whose period shrinks.
//   - mean anomaly M = 2 pi [(t - T0) / Pb - Pb_dot / 2 ((t - T0) / Pb)^2]
//   - longitude of periastron omega = omega0 + k A_e, with A_e the true
//     anomaly counted continuously from T0 and k = omega_dot Pb / (2 pi): the
//     ellipse turns mostly near periastron, at omega_dot on average
//   - semi-major axis from Kepler's third law for the current period
//     Pb + Pb_dot (t - T0)
// The small post-Keplerian corrections to the shape (delta_r, delta_theta, of
// order 1e-6) and the secular drifts of e and of a sin i are left out.
struct PostKeplerianOrbit {
    VisualOrbit elements;   // at T0 = elements.t_peri_tdb; a_km is set from period_s
    double period_s = 0.0;  // Pb at T0
    double omega_dot_rad_s = 0.0;
    double pb_dot = 0.0;    // dimensionless (s/s)
};

// The orbit at one moment: its Keplerian elements then.
struct PostKeplerianPhase {
    double arg_peri_rad = 0.0; // of the companion, like elements.arg_peri_deg
    double period_s = 0.0;
    double a_km = 0.0;
    double mean_anomaly = 0.0; // in [0, 2 pi)
    double ecc_anomaly = 0.0;
    double true_anomaly = 0.0; // in [0, 2 pi)
};

// `gm` is the total G(m1 + m2) (km^3/s^2); sets orbit.elements.a_km.
PostKeplerianOrbit make_post_keplerian(const VisualOrbit& elements, double period_s, double omega_dot_rad_s,
                                       double pb_dot, double gm);

PostKeplerianPhase post_keplerian_phase(const PostKeplerianOrbit& orbit, double gm, double t_tdb);

// The companion relative to the primary (ICRF, km, km/s).
State post_keplerian_state(const PostKeplerianOrbit& orbit, double gm, double t_tdb);

// What Newtonian gravity predicts from the orbit as it was at t_ref: the
// Keplerian ellipse through the same position, with that moment's periastron
// and period kept forever.
VisualOrbit newtonian_continuation(const PostKeplerianOrbit& orbit, double gm, double t_ref_tdb);

class PostKeplerianMotion final : public MotionSource {
public:
    PostKeplerianMotion(const PostKeplerianOrbit& orbit, double gm)
        : m_orbit(orbit)
        , m_gm(gm)
    {
    }
    State eval(double t_tdb) const override { return post_keplerian_state(m_orbit, m_gm, t_tdb); }
    // The period runs out (the stars merge) long after any time the app shows.
    bool valid_at(double t_tdb) const override;
    // The ellipse of the moment, with the current periastron.
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override;
    const PostKeplerianOrbit& orbit() const { return m_orbit; }

private:
    PostKeplerianOrbit m_orbit;
    double m_gm;
};

} // namespace astraxis
