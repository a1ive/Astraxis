#pragma once

#include "ephem/kepler.hpp"
#include "ephem/motion.hpp"

#include <glm/mat3x3.hpp>

namespace astraxis {

// Mean orbital elements in the format of the JPL SSD "Planetary Satellite Mean
// Elements" tables: elements at an epoch, referred to a reference (Laplace)
// plane given by its ICRF pole, with linear apsidal and nodal precession.
struct MeanElements {
    double epoch_tdb = 0.0; // TDB seconds since J2000
    double a_km = 0.0;
    double e = 0.0;
    double arg_peri_deg = 0.0;
    double mean_anomaly_deg = 0.0;
    double i_deg = 0.0;
    double node_deg = 0.0;
    double period_days = 0.0; // period of the mean anomaly
    // Signed precession periods in Julian years: > 0 advancing, < 0 regressing,
    // 0 = no precession. (JPL tabulates magnitudes only.)
    double apsis_period_years = 0.0;
    double node_period_years = 0.0;
    // Quadratic term of the mean anomaly, M += accel * (days since epoch)^2
    // (e.g. Phobos' tidal acceleration).
    double mean_anomaly_accel_deg_per_day2 = 0.0;
    double pole_ra_deg = 0.0; // reference plane pole
    double pole_dec_deg = 90.0;
};

class MeanElementOrbit final : public MotionSource {
public:
    explicit MeanElementOrbit(const MeanElements& el);

    State eval(double t_tdb) const override;
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override;

    // Osculating elements at t_tdb, in the reference plane frame.
    KeplerElements elements_at(double t_tdb) const;
    double mean_motion() const { return m_mean_motion; }

private:
    MeanElements m_el;
    glm::dmat3 m_plane_to_icrf{1.0};
    double m_mean_motion = 0.0; // rad/s
    double m_apsis_rate = 0.0;  // rad/s
    double m_node_rate = 0.0;   // rad/s
};

} // namespace astraxis
