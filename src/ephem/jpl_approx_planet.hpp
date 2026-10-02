#pragma once

#include "ephem/kepler.hpp"
#include "ephem/motion.hpp"

namespace astraxis {

// Keplerian elements and rates in the format of E. M. Standish, "Approximate
// Positions of the Planets" (JPL SSD, https://ssd.jpl.nasa.gov/planets/approx_pos.html),
// referred to the mean ecliptic and equinox of J2000. Rates are per Julian century.
// b, c, s, f are the extra mean-anomaly terms of Table 2b (zero for Table 1).
struct ApproxPlanetElements {
    double a_au = 0.0, a_rate = 0.0;
    double e = 0.0, e_rate = 0.0;
    double i_deg = 0.0, i_rate = 0.0;
    double mean_longitude_deg = 0.0, mean_longitude_rate = 0.0;
    double lon_peri_deg = 0.0, lon_peri_rate = 0.0;
    double node_deg = 0.0, node_rate = 0.0;
    double b = 0.0, c = 0.0, s = 0.0, f = 0.0;
};

// Heliocentric orbit of a planet, returned in ICRF-aligned axes (km, km/s).
class JplApproxPlanetOrbit final : public MotionSource {
public:
    explicit JplApproxPlanetOrbit(const ApproxPlanetElements& el)
        : m_el(el)
    {
    }

    State eval(double t_tdb) const override;
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override;

private:
    KeplerElements ecliptic_elements_at(double t_tdb, double* mean_motion) const;

    ApproxPlanetElements m_el;
};

} // namespace astraxis
