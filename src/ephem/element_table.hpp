#pragma once

#include "ephem/kepler.hpp"
#include "ephem/motion.hpp"

#include <glm/mat3x3.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// A table of slowly varying orbital elements (tools/bake/element_table_bake.py),
// for moons close to an oblate planet whose Kepler-relative ephemerides would
// need a knot every few hours (Janus and Epimetheus). Each knot holds the
// elements of a precessing ellipse fitted to the positions over a window: a, e,
// the longitude of periapsis varpi, i, the node and the mean longitude lambda in
// a reference plane given by its ICRF pole (x axis: the plane's ascending node
// on the ICRF equator), and the mean motion n = dlambda/dt, which the planet's
// oblateness makes faster than Kepler's third law says for a.
//
// Between knots a, e, varpi, i and the node are interpolated linearly and lambda
// as a cubic Hermite curve through lambda and n; outside the table, lambda runs
// on at the end mean motion and varpi and the node at the table's precession
// rates.
//
// File format (little endian):
//   char[8] "AXELM1\0\0"; int32 target NAIF id; int32 center NAIF id;
//   uint32 knot count; uint32 0;
//   float64 reference pole RA, Dec (deg), apsidal and nodal precession rate (rad/s)
//   knots: float64 t (TDB s since J2000), a (km), e, varpi, i, node, lambda (rad;
//          varpi, node and lambda unwrapped), n (rad/s)
class ElementTable {
public:
    struct Knot {
        double t = 0.0;
        double a = 0.0;
        double e = 0.0;
        double varpi = 0.0;
        double i = 0.0;
        double node = 0.0;
        double lambda = 0.0;
        double n = 0.0;
    };

    bool load(const std::filesystem::path& path, std::string* error);

    int target_id() const { return m_target_id; }
    int center_id() const { return m_center_id; }
    double start() const { return m_knots.front().t; }
    double end() const { return m_knots.back().t; }
    const std::vector<Knot>& knots() const { return m_knots; }
    // Reference plane axes in ICRF (columns).
    const glm::dmat3& plane_to_icrf() const { return m_plane_to_icrf; }

    // Interpolated (or extrapolated) elements at t; `lambda_rate` receives dlambda/dt.
    Knot at(double t, double* lambda_rate = nullptr) const;

    double apsis_rate() const { return m_apsis_rate; }
    double node_rate() const { return m_node_rate; }

private:
    int m_target_id = 0;
    int m_center_id = 0;
    glm::dmat3 m_plane_to_icrf{1.0};
    double m_apsis_rate = 0.0; // rad/s
    double m_node_rate = 0.0;  // rad/s
    std::vector<Knot> m_knots;
};

class ElementTableOrbit final : public MotionSource {
public:
    explicit ElementTableOrbit(std::shared_ptr<const ElementTable> table)
        : m_table(std::move(table))
    {
    }

    State eval(double t_tdb) const override;
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override;

    const ElementTable& table() const { return *m_table; }
    // Elements at t_tdb in the reference plane frame (the mean longitude in mean_anomaly
    // minus the longitude of periapsis), and the mean anomaly's rate.
    KeplerElements elements_at(double t_tdb, double* mean_anomaly_rate = nullptr) const;

private:
    std::shared_ptr<const ElementTable> m_table;
};

} // namespace astraxis
