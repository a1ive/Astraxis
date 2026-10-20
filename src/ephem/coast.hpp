#pragma once

#include "ephem/ephemeris.hpp"
#include "ephem/motion.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace astraxis {

// The gravity a coasting spacecraft feels around a center body: the center's
// point mass, its oblateness (zonal harmonics J2, J4 about its pole) and the
// tides of other bodies (e.g. the Earth's on an orbit around the Moon: their
// pull on the spacecraft minus their pull on the center).
struct CoastGravity {
    double mu = 0.0; // km^3/s^2
    // Unnormalized zonal harmonics: zonal[n] is J_n (zonal[0] and [1] unused),
    // about `pole` (unit, ICRF) at reference radius `radius_km`.
    std::vector<double> zonal;
    double radius_km = 0.0;
    glm::dvec3 pole{0.0, 0.0, 1.0};
    struct Tide {
        double gm = 0.0;
        std::function<glm::dvec3(double)> position; // relative to the center at t_tdb
    };
    std::vector<Tide> tides;

    // Only the point mass: two-body (propagate_kepler) is exact.
    bool point_mass() const;
    glm::dvec3 acceleration(const glm::dvec3& r, double t_tdb) const;
};

// Coasts a state from t_tdb by dt seconds (either direction): two-body for a
// point mass, else 4th-order Runge-Kutta in steps of 1/500 of the local
// dynamical time, min(r / v, sqrt(r^3 / mu)).
State coast(const CoastGravity& gravity, const State& s, double t_tdb, double dt);

// The coast of `state` from `epoch_tdb` over [t0, t1] (which contains the
// epoch), one knot per integration step (see coast).
std::vector<EphemerisTable::Knot> integrate_coast(const CoastGravity& gravity, const State& state, double epoch_tdb,
                                                  double t0, double t1);

// Where a body would have gone without thrust from `state` at `epoch_tdb`
// (e.g. a spacecraft's path had it not fired its engine). Valid on [t0, t1].
// A point mass coasts in closed form; anything else is integrated over the
// window when constructed (integrate_coast; the tides' position functions are
// not kept).
class CoastMotion final : public MotionSource {
public:
    CoastMotion(const State& state, const CoastGravity& gravity, double epoch_tdb, double t0, double t1);

    State eval(double t_tdb) const override;
    bool valid_at(double t_tdb) const override { return t_tdb >= m_t0 && t_tdb <= m_t1; }
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override;
    void history_times(double t0, double t1, int max_points, std::vector<double>& out) const override;

private:
    State m_state;
    double m_mu;
    double m_epoch;
    double m_t0;
    double m_t1;
    std::shared_ptr<EphemerisTable> m_table; // integrated; null for a point mass
};

// The velocity change of `motion` (relative to the center of `gravity`) that
// gravity does not explain, from t0 over `duration`: each interval between
// the motion's own sampling times (MotionSource::history_times, an
// ephemeris' knots) coasts the state at its start to its end and adds the
// difference from the motion's velocity there. Over a burn this is the
// thrust's delta-v (the integral of its acceleration, as in the rocket
// equation, not the net change of the velocity vector, which turns meanwhile).
// Near a giant planet the oblateness pulls about as hard as the engine, so
// `gravity` must include it. Knots are exact samples; between them the cubic
// swings across a burn's sudden start or end, or a small seam between
// trajectory files, and finer steps would add up both ways of each swing (13%
// too much for Venus Express). Returns the running total at `steps` + 1 even
// samples (the first 0), linear between the sampling times.
std::vector<double> unexplained_delta_v(const MotionSource& motion, const CoastGravity& gravity, double t0,
                                        double duration, int steps);

} // namespace astraxis
