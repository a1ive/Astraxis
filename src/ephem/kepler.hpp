#pragma once

#include "ephem/motion.hpp"

#include <glm/mat3x3.hpp>

#include <vector>

namespace astraxis {

// Classical elliptic elements. Distances in km, angles in radians.
struct KeplerElements {
    double a = 0.0;
    double e = 0.0;
    double i = 0.0;
    double node = 0.0;     // longitude of ascending node
    double arg_peri = 0.0; // argument of periapsis
    double mean_anomaly = 0.0;
};

// Solves M = E - e sin E for the eccentric anomaly E (elliptic, 0 <= e < 1).
double solve_kepler(double mean_anomaly, double e);

// Rotation from the perifocal frame (x to periapsis, z to orbit normal) to the
// reference frame of the elements: Rz(node) Rx(i) Rz(arg_peri).
glm::dmat3 perifocal_to_reference(const KeplerElements& el);

// State in the reference frame of the elements; `mean_motion` in rad/s.
State kepler_state(const KeplerElements& el, double mean_motion);

// Two-body propagation of a state by dt seconds (either direction) around a
// center with gravitational parameter `mu` (km^3/s^2). Universal-variable
// formulation (H. D. Curtis, Orbital Mechanics for Engineering Students,
// Algorithms 3.3 and 3.4): valid for elliptic, parabolic and hyperbolic orbits.
// `universal_anomaly`, if given and nonzero, starts the Newton iteration (e.g.
// the previous frame's value for the same orbit) and receives the solution.
// mu <= 0 (no net attraction): straight-line motion.
State propagate_kepler(const State& s, double mu, double dt, double* universal_anomaly = nullptr);

// Ellipse points in the reference frame, starting at eccentric anomaly
// `start_E` and going backwards one full revolution (`count` points, closed).
void sample_ellipse(const KeplerElements& el, double start_E, int count, std::vector<glm::dvec3>& out);

// Osculating ellipse through a state vector around a center with gravitational
// parameter `mu` (km^3/s^2), sampled like sample_ellipse (starting at the
// current position). Returns false for unbound or degenerate orbits.
bool sample_osculating_ellipse(const State& s, double mu, int count, std::vector<glm::dvec3>& out);

// Two-body coast from `state` at `epoch_tdb` around a center of `mu`: where a
// body would have gone without thrust or perturbations (e.g. a spacecraft's
// path had it not fired its engine). Valid on [t0, t1].
class CoastMotion final : public MotionSource {
public:
    CoastMotion(const State& state, double mu, double epoch_tdb, double t0, double t1)
        : m_state(state)
        , m_mu(mu)
        , m_epoch(epoch_tdb)
        , m_t0(t0)
        , m_t1(t1)
    {
    }
    State eval(double t_tdb) const override { return propagate_kepler(m_state, m_mu, t_tdb - m_epoch); }
    bool valid_at(double t_tdb) const override { return t_tdb >= m_t0 && t_tdb <= m_t1; }
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override
    {
        return sample_osculating_ellipse(eval(t_tdb), m_mu, count, out);
    }

private:
    State m_state;
    double m_mu;
    double m_epoch;
    double m_t0;
    double m_t1;
};

// The velocity change of `motion` (relative to a center of `mu`) that
// two-body gravity does not explain, summed from t0 over `steps` equal steps
// of `duration`: each step coasts the motion's state to the next sample and
// adds the difference from the motion's velocity there. Over a burn this is
// the thrust's delta-v (the integral of its acceleration, as in the rocket
// equation, not the net change of the velocity vector, which turns meanwhile).
// Returns the running total at each of the steps + 1 samples (the first 0).
std::vector<double> unexplained_delta_v(const MotionSource& motion, double mu, double t0, double duration, int steps);

} // namespace astraxis
