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
State propagate_kepler(const State& s, double mu, double dt);

// Ellipse points in the reference frame, starting at eccentric anomaly
// `start_E` and going backwards one full revolution (`count` points, closed).
void sample_ellipse(const KeplerElements& el, double start_E, int count, std::vector<glm::dvec3>& out);

// Osculating ellipse through a state vector around a center with gravitational
// parameter `mu` (km^3/s^2), sampled like sample_ellipse (starting at the
// current position). Returns false for unbound or degenerate orbits.
bool sample_osculating_ellipse(const State& s, double mu, int count, std::vector<glm::dvec3>& out);

} // namespace astraxis
