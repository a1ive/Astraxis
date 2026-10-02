#pragma once

#include "ephem/ephemeris.hpp"

#include <glm/vec3.hpp>

#include <vector>

namespace astraxis {

struct NBodyParticle {
    double gm = 0.0; // km^3/s^2
    glm::dvec3 position{0.0};
    glm::dvec3 velocity{0.0};
};

struct NBodyOptions {
    // 4: Yoshida's 4th-order scheme. 6: three 4th-order steps composed with
    // weights z1 = 1 / (2 - 2^(1/5)), z0 = -2^(1/5) z1 (Yoshida 1990), three
    // times the work per step but far less phase drift. Tight resonant chains
    // need that: at ~75 steps per orbit the 4th-order scheme shifts
    // TRAPPIST-1 b's transits by an hour within four years.
    int order = 4;
    // Knots are relative to the barycentre of these particles (empty: the
    // integration frame as given).
    std::vector<size_t> origin;
    // Particle i keeps a knot every stride[i] steps (empty: every step). The
    // epoch and both ends are always kept. Fast inner planets need many more
    // knots than slow outer ones, and long spans would not fit in memory.
    std::vector<int> stride;
};

// Newtonian N-body integration with the 4th-order symplectic scheme of
// Yoshida (1990, Phys. Lett. A 150, 262): three drift-kick-drift stages with
// weights w1 = 1 / (2 - 2^(1/3)), w0 = -2^(1/3) w1 (or its 6th-order
// composition, see NBodyOptions). Fixed step; integrates
// from `t_epoch` forward to t_end and backward to t_begin.
// Returns one knot list per particle in time order.
std::vector<std::vector<EphemerisTable::Knot>> integrate_nbody(const std::vector<NBodyParticle>& particles,
                                                               double t_epoch, double t_begin, double t_end,
                                                               double step_s, const NBodyOptions& options = {});

// Total energy (kinetic + potential) per the particles' gm (for tests).
double nbody_energy(const std::vector<NBodyParticle>& particles);

} // namespace astraxis
