#pragma once

#include <glm/vec2.hpp>

namespace astraxis {

// Timelike geodesics in the equatorial plane of a Kerr black hole, in the 3+1
// (ADM) form of Bacchini et al. (2018, ApJS 237, 6). Ported from the reference
// project (ref/.../RelativityMath.cs and GRPhysicsJob.cs).
//
// Geometric units with G = c = 1 and lengths in km: M = GM / c^2, a = spin * M,
// coordinate time t in km (t_seconds = t / c). The state is the Cartesian-like
// position x (from Boyer-Lindquist r, phi) and the covariant spatial 4-velocity
// u_i (dimensionless). Evolution (Bacchini et al. eqs. for dx/dt, du/dt):
//   u^0     = sqrt(gamma^ij u_i u_j + 1) / alpha
//   dx^i/dt = gamma^ij u_j / u^0 - beta^i
//   du_i/dt = -alpha u^0 d_i alpha + u_k d_i beta^k - u_j u_k d_i gamma^jk / (2 u^0)
// The spin axis is +z: a > 0 is prograde for counter-clockwise orbits.
class KerrEquatorial {
public:
    KerrEquatorial(double mass_km, double spin_a_km);

    double mass() const { return m_M; }
    double spin() const { return m_a; }
    // Outer event horizon radius r+ = M + sqrt(M^2 - a^2).
    double horizon_radius() const;

    void derivatives(const glm::dvec2& x, const glm::dvec2& u, glm::dvec2& dxdt, glm::dvec2& dudt) const;

    // One implicit-midpoint step of length dt (km of coordinate time; may be
    // negative). Newton iterations with a central-difference Jacobian, as in the
    // reference integrator. Returns false if Newton did not converge.
    bool step(glm::dvec2& x, glm::dvec2& u, double dt) const;

    // Covariant velocity from a coordinate velocity dx/dt (dimensionless, v/c).
    glm::dvec2 u_from_velocity(const glm::dvec2& x, const glm::dvec2& dxdt) const;

    // Conserved quantities: energy per unit mass E = -u_t and axial angular
    // momentum L = u_phi = x u_y - y u_x.
    double energy(const glm::dvec2& x, const glm::dvec2& u) const;
    double angular_momentum(const glm::dvec2& x, const glm::dvec2& u) const;

private:
    struct Metric;
    Metric metric(const glm::dvec2& x) const;

    double m_M;
    double m_a;
};

} // namespace astraxis
