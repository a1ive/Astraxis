#pragma once

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

namespace astraxis {

// Null geodesics (light rays) around a Kerr black hole in Cartesian
// Kerr-Schild coordinates, in units G = c = M = 1, spin a along +z.
// This is the CPU reference for shaders/black_hole.frag.hlsl (same formulas).
//
// Metric: g_mn = eta_mn + f l_m l_n, inverse g^mn = eta^mn - f l^m l^n, with
//   r^2 = (R^2 - a^2)/2 + sqrt((R^2 - a^2)^2 / 4 + a^2 z^2),   R^2 = x^2 + y^2 + z^2,
//   f   = 2 r^3 / (r^4 + a^2 z^2),
//   l_m = (1, (r x + a y)/(r^2 + a^2), (r y - a x)/(r^2 + a^2), z / r),  l^t = -1.
// Hamiltonian H = (eta^mn p_m p_n - f (l.p)^2) / 2 with l.p = -p_t + l_i p_i:
//   dx^t/dl = -p_t + f (l.p),   dx^i/dl = p_i - f (l.p) l_i,
//   dp_t/dl = 0,                dp_i/dl = [d_i f (l.p)^2 + 2 f (l.p) (d_i l_j) p_j] / 2.
// The spatial derivatives use the implicit derivative of r(x, y, z).
struct KerrSchild {
    double r = 0.0;
    double f = 0.0;
    glm::dvec3 l{0.0};      // spatial l_i
    glm::dvec3 df{0.0};     // d f / d x_j
    glm::dvec3 dl[3];       // dl[i][j] = d l_i / d x_j
};

KerrSchild kerr_schild(const glm::dvec3& x, double a);

// Boyer-Lindquist radius r of a Cartesian Kerr-Schild point.
double kerr_schild_r(const glm::dvec3& x, double a);

// Ray state: position and covariant momentum (p.w = p_t, conserved).
struct Photon {
    glm::dvec3 x{0.0};
    glm::dvec3 p{0.0};
    double pt = 0.0;
};

// Photon with spatial momentum `dir` at x, with p_t chosen so that H = 0. With
// `backward` the ray is traced into the past (from a camera into the scene):
// the momentum is that of the arriving photon with its sign flipped.
Photon make_photon(const glm::dvec3& x, const glm::dvec3& dir, double a);

double hamiltonian(const Photon& ph, double a);

// One RK4 step of affine length h.
void rk4_step(Photon& ph, double a, double h);

// Outer horizon r+ and prograde ISCO radius (Bardeen, Press & Teukolsky 1972).
double horizon_radius(double a);
double isco_radius(double a);

enum class RayOutcome { Captured, Escaped, MaxSteps };

struct RayResult {
    RayOutcome outcome = RayOutcome::MaxSteps;
    glm::dvec3 direction{0.0}; // escape direction (unit) when Escaped
    int disk_crossings = 0;    // z = 0 crossings inside [r_in, r_out]
    double first_crossing_r = 0.0;
    int steps = 0;
};

// Traces a ray from x along dir (both in Kerr-Schild Cartesian, M = 1) into the
// past until it falls into the hole or escapes beyond `escape_radius`.
RayResult trace_ray(const glm::dvec3& x, const glm::dvec3& dir, double a, double escape_radius, double disk_r_in,
                    double disk_r_out, int max_steps = 4000);

// Step length used by trace_ray (also mirrored in the shader).
double ray_step(const glm::dvec3& x, double a);

// Redshift factor g = nu_observed / nu_emitted for a photon (traced backward,
// at the disk point x on z = 0) emitted by matter on a prograde circular
// Keplerian orbit, received by a static observer far away.
double disk_redshift(const Photon& ph, double a);

} // namespace astraxis
