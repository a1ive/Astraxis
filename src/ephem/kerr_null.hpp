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

// g^mn p_m q_n of two covectors at x (w = t component).
double covector_dot(const glm::dvec3& x, double a, const glm::dvec4& p, const glm::dvec4& q);

// Orthonormal frame of a static observer at x: at rest in these coordinates,
// u along the Killing vector d/dt (exists outside the ergosphere). The camera's
// axes (Euclidean unit vectors in the coordinates) are made orthonormal in the
// metric - components along u removed, then Gram-Schmidt in the order forward,
// up, right - so that screen directions are angles in the observer's rest frame.
// Covariant components, w = t; e_t = 0 for all three axes since they are
// orthogonal to u.
struct ObserverFrame {
    glm::dvec4 u{0.0};
    glm::dvec4 right{0.0};
    glm::dvec4 up{0.0};
    glm::dvec4 forward{0.0};
    // Spatial contravariant parts of the axes: coordinate directions of rays.
    glm::dvec3 right_dir{0.0};
    glm::dvec3 up_dir{0.0};
    glm::dvec3 forward_dir{0.0};
    double ut = 1.0; // u^t = 1 / sqrt(1 - f): observed / emitted energy of light from rest at infinity
};

ObserverFrame static_observer_frame(const glm::dvec3& x, double a, const glm::dvec3& forward, const glm::dvec3& up,
                                    const glm::dvec3& right);

// Time-reversed photon arriving at the observer from the local direction
// n = (right, up, forward) components (normalized here), with observed energy 1.
// Linear in the frame vectors, as evaluated per pixel in the shader.
Photon frame_photon(const glm::dvec3& x, const ObserverFrame& frame, const glm::dvec3& n);

// One RK4 step of affine length h.
void rk4_step(Photon& ph, double a, double h);

// Hamilton's equations above (dx^i/dl, dp_i/dl and dt/dl = -p_t + f (l.p)). They
// hold for any H: a massive particle follows the same equations with H = -1/2
// (l = proper time), so the Photon state serves for it too.
void geodesic_derivatives(const Photon& ph, double a, glm::dvec3& dx, glm::dvec3& dp, double& dt);

// A massive particle at x with coordinate velocity v = dx/dt (|v| < 1 well
// outside the hole): future directed, H = -1/2, p_t = -E.
Photon make_particle(const glm::dvec3& x, const glm::dvec3& v, double a);

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
    Photon photon; // state at the end
};

// Traces a ray from x along dir (both in Kerr-Schild Cartesian, M = 1) into the
// past until it falls into the hole or escapes beyond `escape_radius`.
RayResult trace_ray(const glm::dvec3& x, const glm::dvec3& dir, double a, double escape_radius, double disk_r_in,
                    double disk_r_out, int max_steps = 4000);
// The same from a given photon (e.g. frame_photon()).
RayResult trace_ray(Photon ph, double a, double escape_radius, double disk_r_in, double disk_r_out,
                    int max_steps = 4000);

// Weak field beyond the integration radius, exact for Schwarzschild in these
// (areal radius) coordinates to second order in M / r; spin is neglected there.
//
// Angle swept around the hole between infinity and radius r = b / w along a ray
// with impact parameter b = L / E, on a leg where r changes monotonically
// (w <= ~0.8: the expansion breaks down near the periapsis).
double weak_leg_angle(double w, double b);
// For rays passing far from the hole (b >~ 100): how much the direction still
// turns from its tangent at radius r to the outgoing asymptote, moving inward
// (before the periapsis) or outward. Total over a pass: 4/b + 15 pi / (4 b^2).
double weak_pass_deflection(double r, double b, bool inward);

// A camera ray as the shader traces it: the photon from local direction n
// (frame_photon) is integrated inside r_int (entering it on the exact ray when
// the camera is outside) and continued analytically beyond; rays with
// b >= 0.8 r_int are not integrated at all. direction: at infinity.
RayResult trace_camera_ray(const glm::dvec3& cam, const ObserverFrame& frame, const glm::dvec3& n, double a,
                           double r_int, double disk_r_in, double disk_r_out, int max_steps = 4000);

// Step length used by trace_ray (also mirrored in the shader).
double ray_step(const glm::dvec3& x, double a);

// Redshift factor g = nu_observed / nu_emitted for a photon (traced backward,
// at the disk point x on z = 0) emitted by matter on a prograde circular
// Keplerian orbit, received by a static observer with u^t = observer_ut
// (1: far away; ObserverFrame::ut for a camera nearby). Independent of the
// photon's normalization.
double disk_redshift(const Photon& ph, double a, double observer_ut = 1.0);

} // namespace astraxis
