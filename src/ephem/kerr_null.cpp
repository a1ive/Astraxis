#include "ephem/kerr_null.hpp"

#include "core/math.hpp"

#include <glm/geometric.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

double kerr_schild_r(const glm::dvec3& x, double a)
{
    const double a2 = a * a;
    const double w = glm::dot(x, x) - a2;
    const double r2 = 0.5 * (w + std::sqrt(w * w + 4.0 * a2 * x.z * x.z));
    return std::sqrt(std::max(r2, 0.0));
}

KerrSchild kerr_schild(const glm::dvec3& p, double a)
{
    const double x = p.x;
    const double y = p.y;
    const double z = p.z;
    const double a2 = a * a;
    const double w = x * x + y * y + z * z - a2;
    const double root = std::sqrt(w * w + 4.0 * a2 * z * z); // = 2 r^2 - w
    const double r2 = 0.5 * (w + root);
    const double r = std::sqrt(r2);

    KerrSchild k;
    k.r = r;

    // dr/dx_j from r^4 - w r^2 - a^2 z^2 = 0.
    const glm::dvec3 dr(x * r / root, y * r / root, z * (r2 + a2) / (r * root));

    const double num = 2.0 * r2 * r;
    const double den = r2 * r2 + a2 * z * z;
    k.f = num / den;
    for (int j = 0; j < 3; ++j) {
        const double dden = 4.0 * r2 * r * dr[j] + (j == 2 ? 2.0 * a2 * z : 0.0);
        k.df[j] = (6.0 * r2 * dr[j] * den - num * dden) / (den * den);
    }

    const double s = r2 + a2;
    const double nx = r * x + a * y;
    const double ny = r * y - a * x;
    k.l = glm::dvec3(nx / s, ny / s, z / r);
    for (int j = 0; j < 3; ++j) {
        const double dnx = dr[j] * x + (j == 0 ? r : 0.0) + (j == 1 ? a : 0.0);
        const double dny = dr[j] * y + (j == 1 ? r : 0.0) - (j == 0 ? a : 0.0);
        const double ds = 2.0 * r * dr[j];
        k.dl[0][j] = (dnx * s - nx * ds) / (s * s);
        k.dl[1][j] = (dny * s - ny * ds) / (s * s);
        k.dl[2][j] = ((j == 2 ? r : 0.0) - z * dr[j]) / r2;
    }
    return k;
}

Photon make_photon(const glm::dvec3& x, const glm::dvec3& dir, double a)
{
    // Solve H = 0 for p_t: (1 + f) q^2 - 2 f L q - (|p|^2 - f L^2) = 0, L = l.p.
    // The positive root is the past-directed (time-reversed) photon.
    const KerrSchild k = kerr_schild(x, a);
    Photon ph;
    ph.x = x;
    ph.p = glm::normalize(dir);
    const double L = glm::dot(k.l, ph.p);
    const double pp = glm::dot(ph.p, ph.p);
    const double disc = k.f * k.f * L * L + (1.0 + k.f) * (pp - k.f * L * L);
    ph.pt = (k.f * L + std::sqrt(std::max(disc, 0.0))) / (1.0 + k.f);
    return ph;
}

double hamiltonian(const Photon& ph, double a)
{
    const KerrSchild k = kerr_schild(ph.x, a);
    const double L = -ph.pt + glm::dot(k.l, ph.p);
    return 0.5 * (-ph.pt * ph.pt + glm::dot(ph.p, ph.p) - k.f * L * L);
}

double covector_dot(const glm::dvec3& x, double a, const glm::dvec4& p, const glm::dvec4& q)
{
    // g^mn = eta^mn - f l^m l^n with l^t = -1, l^i = l_i.
    const KerrSchild k = kerr_schild(x, a);
    const double lp = -p.w + glm::dot(k.l, glm::dvec3(p));
    const double lq = -q.w + glm::dot(k.l, glm::dvec3(q));
    return -p.w * q.w + glm::dot(glm::dvec3(p), glm::dvec3(q)) - k.f * lp * lq;
}

namespace {

// V_m = eta_mn V^n + f l_m (l_n V^n) for a vector V^m (w = t), with l_t = 1.
glm::dvec4 lower(const KerrSchild& k, const glm::dvec4& v)
{
    const double lv = v.w + glm::dot(k.l, glm::dvec3(v));
    return glm::dvec4(glm::dvec3(v) + k.f * lv * k.l, -v.w + k.f * lv);
}

// g_mn V^m W^n.
double metric_dot(const KerrSchild& k, const glm::dvec4& v, const glm::dvec4& w)
{
    return glm::dot(lower(k, v), w);
}

} // namespace

ObserverFrame static_observer_frame(const glm::dvec3& x, double a, const glm::dvec3& forward, const glm::dvec3& up,
                                    const glm::dvec3& right)
{
    const KerrSchild k = kerr_schild(x, a);
    // g_tt = -1 + f; the camera keeps out of the ergosphere (f < 1), the clamp only avoids NaN.
    const double ut = 1.0 / std::sqrt(std::max(1.0 - k.f, 1e-9));
    const glm::dvec4 u(0.0, 0.0, 0.0, ut);

    auto axis = [&](const glm::dvec3& v, const glm::dvec4* basis, int n) {
        glm::dvec4 e(v, 0.0);
        e += metric_dot(k, e, u) * u; // g(u, u) = -1
        for (int i = 0; i < n; ++i) {
            e -= metric_dot(k, e, basis[i]) * basis[i];
        }
        return e / std::sqrt(metric_dot(k, e, e));
    };
    glm::dvec4 e[3];
    e[0] = axis(forward, e, 0);
    e[1] = axis(up, e, 1);
    e[2] = axis(right, e, 2);

    ObserverFrame frame;
    frame.u = lower(k, u);
    frame.forward = lower(k, e[0]);
    frame.up = lower(k, e[1]);
    frame.right = lower(k, e[2]);
    frame.forward_dir = glm::dvec3(e[0]);
    frame.up_dir = glm::dvec3(e[1]);
    frame.right_dir = glm::dvec3(e[2]);
    frame.ut = ut;
    return frame;
}

Photon frame_photon(const glm::dvec3& x, const ObserverFrame& frame, const glm::dvec3& n)
{
    // Arriving photon k = u - n (observed energy -k.u = 1); traced backward: p = -k.
    const glm::dvec3 m = glm::normalize(n);
    const glm::dvec4 p = m.x * frame.right + m.y * frame.up + m.z * frame.forward - frame.u;
    Photon ph;
    ph.x = x;
    ph.p = glm::dvec3(p);
    ph.pt = p.w;
    return ph;
}

namespace {

void derivatives(const glm::dvec3& x, const glm::dvec3& p, double pt, double a, glm::dvec3& dx, glm::dvec3& dp)
{
    const KerrSchild k = kerr_schild(x, a);
    const double L = -pt + glm::dot(k.l, p);
    dx = p - k.f * L * k.l;
    for (int j = 0; j < 3; ++j) {
        const double dlp = k.dl[0][j] * p.x + k.dl[1][j] * p.y + k.dl[2][j] * p.z;
        dp[j] = 0.5 * (k.df[j] * L * L + 2.0 * k.f * L * dlp);
    }
}

} // namespace

void rk4_step(Photon& ph, double a, double h)
{
    glm::dvec3 k1x, k1p, k2x, k2p, k3x, k3p, k4x, k4p;
    derivatives(ph.x, ph.p, ph.pt, a, k1x, k1p);
    derivatives(ph.x + 0.5 * h * k1x, ph.p + 0.5 * h * k1p, ph.pt, a, k2x, k2p);
    derivatives(ph.x + 0.5 * h * k2x, ph.p + 0.5 * h * k2p, ph.pt, a, k3x, k3p);
    derivatives(ph.x + h * k3x, ph.p + h * k3p, ph.pt, a, k4x, k4p);
    ph.x += h / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    ph.p += h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
}

double horizon_radius(double a)
{
    return 1.0 + std::sqrt(std::max(0.0, 1.0 - a * a));
}

double isco_radius(double a)
{
    const double z1 = 1.0 + std::cbrt(1.0 - a * a) * (std::cbrt(1.0 + a) + std::cbrt(1.0 - a));
    const double z2 = std::sqrt(3.0 * a * a + z1 * z1);
    return 3.0 + z2 - std::sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2));
}

double ray_step(const glm::dvec3& x, double a)
{
    const double r = kerr_schild_r(x, a);
    return 0.07 * std::max(r - horizon_radius(a), 0.0) + 0.015;
}

RayResult trace_ray(const glm::dvec3& x, const glm::dvec3& dir, double a, double escape_radius, double disk_r_in,
                    double disk_r_out, int max_steps)
{
    return trace_ray(make_photon(x, dir, a), a, escape_radius, disk_r_in, disk_r_out, max_steps);
}

RayResult trace_ray(Photon ph, double a, double escape_radius, double disk_r_in, double disk_r_out, int max_steps)
{
    RayResult out;
    const double r_plus = horizon_radius(a);
    for (int i = 0; i < max_steps; ++i) {
        out.steps = i + 1;
        const glm::dvec3 prev = ph.x;
        rk4_step(ph, a, ray_step(ph.x, a));

        if (prev.z * ph.x.z < 0.0) {
            const double t = prev.z / (prev.z - ph.x.z);
            const glm::dvec3 c = prev + (ph.x - prev) * t;
            const double rc = kerr_schild_r(glm::dvec3(c.x, c.y, 0.0), a);
            if (rc >= disk_r_in && rc <= disk_r_out) {
                if (out.disk_crossings == 0) {
                    out.first_crossing_r = rc;
                }
                ++out.disk_crossings;
            }
        }
        const double r = kerr_schild_r(ph.x, a);
        if (r < r_plus + 0.02) {
            out.outcome = RayOutcome::Captured;
            out.photon = ph;
            return out;
        }
        if (glm::length(ph.x) > escape_radius && glm::dot(ph.x, ph.p) > 0.0) {
            glm::dvec3 dx;
            glm::dvec3 dp;
            derivatives(ph.x, ph.p, ph.pt, a, dx, dp);
            out.outcome = RayOutcome::Escaped;
            out.direction = glm::normalize(dx);
            out.photon = ph;
            return out;
        }
    }
    out.photon = ph;
    return out;
}

double weak_leg_angle(double w, double b)
{
    // phi = int_0^w dw' / sqrt(1 - w'^2 + e w'^3), e = 2 / b, expanded to e^2 with
    // w' = sin(th): int cos / sqrt(cos^2 + e sin^3) = th - (e/2) K1 + (3/8) e^2 K2,
    // K1 = int sin^3 / cos^2 = (1 - c)^2 / c, K2 = int sin^6 / cos^4 = tan^3 / 3 - 2 tan + 5/2 th - sin cos / 2.
    const double th = std::asin(w);
    const double c = std::cos(th);
    const double one_minus_c = w * w / (1.0 + c);
    const double k1 = one_minus_c * one_minus_c / c;
    const double t = std::tan(th);
    const double k2 = t * t * t / 3.0 - 2.0 * t + 2.5 * th - 0.5 * w * c;
    return th - k1 / b + 1.5 * k2 / (b * b);
}

double weak_pass_deflection(double r, double b, bool inward)
{
    // First-order orbit u(phi) = sin(phi) / b + (1 + cos^2 phi) / b^2, phi from the
    // incoming asymptote (u = 0 at phi_in ~ -2/b) to the outgoing one (pi - phi_in).
    // sin(phi) at u = 1/r: the small root of s^2 - b s + b^2 / r - 2 = 0.
    const double s_in = -4.0 / (b + std::sqrt(b * b + 8.0));
    const double phi_in = std::asin(s_in);
    const double phi_out = kPi - phi_in;
    double s = 2.0 * (b * b / r - 2.0) / (b + std::sqrt(std::max(b * b * (1.0 - 4.0 / r) + 8.0, 0.0)));
    s = std::clamp(s, -1.0, 1.0);
    const double phi = inward ? std::asin(s) : kPi - std::asin(s);
    // Tangent direction at phi: dx/dphi ~ u e_phi - u' e_r, at atan2(u, -u') from the radius.
    const double u = std::sin(phi) / b + (1.0 + std::cos(phi) * std::cos(phi)) / (b * b);
    const double du = std::cos(phi) * (1.0 / b - 2.0 * std::sin(phi) / (b * b));
    const double first = phi_out - phi - std::atan2(u, -du);
    // Second order: the whole pass bends by 15 pi / (4 b^2) more, shared out by the angle swept.
    return first + (phi_out - phi) / (phi_out - phi_in) * 15.0 * kPi / (4.0 * b * b);
}

namespace {

// v rotated by angle about the unit axis (v perpendicular to the axis).
glm::dvec3 rotate_in_plane(const glm::dvec3& v, const glm::dvec3& axis, double angle)
{
    return v * std::cos(angle) + glm::cross(axis, v) * std::sin(angle);
}

// Asymptotic direction of a ray leaving radius |x| outward (weak field beyond).
glm::dvec3 outgoing_direction(const Photon& ph)
{
    const glm::dvec3 l = glm::cross(ph.x, ph.p);
    const double r = glm::length(ph.x);
    const double b = glm::length(l) / ph.pt;
    const glm::dvec3 axis = l / std::max(glm::length(l), 1e-300);
    return rotate_in_plane(ph.x / r, axis, weak_leg_angle(std::min(b / r, 0.95), b));
}

} // namespace

RayResult trace_camera_ray(const glm::dvec3& cam, const ObserverFrame& frame, const glm::dvec3& n, double a,
                           double r_int, double disk_r_in, double disk_r_out, int max_steps)
{
    const Photon pc = frame_photon(cam, frame, n);
    const glm::dvec3 d = glm::normalize(n.x * frame.right_dir + n.y * frame.up_dir + n.z * frame.forward_dir);
    const double dist = glm::length(cam);
    const glm::dvec3 l = glm::cross(cam, pc.p);
    const double b = glm::length(l) / pc.pt; // L / E
    const glm::dvec3 axis = l / std::max(glm::length(l), 1e-300);
    const bool inward = glm::dot(cam, d) < 0.0;
    const double b_int = 0.8 * r_int; // beyond: legs with b / r > 0.8, where the expansion fails

    RayResult out;
    if (b >= b_int) {
        out.outcome = RayOutcome::Escaped;
        out.direction = glm::normalize(rotate_in_plane(d, axis, weak_pass_deflection(dist, b, inward)));
        return out;
    }
    Photon ph = pc;
    if (dist > r_int) {
        if (!inward) {
            out.outcome = RayOutcome::Escaped;
            out.direction = glm::normalize(rotate_in_plane(cam / dist, axis, weak_leg_angle(b / dist, b)));
            return out;
        }
        // Enter the integration sphere where the exact (Schwarzschild) ray crosses it:
        // swept angle from the leg formula, tangential p = L / r, radial p from H = 0.
        const double shift = weak_leg_angle(b / r_int, b) - weak_leg_angle(b / dist, b);
        const glm::dvec3 xr = rotate_in_plane(cam / dist, axis, shift);
        const double f = 2.0 / r_int;
        const double e = pc.pt;
        const double p_t = b * e / r_int;
        const double p_r = (-f * e - std::sqrt(std::max(e * e - (1.0 - f) * p_t * p_t, 0.0))) / (1.0 - f);
        ph = make_photon(r_int * xr, p_r * xr + p_t * glm::cross(axis, xr), a);
    }
    out = trace_ray(ph, a, r_int, disk_r_in, disk_r_out, max_steps);
    if (out.outcome == RayOutcome::Escaped) {
        out.direction = glm::normalize(outgoing_direction(out.photon));
    }
    return out;
}

double disk_redshift(const Photon& ph, double a, double observer_ut)
{
    // nu_obs / nu_emit = (p.u_observer) / (p.u_emitter); a static observer has p.u = p_t u^t.
    const glm::dvec3 x(ph.x.x, ph.x.y, 0.0);
    const KerrSchild k = kerr_schild(x, a);
    const double r = k.r;
    const double omega = 1.0 / (r * std::sqrt(r) + a); // prograde Keplerian (M = 1)
    // v = (1, -omega y, omega x, 0); g(v, v) = eta(v, v) + f (l_m v^m)^2 with l_t = 1.
    const double lv = 1.0 + omega * (-x.y * k.l.x + x.x * k.l.y);
    const double gvv = -1.0 + omega * omega * (x.x * x.x + x.y * x.y) + k.f * lv * lv;
    const double ut = 1.0 / std::sqrt(-gvv);
    const double pu_emitter = ut * (ph.pt - omega * x.y * ph.p.x + omega * x.x * ph.p.y);
    return observer_ut * ph.pt / pu_emitter;
}

} // namespace astraxis
