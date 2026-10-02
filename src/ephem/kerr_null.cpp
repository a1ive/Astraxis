#include "ephem/kerr_null.hpp"

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
    RayResult out;
    Photon ph = make_photon(x, dir, a);
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
            return out;
        }
        if (glm::length(ph.x) > escape_radius && glm::dot(ph.x, ph.p) > 0.0) {
            glm::dvec3 dx;
            glm::dvec3 dp;
            derivatives(ph.x, ph.p, ph.pt, a, dx, dp);
            out.outcome = RayOutcome::Escaped;
            out.direction = glm::normalize(dx);
            return out;
        }
    }
    return out;
}

double disk_redshift(const Photon& ph, double a)
{
    const glm::dvec3 x(ph.x.x, ph.x.y, 0.0);
    const KerrSchild k = kerr_schild(x, a);
    const double r = k.r;
    const double omega = 1.0 / (r * std::sqrt(r) + a); // prograde Keplerian (M = 1)
    // v = (1, -omega y, omega x, 0); g(v, v) = eta(v, v) + f (l_m v^m)^2 with l_t = 1.
    const double lv = 1.0 + omega * (-x.y * k.l.x + x.x * k.l.y);
    const double gvv = -1.0 + omega * omega * (x.x * x.x + x.y * x.y) + k.f * lv * lv;
    const double ut = 1.0 / std::sqrt(-gvv);
    const double pu_emitter = ut * (ph.pt - omega * x.y * ph.p.x + omega * x.x * ph.p.y);
    return ph.pt / pu_emitter;
}

} // namespace astraxis
