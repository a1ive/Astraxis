#include "ephem/kerr_geodesic.hpp"

#include <glm/geometric.hpp>
#include <glm/mat2x2.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

// Equatorial Kerr metric in 3+1 form and its spatial derivatives, written with
// the auxiliary functions of the reference project (RelativityMath.cs):
//   D = r^4 + a^2 r^2 + 2 M a^2 r,  Delta = r^2 - 2 M r + a^2,
//   G = 2 M a r / D (frame-dragging rate), H = 1 / Delta, K = D / r^6,
//   alpha = sqrt(Delta r^2 / D),  beta^i = (y G, -x G),
//   gamma_ij = H x_i x_j + K (r^2 delta_ij - x_i x_j)   (in 2D).
struct KerrEquatorial::Metric {
    double alpha = 1.0;
    glm::dvec2 beta{0.0};
    glm::dmat2 gamma_inv{1.0};
    glm::dvec2 grad_alpha{0.0};
    glm::dvec2 dbeta_dx{0.0}; // d(beta^k)/dx
    glm::dvec2 dbeta_dy{0.0}; // d(beta^k)/dy
    glm::dmat2 dgamma_inv_dx{0.0};
    glm::dmat2 dgamma_inv_dy{0.0};
};

KerrEquatorial::KerrEquatorial(double mass_km, double spin_a_km)
    : m_M(mass_km)
    , m_a(spin_a_km)
{
}

double KerrEquatorial::horizon_radius() const
{
    return m_M + std::sqrt(std::max(0.0, m_M * m_M - m_a * m_a));
}

KerrEquatorial::Metric KerrEquatorial::metric(const glm::dvec2& p) const
{
    const double M = m_M;
    const double a = m_a;
    const double a2 = a * a;
    const double x = p.x;
    const double y = p.y;
    const double r2 = x * x + y * y;
    const double r = std::sqrt(r2);
    const double r6 = r2 * r2 * r2;

    const double D = r2 * r2 + a2 * r2 + 2.0 * M * a2 * r;
    const double dD = 4.0 * r2 * r + 2.0 * a2 * r + 2.0 * M * a2;
    const double delta = r2 - 2.0 * M * r + a2;
    const double d_delta = 2.0 * (r - M);
    const double G = a * 2.0 * M * r / D;
    const double dG = a * (2.0 * M / D - 2.0 * M * r * dD / (D * D));
    const double H = 1.0 / delta;
    const double dH = -d_delta / (delta * delta);
    const double K = D / r6;
    const double dK = dD / r6 - 6.0 * D / (r6 * r);

    Metric m;
    m.alpha = std::sqrt(delta * r2 / D);
    m.beta = glm::dvec2(y * G, -x * G);

    // glm matrices are column-major; all of these are symmetric.
    const glm::dmat2 gamma(x * x * H + y * y * K, x * y * (H - K), x * y * (H - K), y * y * H + x * x * K);
    m.gamma_inv = glm::inverse(gamma);

    const double common = (d_delta / delta + 2.0 / r - dD / D) / (2.0 * r);
    m.grad_alpha = glm::dvec2(m.alpha * x * common, m.alpha * y * common);

    m.dbeta_dx = glm::dvec2(y * dG * x / r, -G - x * dG * x / r);
    m.dbeta_dy = glm::dvec2(G + y * dG * y / r, -x * dG * y / r);

    const double rx = x / r;
    const double ry = y / r;
    const glm::dmat2 dgamma_dx(2.0 * x * H + (x * x * dH + y * y * dK) * rx,
                               y * (H - K) + x * y * (dH - dK) * rx,
                               y * (H - K) + x * y * (dH - dK) * rx,
                               2.0 * x * K + (x * x * dK + y * y * dH) * rx);
    const glm::dmat2 dgamma_dy(2.0 * y * K + (x * x * dH + y * y * dK) * ry,
                               x * (H - K) + x * y * (dH - dK) * ry,
                               x * (H - K) + x * y * (dH - dK) * ry,
                               2.0 * y * H + (x * x * dK + y * y * dH) * ry);
    // d(gamma^-1) = -gamma^-1 d(gamma) gamma^-1
    m.dgamma_inv_dx = -(m.gamma_inv * dgamma_dx * m.gamma_inv);
    m.dgamma_inv_dy = -(m.gamma_inv * dgamma_dy * m.gamma_inv);
    return m;
}

void KerrEquatorial::derivatives(const glm::dvec2& x, const glm::dvec2& u, glm::dvec2& dxdt, glm::dvec2& dudt) const
{
    const Metric m = metric(x);
    const glm::dvec2 gu = m.gamma_inv * u;
    const double u0 = std::sqrt(glm::dot(u, gu) + 1.0) / m.alpha;

    dxdt = gu / u0 - m.beta;
    dudt = -m.alpha * u0 * m.grad_alpha + glm::dvec2(glm::dot(u, m.dbeta_dx), glm::dot(u, m.dbeta_dy)) -
           glm::dvec2(glm::dot(u, m.dgamma_inv_dx * u), glm::dot(u, m.dgamma_inv_dy * u)) / (2.0 * u0);
}

bool KerrEquatorial::step(glm::dvec2& x, glm::dvec2& u, double dt) const
{
    auto f = [&](const glm::dvec4& y) {
        glm::dvec2 dx;
        glm::dvec2 du;
        derivatives(glm::dvec2(y.x, y.y), glm::dvec2(y.z, y.w), dx, du);
        return glm::dvec4(dx, du);
    };

    const glm::dvec4 y0(x, u);
    const double r = glm::length(x);
    const double un = glm::length(u) + 1e-3;
    const glm::dvec4 scale(r, r, un, un);

    // Explicit Euler predictor, then Newton on F(y1) = y1 - y0 - dt f((y0 + y1) / 2).
    glm::dvec4 y1 = y0 + f(y0) * dt;
    bool converged = false;
    for (int iter = 0; iter < 12; ++iter) {
        const glm::dvec4 mid = (y0 + y1) * 0.5;
        const glm::dvec4 F = y1 - y0 - f(mid) * dt;

        double err = 0.0;
        for (int i = 0; i < 4; ++i) {
            err = std::max(err, std::abs(F[i]) / scale[i]);
        }
        if (err < 1e-15) {
            converged = true;
            break;
        }

        // J = I - dt/2 df/dy at the midpoint (central differences, relative steps).
        glm::dmat4 J(1.0);
        for (int i = 0; i < 4; ++i) {
            const double h = 1e-7 * scale[i];
            glm::dvec4 yp = mid;
            glm::dvec4 ym = mid;
            yp[i] += h;
            ym[i] -= h;
            J[i] -= (f(yp) - f(ym)) / (2.0 * h) * (dt * 0.5);
        }
        y1 -= glm::inverse(J) * F;
    }

    x = glm::dvec2(y1.x, y1.y);
    u = glm::dvec2(y1.z, y1.w);
    return converged;
}

glm::dvec2 KerrEquatorial::u_from_velocity(const glm::dvec2& x, const glm::dvec2& dxdt) const
{
    const Metric m = metric(x);
    const glm::dmat2 gamma = glm::inverse(m.gamma_inv);
    const glm::dvec2 w = dxdt + m.beta;
    const double u0 = 1.0 / std::sqrt(m.alpha * m.alpha - glm::dot(w, gamma * w));
    return u0 * (gamma * w);
}

double KerrEquatorial::energy(const glm::dvec2& x, const glm::dvec2& u) const
{
    const Metric m = metric(x);
    const double u0 = std::sqrt(glm::dot(u, m.gamma_inv * u) + 1.0) / m.alpha;
    return m.alpha * m.alpha * u0 - glm::dot(m.beta, u);
}

double KerrEquatorial::angular_momentum(const glm::dvec2& x, const glm::dvec2& u) const
{
    return x.x * u.y - x.y * u.x;
}

} // namespace astraxis
