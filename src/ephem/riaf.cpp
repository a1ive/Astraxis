#include "ephem/riaf.hpp"

#include "core/math.hpp"

#include <glm/geometric.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

// CGS constants (CODATA 2018).
constexpr double kElectronCharge = 4.80320471e-10; // esu
constexpr double kElectronMass = 9.1093837015e-28; // g
constexpr double kProtonMass = 1.67262192369e-24;  // g
constexpr double kLightCm = 2.99792458e10;         // cm/s
constexpr double kBoltzmann = 1.380649e-16;        // erg/K

double bessel_i0_small(double x) // |x| < 3.75
{
    const double y = (x / 3.75) * (x / 3.75);
    return 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
}

double bessel_i1_small(double x) // |x| < 3.75
{
    const double y = (x / 3.75) * (x / 3.75);
    return x * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
}

double bessel_k0(double x)
{
    if (x <= 2.0) {
        const double y = 0.25 * x * x;
        return -std::log(0.5 * x) * bessel_i0_small(x) +
               (-0.57721566 + y * (0.42278420 + y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2 + y * (0.10750e-3 + y * 0.74e-5))))));
    }
    const double y = 2.0 / x;
    return std::exp(-x) / std::sqrt(x) *
           (1.25331414 + y * (-0.7832358e-1 + y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2 + y * (-0.251540e-2 + y * 0.53208e-3))))));
}

double bessel_k1(double x)
{
    if (x <= 2.0) {
        const double y = 0.25 * x * x;
        return std::log(0.5 * x) * bessel_i1_small(x) +
               (1.0 / x) * (1.0 + y * (0.15443144 + y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
    }
    const double y = 2.0 / x;
    return std::exp(-x) / std::sqrt(x) *
           (1.25331414 + y * (0.23498619 + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3)))))));
}

} // namespace

double bessel_k2(double x)
{
    return bessel_k0(x) + 2.0 / x * bessel_k1(x);
}

double riaf_far_view_fade(double ring_pixels)
{
    const double t = std::clamp((ring_pixels - 1.5) / 2.5, 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

RiafPlasma riaf_plasma(const glm::dvec3& x, double a, const RiafParams& p)
{
    RiafPlasma out;
    const double r = kerr_schild_r(x, a);
    if (r > p.outer_m || r <= horizon_radius(a)) {
        return out;
    }
    const double s2 = std::max(r * r - x.z * x.z, 1e-300);
    const double half_r = 0.5 * r; // r / r_S
    out.n_e_cm3 = p.n_e0_cm3 * std::pow(half_r, -1.1) * std::exp(-x.z * x.z / (2.0 * s2));
    out.t_e_k = p.t_e0_k * std::pow(half_r, -0.84);
    // B^2 / 8 pi = n m_p c^2 r_S / (12 r) / beta, with r_S / r = 1 / half_r.
    out.b_gauss = std::sqrt(8.0 * kPi * out.n_e_cm3 * kProtonMass * kLightCm * kLightCm / (12.0 * half_r * p.beta));
    return out;
}

double thermal_synchrotron_emissivity(double nu_hz, double n_e_cm3, double t_e_k, double b_gauss, double sin_pitch)
{
    const double theta = kBoltzmann * t_e_k / (kElectronMass * kLightCm * kLightCm);
    if (n_e_cm3 <= 0.0 || theta < 0.02 || b_gauss <= 0.0) {
        // Non-relativistic electrons: no synchrotron light here (and below ~0.0116 the
        // shader's float K_2 and emission underflow to 0 / 0).
        return 0.0;
    }
    const double nu_c = kElectronCharge * b_gauss / (2.0 * kPi * kElectronMass * kLightCm);
    const double nu_s = 2.0 / 9.0 * nu_c * theta * theta * sin_pitch;
    const double x = nu_hz / nu_s;
    const double shape = std::sqrt(x) + std::pow(2.0, 11.0 / 12.0) * std::cbrt(std::sqrt(x));
    return n_e_cm3 * std::sqrt(2.0) * kPi * kElectronCharge * kElectronCharge * nu_s /
           (3.0 * bessel_k2(1.0 / theta) * kLightCm) * shape * shape * std::exp(-std::cbrt(x));
}

glm::dvec4 riaf_four_velocity(const glm::dvec3& x, double a)
{
    const double r = kerr_schild_r(x, a);
    const double a2 = a * a;
    const double cos_t = std::clamp(x.z / r, -1.0, 1.0);
    const double cos2 = cos_t * cos_t;
    const double sin2 = 1.0 - cos2;
    const double r_isco = isco_radius(a);

    // Keplerian shell at rs = max(r, r_isco): u = u^t (d_t + Omega d_phi) in Boyer-Lindquist.
    const double rs = std::max(r, r_isco);
    const double omega = 1.0 / (rs * std::sqrt(rs) + a);
    const double sigma_s = rs * rs + a2 * cos2;
    const double delta_s = rs * rs - 2.0 * rs + a2;
    const double big_a_s = (rs * rs + a2) * (rs * rs + a2) - a2 * delta_s * sin2;
    const double g_tt = -(1.0 - 2.0 * rs / sigma_s);
    const double g_tp_s = -2.0 * a * rs / sigma_s; // g_tphi / sin^2
    const double g_pp_s = big_a_s / sigma_s;       // g_phiphi / sin^2
    const double ut_s = 1.0 / std::sqrt(-(g_tt + sin2 * (2.0 * omega * g_tp_s + omega * omega * g_pp_s)));

    if (r >= r_isco) {
        // In Kerr-Schild, d_t and d_phi are the same Killing vectors: d_phi = (-y, x, 0).
        return glm::dvec4(-omega * ut_s * x.y, omega * ut_s * x.x, 0.0, ut_s);
    }

    // Plunge at fixed theta with the ISCO shell's u_t and u_phi = sin^2 L_hat.
    const double u_t = ut_s * (g_tt + sin2 * omega * g_tp_s);
    const double l_hat = ut_s * (g_tp_s + omega * g_pp_s);
    const double sigma = r * r + a2 * cos2;
    const double delta = r * r - 2.0 * r + a2;
    const double big_a = (r * r + a2) * (r * r + a2) - a2 * delta * sin2;
    const double sd = sigma * delta;
    const double gi_tt = -big_a / sd;
    const double gi_tp = -2.0 * a * r / sd;
    const double gi_pp_s = (delta - a2 * sin2) / sd; // g^phiphi sin^2
    const double gi_rr = delta / sigma;
    const double up_bl = gi_tp * u_t + gi_pp_s * l_hat;
    const double ut_bl = gi_tt * u_t + gi_tp * sin2 * l_hat;
    const double rest = -1.0 - (gi_tt * u_t * u_t + 2.0 * gi_tp * sin2 * u_t * l_hat + gi_pp_s * sin2 * l_hat * l_hat);
    const double ur = -std::sqrt(std::max(gi_rr * rest, 0.0));

    // Boyer-Lindquist -> Kerr-Schild: dt += 2 r / Delta dr, dphi += a / Delta dr; the
    // Cartesian position moves by l dr (l_x = sin cos phi, ...) and (-y, x, 0) dphi.
    const double ut_ks = ut_bl + 2.0 * r / delta * ur;
    const double up_ks = up_bl + a / delta * ur;
    const double s = r * r + a2;
    const glm::dvec3 l((r * x.x + a * x.y) / s, (r * x.y - a * x.x) / s, x.z / r);
    return glm::dvec4(l * ur + up_ks * glm::dvec3(-x.y, x.x, 0.0), ut_ks);
}

double hot_spot_period_s(double a, double radius_m, double mass_km)
{
    return kTwoPi * (radius_m * std::sqrt(radius_m) + a) * mass_km / (kLightCm * 1e-5);
}

glm::dvec3 hot_spot_center(const HotSpot& s, double a, double t_m)
{
    const double phi = s.phase_rad + t_m / (s.radius_m * std::sqrt(s.radius_m) + a);
    const double c = std::cos(phi);
    const double sn = std::sin(phi);
    return glm::dvec3(s.radius_m * c - a * sn, s.radius_m * sn + a * c, 0.0);
}

namespace {

double coordinate_time_rate(const Photon& ph, double a)
{
    glm::dvec3 dx;
    glm::dvec3 dp;
    double dt = 0.0;
    geodesic_derivatives(ph, a, dx, dp, dt);
    return dt;
}

} // namespace

FlowRay trace_flow_ray(Photon ph, double a, double mass_cm, double nu_obs_hz, const RiafParams& p,
                       double escape_radius, int max_steps, const HotSpot* spot)
{
    FlowRay out;
    const double r_plus = horizon_radius(a);
    double transmit = 1.0;
    const bool with_spot = spot && spot->strength > 0.0;
    const double camera_distance = glm::length(ph.x);
    double ray_time = 0.0; // coordinate time since leaving the camera (negative: into the past)
    double rate = with_spot ? coordinate_time_rate(ph, a) : 0.0;
    for (int i = 0; i < max_steps; ++i) {
        const double h = ray_step(ph.x, a);
        rk4_step(ph, a, h);
        if (ray_runaway(ph, a)) {
            out.outcome = RayOutcome::Captured; // the step is garbage: add no light from it
            return out;
        }
        if (with_spot) {
            const double next_rate = coordinate_time_rate(ph, a);
            ray_time += 0.5 * h * (rate + next_rate);
            rate = next_rate;
        }

        const RiafPlasma gas = riaf_plasma(ph.x, a, p);
        // The spot: within 5 sigma of where its centre was when this light left it.
        double spot_d2 = 0.0;
        bool near_spot = false;
        if (with_spot) {
            const glm::dvec3 d = ph.x - hot_spot_center(*spot, a, ray_time + camera_distance);
            spot_d2 = glm::dot(d, d);
            near_spot = spot_d2 < 25.0 * spot->sigma_m * spot->sigma_m;
        }
        if (near_spot) {
            const glm::dvec4 u = riaf_four_velocity(ph.x, a);
            const double pu = ph.pt * u.w + glm::dot(ph.p, glm::dvec3(u));
            const double sigma = spot->sigma_m;
            const double dtau = spot->strength * spot->tau0 * std::exp(-0.5 * spot_d2 / (sigma * sigma)) /
                                (sigma * std::sqrt(kTwoPi)) * pu * h;
            out.tb_k += transmit * spot->t_s_k / pu * -std::expm1(-dtau);
            transmit *= std::exp(-dtau);
        }
        if (gas.n_e_cm3 > 0.0) {
            const glm::dvec4 u = riaf_four_velocity(ph.x, a);
            const double pu = ph.pt * u.w + glm::dot(ph.p, glm::dvec3(u)); // nu_emit / nu_obs
            const double g = 1.0 / pu;
            const double nu = nu_obs_hz * pu;
            const double j = thermal_synchrotron_emissivity(nu, gas.n_e_cm3, gas.t_e_k, gas.b_gauss, kRiafSinPitch);
            const double alpha = j * kLightCm * kLightCm / (2.0 * kBoltzmann * nu * nu * gas.t_e_k);
            const double dtau = alpha * pu * h * mass_cm;
            const double source = g * gas.t_e_k;
            out.tb_k += transmit * source * -std::expm1(-dtau);
            transmit *= std::exp(-dtau);
            if (dtau > 0.0) {
                out.max_source_k = std::max(out.max_source_k, source);
                if (out.tau < 1.0 && out.tau + dtau >= 1.0) {
                    out.tau1_source_k = source;
                }
            }
            out.tau += dtau;
        }

        if (kerr_schild_r(ph.x, a) < r_plus + 0.02) {
            out.outcome = RayOutcome::Captured;
            return out;
        }
        if (glm::length(ph.x) > escape_radius && glm::dot(ph.x, ph.p) > 0.0) {
            out.outcome = RayOutcome::Escaped;
            return out;
        }
    }
    return out;
}

} // namespace astraxis
