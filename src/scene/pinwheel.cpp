#include "scene/pinwheel.hpp"

#include "core/color.hpp"
#include "core/math.hpp"
#include "scene/scene.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

// Mean anomaly at true anomaly nu (radians) on an orbit of eccentricity e.
double mean_anomaly_at(double nu, double e)
{
    const double ecc = 2.0 * std::atan2(std::sqrt(1.0 - e) * std::sin(0.5 * nu), std::sqrt(1.0 + e) * std::cos(0.5 * nu));
    return ecc - e * std::sin(ecc);
}

// The secondary's orbit around the primary, from their states at t_tdb and
// the pair's G(m1 + m2): eccentricity, periastron direction and its normal in
// the orbit plane (90 deg ahead), mean motion and the time of periastron.
struct PairOrbit {
    double e = 0.0;
    glm::dvec3 p{1.0, 0.0, 0.0};
    glm::dvec3 q{0.0, 1.0, 0.0};
    double n = 0.0;
    double t_peri = 0.0;
};

PairOrbit pair_orbit(const Scene& scene, const Pinwheel& pinwheel, double t_tdb)
{
    const double gm = scene.bodies[static_cast<size_t>(pinwheel.primary)].body_gm_km3_s2 +
                      scene.bodies[static_cast<size_t>(pinwheel.secondary)].body_gm_km3_s2;
    const State a = scene.icrf_state_at(pinwheel.primary, t_tdb);
    const State b = scene.icrf_state_at(pinwheel.secondary, t_tdb);
    const glm::dvec3 r = b.position - a.position;
    const glm::dvec3 v = b.velocity - a.velocity;
    const glm::dvec3 h = glm::cross(r, v);
    const glm::dvec3 e_vec = glm::cross(v, h) / gm - glm::normalize(r);
    PairOrbit o;
    o.e = glm::length(e_vec);
    o.p = o.e > 1e-12 ? e_vec / o.e : glm::normalize(r);
    o.q = glm::cross(glm::normalize(h), o.p);
    const double a_km = 1.0 / (2.0 / glm::length(r) - glm::dot(v, v) / gm);
    o.n = std::sqrt(gm / (a_km * a_km * a_km));
    const double nu = std::atan2(glm::dot(r, o.q), glm::dot(r, o.p));
    o.t_peri = t_tdb - mean_anomaly_at(nu, o.e) / o.n;
    return o;
}

// Han et al. 2022, eq. 1, beyond the end of the constant acceleration
// (radius r_t, speed v_t): a = a_max (r_t / r)^2, so v^2 = A - B / r with
// A = v_t^2 + 2 a_max r_t and B = 2 a_max r_t^2, and the time to reach r from
// r_t is F(r) - F(r_t) with
//   F(r) = sqrt(r (A r - B)) / A + B / A^1.5 ln(sqrt(A r) + sqrt(A r - B)).
double coast_time(double r, double a, double b)
{
    return std::sqrt(r * (a * r - b)) / a + b / std::pow(a, 1.5) * std::log(std::sqrt(a * r) + std::sqrt(a * r - b));
}

} // namespace

double Pinwheel::step_km() const
{
    return pinwheel_radius(*this, turns * period_s) - pinwheel_radius(*this, (turns - 1) * period_s);
}

double Pinwheel::window_rad() const
{
    return windowed ? (window_off_deg - window_on_deg) * kDegToRad : kTwoPi;
}

void build_pinwheel_grains(const Scene& scene, Pinwheel& pinwheel, int rings_per_turn, int grains_per_ring)
{
    pinwheel.rings_per_turn = rings_per_turn;
    pinwheel.grains_per_ring = grains_per_ring;
    pinwheel.grains.clear();
    pinwheel.grains.reserve(static_cast<size_t>(pinwheel.turns) * static_cast<size_t>(rings_per_turn) *
                            static_cast<size_t>(grains_per_ring));
    const double half = 0.5 * pinwheel.opening_angle_deg * kDegToRad;

    // When each ring leaves (as a fraction of the period after t_ref) and how
    // much dust it carries: evenly in time all around the orbit, or evenly in
    // true anomaly over the window, each ring then standing for the time the
    // orbit takes to sweep its share, dt/dnu ~ (1 + e cos nu)^-2.
    std::vector<double> phases(static_cast<size_t>(rings_per_turn));
    std::vector<double> weights(static_cast<size_t>(rings_per_turn), 1.0);
    if (pinwheel.windowed) {
        const PairOrbit orbit = pair_orbit(scene, pinwheel, pinwheel.t_ref_tdb);
        double heaviest = 0.0;
        for (int k = 0; k < rings_per_turn; ++k) {
            const double nu_deg = pinwheel.window_on_deg +
                                  (k + 0.5) / rings_per_turn * (pinwheel.window_off_deg - pinwheel.window_on_deg);
            const double nu = nu_deg * kDegToRad;
            const double t = orbit.t_peri + mean_anomaly_at(nu, orbit.e) / orbit.n;
            const double x = (t - pinwheel.t_ref_tdb) / pinwheel.period_s;
            phases[static_cast<size_t>(k)] = x - std::floor(x);
            const double sweep = 1.0 / ((1.0 + orbit.e * std::cos(nu)) * (1.0 + orbit.e * std::cos(nu)));
            weights[static_cast<size_t>(k)] = pinwheel_orbit_rate(pinwheel, nu_deg) * sweep;
            heaviest = std::max(heaviest, weights[static_cast<size_t>(k)]);
        }
        for (double& w : weights) {
            w /= heaviest;
        }
    } else {
        for (int k = 0; k < rings_per_turn; ++k) {
            phases[static_cast<size_t>(k)] = (k + 0.5) / rings_per_turn;
        }
    }

    pinwheel.rings.clear();
    for (int k = 0; k < rings_per_turn; ++k) {
        pinwheel.rings.push_back({phases[static_cast<size_t>(k)], weights[static_cast<size_t>(k)]});
    }

    // Golden-ratio offsets of the grains around each ring, so that grains of
    // neighboring rings and turns do not line up into streaks.
    constexpr double kGolden = 0.6180339887498949;
    for (int k = 0; k < rings_per_turn; ++k) {
        const double phase = phases[static_cast<size_t>(k)];
        const double t = pinwheel.t_ref_tdb + phase * pinwheel.period_s;
        const State p = scene.icrf_state_at(pinwheel.primary, t);
        const State s = scene.icrf_state_at(pinwheel.secondary, t);
        // The cone (White et al. 2025, eqs. B9-B12): its axis from the primary
        // to the secondary, the ring around it in the plane of the axis' normal
        // in the orbit plane (`side`, toward the orbital motion: the leading
        // edge) and the orbit normal.
        const glm::dvec3 rel = s.position - p.position;
        const glm::dvec3 axis = glm::normalize(rel);
        const glm::dvec3 normal = glm::normalize(glm::cross(rel, s.velocity - p.velocity));
        const glm::dvec3 side = glm::cross(normal, axis);
        for (int n = 0; n < pinwheel.turns; ++n) {
            const double offset = std::fmod((k + n * rings_per_turn) * kGolden, 1.0);
            for (int j = 0; j < grains_per_ring; ++j) {
                const double psi = kTwoPi * (j + offset) / grains_per_ring; // the leading edge at 90 deg
                const glm::dvec3 d = axis * std::cos(half) +
                                     (side * std::sin(psi) + normal * std::cos(psi)) * std::sin(half);
                PinwheelGrain g;
                g.direction = glm::vec3(d);
                g.tangent = glm::vec3(side * std::cos(psi) - normal * std::sin(psi));
                g.phase = static_cast<float>(phase);
                g.turn = static_cast<float>(n);
                g.weight = static_cast<float>(weights[static_cast<size_t>(k)] *
                                              pinwheel_azimuth_rate(pinwheel, psi * kRadToDeg - 90.0));
                pinwheel.grains.push_back(g);
            }
        }
    }
}

double pinwheel_radius(const Pinwheel& pinwheel, double age_s)
{
    const double v0 = pinwheel.speed_km_s;
    const double a = pinwheel.acceleration_km_s2;
    const double r_nuc = pinwheel.nucleation_radius_km;
    const double t_nuc = r_nuc / v0;
    if (a <= 0.0 || age_s <= t_nuc) {
        return v0 * age_s;
    }
    // Constant acceleration from r_nuc to r_t.
    const double r_t = std::max(pinwheel.acceleration_end_km, r_nuc);
    const double v_t = std::sqrt(v0 * v0 + 2.0 * a * (r_t - r_nuc));
    const double tau = age_s - t_nuc;
    const double tau_t = (v_t - v0) / a;
    if (tau <= tau_t) {
        return r_nuc + v0 * tau + 0.5 * a * tau * tau;
    }
    // Beyond r_t: solve F(r) - F(r_t) = tau - tau_t (see coast_time) by Newton's
    // method, dt/dr = 1 / v, kept within the bracket of the slowest and
    // fastest speeds (v_t, sqrt(A)).
    const double big_a = v_t * v_t + 2.0 * a * r_t;
    const double big_b = 2.0 * a * r_t * r_t;
    const double dt = tau - tau_t;
    const double target = coast_time(r_t, big_a, big_b) + dt;
    double lo = r_t + v_t * dt;
    double hi = r_t + std::sqrt(big_a) * dt;
    double r = 0.5 * (lo + hi);
    for (int i = 0; i < 60; ++i) {
        const double f = coast_time(r, big_a, big_b) - target;
        (f < 0.0 ? lo : hi) = r;
        double next = r - f * std::sqrt(big_a - big_b / r);
        if (next <= lo || next >= hi) {
            next = 0.5 * (lo + hi);
        }
        if (std::abs(next - r) <= 1e-12 * r) {
            return next;
        }
        r = next;
    }
    return r;
}

double pinwheel_grain_age(const Pinwheel& pinwheel, const PinwheelGrain& grain, double t_tdb)
{
    const double n = static_cast<double>(pinwheel.turns);
    const double y = ((t_tdb - pinwheel.t_ref_tdb) / pinwheel.period_s - static_cast<double>(grain.phase) -
                      static_cast<double>(grain.turn)) /
                     n;
    return n * pinwheel.period_s * (y - std::floor(y));
}

glm::dvec3 pinwheel_grain_offset(const Pinwheel& pinwheel, const PinwheelGrain& grain, double t_tdb)
{
    return glm::dvec3(grain.direction) * pinwheel_radius(pinwheel, pinwheel_grain_age(pinwheel, grain, t_tdb));
}

double pinwheel_orbit_rate(const Pinwheel& pinwheel, double nu_deg)
{
    if (!pinwheel.windowed) {
        return 1.0;
    }
    if (nu_deg < pinwheel.window_on_deg || nu_deg > pinwheel.window_off_deg) {
        return 0.0;
    }
    const double s = pinwheel.modulation_sigma_deg;
    if (s <= 0.0) {
        return 1.0;
    }
    // Han et al. 2022, Methods: two half-Gaussians peaking where production
    // turns on and off, with a minimum at periastron in between.
    const double on = (nu_deg - pinwheel.window_on_deg) / s;
    const double off = (pinwheel.window_off_deg - nu_deg) / s;
    return std::exp(-0.5 * on * on) + std::exp(-0.5 * off * off);
}

double pinwheel_azimuth_rate(const Pinwheel& pinwheel, double psi_deg)
{
    const double s = pinwheel.trailing_sigma_deg;
    if (s <= 0.0) {
        return 1.0;
    }
    const double d = wrap_pi((psi_deg - 180.0) * kDegToRad) * kRadToDeg / s;
    return std::exp(-0.5 * d * d);
}

double pinwheel_peak_brightness(const Pinwheel& pinwheel, double t_tdb)
{
    const double reference =
        pinwheel.display_reference_km > 0.0 ? pinwheel.display_reference_km : pinwheel.nucleation_radius_km;
    const double x = (t_tdb - pinwheel.t_ref_tdb) / pinwheel.period_s;
    double peak = 0.0;
    for (const Pinwheel::Ring& ring : pinwheel.rings) {
        const double y = x - ring.phase;
        const double age = (y - std::floor(y)) * pinwheel.period_s;
        const double r = std::max(pinwheel_radius(pinwheel, age), pinwheel.nucleation_radius_km);
        const double q = reference / r;
        peak = std::max(peak, ring.weight * std::pow(q * q * q * q, pinwheel.display_gamma));
    }
    return peak;
}

double pinwheel_exposure(const Pinwheel& pinwheel, double t_tdb, double camera_km)
{
    // WR 140's newest shell peaks at thousands times the reference: dimming
    // it all the way to kPinwheelPeak would leave the older shells ~1500 times
    // darker, gone until the next periastron.
    constexpr double kPinwheelPeak = 2.0;
    constexpr double kPinwheelMaxDimming = 20.0;
    const double peak = pinwheel_peak_brightness(pinwheel, t_tdb);
    const double inside = std::max(pinwheel.step_km() / std::max(camera_km, 1.0), 1.0);
    return std::clamp(kPinwheelPeak / peak, 1.0 / (kPinwheelMaxDimming * inside * inside), 1.0);
}

double pinwheel_dust_temperature(const Pinwheel& pinwheel, double r_km)
{
    const double r = std::max(r_km, pinwheel.nucleation_radius_km);
    return pinwheel.nucleation_temperature_k * std::pow(r / pinwheel.nucleation_radius_km, -pinwheel.temperature_index);
}

double pinwheel_grain_brightness(const Pinwheel& pinwheel, double r_km)
{
    const double r_nuc = pinwheel.nucleation_radius_km;
    if (r_km <= r_nuc) {
        return 0.0;
    }
    const double fade_in = std::min((r_km - r_nuc) / (0.25 * r_nuc), 1.0);
    const double end = pinwheel_radius(pinwheel, pinwheel.turns * pinwheel.period_s);
    const double fade_out = std::clamp((end - r_km) / pinwheel.step_km(), 0.0, 1.0);
    const double q = r_nuc / r_km;
    return q * q * fade_in * fade_out;
}

PinwheelVisibleLight pinwheel_visible_light(const Pinwheel& pinwheel, double r_km, double star_visible_share)
{
    const double t = pinwheel_dust_temperature(pinwheel, r_km);
    PinwheelVisibleLight light;
    light.thermal = (1.0 - pinwheel.visible_albedo) * blackbody_band_share(t, kVisibleLoNm, kVisibleHiNm);
    light.scattered = pinwheel.visible_albedo * star_visible_share;
    return light;
}

double pinwheel_star_visible_share(const Scene& scene, const Pinwheel& pinwheel)
{
    double visible = 0.0;
    double total = 0.0;
    for (int star : {pinwheel.primary, pinwheel.secondary}) {
        const Body& b = scene.bodies[static_cast<size_t>(star)];
        visible += b.luminosity_solar * blackbody_band_share(b.temperature_k, kVisibleLoNm, kVisibleHiNm);
        total += b.luminosity_solar;
    }
    return total > 0.0 ? visible / total : 0.0;
}

} // namespace astraxis
