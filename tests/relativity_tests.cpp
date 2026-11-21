// Tests for relativistic motion: Kerr geodesics, binary pulsars, Sgr A*, Gargantua.

#include "test_util.hpp"

#include "scene/flares.hpp"
#include "scene/lensing.hpp"
#include "scene/pulsar.hpp"
#include "scene/simulation.hpp"

#include "core/color.hpp"
#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/jet.hpp"
#include "ephem/kerr_lens.hpp"
#include "ephem/kerr_null.hpp"
#include "ephem/kerr_orbit.hpp"
#include "ephem/post_keplerian.hpp"
#include "ephem/riaf.hpp"
#include "ephem/schwarzschild_light.hpp"
#include "ephem/visual_orbit.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

using namespace astraxis;

namespace {

// Periastron advance and orbital decay of general relativity at first
// post-Newtonian order for two point masses (M_sun), e.g. Lorimer & Kramer,
// Handbook of Pulsar Astronomy (2005), eqs. 8.53 and 8.56.
double gr_omega_dot_deg_per_year(double m1, double m2, double pb_days, double e)
{
    const double t_sun = kSunGmKm3S2 / (kSpeedOfLightKmS * kSpeedOfLightKmS * kSpeedOfLightKmS); // s
    const double n = kTwoPi / (pb_days * kSecondsPerDay);
    const double rate = 3.0 * std::pow(t_sun * (m1 + m2) * n, 2.0 / 3.0) * n / (1.0 - e * e);
    return rate * kRadToDeg * kDaysPerJulianYear * kSecondsPerDay;
}

double gr_pb_dot(double m1, double m2, double pb_days, double e)
{
    const double t_sun = kSunGmKm3S2 / (kSpeedOfLightKmS * kSpeedOfLightKmS * kSpeedOfLightKmS);
    const double n = kTwoPi / (pb_days * kSecondsPerDay);
    const double e2 = e * e;
    const double f = (1.0 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e2 * e2) / std::pow(1.0 - e2, 3.5);
    return -192.0 * kPi / 5.0 * std::pow(t_sun * n, 5.0 / 3.0) * f * m1 * m2 / std::cbrt(m1 + m2);
}

// The post-Keplerian orbit inside a split pair (ScaledMotion), or a ghost's ellipse.
const PostKeplerianMotion* post_keplerian_of(const Scene& scene, int body)
{
    const auto* scaled = dynamic_cast<const ScaledMotion*>(scene.bodies[static_cast<size_t>(body)].motion.get());
    return scaled ? dynamic_cast<const PostKeplerianMotion*>(&scaled->inner()) : nullptr;
}

const VisualOrbitMotion* ghost_orbit_of(const Scene& scene, int body)
{
    const auto* scaled = dynamic_cast<const ScaledMotion*>(scene.bodies[static_cast<size_t>(body)].motion.get());
    return scaled ? dynamic_cast<const VisualOrbitMotion*>(&scaled->inner()) : nullptr;
}

// A pulsar's distance behind the system's barycentre (the root) along the line of sight, in lt-s.
double pulsar_depth_s(const Scene& scene, int pulsar, const glm::dvec3& away, double t)
{
    return glm::dot(scene.icrf_state_at(pulsar, t).position - scene.icrf_state_at(0, t).position, away) /
           kSpeedOfLightKmS;
}

// Exact Schwarzschild periapsis advance per orbit for turning points r_p, r_a
// (units of M): (du/dphi)^2 = 2 (u1 - u)(u - u2)(u3 - u), u3 = 1/2 - u1 - u2;
// with u = u2 + (u1 - u2) sin^2(theta) the period in phi is
// 4 * integral_0^{pi/2} dtheta / sqrt(2 (u3 - u2 - (u1 - u2) sin^2 theta)).
double schwarzschild_advance(double r_p, double r_a)
{
    const double u1 = 1.0 / r_p;
    const double u2 = 1.0 / r_a;
    const double u3 = 0.5 - u1 - u2;
    const int n = 20000;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        const double th = (i + 0.5) / n * 0.5 * kPi;
        const double s = std::sin(th);
        sum += 1.0 / std::sqrt(2.0 * (u3 - u2 - (u1 - u2) * s * s));
    }
    return 4.0 * sum * (0.5 * kPi / n) - kTwoPi;
}

void test_kerr_conservation()
{
    // Strong field, rapidly spinning hole: energy and angular momentum stay constant.
    const double m = 10.0;
    const KerrEquatorial kerr(m, 0.9 * m);
    glm::dvec2 x(40.0 * m, 0.0);
    glm::dvec2 u = kerr.u_from_velocity(x, glm::dvec2(0.0, 0.17));
    const double e0 = kerr.energy(x, u);
    const double l0 = kerr.angular_momentum(x, u);
    // IMR is symplectic: the energy error is O(dt^2) and bounded (no drift).
    double worst_e[2] = {0.0, 0.0};
    double worst_l = 0.0;
    const int steps = 40000;
    for (int i = 0; i < steps; ++i) {
        const double r = glm::length(x);
        check(kerr.step(x, u, 0.01 * std::sqrt(r * r * r / m)), "Kerr IMR Newton converges");
        double& w = worst_e[i < steps / 2 ? 0 : 1];
        w = std::max(w, std::abs(kerr.energy(x, u) / e0 - 1.0));
        worst_l = std::max(worst_l, std::abs(kerr.angular_momentum(x, u) / l0 - 1.0));
    }
    check(worst_e[1] < 1e-6, "Kerr energy error small", worst_e[1]);
    check(worst_e[1] < 1.5 * worst_e[0] + 1e-12, "Kerr energy error bounded (no drift)", worst_e[1] / worst_e[0]);
    check(worst_l < 1e-6, "Kerr angular momentum conserved", worst_l);
    check(kerr.horizon_radius() > m && kerr.horizon_radius() < 2.0 * m, "Kerr horizon radius", kerr.horizon_radius());
}

// A tilted orbit around a fast-spinning hole keeps its energy, axial angular
// momentum and Carter constant (particle_step, the integrator of tilted orbits).
void test_kerr_tilted_conservation()
{
    const double a = 0.94;
    Photon s = make_particle(glm::dvec3(30.0, 0.0, 0.0), glm::dvec3(0.0, 0.12, 0.12), a);
    const double e0 = -s.pt;
    const double l0 = s.x.x * s.p.y - s.x.y * s.p.x;
    const double q0 = carter_constant(s, a);
    check(q0 > 1.0, "tilted orbit has a Carter constant", q0);
    double t = 0.0;
    double worst_e = 0.0;
    double worst_l = 0.0;
    double worst_q = 0.0;
    double z_max = 0.0;
    for (int i = 0; i < 200000 && t < 20.0 * kTwoPi * std::pow(30.0, 1.5); ++i) {
        const double r = kerr_schild_r(s.x, a);
        particle_step(s, t, a, 0.003 * r * std::sqrt(r));
        worst_e = std::max(worst_e, std::abs(-s.pt / e0 - 1.0));
        worst_l = std::max(worst_l, std::abs((s.x.x * s.p.y - s.x.y * s.p.x) / l0 - 1.0));
        if (std::abs(s.x.z) < 0.9 * kerr_schild_r(s.x, a)) { // Q's formula divides by sin^2 theta
            worst_q = std::max(worst_q, std::abs(carter_constant(s, a) / q0 - 1.0));
        }
        z_max = std::max(z_max, std::abs(s.x.z));
    }
    check(z_max > 10.0, "tilted orbit leaves the equatorial plane", z_max);
    check(worst_e < 1e-9, "tilted Kerr orbit: energy conserved", worst_e);
    check(worst_l < 1e-7, "tilted Kerr orbit: axial angular momentum conserved", worst_l);
    check(worst_q < 1e-6, "tilted Kerr orbit: Carter constant conserved", worst_q);
}

// Without spin the plane does not matter: the 3D integrator about any axis follows
// the equatorial integrator (independent code, Boyer-Lindquist-like coordinates).
void test_tilted_orbit_matches_equatorial()
{
    KerrOrbitSetup setup;
    setup.gm_km3_s2 = 4.261e6 * kSunGmKm3S2;
    setup.elements.a_km = 0.125058 * 8246.7 * kAuKm;
    setup.elements.e = 0.884649;
    setup.elements.i_deg = 134.567;
    setup.elements.node_deg = 228.171;
    setup.elements.arg_peri_deg = 66.263;
    setup.elements.t_peri_tdb = tdb_from_julian_year(2018.379);
    setup.elements.ra_deg = 266.41683;
    setup.elements.dec_deg = -29.00782;
    const double period = kTwoPi * std::sqrt(std::pow(setup.elements.a_km, 3) / setup.gm_km3_s2);
    setup.elements_epoch_tdb = setup.elements.t_peri_tdb - 0.5 * period;
    const double t_end = setup.elements_epoch_tdb + 1.2 * period;
    const KerrOrbitMotion flat(setup, setup.elements_epoch_tdb, t_end);
    KerrOrbitSetup tilted = setup;
    tilted.tilted = true;
    tilted.pole = glm::normalize(glm::dvec3(0.3, 0.4, 0.866));
    const KerrOrbitMotion three_d(tilted, setup.elements_epoch_tdb, t_end);
    check(std::abs(three_d.table().end() - t_end) < 0.01 * period, "tilted orbit covers the span");
    double worst = 0.0;
    for (double t = setup.elements_epoch_tdb; t < t_end; t += 0.01 * period) {
        const glm::dvec3 a = flat.eval(t).position;
        const glm::dvec3 b = three_d.eval(t).position;
        worst = std::max(worst, glm::length(a - b) / glm::length(a));
    }
    check(worst < 2e-4, "tilted integrator = equatorial integrator at a = 0", worst);
}

void test_kerr_schwarzschild_precession()
{
    // S2-like orbit (GRAVITY 2020 elements) around a non-spinning hole.
    KerrOrbitSetup setup;
    setup.gm_km3_s2 = 4.261e6 * kSunGmKm3S2;
    setup.spin = 0.0;
    setup.elements.a_km = 0.125058 * 8246.7 * kAuKm;
    setup.elements.e = 0.884649;
    setup.elements.i_deg = 134.567;
    setup.elements.node_deg = 228.171;
    setup.elements.arg_peri_deg = 66.263;
    setup.elements.t_peri_tdb = tdb_from_julian_year(2018.379);
    setup.elements.ra_deg = 266.41683;
    setup.elements.dec_deg = -29.00782;
    const double period = kTwoPi * std::sqrt(std::pow(setup.elements.a_km, 3) / setup.gm_km3_s2);
    setup.elements_epoch_tdb = setup.elements.t_peri_tdb - 0.5 * period; // 2010 apocentre

    const double t_end = setup.elements_epoch_tdb + 4.2 * period;
    const KerrOrbitMotion orbit(setup, setup.elements_epoch_tdb, t_end);
    const auto& knots = orbit.table().knots();
    // A failed Newton solve ends the trajectory: the whole span must be covered.
    check(std::abs(knots.back().t - t_end) < 1.0, "S2 integration covers the span",
          (t_end - knots.back().t) / kSecondsPerDay);

    // Every Newton solve converges at the scene's step size (the residual floor
    // is ~1e-13 here; a tolerance below it would fail ~3% of steps).
    const KerrEquatorial& kerr = orbit.metric();
    const glm::dvec3 P = orbit.plane()[0];
    const glm::dvec3 Q = orbit.plane()[1];
    {
        const double m = kerr.mass();
        glm::dvec2 x(glm::dot(knots.front().position, P), glm::dot(knots.front().position, Q));
        glm::dvec2 u = kerr.u_from_velocity(
            x, glm::dvec2(glm::dot(knots.front().velocity, P), glm::dot(knots.front().velocity, Q)) /
                   kSpeedOfLightKmS);
        int failures = 0;
        for (double t = 0.0; t < period * kSpeedOfLightKmS;) {
            const double r = glm::length(x);
            const double dt = 0.003 * std::sqrt(r * r * r / m);
            failures += kerr.step(x, u, dt) ? 0 : 1;
            t += dt;
        }
        check(failures == 0, "S2 Kerr IMR Newton converges over one orbit", failures);
    }

    // Pericentre longitudes in the orbit plane: the knot with minimal r, refined
    // by a parabola through r(phi) of its neighbours (knots are ~14' apart in
    // angle at pericentre, comparable to the 12' precession being measured).
    auto phi_of = [&](const glm::dvec3& p) { return std::atan2(glm::dot(p, Q), glm::dot(p, P)); };
    std::vector<double> peri;
    for (size_t i = 1; i + 1 < knots.size(); ++i) {
        const double r0 = glm::length(knots[i - 1].position);
        const double r1 = glm::length(knots[i].position);
        const double r2 = glm::length(knots[i + 1].position);
        if (r1 < r0 && r1 <= r2) {
            const double p1 = phi_of(knots[i].position);
            const double a0 = wrap_pi(phi_of(knots[i - 1].position) - p1);
            const double a2 = wrap_pi(phi_of(knots[i + 1].position) - p1);
            // Vertex of the parabola through (a0, r0), (0, r1), (a2, r2).
            const double num = a0 * a0 * (r2 - r1) - a2 * a2 * (r0 - r1);
            const double den = a0 * (r2 - r1) - a2 * (r0 - r1);
            peri.push_back(p1 + 0.5 * num / den);
        }
    }
    check(peri.size() == 4, "S2: 4 pericentre passages", static_cast<double>(peri.size()));
    if (peri.size() >= 2) {
        const double advance = wrap_pi(peri.back() - peri.front()) / static_cast<double>(peri.size() - 1);
        const double m = setup.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
        const double expected =
            6.0 * kPi * m / (setup.elements.a_km * (1.0 - setup.elements.e * setup.elements.e));
        // ~12 arcmin per orbit; 1PN is accurate to O(M/r) ~ 1e-3 here.
        check(std::abs(advance / expected - 1.0) < 0.02, "S2 Schwarzschild precession / 1PN", advance / expected);
        std::printf("info: S2 periapsis advance %.2f arcmin/orbit (1PN %.2f)\n", advance * kRadToDeg * 60.0,
                    expected * kRadToDeg * 60.0);
    }
    // Radial period close to Keplerian (GR correction ~1e-3).
    check(knots.back().t > knots.front().t, "S2 integration span");
}

// The post-Keplerian (DD) orbit: Keplerian without the relativistic terms,
// velocity = d position / dt, the periastron turning by k 2 pi per orbit, the
// quadratic phase of a changing period, and the Newtonian continuation.
void test_post_keplerian()
{
    VisualOrbit el;
    el.e = 0.6;
    el.i_deg = 50.0;
    el.node_deg = 30.0;
    el.arg_peri_deg = 100.0;
    el.t_peri_tdb = 1.0e8;
    el.ra_deg = 100.0;
    el.dec_deg = 20.0;
    const double gm = 2.8 * kSunGmKm3S2;
    const double pb = 0.3 * kSecondsPerDay;

    const PostKeplerianOrbit kepler = make_post_keplerian(el, pb, 0.0, 0.0, gm);
    const double t = el.t_peri_tdb + 12345.6;
    const State a = post_keplerian_state(kepler, gm, t);
    const State b = visual_orbit_state(kepler.elements, gm, t);
    check(glm::length(a.position - b.position) < 1e-6 && glm::length(a.velocity - b.velocity) < 1e-9,
          "post-Keplerian without relativistic terms is Keplerian", glm::length(a.position - b.position));

    // Exaggerated: the periastron turns 0.05 rad per radian of true anomaly.
    const double k = 0.05;
    const PostKeplerianOrbit fast = make_post_keplerian(el, pb, k * kTwoPi / pb, 0.0, gm);
    double worst = 0.0;
    for (double dt : {0.0, 3.0, 1234.5, 0.5 * pb, 7.9 * pb, -2.3 * pb}) {
        const double h = 0.5;
        const State s = post_keplerian_state(fast, gm, el.t_peri_tdb + dt);
        const glm::dvec3 fd = (post_keplerian_state(fast, gm, el.t_peri_tdb + dt + h).position -
                               post_keplerian_state(fast, gm, el.t_peri_tdb + dt - h).position) /
                              (2.0 * h);
        worst = std::max(worst, glm::length(fd - s.velocity) / glm::length(s.velocity));
    }
    check(worst < 1e-6, "post-Keplerian velocity is the derivative of position", worst);

    VisualOrbit turned = fast.elements;
    turned.arg_peri_deg += k * kTwoPi * 7.0 * kRadToDeg;
    const glm::dvec3 peri = glm::normalize(post_keplerian_state(fast, gm, el.t_peri_tdb + 7.0 * pb).position);
    const double off = std::acos(std::min(1.0, glm::dot(peri, visual_orbit_frame(turned)[0])));
    check(off < 1e-7, "periastron turns by 2 pi k per orbit", off);

    // A shrinking period: the N-th periastron comes ~Pb_dot N^2 / 2 periods early
    // (solving N = o - Pb_dot o^2 / 2 for the orbits o to second order in Pb_dot).
    const double pb_dot = -1e-6;
    const PostKeplerianOrbit decaying = make_post_keplerian(el, pb, 0.0, pb_dot, gm);
    const double n = 100.0;
    const double t_n = el.t_peri_tdb + pb * (n + 0.5 * pb_dot * n * n + 0.5 * pb_dot * pb_dot * n * n * n);
    const PostKeplerianPhase at_n = post_keplerian_phase(decaying, gm, t_n);
    check(std::abs(wrap_pi(at_n.mean_anomaly)) < 1e-6, "period decay: quadratic orbital phase", at_n.mean_anomaly);
    const double shrink = std::pow(1.0 + pb_dot * (t_n - el.t_peri_tdb) / pb, 2.0 / 3.0);
    check(std::abs(at_n.a_km / decaying.elements.a_km - shrink) < 1e-12, "semi-major axis follows the period",
          at_n.a_km);

    // Newton from t_ref: through the same point, with that moment's periastron, forever.
    const double t_ref = el.t_peri_tdb + 3.3 * pb;
    const VisualOrbit newton = newtonian_continuation(fast, gm, t_ref);
    const double moved = glm::length(visual_orbit_state(newton, gm, t_ref).position -
                                     post_keplerian_state(fast, gm, t_ref).position);
    check(moved < 1e-5, "Newtonian continuation starts at the same point", moved);
    check(std::abs(newton.arg_peri_deg * kDegToRad - post_keplerian_phase(fast, gm, t_ref).arg_peri_rad) < 1e-12,
          "Newtonian continuation keeps that moment's periastron");

    // The drawn ellipse starts at the body and closes.
    const PostKeplerianMotion motion(fast, gm);
    std::vector<glm::dvec3> ring;
    motion.sample_orbit(t, 64, ring);
    check(ring.size() == 64 && glm::length(ring.front() - motion.eval(t).position) < 1e-6 &&
              glm::length(ring.back() - ring.front()) < 1e-6,
          "post-Keplerian orbit line starts at the body and closes");
}

void test_psr_b1913_scene()
{
    Scene scene = load_scene_or_die("psr_b1913.toml");
    check_model_sky(scene, "PSR B1913+16 model sky");
    const int pulsar = scene.find("PSR B1913+16");
    const int companion = scene.find("B1913+16 companion");
    const int ghost_p = scene.find("PSR B1913+16 (Newton)");
    const int ghost_c = scene.find("Companion (Newton)");
    check(pulsar > 0 && companion > 0 && ghost_p > 0 && ghost_c > 0, "PSR B1913+16 bodies");
    if (pulsar <= 0 || companion <= 0 || ghost_p <= 0 || ghost_c <= 0) {
        return;
    }
    const PostKeplerianMotion* orbit = post_keplerian_of(scene, pulsar);
    const VisualOrbitMotion* newton = ghost_orbit_of(scene, ghost_p);
    check(orbit && newton, "PSR B1913+16 post-Keplerian orbit and Newtonian ghost");
    if (!orbit || !newton) {
        return;
    }

    // Roemer delay at periastron T0: x (1 - e) sin(omega0) [WH16 Table 2].
    const glm::dvec3 away = unit_toward(288.86666425, 16.10760744);
    const double x = 2.341776;
    const double e = 0.6171340;
    const double w = 292.54450 * kDegToRad;
    const double t0 = tdb_from_jd_tdb(2452145.40097849);
    const double depth = pulsar_depth_s(scene, pulsar, away, t0);
    check(std::abs(depth - x * (1.0 - e) * std::sin(w)) < 2e-5, "PSR B1913+16 Roemer delay at T0 (s)", depth);
    // Over one orbit the line-of-sight depth spans 2 x sqrt(1 - e^2 cos^2 omega).
    double z_min = 1e30;
    double z_max = -1e30;
    for (int i = 0; i < 4000; ++i) {
        const double z = pulsar_depth_s(scene, pulsar, away, t0 + 0.322997448918 * kSecondsPerDay * i / 4000.0);
        z_min = std::min(z_min, z);
        z_max = std::max(z_max, z);
    }
    const double span = 2.0 * x * std::sqrt(1.0 - e * e * std::cos(w) * std::cos(w));
    check(std::abs((z_max - z_min) / span - 1.0) < 1e-4, "PSR B1913+16 projected semi-major axis", z_max - z_min);

    // The measured advance and decay against GR for the [WH16] masses.
    const double omdot_gr = gr_omega_dot_deg_per_year(1.438, 1.390, 0.322997448918, e);
    check(std::abs(omdot_gr / 4.226585 - 1.0) < 2e-4, "PSR B1913+16 omega_dot is GR's for its masses", omdot_gr);
    const double pbdot_gr = gr_pb_dot(1.438, 1.390, 0.322997448918, e);
    check(std::abs(pbdot_gr / -2.40263e-12 - 1.0) < 1e-3, "PSR B1913+16 Pb_dot(GR) [WH16 Sec. 7]", pbdot_gr);
    std::printf("info: PSR B1913+16 GR: omega_dot %.5f deg/yr, Pb_dot %.5e\n", omdot_gr, pbdot_gr);

    // The Newtonian ghost leaves the real orbit on the day of discovery; by
    // 2026-10 the real periastron has turned ~221 deg away from it.
    double t_disc = 0.0;
    check(parse_utc("1974-07-02", &t_disc), "parse the discovery date");
    const double split = glm::length(scene.icrf_state_at(pulsar, t_disc).position -
                                     scene.icrf_state_at(ghost_p, t_disc).position);
    check(split < 1e-3, "PSR B1913+16 ghost starts on the real orbit (km)", split);
    const double t_now = tdb_from_julian_year(2026.76);
    const double turned = (post_keplerian_phase(orbit->orbit(), 3.7531078044e11, t_now).arg_peri_rad -
                           newton->orbit().arg_peri_deg * kDegToRad) *
                          kRadToDeg;
    const double expected = 4.226585 * (2026.76 - (2000.0 + t_disc / (kDaysPerJulianYear * kSecondsPerDay)));
    check(std::abs(turned - expected) < 0.01, "PSR B1913+16 periastron advance since 1974 (deg)", turned);
    std::printf("info: PSR B1913+16 periastron advance since discovery: %.2f deg\n", turned);

    // Ghosts pair among themselves and carry no mass; the pulsar's beams turn
    // about a tilted spin axis precessing around the real orbit normal.
    check(scene.satellite_host(companion) == pulsar && scene.satellite_host(ghost_c) == ghost_p &&
              scene.satellite_host(ghost_p) == -1,
          "PSR B1913+16 hosts (ghosts pair with ghosts)");
    check(scene.bodies[static_cast<size_t>(ghost_p)].body_gm_km3_s2 == 0.0 &&
              std::abs(scene.bodies[static_cast<size_t>(pulsar)].body_gm_km3_s2 / kSunGmKm3S2 - 1.438) < 1e-6,
          "PSR B1913+16 masses (ghosts massless)");
    check(scene.bodies[static_cast<size_t>(pulsar)].pulsar.spin_axis_orbit_of == companion,
          "PSR B1913+16 precession uses the companion orbit");

    Body& body = scene.bodies[static_cast<size_t>(pulsar)];
    const auto& prec = body.pulsar.precession;
    const double year_s = kDaysPerJulianYear * kSecondsPerDay;
    const double rate = prec.rate_rad_s * year_s * kRadToDeg;
    check(prec.enabled && std::abs(rate - 1.213) < 0.003,
          "B1913 spin precession against WT02 1.213 deg/yr (older masses)", rate);
    // Independent relation to the GR apsidal advance: Omega_spin / omega_dot
    // = m_c (4 m_p + 3 m_c) / [6 (m_p + m_c)^2].
    const double rate_from_apsides = omdot_gr * 1.390 * (4.0 * 1.438 + 3.0 * 1.390) /
                                      (6.0 * 2.828 * 2.828);
    check(std::abs(rate - rate_from_apsides) < 1e-9, "B1913 precession rate from the scene masses", rate);
    const double inclination = 47.1254 * kDegToRad;
    const double tilt = 21.1 * kDegToRad;
    check(std::abs(glm::dot(prec.axis_icrf, away) - std::cos(inclination)) < 1e-12,
          "visual inclination converted to WT02 signed orbit normal");
    check(prec.misalignment_deg == 21.1 && body.pulsar.magnetic_inclination_deg == 157.5 &&
              prec.epoch_tdb == tdb_from_julian_year(1981.0), "WT02 Table 1 solution 1 parameters");
    const double cycle = kTwoPi / prec.rate_rad_s;
    for (double turns : {0.0, 0.125, 0.25, 0.5, 0.75, 1.0, -0.25, 100.0}) {
        const double t = prec.epoch_tdb + turns * cycle;
        scene.update(t);
        const glm::dvec3 spin = pulsar_spin_axis(scene, body, t);
        // WT02 Eq. 5, with the scene's updated timing inclination.
        const double cos_zeta = -std::cos(kTwoPi * turns) * std::sin(tilt) * std::sin(inclination) -
                                 std::cos(tilt) * std::cos(inclination);
        check(std::abs(glm::dot(spin, -away) - cos_zeta) < 1e-11,
              "B1913 line-of-sight colatitude follows WT02 Eq. 5");
        check(std::abs(glm::length(spin) - 1.0) < 1e-12 &&
                  std::abs(glm::dot(spin, scene.orbit_normal(companion)) - std::cos(tilt)) < 1e-12,
              "B1913 spin stays on the 21.1 degree cone about orbital angular momentum");
    }
    const glm::dvec3 spin0 = pulsar_spin_axis(scene, body, prec.epoch_tdb);
    const glm::dvec3 spin1 = pulsar_spin_axis(scene, body, prec.epoch_tdb + 0.01 * year_s);
    check(glm::dot(spin1 - spin0, glm::cross(prec.axis_icrf, spin0)) > 0.0,
          "B1913 precesses in the angular-momentum cross spin direction");
    // The unknown sky node rotates both basis vectors around the line of sight.
    const auto saved = body.pulsar.precession;
    auto rotate_node = [&](const glm::dvec3& v) {
        return std::cos(0.7) * v + std::sin(0.7) * glm::cross(away, v) +
               (1.0 - std::cos(0.7)) * glm::dot(away, v) * away;
    };
    body.pulsar.precession.axis_icrf = rotate_node(saved.axis_icrf);
    body.pulsar.precession.x_icrf = rotate_node(saved.x_icrf);
    check(std::abs(glm::dot(pulsar_spin_axis(scene, body, prec.epoch_tdb + 0.01 * year_s), away) -
                       glm::dot(spin1, away)) < 1e-12, "B1913 zeta is independent of the sky node");
    body.pulsar.precession = saved;

    auto max_flash = [&](double year) {
        const double t = tdb_from_julian_year(year);
        scene.update(t);
        double peak = 0.0;
        for (int j = 0; j < 2048; ++j) {
            glm::dvec3 beams[2];
            check(pulsar_beams(scene, body, t, body.pulsar.display_period_s * j / 2048.0, beams),
                  "B1913 renders pulsar beams");
            const glm::dvec3 spin = pulsar_spin_axis(scene, body, t);
            check(glm::length(beams[0] + beams[1]) < 1e-12 &&
                      std::abs(glm::dot(spin, beams[0]) - std::cos(157.5 * kDegToRad)) < 1e-12,
                  "B1913 antipodal beams preserve magnetic inclination");
            peak = std::max(peak, pulsar_beam_flash(body.pulsar, beams[0], -away) +
                                  pulsar_beam_flash(body.pulsar, beams[1], -away));
        }
        return peak;
    };
    const double peak1974 = max_flash(1974.5);
    const double peak2019 = max_flash(2019.0);
    const double peak2026 = max_flash(2026.0);
    check(peak1974 > 0.5 && peak2019 > 0.0 && peak2019 < peak1974 && peak2026 == 0.0 &&
              max_flash(2035.0) == 0.0 && max_flash(1974.5) == peak1974,
          "B1913 finite beam fades, misses Earth after 2025, and restores on rewind");
    // Solve the edge crossing from WT02's analytic colatitude, independently of beam sampling.
    const double zeta_edge = (157.5 - body.pulsar.beam_half_angle_deg) * kDegToRad;
    const double phase_edge = std::acos((-std::cos(zeta_edge) - std::cos(tilt) * std::cos(inclination)) /
                                        (std::sin(tilt) * std::sin(inclination)));
    const double exit_year = 1981.0 + phase_edge / prec.rate_rad_s / year_s;
    check(exit_year > 2024.5 && exit_year < 2025.5, "B1913 illustrative beam exits near 2025", exit_year);
    // Changing display frames must rotate the same beams, preserving pulse phase.
    const double render_time = tdb_from_julian_year(1998.0);
    scene.update(render_time);
    glm::dvec3 inertial_beams[2], rotating_beams[2];
    pulsar_beams(scene, body, render_time, 0.37, inertial_beams);
    scene.frames.push_back({"test rotating", DisplayFrame::Type::Rotating, 0, pulsar, companion});
    scene.set_active_frame(static_cast<int>(scene.frames.size()) - 1);
    scene.update(render_time);
    pulsar_beams(scene, body, render_time, 0.37, rotating_beams);
    check(glm::length(scene.current_transform().axes * rotating_beams[0] - inertial_beams[0]) < 1e-12,
          "pulsar beam phase is independent of the display frame");
    std::printf("info: B1913 spin precession %.6f deg/yr, illustrative beam exit %.3f\n", rate, exit_year);

}

// The corotating dipole inside the light cylinder (Goldreich & Julian 1969).
void test_pulsar_field_lines()
{
    // Aligned: the last closed line crosses the spin equator at the cylinder, so L = c / Omega,
    // and the polar cap's half-angle is asin(sqrt(R / L)).
    check(std::abs(light_cylinder_km(0.0590301) - 2816.6) < 0.1, "B1913's light cylinder from 59 ms (km)",
          light_cylinder_km(0.0590301));
    check(std::abs(pulsar_last_closed_l(0.0, 1.3, 1000.0) - 1000.0) < 1e-6, "aligned: the last closed line has L = R_LC");
    // Orthogonal, in the plane of both axes: the line crosses the spin axis at its magnetic
    // equator, its distance from the axis is L sin^2(theta) |cos(theta)|, at most 2 / (3 sqrt 3) L.
    // Across that plane (phi = 90 deg) the whole line lies in the spin equator.
    check(std::abs(pulsar_last_closed_l(90.0, 0.0, 1000.0) / 1000.0 - 1.5 * std::sqrt(3.0)) < 1e-6,
          "orthogonal, in the plane of both axes: L = (3 sqrt 3 / 2) R_LC", pulsar_last_closed_l(90.0, 0.0, 1000.0));
    check(std::abs(pulsar_last_closed_l(90.0, 0.5 * kPi, 1000.0) - 1000.0) < 1e-6,
          "orthogonal, across it: L = R_LC");

    struct Case { const char* file; const char* name; double period, cap_low, cap_high; };
    // Cap diameters as in the events' captions (12 km star): B1913 about 1.5 km, B1620 some 3.5 km.
    for (const Case& c : {Case{"psr_b1913.toml", "PSR B1913+16", 0.0590301, 1.3, 1.8},
                          Case{"psr_j1141.toml", "PSR J1141-6545", 0.393902, 0.3, 1.0},
                          Case{"psr_b1620.toml", "PSR B1620-26 A", 0.0110757509142, 3.0, 4.0}}) {
        Scene scene = load_scene_or_die(c.file);
        const int index = scene.find(c.name);
        check(index >= 0, c.name);
        if (index < 0) continue;
        const Body& body = scene.bodies[static_cast<size_t>(index)];
        const Body::Pulsar& p = body.pulsar;
        const double rlc = p.light_cylinder_km, radius = body.equatorial_radius_km;
        check(p.spin_period_s == c.period && std::abs(rlc - light_cylinder_km(c.period)) < 1e-9 && !p.field_lines.empty(),
              "the pulsar's lines come from its true period");
        const double alpha = p.magnetic_inclination_deg * kDegToRad;
        const glm::dvec3 m(std::sin(alpha), 0.0, std::cos(alpha));
        double worst_surface = 0.0, worst_outside = 0.0, worst_end = 0.0, worst_tangent = 0.0, widest_closed = 0.0;
        int closed = 0, open = 0;
        for (const Body::Pulsar::FieldLine& line : p.field_lines) {
            worst_surface = std::max(worst_surface, std::abs(glm::length(glm::dvec3(line.points.front())) - radius));
            for (size_t k = 0; k < line.points.size(); ++k) {
                const glm::dvec3 q(line.points[k]);
                const double rho = std::hypot(q.x, q.y);
                worst_outside = std::max(worst_outside, rho / rlc - 1.0);
                if (!line.open) widest_closed = std::max(widest_closed, rho / rlc);
                if (k == 0) continue;
                // Along the dipole field 3 (m.n) n - m at the segment's middle.
                const glm::dvec3 a(line.points[k - 1]);
                const glm::dvec3 mid = 0.5 * (a + q), n = glm::normalize(mid);
                const glm::dvec3 b = 3.0 * glm::dot(m, n) * n - m;
                worst_tangent = std::max(worst_tangent, std::asin(std::min(glm::length(glm::cross(glm::normalize(q - a), glm::normalize(b))), 1.0)));
            }
            if (line.open) {
                ++open;
                const glm::dvec3 end(line.points.back());
                // The magnetic axis lines stop at the cylinder or 2 R_LC out.
                if (std::abs(glm::dot(glm::normalize(end), m)) < 0.9999) {
                    worst_end = std::max(worst_end, std::abs(std::hypot(end.x, end.y) / rlc - 1.0));
                }
            } else {
                ++closed;
                worst_surface = std::max(worst_surface, std::abs(glm::length(glm::dvec3(line.points.back())) - radius));
            }
        }
        check(worst_surface < 1e-6 * radius, "every line starts (and a closed one ends) on the star (km)", worst_surface);
        check(worst_outside < 1e-9 && widest_closed > 0.995, "nothing beyond the cylinder; the last closed lines touch it",
              widest_closed);
        check(worst_end < 1e-6, "open lines end on the cylinder", worst_end);
        check(worst_tangent < 0.03, "the lines follow the dipole (rad)", worst_tangent);
        // The polar cap: lines rooted within asin(sqrt(R / L_last)) of the magnetic pole open.
        double cap = 0.0;
        for (int k = 0; k < 64; ++k) {
            cap += std::asin(std::sqrt(radius / pulsar_last_closed_l(p.magnetic_inclination_deg, kTwoPi * k / 64.0, rlc))) / 64.0;
        }
        const double cap_km = 2.0 * radius * cap;
        check(cap_km > c.cap_low && cap_km < c.cap_high, "the polar cap is as wide as the caption says (km)", cap_km);
        // The lines turn with the beams: the corotating frame's tilted axis is beam 0.
        scene.update(tdb_from_julian_year(2026.0));
        glm::dmat3 axes;
        glm::dvec3 beams[2];
        check(pulsar_corotating_axes(scene, body, tdb_from_julian_year(2026.0), 0.61, axes) &&
                  pulsar_beams(scene, body, tdb_from_julian_year(2026.0), 0.61, beams) &&
                  glm::length(axes * m - beams[0]) < 1e-12 && std::abs(glm::determinant(axes) - 1.0) < 1e-12,
              "the field turns with the beams");
        std::printf("info: %s: light cylinder %.0f km, %d closed and %d open lines, polar cap %.2f km across\n", c.name,
                    rlc, closed, open, cap_km);
    }
}

void test_psr_j1141_scene()
{
    Scene scene = load_scene_or_die("psr_j1141.toml");
    check_model_sky(scene, "PSR J1141-6545 model sky");
    const int pulsar = scene.find("PSR J1141-6545");
    const int wd = scene.find("J1141-6545 white dwarf");
    const int ghost_p = scene.find("PSR J1141-6545 (Newton)");
    check(pulsar > 0 && wd > 0 && ghost_p > 0, "PSR J1141-6545 bodies");
    if (pulsar <= 0 || wd <= 0 || ghost_p <= 0) {
        return;
    }
    const PostKeplerianMotion* orbit = post_keplerian_of(scene, pulsar);
    check(orbit != nullptr, "PSR J1141-6545 post-Keplerian orbit");
    if (!orbit) {
        return;
    }
    const double gm = 3.03867692589e11;
    const glm::dvec3 away = unit_toward(175.27919583, -65.75531667);
    const double x = 1.858915;
    const double e = 0.171876;
    const double t0 = tdb_from_jd_tdb(2454000.4960283);
    const double depth = pulsar_depth_s(scene, pulsar, away, t0);
    check(std::abs(depth - x * (1.0 - e) * std::sin(80.6911 * kDegToRad)) < 2e-5,
          "PSR J1141-6545 Roemer delay at T0 (s)", depth);

    // [BBV08]'s measured advance is GR's for the [VK20] masses...
    const double omdot_gr = gr_omega_dot_deg_per_year(2.28967 - 1.02, 1.02, 0.19765096149, e);
    check(std::abs(omdot_gr / 5.3096 - 1.0) < 3e-4, "PSR J1141-6545 omega_dot is GR's for its masses", omdot_gr);
    const double pbdot_gr = gr_pb_dot(2.28967 - 1.02, 1.02, 0.19765096149, e);
    check(std::abs(pbdot_gr - -0.403e-12) < 0.025e-12, "PSR J1141-6545 Pb_dot within 1 sigma of GR", pbdot_gr);
    std::printf("info: PSR J1141-6545 GR: omega_dot %.5f deg/yr, Pb_dot %.4e\n", omdot_gr, pbdot_gr);

    // ...and the model, run back seven years from [VK20]'s T0, finds [BBV08]'s 1999
    // periastron (an independent timing solution) and its omega there. [VK20] used
    // tempo2, whose default time scale is TCB (TCB - TDB ~ 11-15 s then), [BBV08]
    // TEMPO (TDB): seconds of difference are expected.
    const double t_bbv = tdb_from_jd_tdb(2451370.3545515);
    const PostKeplerianPhase then = post_keplerian_phase(orbit->orbit(), gm, t_bbv);
    const double late_s = wrap_pi(then.mean_anomaly) / kTwoPi * 0.19765096149 * kSecondsPerDay;
    check(std::abs(late_s) < 30.0, "PSR J1141-6545 periastron of 1999 [BBV08] (s off)", late_s);
    const double w_pulsar = wrap_two_pi(then.arg_peri_rad - kPi) * kRadToDeg;
    check(std::abs(w_pulsar - 42.4561) < 0.003, "PSR J1141-6545 omega in 1999 [BBV08] (deg)", w_pulsar);
    std::printf("info: PSR J1141-6545 in 1999: periastron %.2f s off, omega %.4f deg\n", late_s, w_pulsar);

    double t_ref = 0.0;
    check(parse_utc("1999-07-10 20:30", &t_ref), "parse the J1141 reference date");
    const double split =
        glm::length(scene.icrf_state_at(pulsar, t_ref).position - scene.icrf_state_at(ghost_p, t_ref).position);
    check(split < 1e-3, "PSR J1141-6545 ghost starts on the real orbit (km)", split);
    check(std::abs(scene.bodies[static_cast<size_t>(wd)].body_gm_km3_s2 / kSunGmKm3S2 - 1.02) < 1e-6,
          "J1141-6545 white dwarf mass");
}

void test_schwarzschild_star_light()
{
    const SchwarzschildStarLight flat{0.0, 12.0, 144.0};
    check(std::abs(flat.limb_angle() - std::asin(12.0 / 144.0)) < 1e-12, "flat star limb");
    check(std::abs(flat.surface_angle(0.0) - std::acos(12.0 / 144.0)) < 1e-7, "flat star limb surface point");
    for (double angle : {flat.limb_angle(), 0.3, 0.5 * kPi, 2.0, kPi}) {
        check(std::abs(flat.sky_deflection(angle)) < 1e-7, "flat exterior ray stays straight", angle);
    }

    const SchwarzschildStarLight light{1.438 * kSunGmKm3S2 / (kSpeedOfLightKmS * kSpeedOfLightKmS), 12.0, 144.0};
    // Independent high-resolution quadrature of the radial integral,
    // regularized with t = sqrt(1 - r / r'), for the B1913+16 close view.
    check(std::abs(light.limb_angle() - 0.102311809) < 1e-8, "B1913+16 apparent limb");
    check(std::abs(light.surface_angle(0.0) - 2.025222767) < 1e-7, "B1913+16 far side visible at limb");
    check(std::abs(light.surface_angle(0.5) - 1.251830223) < 1e-7, "B1913+16 surface ray inside limb");
    check(std::abs(light.sky_deflection(light.limb_angle()) - 1.115020751) < 1e-7,
          "B1913+16 sky ray grazing limb");
    check(std::abs(light.sky_deflection(0.5 * kPi) - 0.030365823) < 1e-7, "B1913+16 tangential observer ray");
    check(std::abs(light.sky_deflection(0.5 * kPi - 1e-5) - light.sky_deflection(0.5 * kPi + 1e-5)) < 1e-6,
          "B1913+16 sky map is continuous across the observer tangent");
}

void test_kerr_null_geodesics()
{
    // Analytic metric derivatives (also used by the shader) vs finite differences.
    const double a = 0.9;
    double worst = 0.0;
    for (const glm::dvec3 x : {glm::dvec3(3.1, -2.2, 1.7), glm::dvec3(-0.4, 5.0, -3.3), glm::dvec3(12.0, 1.0, 0.2)}) {
        const KerrSchild k = kerr_schild(x, a);
        for (int j = 0; j < 3; ++j) {
            glm::dvec3 d(0.0);
            d[j] = 1e-6;
            const KerrSchild kp = kerr_schild(x + d, a);
            const KerrSchild km = kerr_schild(x - d, a);
            worst = std::max(worst, std::abs((kp.f - km.f) / 2e-6 - k.df[j]));
            for (int i = 0; i < 3; ++i) {
                worst = std::max(worst, std::abs((kp.l[i] - km.l[i]) / 2e-6 - k.dl[i][j]));
            }
        }
    }
    check(worst < 1e-6, "Kerr-Schild analytic derivatives", worst);

    // The ray stays null (H = 0) through a close passage.
    Photon ph = make_photon(glm::dvec3(-50.0, 6.0, 3.0), glm::dvec3(1.0, 0.0, 0.0), a);
    double worst_h = 0.0;
    for (int i = 0; i < 3000; ++i) {
        rk4_step(ph, a, ray_step(ph.x, a));
        worst_h = std::max(worst_h, std::abs(hamiltonian(ph, a)));
        if (glm::length(ph.x) > 60.0) {
            break;
        }
    }
    check(worst_h < 1e-6, "null geodesic stays null", worst_h);

    // Schwarzschild: capture iff b < 3 sqrt(3) = 5.196; weak-field deflection 4/b + 15 pi / (4 b^2).
    auto shoot = [](double spin, double b) {
        return trace_ray(glm::dvec3(-2000.0, b, 0.0), glm::dvec3(1.0, 0.0, 0.0), spin, 2100.0, 0.0, 0.0, 20000);
    };
    check(shoot(0.0, 5.15).outcome == RayOutcome::Captured, "b = 5.15 captured");
    check(shoot(0.0, 5.25).outcome == RayOutcome::Escaped, "b = 5.25 escapes");
    const RayResult far = shoot(0.0, 100.0);
    const double deflection = std::acos(std::clamp(far.direction.x, -1.0, 1.0));
    const double expected = 4.0 / 100.0 + 15.0 * kPi / (4.0 * 100.0 * 100.0);
    check(std::abs(deflection - expected) < 2e-3 * expected * 10.0, "weak-field deflection", deflection / expected);

    // Spin: the traced ray is time-reversed, so the real photon moves along -x and
    // has L_z = y. Prograde rays (y > 0) get closer before capture than retrograde ones.
    check(shoot(0.9, 4.0).outcome == RayOutcome::Escaped && shoot(0.9, -4.0).outcome == RayOutcome::Captured,
          "Kerr prograde/retrograde capture asymmetry");

    check(std::abs(isco_radius(0.0) - 6.0) < 1e-12, "ISCO a = 0", isco_radius(0.0));
    check(std::abs(isco_radius(0.9) - 2.3209) < 1e-3, "ISCO a = 0.9", isco_radius(0.9));
    check(std::abs(horizon_radius(0.0) - 2.0) < 1e-12, "horizon a = 0");

    // Doppler: a disk viewed edge-on from -x; the side moving toward the observer is blueshifted.
    const Photon toward = make_photon(glm::dvec3(0.0, 10.0, 0.0), glm::dvec3(1.0, 0.0, 0.0), 0.0);
    const Photon away = make_photon(glm::dvec3(0.0, -10.0, 0.0), glm::dvec3(1.0, 0.0, 0.0), 0.0);
    const double g_toward = disk_redshift(toward, 0.0);
    const double g_away = disk_redshift(away, 0.0);
    check(g_toward > 1.0 && g_away < 1.0, "disk Doppler shift sides", g_toward);
    // Face-on emission: only gravitational + transverse redshift, g = sqrt(1 - 3/r).
    const Photon face_on = make_photon(glm::dvec3(10.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 1.0), 0.0);
    check(std::abs(disk_redshift(face_on, 0.0) - std::sqrt(1.0 - 3.0 / 10.0)) < 1e-9, "face-on redshift sqrt(1-3/r)",
          disk_redshift(face_on, 0.0));
}

// Local angle from the camera's forward axis at which rays stop being captured
// (bisection), for rays in the plane of forward and `side` (a unit 2-vector in
// the (right, up) screen components).
double shadow_edge_angle(const glm::dvec3& cam, double a, const ObserverFrame& frame, const glm::dvec2& side)
{
    double lo = 0.0;      // captured
    double hi = kPi / 2.0; // escapes
    for (int i = 0; i < 40; ++i) {
        const double t = 0.5 * (lo + hi);
        const glm::dvec3 n(std::sin(t) * side.x, std::sin(t) * side.y, std::cos(t));
        const RayResult r = trace_ray(frame_photon(cam, frame, n), a, 1000.0, 0.0, 0.0, 20000);
        (r.outcome == RayOutcome::Escaped ? hi : lo) = t;
    }
    return 0.5 * (lo + hi);
}

void test_observer_frame()
{
    // The frame is orthonormal in the metric, and its photons are null with observed energy 1.
    const double a = 0.9;
    const glm::dvec3 x(7.0, -3.0, 2.0);
    const glm::dvec3 fwd = glm::normalize(-x);
    const glm::dvec3 up = glm::normalize(glm::cross(glm::cross(fwd, glm::dvec3(0.0, 0.0, 1.0)), fwd));
    const ObserverFrame frame = static_observer_frame(x, a, fwd, up, glm::cross(fwd, up));
    const glm::dvec4 basis[4] = {frame.u, frame.right, frame.up, frame.forward};
    double worst = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const double expected = i == j ? (i == 0 ? -1.0 : 1.0) : 0.0;
            worst = std::max(worst, std::abs(covector_dot(x, a, basis[i], basis[j]) - expected));
        }
    }
    check(worst < 1e-12, "static observer frame orthonormal", worst);
    const Photon ph = frame_photon(x, frame, glm::dvec3(0.3, -0.2, 0.9));
    check(std::abs(hamiltonian(ph, a)) < 1e-12, "frame photon is null", hamiltonian(ph, a));
    check(std::abs(ph.pt * frame.ut - 1.0) < 1e-12, "frame photon observed energy 1", ph.pt * frame.ut);

    // Far away the frame is the coordinate one.
    const glm::dvec3 far_x(0.0, 3e6, 1e6);
    const glm::dvec3 far_fwd = glm::normalize(-far_x);
    const glm::dvec3 far_up = glm::normalize(glm::cross(glm::cross(far_fwd, glm::dvec3(1.0, 0.0, 0.0)), far_fwd));
    const ObserverFrame far = static_observer_frame(far_x, a, far_fwd, far_up, glm::cross(far_fwd, far_up));
    const glm::dvec3 n(0.1, 0.05, 1.0);
    const glm::dvec3 coordinate = glm::normalize(n.x * glm::cross(far_fwd, far_up) + n.y * far_up + n.z * far_fwd);
    const double far_diff = glm::length(glm::normalize(frame_photon(far_x, far, n).p) - coordinate);
    check(far_diff < 1e-6, "observer frame far away = coordinate directions", far_diff);

    // Schwarzschild shadow seen by a static observer at r: sin(alpha) = 3 sqrt(3) sqrt(1 - 2/r) / r
    // (Synge 1966); the coordinate directions used before gave 18.4 deg instead of 23.3 at r = 12.
    for (const double r : {4.0, 12.0, 30.0}) {
        const glm::dvec3 cam = r * glm::normalize(glm::dvec3(1.0, 2.0, 0.5));
        const glm::dvec3 f = glm::normalize(-cam);
        const glm::dvec3 u = glm::normalize(glm::cross(glm::cross(f, glm::dvec3(0.0, 0.0, 1.0)), f));
        const ObserverFrame fr = static_observer_frame(cam, 0.0, f, u, glm::cross(f, u));
        const double expected = std::asin(3.0 * std::sqrt(3.0) * std::sqrt(1.0 - 2.0 / r) / r);
        const double edge = shadow_edge_angle(cam, 0.0, fr, glm::normalize(glm::dvec2(1.0, 0.6)));
        check(std::abs(edge - expected) < 1e-4, "Schwarzschild shadow radius for a static observer (rad)",
              edge - expected);
    }

    // Kerr, observer on the spin axis at r = 12: photons asymptotic to the spherical photon
    // orbit with L_z = 0, where r^3 - 3 r^2 + a^2 r + a^2 = 0 and the Carter constant is
    // Q = -r^3 (r^3 - 6 r^2 + 9 r - 4 a^2) / (a^2 (r - 1)^2) (Teo 2003, GRG 35, 1909: Phi(r) = 0
    // and Q(r) of the spherical photon orbits). For L_z = 0, R(r) = (r^2 + a^2)^2 - Delta (Q + a^2); on the axis the
    // static observer measures cos(alpha) = sqrt(R) / (r^2 + a^2), so
    // sin(alpha) = sqrt(Delta (Q + a^2)) / (r^2 + a^2). The shadow is a circle.
    double rp = 3.0;
    for (int i = 0; i < 50; ++i) {
        rp -= (rp * rp * rp - 3.0 * rp * rp + a * a * rp + a * a) / (3.0 * rp * rp - 6.0 * rp + a * a);
    }
    const double q = -rp * rp * rp * (rp * rp * rp - 6.0 * rp * rp + 9.0 * rp - 4.0 * a * a) /
                     (a * a * (rp - 1.0) * (rp - 1.0));
    const double r_obs = 12.0;
    const double delta = r_obs * r_obs - 2.0 * r_obs + a * a;
    const double expected_pole = std::asin(std::sqrt(delta * (q + a * a)) / (r_obs * r_obs + a * a));
    const glm::dvec3 pole_cam(0.0, 0.0, r_obs);
    const ObserverFrame pole = static_observer_frame(pole_cam, a, glm::dvec3(0.0, 0.0, -1.0),
                                                     glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(0.0, -1.0, 0.0));
    for (const glm::dvec2 side : {glm::dvec2(1.0, 0.0), glm::dvec2(-0.6, 0.8)}) {
        const double edge = shadow_edge_angle(pole_cam, a, pole, side);
        check(std::abs(edge - expected_pole) < 1e-4, "Kerr a = 0.9 shadow radius on the axis (rad)",
              edge - expected_pole);
    }

    // Weak-field legs vs. quadrature of dphi = dw / sqrt(1 - w^2 + 2 w^3 / b) (Schwarzschild):
    // from infinity to r = 150, b = 6 sweeps 0.040010568 rad and b = 100 0.728861494 rad
    // (Gauss-Legendre); a whole pass at b = 300 bends by 0.0134658 rad.
    check(std::abs(weak_leg_angle(6.0 / 150.0, 6.0) - 0.040010568) < 1e-8, "weak leg b = 6",
          weak_leg_angle(6.0 / 150.0, 6.0));
    check(std::abs(weak_leg_angle(100.0 / 150.0, 100.0) - 0.728861494) < 1e-7, "weak leg b = 100",
          weak_leg_angle(100.0 / 150.0, 100.0));
    const double pass = weak_pass_deflection(1e9, 300.0, true);
    check(std::abs(pass - 0.0134658) < 3e-6, "weak pass deflection b = 300", pass);

    // Camera rays as the shader traces them (integration inside 150 M, exact entry and
    // analytic legs outside) vs. integrating all the way from the camera to 1e7 M.
    for (const double spin : {0.0, 0.6}) {
        for (const double dist : {40.0, 400.0, 3000.0}) {
            const glm::dvec3 cam = dist * glm::normalize(glm::dvec3(0.3, -1.0, 0.25));
            const glm::dvec3 f = glm::normalize(-cam);
            const glm::dvec3 u = glm::normalize(glm::cross(glm::cross(f, glm::dvec3(0.0, 0.0, 1.0)), f));
            const ObserverFrame fr = static_observer_frame(cam, spin, f, u, glm::cross(f, u));
            double worst_dir = 0.0;
            int mismatched = 0;
            int compared = 0;
            const double max_angle = std::min(1.4, 400.0 / dist); // out to b ~ 400 M
            for (int i = 0; i <= 40; ++i) {
                for (int k = 0; k < 3; ++k) {
                    const double t = max_angle * i / 40.0;
                    const double az = 0.4 + 2.1 * k;
                    const glm::dvec3 dir(std::sin(t) * std::cos(az), std::sin(t) * std::sin(az), std::cos(t));
                    const RayResult fast = trace_camera_ray(cam, fr, dir, spin, 150.0, 0.0, 0.0, 20000);
                    const RayResult full = trace_ray(frame_photon(cam, fr, dir), spin, 1e7, 0.0, 0.0, 20000);
                    if (fast.outcome != full.outcome) {
                        ++mismatched;
                    } else if (fast.outcome == RayOutcome::Escaped) {
                        ++compared;
                        worst_dir = std::max(worst_dir, std::atan2(glm::length(glm::cross(fast.direction, full.direction)),
                                                                   glm::dot(fast.direction, full.direction)));
                    }
                }
            }
            std::printf("info: camera ray a = %.1f at %.0f M: %d escaped, worst %.2e rad, %d outcomes differ\n", spin,
                        dist, compared, worst_dir, mismatched);
            check(mismatched == 0 && worst_dir < (spin == 0.0 ? 1.5e-4 : 4e-4),
                  "camera rays: analytic far field matches full integration (rad)", worst_dir);
        }
    }

    // Disk seen by a nearby static observer: blueshifted by u^t relative to infinity.
    const Photon face_on = make_photon(glm::dvec3(10.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 1.0), 0.0);
    check(std::abs(disk_redshift(face_on, 0.0, 1.25) / disk_redshift(face_on, 0.0) - 1.25) < 1e-12,
          "disk redshift scales with the observer's u^t");
}

void test_sgr_a_scene()
{
    Scene scene = load_scene_or_die("sgr_a.toml");
    const int s2 = scene.find("S2");
    check(scene.bodies[0].kind == BodyKind::BlackHole, "Sgr A* is a black hole");
    // Horizon of a 4.261e6 M_sun hole with a = 0.94: (1 + sqrt(1 - a^2)) GM / c^2.
    check(std::abs(scene.bodies[0].equatorial_radius_km / (6.2921e6 * (1.0 + std::sqrt(1.0 - 0.94 * 0.94))) - 1.0) < 1e-3,
          "Sgr A* horizon", scene.bodies[0].equatorial_radius_km);

    // Pericentre of the GR-integrated orbit (started at the 2010 apocentre).
    double best_t = 0.0;
    double best_r = 1e300;
    for (double y = 2017.5; y < 2019.5; y += 1e-4) {
        const double r = glm::length(scene.icrf_state_at(s2, tdb_from_julian_year(y)).position);
        if (r < best_r) {
            best_r = r;
            best_t = y;
        }
    }
    check(std::abs(best_t - 2018.379) < 0.01, "S2 pericentre epoch", best_t);
    check(std::abs(best_r / kAuKm - 119.0) < 3.0, "S2 pericentre distance (au)", best_r / kAuKm);
    std::printf("info: S2 pericentre %.4f at %.1f au\n", best_t, best_r / kAuKm);

    const int s301 = scene.find("S301");
    check(s301 > 0, "S301 body");
    if (s301 > 0) {
        const std::vector<double> peri = scene.periapsis_times(s301, tdb_from_julian_year(2023.0),
                                                               tdb_from_julian_year(2032.0));
        check(peri.size() == 2, "S301 has two successive pericentres from 2023 to 2032");
        if (peri.size() == 2) {
            // [GR26] Extended Data Table 2 and Results (first orientation).
            const State first = scene.icrf_state_at(s301, peri[0]);
            const State second = scene.icrf_state_at(s301, peri[1]);
            const double epoch = 2000.0 + peri[0] / (365.25 * kSecondsPerDay);
            const double period = (peri[1] - peri[0]) / (365.25 * kSecondsPerDay);
            const double radius = glm::length(first.position) / kAuKm;
            const double speed = glm::length(first.velocity);
            const double advance = std::acos(glm::dot(glm::normalize(first.position),
                                                      glm::normalize(second.position))) / kDegToRad;
            check(std::abs(epoch - 2023.126) < 0.010, "S301 pericentre epoch", epoch);
            check(std::abs(period - 8.68) < 0.11, "S301 radial period (yr)", period);
            check(std::abs(radius - 11.5) < 0.8, "S301 pericentre distance (au)", radius);
            check(std::abs(speed - 25600.0) < 1000.0, "S301 pericentre speed (km/s)", speed);
            check(std::abs(advance - 2.0) < 0.1, "S301 apsidal advance (deg)", advance);
            std::printf("info: S301 pericentre %.6f at %.2f au, %.0f km/s, period %.3f yr, advance %.2f deg\n",
                        epoch, radius, speed, period, advance);
        }

        // [GR26] places the April 2023 detection about 15 mas northwest of Sgr A*.
        const glm::dvec3 position = scene.icrf_state_at(s301, tdb_from_julian_year(2023.2604)).position;
        const glm::dvec3 east = unit_toward(266.41683333 + 90.0, 0.0);
        const glm::dvec3 north = glm::cross(glm::normalize(scene.sky.viewer_pc), east);
        const double km_per_mas = glm::length(scene.sky.viewer_pc) * kAuKm / 1000.0;
        const double ra_mas = glm::dot(position, east) / km_per_mas;
        const double dec_mas = glm::dot(position, north) / km_per_mas;
        check(ra_mas < 0.0 && dec_mas > 0.0 && std::abs(std::hypot(ra_mas, dec_mas) - 15.0) < 3.0,
              "S301 sky position at discovery (mas)", std::hypot(ra_mas, dec_mas));
    }

    // The Death Star's rolling window follows the clock, including big jumps.
    const int hyp = scene.find("Death Star");
    check(hyp >= 0 && scene.bodies[static_cast<size_t>(hyp)].style == SurfaceStyle::DeathStar &&
              scene.bodies[static_cast<size_t>(hyp)].mark_periapsides,
          "Death Star body");
    // The sky is the model seen from Sgr A* itself.
    check_model_sky(scene, "Sgr A* model sky");
    const glm::dvec3 hole_pc = unit_toward(266.41683333, -29.00781611) * 8246.7;
    check(glm::length(scene.sky.viewer_pc - hole_pc) < 1e-6, "Sgr A* sky seen from Sgr A*",
          glm::length(scene.sky.viewer_pc - hole_pc));

    // The other S-stars: Keplerian orbits with [GI17]'s periods in this scene's
    // potential have [GI17]'s angular semi-major axes (they used 4.28e6 M_sun at
    // 8.32 kpc; here 4.261e6 at 8.2467: ~0.7% larger) - Table 3.
    struct SStar {
        const char* name;
        double a_arcsec;
    };
    for (const SStar& s : {SStar{"S1", 0.595}, SStar{"S4", 0.3570}, SStar{"S6", 0.6574}, SStar{"S8", 0.4047},
                           SStar{"S12", 0.2987}, SStar{"S13", 0.2641}, SStar{"S14", 0.2863}, SStar{"S31", 0.449}}) {
        const int index = scene.find(s.name);
        check(index > 0 && !scene.bodies[static_cast<size_t>(index)].label, s.name);
        if (index <= 0) {
            continue;
        }
        const State st = scene.bodies[static_cast<size_t>(index)].motion->eval(0.0);
        const double r = glm::length(st.position);
        const double a = 1.0 / (2.0 / r - glm::dot(st.velocity, st.velocity) / scene.bodies[0].gm_km3_s2);
        const double a_arcsec = a / kAuKm / 8246.7;
        check(std::abs(a_arcsec / s.a_arcsec - 1.0) < 0.02, "S-star semi-major axis vs GI17", a_arcsec);
    }
    const double m = scene.bodies[0].gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    double t = tdb_from_julian_year(2026.5);
    double r_min = 1e300;
    for (int i = 0; i < 2000; ++i) {
        t += 3600.0;
        r_min = std::min(r_min, glm::length(scene.icrf_state_at(hyp, t).position));
    }
    check(std::abs(r_min / m - 40.0) < 2.0, "Death Star pericentre ~40 M", r_min / m);
    const glm::dvec3 far = scene.icrf_state_at(hyp, tdb_from_julian_year(1900.0)).position;
    check(glm::length(far) > 30.0 * m && glm::length(far) < 400.0 * m, "rolling orbit after a jump", glm::length(far) / m);
}

// The Death Star's orbit (as marked in the app) around the spinning hole: its
// periapsis advances in the orbit plane by the exact Schwarzschild amount for its
// turning points (~18 deg per orbit, vs 16.6 at 1PN) plus the spin's Lense-Thirring
// and quadrupole terms, and its node turns about the spin axis; the event "at
// periapsis" really is at one.
void test_death_star_precession()
{
    Scene scene = load_scene_or_die("sgr_a.toml");
    const int ds = scene.find("Death Star");
    const double m = scene.bodies[0].gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);

    const double t0 = tdb_from_julian_year(2026.5);
    scene.update(t0);
    const double span = 10.0 * kSecondsPerDay; // inside the 30-day rolling window
    const std::vector<double> peri = scene.periapsis_times(ds, t0, t0 + span);
    check(peri.size() == 4 || peri.size() == 5, "periapsides in 10 days", static_cast<double>(peri.size()));
    if (peri.size() < 3) {
        return;
    }

    double r_p = 0.0;
    double r_a = 0.0;
    for (double tp : peri) {
        r_p += glm::length(scene.icrf_state_at(ds, tp).position) / m / static_cast<double>(peri.size());
    }
    for (double t = peri[0]; t < peri[1]; t += 60.0) {
        r_a = std::max(r_a, glm::length(scene.icrf_state_at(ds, t).position) / m);
    }
    // Orbit angles about the spin axis z: the node n on the hole's equator and the
    // periapsis's angle omega from it in the orbit plane.
    const glm::dvec3 z = scene.bodies[0].pole;
    const glm::dvec3 x_axis = glm::normalize(glm::cross(glm::dvec3(0.0, 0.0, 1.0), z));
    const glm::dvec3 y_axis = glm::cross(z, x_axis);
    struct Angles {
        double node = 0.0;
        double omega = 0.0;
        double inclination = 0.0;
    };
    auto angles_at = [&](double t) {
        const State s = scene.icrf_state_at(ds, t);
        const glm::dvec3 l = glm::normalize(glm::cross(s.position, s.velocity));
        const glm::dvec3 n = glm::normalize(glm::cross(z, l)); // ascending node
        Angles out;
        out.node = std::atan2(glm::dot(n, y_axis), glm::dot(n, x_axis));
        out.omega = std::atan2(glm::dot(s.position, glm::cross(l, n)), glm::dot(s.position, n));
        out.inclination = std::acos(std::clamp(glm::dot(l, z), -1.0, 1.0));
        return out;
    };
    double node_sum = 0.0;
    double omega_sum = 0.0;
    for (size_t k = 1; k < peri.size(); ++k) {
        const Angles before = angles_at(peri[k - 1]);
        const Angles after = angles_at(peri[k]);
        node_sum += wrap_pi(after.node - before.node);
        omega_sum += wrap_pi(after.omega - before.omega);
    }
    const double orbits = static_cast<double>(peri.size() - 1);
    const double node_rate = node_sum / orbits;
    const double omega_rate = omega_sum / orbits;

    // Predictions per orbit for a = 0.94: Lense-Thirring (first post-Newtonian order)
    // dOmega = 4 pi a / (A^3/2 (1 - e^2)^3/2), domega = -3 cos i dOmega; the quadrupole
    // Q = -a^2 M as an oblate J2 R^2 = a^2: dOmega = -3 pi a^2 cos i / p^2,
    // domega = 3 pi a^2 (2 - 5/2 sin^2 i) / p^2, p = A (1 - e^2).
    const double a_spin = 0.94;
    const double inc = angles_at(peri[0]).inclination;
    const double big_a = 0.5 * (r_p + r_a);
    const double ecc = (r_a - r_p) / (r_a + r_p);
    const double semi_latus = big_a * (1.0 - ecc * ecc);
    const double lt = 4.0 * kPi * a_spin / (std::pow(big_a, 1.5) * std::pow(1.0 - ecc * ecc, 1.5));
    const double quad = 3.0 * kPi * a_spin * a_spin / (semi_latus * semi_latus);
    const double node_expected = lt - quad * std::cos(inc);
    const double schwarzschild = schwarzschild_advance(r_p, r_a);
    const double omega_expected =
        schwarzschild - 3.0 * std::cos(inc) * lt + quad * (2.0 - 2.5 * std::sin(inc) * std::sin(inc));
    std::printf("info: Death Star r_p %.2f M, r_a %.1f M, i %.1f deg: node %.3f deg/orbit (expected %.3f), "
                "periapsis %.3f (expected %.3f, Schwarzschild alone %.3f)\n",
                r_p, r_a, inc / kDegToRad, node_rate / kDegToRad, node_expected / kDegToRad, omega_rate / kDegToRad,
                omega_expected / kDegToRad, schwarzschild / kDegToRad);
    check(node_rate > 0.0 && std::abs(node_rate / node_expected - 1.0) < 0.03,
          "Death Star node: Lense-Thirring + quadrupole precession", node_rate / node_expected);
    check(omega_rate > 0.0 && std::abs(omega_rate / omega_expected - 1.0) < 0.03,
          "Death Star periapsis: Schwarzschild + Lense-Thirring + quadrupole advance", omega_rate / omega_expected);

    // Event "Death Star at periapsis": a periapsis within the following 3 hours.
    for (const SceneEvent& e : scene.events) {
        if (e.focus == ds && e.name.find("periapsis") != std::string::npos) {
            scene.update(e.t_tdb);
            const std::vector<double> next = scene.periapsis_times(ds, e.t_tdb, e.t_tdb + 3.0 * 3600.0);
            check(next.size() == 1, "periapsis event is 0-3 h before a periapsis",
                  next.empty() ? -1.0 : (next[0] - e.t_tdb) / 3600.0);
        }
        if (e.name == "Death Star: apsidal advance (top-down)") {
            check(e.from_orbit_normal == ds, "face-on event looks along the Death Star's orbit normal");
        }
    }
}

// The Newtonian ghosts of S2 and the Death Star start where the integrations do.
void test_sgr_a_ghosts()
{
    Scene scene = load_scene_or_die("sgr_a.toml");
    const int s2 = scene.find("S2");
    const int ghost = scene.find("S2 (Newton)");
    const int ds = scene.find("Death Star");
    const int ds_ghost = scene.find("Death Star (Newton)");
    check(s2 > 0 && ghost > 0 && ds > 0 && ds_ghost > 0, "Sgr A* ghost bodies");
    if (s2 <= 0 || ghost <= 0 || ds <= 0 || ds_ghost <= 0) {
        return;
    }
    // The integration starts from the elements at their 2010 apocentre (the
    // loader's "previous_apocentre": half the elements' Keplerian period before t_peri).
    const double a = 125.058 / 1000.0 * 8246.7 * kAuKm;
    const double apo = tdb_from_julian_year(2018.37900) - kPi * std::sqrt(a * a * a / 5.65487707e17);
    const double met = glm::length(scene.icrf_state_at(s2, apo).position - scene.icrf_state_at(ghost, apo).position);
    check(met < 1e3, "S2 and its ghost meet at the 2010 apocentre (km)", met);
    const double later = tdb_from_julian_year(2034.4);
    const double apart =
        glm::length(scene.icrf_state_at(s2, later).position - scene.icrf_state_at(ghost, later).position) / kAuKm;
    check(apart > 1.0, "S2 leaves its ghost by the next pericentre (au)", apart);
    std::printf("info: S2 vs Newton at the 2034 pericentre: %.1f au apart\n", apart);

    // The rolling Death Star restarts from its elements at the requested time.
    const double t = tdb_from_julian_year(2026.0);
    const double start = glm::length(scene.icrf_state_at(ds, t).position - scene.icrf_state_at(ds_ghost, t).position);
    check(start < 1e3, "Death Star starts on its ghost's ellipse (km)", start);
}


// u^t = dt/dtau of a prograde circular equatorial orbit at Boyer-Lindquist r
// (units of M): the rate at which a distant clock runs against the orbiting
// one. Omega from Bardeen, Press & Teukolsky 1972 (ApJ 178, 347), eq. 2.16,
// with the equatorial Kerr metric.
double circular_orbit_ut(double r, double a)
{
    const double omega = 1.0 / (r * std::sqrt(r) + a);
    const double g_tt = -(1.0 - 2.0 / r);
    const double g_tphi = -2.0 * a / r;
    const double g_phiphi = r * r + a * a + 2.0 * a * a / r;
    return 1.0 / std::sqrt(-(g_tt + 2.0 * g_tphi * omega + g_phiphi * omega * omega));
}

// Gargantua (Interstellar): the hole, disk and camera of James et al. 2015, and
// the numbers in the captions.
void test_gargantua_scene()
{
    Scene scene = load_scene_or_die("gargantua.toml");
    const Body& hole = scene.bodies[0];
    const double m = hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    check(hole.kind == BodyKind::BlackHole && hole.spin == 0.6, "Gargantua spins at a/M = 0.6 as filmed");
    check(std::abs(hole.gm_km3_s2 / kSunGmKm3S2 / 1e8 - 1.0) < 1e-9, "Gargantua weighs 1e8 M_sun");
    check(hole.disk_inner_m == 9.26 && hole.disk_outer_m == 18.70 && hole.disk_uniform &&
              hole.disk_temperature_k == 4500.0,
          "Gargantua's disk: 9.26-18.70 M at a uniform 4500 K");
    check(hole.disk_inner_m > isco_radius(hole.spin), "Gargantua's disk is outside the ISCO");

    // The camera of JA15's Figure 13: r_c = 74.1 M, theta_c = 1.511 (86.56 deg from the axis).
    check(scene.events.size() == 8, "Gargantua events", static_cast<double>(scene.events.size()));
    if (scene.events.size() < 3) {
        return;
    }
    for (size_t i = 0; i < 2; ++i) {
        const SceneEvent& e = scene.events[i];
        check(std::abs(e.distance_km / m - 74.1) < 0.05, "Gargantua camera at 74.1 M", e.distance_km / m);
        // 1.511 rad is given to +-0.0005 rad (0.03 deg); JA15 also write 86.56 deg.
        check(std::abs(e.elevation_deg - (90.0 - 1.511 / kDegToRad)) < 0.03, "Gargantua camera 86.56 deg off the axis",
              e.elevation_deg);
    }
    check(!scene.events[0].disk_shifts && scene.events[1].disk_shifts, "Gargantua as filmed, then with the shifts");

    // Caption 2: the gas at the inner edge moves at "about a third of the speed of
    // light", measured by the locally non-rotating observer (Bardeen et al. 1972, eq. 3.8).
    {
        const double r = hole.disk_inner_m;
        const double a = hole.spin;
        const double delta = r * r - 2.0 * r + a * a;
        const double big_a = (r * r + a * a) * (r * r + a * a) - a * a * delta;
        const double omega = 1.0 / (r * std::sqrt(r) + a);
        const double v = (omega - 2.0 * a * r / big_a) * big_a / (r * r * std::sqrt(delta));
        check(v > 0.3 && v < 0.37, "Gargantua's disk edge at about c/3", v);
    }

    // Caption 3: the ISCO at a/M = 0.6 (3.83 M) slows time 1.76-fold.
    const double ut_06 = circular_orbit_ut(isco_radius(0.6), 0.6);
    check(std::abs(ut_06 - 1.76) < 0.005, "time dilation at the a/M = 0.6 ISCO", ut_06);

    // The example planet on the ISCO: prograde in the disk plane at the Kerr-Schild
    // radius of r = 3.83 M, once every 6.96 hours.
    const int planet = scene.find("Planet at the ISCO");
    check(planet > 0, "Gargantua's example planet");
    if (planet > 0) {
        const double r_isco = isco_radius(hole.spin);
        const double period = 2.0 * kPi * (r_isco * std::sqrt(r_isco) + hole.spin) * m / kSpeedOfLightKmS;
        check(std::abs(period / 3600.0 - 6.96) < 0.01, "ISCO period 6.96 h", period / 3600.0);
        const double t0 = tdb_from_julian_year(2026.8);
        const State s0 = scene.icrf_state_at(planet, t0);
        const State s_half = scene.icrf_state_at(planet, t0 + 0.5 * period);
        const State s_full = scene.icrf_state_at(planet, t0 + period);
        const double radius = std::sqrt(r_isco * r_isco + hole.spin * hole.spin) * m;
        check(std::abs(glm::length(s0.position) / radius - 1.0) < 1e-9, "example planet on the ISCO");
        check(glm::length(s_full.position - s0.position) < 1e-6 * radius, "example planet back after a period");
        check(glm::length(s_half.position + s0.position) < 1e-6 * radius, "example planet opposite after half");
        check(std::abs(glm::dot(glm::normalize(s0.position), hole.pole)) < 1e-12, "example planet in the disk plane");
        check(glm::dot(glm::cross(s0.position, s0.velocity), hole.pole) > 0.0, "example planet is prograde");
    }

    // Endurance on its assumed orbit: prograde in the disk plane at 30 M, beyond the
    // disk, once every 5.9 days; time aboard 5 percent slower (the last caption), and
    // 4.4 billion km from the center.
    const int craft = scene.find("Endurance");
    check(craft > 0, "Endurance in the Gargantua scene");
    if (craft > 0) {
        const double r = 30.0;
        check(r > hole.disk_outer_m, "Endurance orbits beyond the disk");
        const double period = 2.0 * kPi * (r * std::sqrt(r) + hole.spin) * m / kSpeedOfLightKmS;
        check(std::abs(period / kSecondsPerDay - 5.9) < 0.05, "Endurance's orbit in 5.9 days",
              period / kSecondsPerDay);
        const double ut = circular_orbit_ut(r, hole.spin);
        check(std::abs(ut - 1.05) < 0.005, "time aboard Endurance 5 percent slower", ut);
        const double t0 = tdb_from_julian_year(2026.8);
        const State s0 = scene.icrf_state_at(craft, t0);
        const State s_full = scene.icrf_state_at(craft, t0 + period);
        const double radius = std::sqrt(r * r + hole.spin * hole.spin) * m;
        check(std::abs(glm::length(s0.position) / radius - 1.0) < 1e-9, "Endurance at 30 M");
        check(std::abs(glm::length(s0.position) / 4.4e9 - 1.0) < 0.01, "Endurance 4.4 billion km out",
              glm::length(s0.position));
        check(glm::length(s_full.position - s0.position) < 1e-6 * radius, "Endurance back after a period");
        check(std::abs(glm::dot(glm::normalize(s0.position), hole.pole)) < 1e-12, "Endurance in the disk plane");
        check(glm::dot(glm::cross(s0.position, s0.velocity), hole.pole) > 0.0, "Endurance is prograde");
        const SceneEvent& e = scene.events.back();
        check(e.focus == craft && e.warp == 1.0, "the Endurance event shows its ring at the true rate");
        std::printf("info: Endurance: period %.3f d, u^t %.4f\n", period / kSecondsPerDay, ut);
    }

    // Miller's planet [DH26]: with 1 - a/M = 1.4e-14 the ISCO is at (1 + 3.8e-5) M,
    // time runs ~60,000 times slower there (1 hour : 7 years), and it is circled
    // in 1.7 hours (seen from afar) for a 1e8 M_sun hole.
    const double a_miller = 1.0 - 1.4e-14;
    const double r_miller = isco_radius(a_miller);
    check(std::abs((r_miller - 1.0) / 3.8e-5 - 1.0) < 0.02, "Miller's ISCO at (1 + 3.8e-5) M", r_miller - 1.0);
    const double ut_miller = circular_orbit_ut(r_miller, a_miller);
    const double hour_to_years = ut_miller / (365.25 * 24.0);
    check(std::abs(hour_to_years - 7.0) < 0.2, "an hour on Miller's planet is 7 years", hour_to_years);
    const double period_h = 2.0 * kPi * (r_miller * std::sqrt(r_miller) + a_miller) * m / kSpeedOfLightKmS / 3600.0;
    check(std::abs(period_h - 1.72) < 0.01, "Miller's planet circles in 1.7 hours", period_h);
    std::printf("info: Gargantua: ISCO(0.6) u^t %.3f; Miller u^t %.0f (1 h = %.2f yr), period %.3f h\n", ut_06,
                ut_miller, hour_to_years, period_h);
}

// Passes of Mann's planet by Gargantua: periapsis time and radius (units of M),
// turns swept while within 30 M, Boyer-Lindquist radii where it crosses the
// hole's equatorial plane, and the time from the previous apoapsis.
struct MannPass {
    double t_peri = 0.0;
    double r_peri = 0.0;
    double turns = 0.0;
    std::vector<double> crossings;
    double fall_days = 0.0;
};

std::vector<MannPass> mann_passes(const Scene& scene, int planet, double t0, double t1)
{
    const Body& hole = scene.bodies[0];
    const double m = hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    std::vector<MannPass> passes;
    MannPass pass;
    bool inside = false;
    double r_max = 0.0;
    double t_apo = t0;
    glm::dvec3 prev = scene.icrf_state_at(planet, t0).position / m;
    for (double t = t0; t < t1;) {
        const glm::dvec3 x = scene.icrf_state_at(planet, t).position / m;
        t += glm::length(x) < 40.0 ? 60.0 : 3600.0;
        const double z = glm::dot(x, hole.pole);
        const double r = kerr_schild_r(glm::dvec3(glm::length(x - z * hole.pole), 0.0, z), hole.spin);
        const double rc = glm::length(x);
        if (rc > r_max) {
            r_max = rc;
            t_apo = t;
        }
        if (rc < 30.0) {
            if (!inside) {
                pass = MannPass{};
                pass.r_peri = 1e9;
                pass.fall_days = 0.0;
                inside = true;
            }
            pass.turns += std::acos(std::clamp(glm::dot(glm::normalize(prev), glm::normalize(x)), -1.0, 1.0)) / kTwoPi;
            if (glm::dot(prev, hole.pole) * glm::dot(x, hole.pole) < 0.0) {
                pass.crossings.push_back(r);
            }
            if (rc < pass.r_peri) {
                pass.r_peri = rc;
                pass.t_peri = t;
            }
        } else if (inside) {
            inside = false;
            pass.fall_days = (pass.t_peri - t_apo) / kSecondsPerDay;
            passes.push_back(pass);
            r_max = 0.0;
        }
        prev = x;
    }
    return passes;
}

void test_mann_orbit()
{
    Scene scene = load_scene_or_die("gargantua.toml");
    const int planet = scene.find("Mann's planet");
    check(planet > 0, "Mann's planet");
    if (planet < 0) {
        return;
    }
    const Body& hole = scene.bodies[0];
    const double m = hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    const std::vector<MannPass> passes =
        mann_passes(scene, planet, tdb_from_julian_year(2025.0), tdb_from_julian_year(2028.0));
    check(passes.size() == 6, "Mann's planet passes Gargantua twice a year", static_cast<double>(passes.size()));
    for (size_t i = 0; i < passes.size(); ++i) {
        const MannPass& p = passes[i];
        std::printf("info: Mann pass %.4f r_p %.3f turns %.2f fall %.1f d, crosses the disk plane at",
                    2000.0 + p.t_peri / (kDaysPerJulianYear * kSecondsPerDay), p.r_peri, p.turns, p.fall_days);
        for (double c : p.crossings) {
            std::printf(" %.2f", c);
            // [TH14] Ch. 19: the orbit keeps off the disk.
            check(c < hole.disk_inner_m || c > hole.disk_outer_m, "Mann's planet misses the disk", c);
        }
        std::printf(" M\n");
        check(std::abs(p.r_peri - 5.0) < 0.05, "Mann's periapsis at 5 M", p.r_peri);
        // [TH14] Ch. 19: it flies around Gargantua "once or twice on each swing by".
        check(p.turns > 1.0 && p.turns < 2.0, "Mann's planet whirls once or twice", p.turns);
        // [TH14] note 33: the fall from ~600 radii takes at least forty days.
        if (i > 0) {
            check(p.fall_days > 40.0, "Mann's fall takes over 40 days", p.fall_days);
        }
    }
    // [TH14] Ch. 19: the orbit "precess[es] through a large angle" from one excursion to the next.
    if (passes.size() >= 2) {
        const glm::dvec3 a = scene.icrf_state_at(planet, passes[0].t_peri).position;
        const glm::dvec3 b = scene.icrf_state_at(planet, passes[1].t_peri).position;
        const double turn = std::acos(glm::dot(glm::normalize(a), glm::normalize(b))) / kDegToRad;
        check(turn > 30.0, "Mann's orbit precesses by a large angle", turn);
        std::printf("info: Mann's periapsis turns %.1f deg per orbit\n", turn);
    }

    // [TH14] note 33: from Mann's planet Gargantua's shadow is ~0.9 deg across, which
    // places it ~600 radii out. Traced for this scene's spin, from 600 M in the
    // equatorial plane: the shadow's edges left and right of the center.
    const glm::dvec3 cam(600.0, 0.0, 0.0);
    double width = 0.0;
    for (double side : {-1.0, 1.0}) {
        double lo = 0.0;
        double hi = 0.02; // rad
        for (int k = 0; k < 40; ++k) {
            const double mid = 0.5 * (lo + hi);
            const glm::dvec3 dir(-std::cos(mid), side * std::sin(mid), 0.0);
            const RayResult ray = trace_ray(cam, dir, hole.spin, 1000.0, 0.0, 0.0, 20000);
            (ray.outcome == RayOutcome::Captured ? lo : hi) = mid;
        }
        width += lo;
    }
    check(std::abs(width / kDegToRad - 0.9) < 0.1, "Gargantua is 0.9 deg wide from 600 M", width / kDegToRad);
    std::printf("info: Gargantua's shadow from 600 M: %.3f deg (%.0f km = %.0f au away)\n", width / kDegToRad,
                600.0 * m, 600.0 * m / kAuKm);
}

// Full angular width (rad) of a black hole's shadow across the spin axis's
// normal direction, seen from cam (Kerr-Schild Cartesian, units of M).
double shadow_width(const glm::dvec3& cam, double a)
{
    const glm::dvec3 in = -glm::normalize(cam);
    const glm::dvec3 side = glm::normalize(glm::cross(in, glm::dvec3(0.0, 0.0, 1.0)));
    double width = 0.0;
    for (double s : {-1.0, 1.0}) {
        double lo = 0.0;
        double hi = 0.05;
        for (int k = 0; k < 40; ++k) {
            const double mid = 0.5 * (lo + hi);
            const RayResult ray = trace_ray(cam, std::cos(mid) * in + s * std::sin(mid) * side, a,
                                            1.5 * glm::length(cam), 0.0, 0.0, 40000);
            (ray.outcome == RayOutcome::Captured ? lo : hi) = mid;
        }
        width += lo;
    }
    return width;
}

// Bolometric flux (W/m^2) on a static observer at cam (units of M) from a thin
// disk at uniform temperature t_k between r_in and r_out, traced through the hole's
// lensing: each pixel's ray is followed back to its first crossing of the disk,
// where the emitter's bolometric intensity sigma T^4 / pi arrives scaled by g^4
// (shifts) or unchanged. Pixels cover a box of +-half_width rad around the hole.
double traced_disk_flux(const glm::dvec3& cam, double a, double r_in, double r_out, double t_k, bool shifts,
                        double half_width, int n)
{
    constexpr double kSigma = 5.670374419e-8; // W m^-2 K^-4, CODATA 2018 (exact)
    const glm::dvec3 forward = -glm::normalize(cam);
    const glm::dvec3 right = glm::normalize(glm::cross(forward, glm::dvec3(0.0, 0.0, 1.0)));
    const glm::dvec3 up = glm::cross(right, forward);
    const ObserverFrame frame = static_observer_frame(cam, a, forward, up, right);
    const double r_plus = horizon_radius(a);
    const double cell = 2.0 * half_width / n;
    double flux = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double u = -half_width + (i + 0.5) * cell;
            const double v = -half_width + (j + 0.5) * cell;
            Photon ph = frame_photon(cam, frame, glm::dvec3(u, v, 1.0));
            for (int step = 0; step < 4000; ++step) {
                const Photon prev = ph;
                rk4_step(ph, a, ray_step(ph.x, a));
                if (prev.x.z * ph.x.z < 0.0) {
                    Photon hit = prev;
                    hit.x = prev.x + (ph.x - prev.x) * (prev.x.z / (prev.x.z - ph.x.z));
                    const double rc = kerr_schild_r(glm::dvec3(hit.x.x, hit.x.y, 0.0), a);
                    if (rc >= r_in && rc <= r_out) {
                        const double g = shifts ? disk_redshift(hit, a, frame.ut) : 1.0;
                        flux += kSigma * std::pow(t_k, 4.0) / kPi * std::pow(g, 4.0) * cell * cell;
                        break;
                    }
                }
                const double r = kerr_schild_r(ph.x, a);
                if (r < r_plus + 0.02 || glm::length(ph.x) > 1.5 * glm::length(cam)) {
                    break;
                }
            }
        }
    }
    return flux;
}

// Edmunds' planet on the scene's assumed orbit: Gargantua's shadow about as wide
// as the Sun from Earth, with the film's disk supplying about twice Earth's sunlight at apoapsis.
void test_edmunds_planet()
{
    Scene scene = load_scene_or_die("gargantua.toml");
    const Body& hole = scene.bodies[0];
    const int planet = scene.find("Edmunds' planet");
    check(planet > 0, "Edmunds' planet");
    if (planet < 0) {
        return;
    }
    const double m = hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    const double year = kDaysPerJulianYear * kSecondsPerDay;

    // The orbit as assumed: 980-1020 au, 30 deg to the disk, 3.16 years.
    const double t0 = tdb_from_julian_year(2026.5);
    double r_min = 1e300;
    double r_max = 0.0;
    double elevation_max = 0.0;
    double t_next_apo = 0.0;
    double r_next_apo = 0.0;
    for (double t = t0; t < t0 + 3.5 * year; t += kSecondsPerDay) {
        const glm::dvec3 x = scene.icrf_state_at(planet, t).position;
        const double r = glm::length(x);
        r_min = std::min(r_min, r);
        r_max = std::max(r_max, r);
        elevation_max = std::max(elevation_max, std::asin(glm::dot(x, hole.pole) / r));
        if (t > t0 + 2.5 * year && r > r_next_apo) {
            r_next_apo = r;
            t_next_apo = t;
        }
    }
    check(std::abs(r_max / kAuKm - 1020.0) < 0.5 && std::abs(r_min / kAuKm - 980.0) < 0.5,
          "Edmunds' orbit 980-1020 au", r_min / kAuKm);
    check(std::abs(elevation_max / kDegToRad - 30.0) < 0.2, "Edmunds' orbit 30 deg to the disk",
          elevation_max / kDegToRad);
    double t_apo = t0;
    double r_apo = 0.0;
    for (double t = t0; t < t0 + 0.5 * year; t += kSecondsPerDay) {
        const double r = glm::length(scene.icrf_state_at(planet, t).position);
        if (r > r_apo) {
            r_apo = r;
            t_apo = t;
        }
    }
    const double period = (t_next_apo - t_apo) / year;
    check(std::abs(period - 3.16) < 0.03, "Edmunds' year is 3.2 Earth years", period);

    // At the event (the apoapsis, 30 deg above the disk).
    const glm::dvec3 x = scene.icrf_state_at(planet, t_apo).position;
    const double d = glm::length(x);
    const double beta = std::asin(glm::dot(x, hole.pole) / d);
    const glm::dvec3 cam = glm::dvec3(std::cos(beta), 0.0, std::sin(beta)) * (d / m);
    const double shadow = shadow_width(cam, hole.spin);
    const double sun = 2.0 * kSunRadiusKm / kAuKm;
    check(std::abs(shadow / sun - 1.0) < 0.15, "Gargantua's shadow as wide as the Sun from Edmunds'",
          shadow / kDegToRad);

    // Flux from the uniform disk, first flat: both faces radiate sigma T^4; seen from
    // elevation beta at distance d it is sigma T^4 A sin(beta) / (pi d^2), A = one face.
    constexpr double kSigma = 5.670374419e-8;   // W m^-2 K^-4, CODATA 2018 (exact)
    constexpr double kSunLuminosityW = 3.828e26; // IAU 2015 Resolution B3, nominal
    const double r_in = hole.disk_inner_m * m * 1e3;
    const double r_out = hole.disk_outer_m * m * 1e3;
    const double area = kPi * (r_out * r_out - r_in * r_in);
    const double d_m = d * 1e3;
    const double flat = kSigma * std::pow(hole.disk_temperature_k, 4.0) * area * std::sin(beta) / (kPi * d_m * d_m);
    const double au_m = kAuKm * 1e3;
    const double earth = kSunLuminosityW / (4.0 * kPi * au_m * au_m); // 1361 W/m^2

    // Then traced through the lensing (which adds the far side's bent image), without
    // and with the Doppler and gravitational shifts. 80 x 80 rays are within ~1.5% of
    // 240 x 240.
    const double half_width = 1.2 * hole.disk_outer_m / (d / m);
    auto traced = [&](const glm::dvec3& c, bool shifts) {
        return traced_disk_flux(c, hole.spin, hole.disk_inner_m, hole.disk_outer_m, hole.disk_temperature_k, shifts,
                                half_width, 80);
    };
    const double lensed = traced(cam, false);
    const double shifted = traced(cam, true);
    check(lensed / flat > 1.1 && lensed / flat < 1.35, "lensing adds to the disk's light", lensed / flat);
    // Caption: "twice Earth's sunlight" 30 deg above the disk.
    check(shifted / earth > 2.0 && shifted / earth < 2.4, "Edmunds' gets twice Earth's sunlight", shifted / earth);
    // Caption: in the disk's plane, bent light "still brings half".
    const double in_plane = traced(glm::dvec3(1000.0 * kAuKm / m, 0.0, 0.0), true);
    check(in_plane / earth > 0.45 && in_plane / earth < 0.56, "Edmunds' gets half in the disk plane",
          in_plane / earth);
    std::printf("info: Edmunds': shadow %.3f deg (Sun %.3f); 30 deg above the disk: flat %.0f, lensed %.0f, "
                "shifted %.0f W/m^2 (%.2f x Earth's %.0f); in its plane %.0f (%.2f x)\n",
                shadow / kDegToRad, sun / kDegToRad, flat, lensed, shifted, shifted / earth, earth, in_plane,
                in_plane / earth);
}

// Point lens with the deflection to second order, alpha(b) = 4/b + 15 pi/(4 b^2)
// (Keeton & Petters 2005, Phys. Rev. D 72, 104006, with a = 0; also the total of
// weak_pass_deflection): beta = theta - (D_ls / D_s) alpha(theta D_l). The image
// angle from the hole on the source's side (side = 1), or the other image's
// (side = -1, as a positive angle).
double point_lens_angle(double beta, double d_l, double d_ls, double side)
{
    const double d_s = d_l + d_ls;
    auto lens = [&](double theta) {
        const double b = theta * d_l;
        const double alpha = 4.0 / b + 15.0 * kPi / (4.0 * b * b);
        return theta - d_ls / d_s * alpha - side * beta;
    };
    double lo = 1e-6;
    double hi = 1.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        (lens(mid) < 0.0 ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

double angle_between(const glm::dvec3& a, const glm::dvec3& b)
{
    return std::atan2(glm::length(glm::cross(a, b)), glm::dot(a, b));
}

// Images of a point bent around a black hole (ephem/kerr_lens.hpp).
void test_kerr_lens()
{
    // Far from the hole (b ~ 140 M): the point-lens images, as angles from the hole seen
    // from the camera, and their magnifications (Paczynski's formula).
    const double d_l = 1e4;
    const double d_ls = 1e4;
    const double d_s = d_l + d_ls;
    const double theta_e = std::sqrt(4.0 * d_ls / (d_l * d_s));
    const glm::dvec3 cam(-d_l, 0.0, 0.0);
    const glm::dvec3 toward_hole(1.0, 0.0, 0.0);
    for (double u : {0.3, 1.0}) {
        const double beta = u * theta_e;
        const glm::dvec3 target(d_ls, d_s * std::tan(beta), 0.0);
        LensSettings settings;
        settings.tolerance_rad = 1e-9;
        LensImage images[2];
        const int found = find_lens_images(cam, target, settings, nullptr, 0, images);
        check(found == 3, "point lens: both images found", found);
        if (found != 3) {
            continue;
        }
        const double theta_plus = angle_between(images[0].direction, toward_hole);
        const double theta_minus = angle_between(images[1].direction, toward_hole);
        const double expect_plus = point_lens_angle(beta, d_l, d_ls, 1.0);
        const double expect_minus = point_lens_angle(beta, d_l, d_ls, -1.0);
        std::printf("info: point lens u = %.1f: theta+ %.6f (lens equation %.6f), theta- %.6f (%.6f); "
                    "mu %.4f, %.4f\n",
                    u, theta_plus, expect_plus, theta_minus, expect_minus, images[0].magnification,
                    images[1].magnification);
        check(std::abs(theta_plus / expect_plus - 1.0) < 2e-3, "point lens: primary image angle",
              theta_plus / expect_plus);
        check(std::abs(theta_minus / expect_minus - 1.0) < 2e-3, "point lens: secondary image angle",
              theta_minus / expect_minus);
        check(images[0].direction.y > 0.0 && images[1].direction.y < 0.0,
              "point lens: the primary on the source's side, the secondary opposite");
        const double total = (u * u + 2.0) / (u * std::sqrt(u * u + 4.0));
        const double mu_sum = std::abs(images[0].magnification) + std::abs(images[1].magnification);
        const double mu_diff = images[0].magnification + images[1].magnification; // the secondary is mirrored
        check(std::abs(mu_sum / total - 1.0) < 0.03, "point lens: total magnification", mu_sum / total);
        check(std::abs(mu_diff - 1.0) < 0.03, "point lens: mu+ - |mu-| = 1", mu_diff);
        check(images[1].magnification < 0.0 && images[0].magnification > 1.0,
              "point lens: the secondary image is mirrored, the primary brightened");
        // The image of a small sphere is drawn out along the circle around the hole.
        const glm::dvec3 radial =
            glm::normalize(toward_hole - glm::dot(toward_hole, images[0].direction) * images[0].direction);
        check(std::abs(glm::dot(images[0].axis[1], radial)) < 1e-3 &&
                  images[0].angle_per_m[1] > images[0].angle_per_m[0],
              "point lens: the primary image is longest around the hole");

        // From the last solution (the target moved a little), the same images again.
        const glm::dvec3 moved = target + glm::dvec3(0.0, 0.2, 0.1);
        LensImage cold[2];
        LensImage warm[2];
        const int cold_found = find_lens_images(cam, moved, settings, nullptr, 0, cold);
        const int warm_found = find_lens_images(cam, moved, settings, images, found, warm);
        check(cold_found == 3 && warm_found == 3 && angle_between(cold[0].direction, warm[0].direction) < 1e-8 &&
                  angle_between(cold[1].direction, warm[1].direction) < 1e-8,
              "point lens: refining the last solution finds the same images");
    }

    // A target right behind the hole: an Einstein ring, both images at theta_E.
    {
        LensSettings settings;
        settings.tolerance_rad = 1e-9;
        LensImage images[2];
        const int found = find_lens_images(cam, glm::dvec3(d_ls, 1e-9, 0.0), settings, nullptr, 0, images);
        const double ring = point_lens_angle(0.0, d_l, d_ls, 1.0);
        check(found == 3 && std::abs(angle_between(images[0].direction, toward_hole) / ring - 1.0) < 2e-3 &&
                  std::abs(angle_between(images[1].direction, toward_hole) / ring - 1.0) < 2e-3,
              "point lens: a source behind the hole makes a ring");
    }

    // A target near the camera and off to the side: hardly bent, nor brightened.
    {
        LensSettings settings;
        LensImage images[2];
        const glm::dvec3 target = cam + glm::dvec3(5.0, 30.0, -2.0);
        const int found = find_lens_images(cam, target, settings, nullptr, 0, images);
        check((found & 1) && angle_between(images[0].direction, glm::normalize(target - cam)) < 1e-5 &&
                  std::abs(images[0].magnification - 1.0) < 1e-3,
              "lens: a nearby target is seen where it is", images[0].magnification);
    }

    // Gargantua (a/M = 0.6) from the camera of [JA15] Figure 13 (74.1 M, 3.44 deg above
    // the disk): the example planet on the ISCO, straight behind the hole. Unlensed it
    // would be hidden by the shadow; its light comes over the top of the shadow, in
    // front of the disk (whose far side lies farther out), and under the bottom,
    // through the disk.
    {
        const double a = 0.6;
        const double elevation = 3.44 * kDegToRad;
        const glm::dvec3 camera(-74.1 * std::cos(elevation), 0.0, 74.1 * std::sin(elevation));
        const glm::dvec3 planet(isco_radius(a), 0.0, 0.0);
        LensSettings settings;
        settings.spin = a;
        settings.disk_r_in = 9.26;
        settings.disk_r_out = 18.70;
        LensImage images[2];
        const int found = find_lens_images(camera, planet, settings, nullptr, 0, images);
        check(found == 3, "Gargantua: both images of the ISCO planet behind the hole", found);
        const glm::dvec3 hole = -glm::normalize(camera);
        const glm::dvec3 up = glm::normalize(glm::dvec3(0.0, 0.0, 1.0) - hole.z * hole);
        const double straight = angle_between(glm::normalize(planet - camera), hole);
        const double first = angle_between(images[0].direction, hole);
        const double second = angle_between(images[1].direction, hole);
        std::printf("info: Gargantua ISCO planet behind the hole: straight %.4f rad, images %.4f (up %+.4f) and "
                    "%.4f (up %+.4f) rad from the center, mu %.3f, %.3f, behind the disk %d %d\n",
                    straight, first, glm::dot(images[0].direction, up), second,
                    glm::dot(images[1].direction, up), images[0].magnification, images[1].magnification,
                    images[0].disk_first, images[1].disk_first);
        // The shadow's edge is at an impact parameter of ~4.5-6 M at this spin: 0.06-0.08 rad.
        check(straight < 0.05 && std::min(first, second) > 0.06,
              "Gargantua: the planet behind the hole shows outside the shadow", std::min(first, second));
        const bool first_on_top = glm::dot(images[0].direction, up) > glm::dot(images[1].direction, up);
        const LensImage& top = first_on_top ? images[0] : images[1];
        const LensImage& bottom = first_on_top ? images[1] : images[0];
        check(glm::dot(top.direction, up) > 0.0 && glm::dot(bottom.direction, up) < 0.0,
              "Gargantua: one image above the shadow, one below");
        check(!top.disk_first && bottom.disk_first,
              "Gargantua: the image over the top is in front of the disk, the one below behind it");
    }
}

// The scene's bodies are drawn at their images (Simulation::compute_view, lens_bodies).
void test_gargantua_lensed_bodies()
{
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/gargantua.toml", &error), "gargantua.toml loads");
    sim.jump_to_event(0);
    sim.update(0.1);
    OutputView view;
    sim.compute_view(1280, 720, 1.0, view);
    const Scene& scene = sim.scene();
    int lensed = 0;
    for (size_t i = 0; i < scene.bodies.size() && i < view.lens.size(); ++i) {
        const Body& body = scene.bodies[i];
        const bool child =
            body.parent >= 0 && scene.bodies[static_cast<size_t>(body.parent)].kind == BodyKind::BlackHole;
        check(!view.lens[i].lensed || child, "Gargantua: only bodies around the hole are lensed");
        if (view.lens[i].lensed) {
            ++lensed;
            const LensedView& image = view.lens[i].images[0];
            check(image.angular_radius_per_km > 0.0 && glm::length(image.offset) > 0.0 &&
                      std::abs(glm::determinant(image.turn) - 1.0) < 1e-9,
                  "Gargantua: a lensed body has a place, a size and a rotation");
        }
    }
    check(view.lens.size() == scene.bodies.size() && lensed == 4,
          "Gargantua: the four bodies around the hole are lensed", lensed);
    // Again from the last frame's solution, a little later.
    sim.update(0.05);
    sim.compute_view(1280, 720, 1.0, view);
    int still = 0;
    for (const BodyLens& lens : view.lens) {
        still += lens.lensed ? 1 : 0;
    }
    check(still == 4, "Gargantua: refined images next frame", still);

    // Endurance's own event: the camera 160 m from it, which is then drawn where it is.
    sim.jump_to_event(scene.events.size() - 1);
    for (int i = 0; i < 100; ++i) {
        sim.update(0.1); // the camera's flight
    }
    sim.compute_view(1280, 720, 1.0, view);
    const int craft = scene.find("Endurance");
    check(craft >= 0 && !view.lens[static_cast<size_t>(craft)].lensed, "Gargantua: Endurance seen from up close");
}

// Hot accretion flow (RIAF, src/ephem/riaf.cpp).

RiafParams test_riaf(double n_e0)
{
    RiafParams p;
    p.n_e0_cm3 = n_e0;
    p.t_e0_k = 1e11;
    p.beta = 10.0;
    p.outer_m = 50.0;
    return p;
}

// Sgr A*'s hot flow as seen from the Earth at 230 GHz, traced from a camera at
// 1000 M along the line of sight (the flow is axisymmetric: only the angle between
// the spin and the line of sight matters).
struct RiafImage {
    double flux_jy = 0.0;
    double ring_radius_m = 0.0; // peak of the azimuthally averaged brightness, after the blur
};

// blur_fwhm_m: Gaussian blur (FWHM, units of M) before the ring is measured, as an
// interferometer's resolution blurs it.
RiafImage riaf_image(double a, double gm_km3_s2, const RiafParams& flow, double inclination, double distance_pc, int grid,
                     double half_width_m, double blur_fwhm_m)
{
    const double m_km = gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    const double r_cam = 1000.0;
    const glm::dvec3 cam = r_cam * glm::dvec3(std::sin(inclination), 0.0, std::cos(inclination));
    const glm::dvec3 forward = -glm::normalize(cam);
    const glm::dvec3 right = glm::normalize(glm::cross(forward, glm::dvec3(0.0, 0.0, 1.0)));
    const glm::dvec3 up = glm::cross(right, forward);
    const ObserverFrame frame = static_observer_frame(cam, a, forward, up, right);
    const double pixel = 2.0 * half_width_m / grid / r_cam; // radians at the camera
    const auto n = static_cast<size_t>(grid);
    std::vector<double> tb(n * n);
    double tb_sum = 0.0;
    for (int iy = 0; iy < grid; ++iy) {
        for (int ix = 0; ix < grid; ++ix) {
            const double sx = (ix + 0.5 - 0.5 * grid) * pixel;
            const double sy = (iy + 0.5 - 0.5 * grid) * pixel;
            const FlowRay ray = trace_flow_ray(frame_photon(cam, frame, glm::dvec3(sx, sy, 1.0)), a, m_km * 1e5, 230e9,
                                               flow, 1.01 * r_cam);
            tb[static_cast<size_t>(iy) * n + static_cast<size_t>(ix)] = ray.tb_k;
            tb_sum += ray.tb_k;
        }
    }
    // F = (2 k nu^2 / c^2) sum T_b dOmega, the solid angles shrunk from the camera's distance to ours.
    constexpr double kBoltzmann = 1.380649e-16;
    constexpr double kLightCm = 2.99792458e10;
    const double scale = r_cam * m_km / (distance_pc * kParsecKm);
    RiafImage out;
    out.flux_jy = 2.0 * kBoltzmann * 230e9 * 230e9 / (kLightCm * kLightCm) * tb_sum * pixel * pixel * scale * scale / 1e-23;

    // Separable Gaussian blur.
    const double pixel_m = 2.0 * half_width_m / grid;
    const double sigma = blur_fwhm_m / 2.3548200450309493 / pixel_m;
    if (sigma > 0.0) {
        const int reach = static_cast<int>(std::ceil(3.0 * sigma));
        std::vector<double> kernel(static_cast<size_t>(2 * reach + 1));
        for (int k = -reach; k <= reach; ++k) {
            kernel[static_cast<size_t>(k + reach)] = std::exp(-0.5 * k * k / (sigma * sigma));
        }
        std::vector<double> tmp(n * n, 0.0);
        for (int pass = 0; pass < 2; ++pass) {
            std::fill(tmp.begin(), tmp.end(), 0.0);
            for (int iy = 0; iy < grid; ++iy) {
                for (int ix = 0; ix < grid; ++ix) {
                    double sum = 0.0;
                    for (int k = -reach; k <= reach; ++k) {
                        const int jx = pass == 0 ? ix + k : ix;
                        const int jy = pass == 0 ? iy : iy + k;
                        if (jx >= 0 && jx < grid && jy >= 0 && jy < grid) {
                            sum += kernel[static_cast<size_t>(k + reach)] * tb[static_cast<size_t>(jy) * n + static_cast<size_t>(jx)];
                        }
                    }
                    tmp[static_cast<size_t>(iy) * n + static_cast<size_t>(ix)] = sum;
                }
            }
            tb.swap(tmp);
        }
    }

    constexpr int kBins = 60;
    double ring_sum[kBins] = {};
    int ring_count[kBins] = {};
    for (int iy = 0; iy < grid; ++iy) {
        for (int ix = 0; ix < grid; ++ix) {
            const double radius = std::hypot(ix + 0.5 - 0.5 * grid, iy + 0.5 - 0.5 * grid) * pixel_m;
            const int bin = static_cast<int>(radius / half_width_m * kBins);
            if (bin < kBins) {
                ring_sum[bin] += tb[static_cast<size_t>(iy) * n + static_cast<size_t>(ix)];
                ++ring_count[bin];
            }
        }
    }
    double best = 0.0;
    for (int k = 0; k < kBins; ++k) {
        const double mean = ring_count[k] > 0 ? ring_sum[k] / ring_count[k] : 0.0;
        if (mean > best) {
            best = mean;
            out.ring_radius_m = (k + 0.5) * half_width_m / kBins;
        }
    }
    return out;
}

// A ray that spirals into a fast-spinning hole (a = 0.94) gets numerically unstable near
// the horizon before reaching it (H grew to 1e5 and the ray was thrown back out). It
// is captured, and stays finite.
void test_near_horizon_spiral_captured()
{
    const double a = 0.94;
    const double r_cam = 1000.0;
    const double inclination = 150.0 * kDegToRad;
    const glm::dvec3 cam = r_cam * glm::dvec3(std::sin(inclination), 0.0, std::cos(inclination));
    const glm::dvec3 forward = -glm::normalize(cam);
    const glm::dvec3 right = glm::normalize(glm::cross(forward, glm::dvec3(0.0, 0.0, 1.0)));
    const ObserverFrame frame = static_observer_frame(cam, a, forward, glm::cross(right, forward), right);
    const double pixel = 30.0 / 96.0 / r_cam;
    const glm::dvec3 n((41.5 - 48.0) * pixel, (51.5 - 48.0) * pixel, 1.0);
    const RayResult ray = trace_ray(frame_photon(cam, frame, n), a, 1.01 * r_cam, 0.0, 0.0);
    check(ray.outcome == RayOutcome::Captured, "near-horizon spiral captured (trace_ray)",
          static_cast<double>(ray.outcome));
    const FlowRay flow = trace_flow_ray(frame_photon(cam, frame, n), a, 6.2921e11, 230e9, test_riaf(1e7), 1.01 * r_cam);
    check(flow.outcome == RayOutcome::Captured && std::isfinite(flow.tb_k) && std::isfinite(flow.tau),
          "near-horizon spiral captured and finite (trace_flow_ray)", flow.tb_k);
}

// Sgr A* at 230 GHz as seen from the Earth, north up and east left: the camera at 1000 M
// toward us in the hole frame of hole_frame_from_pole(). Returns the flux (Jy) and the
// brightness centroid (west, north; units of M) of the image.
struct SkyImage {
    double flux_jy = 0.0;
    glm::dvec2 centroid_m{0.0};
};

SkyImage riaf_sky_image(const Body& hole, const HotSpot* spot, int grid, double half_width_m)
{
    const double a = hole.spin;
    const double m_km = hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    const double distance_pc = 8246.7;
    const double ra = 266.41683333 * kDegToRad;
    const double dec = -29.00781611 * kDegToRad;
    const glm::dmat3 to_hole = glm::transpose(hole_frame_from_pole(hole.pole));
    const glm::dvec3 toward_earth = -unit_toward(266.41683333, -29.00781611);
    const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
    const glm::dvec3 west(std::sin(ra), -std::cos(ra), 0.0);
    const double r_cam = 1000.0;
    const glm::dvec3 cam = r_cam * (to_hole * toward_earth);
    const ObserverFrame frame = static_observer_frame(cam, a, -glm::normalize(cam), to_hole * north, to_hole * west);
    const double pixel = 2.0 * half_width_m / grid / r_cam;
    SkyImage out;
    double sum = 0.0;
    for (int iy = 0; iy < grid; ++iy) {
        for (int ix = 0; ix < grid; ++ix) {
            const double sx = (ix + 0.5 - 0.5 * grid) * pixel;
            const double sy = (iy + 0.5 - 0.5 * grid) * pixel;
            const double tb = trace_flow_ray(frame_photon(cam, frame, glm::dvec3(sx, sy, 1.0)), a, m_km * 1e5, 230e9,
                                             *hole.hot_flow, 1.01 * r_cam, 4000, spot)
                                  .tb_k;
            sum += tb;
            out.centroid_m += tb * glm::dvec2(sx, sy) * r_cam;
        }
    }
    constexpr double kBoltzmann = 1.380649e-16;
    constexpr double kLightCm = 2.99792458e10;
    const double scale = r_cam * m_km / (distance_pc * kParsecKm);
    out.flux_jy = 2.0 * kBoltzmann * 230e9 * 230e9 / (kLightCm * kLightCm) * sum * pixel * pixel * scale * scale / 1e-23;
    out.centroid_m = sum > 0.0 ? out.centroid_m / sum : glm::dvec2(0.0);
    return out;
}

// Sgr A*'s flares: the recorded ones and the random rate, the hot spot's 230 GHz flux at
// the ALMA flare's peak (Wielgus et al. 2022, sec. 3.1: I_hotspot <~ 0.3 Jy; flare_t_s_k is
// set to it), and its sense of motion on the sky - clockwise, as GRAVITY 2018, 2023 and
// ALMA all saw.
void test_sgr_a_flares()
{
    Scene scene = load_scene_or_die("sgr_a.toml");
    const Body& hole = scene.bodies[0];
    const FlareSettings& flares = hole.hot_flow_flares;
    check(flares.recorded.size() == 2 && flares.random_per_day == 1.5 && hole.flare_t_s_k > 0.0,
          "Sgr A* has two recorded flares and random ones");
    if (flares.recorded.size() != 2 || !hole.hot_flow) {
        return;
    }
    const FlareRecord& gravity = flares.recorded[0];
    const FlareRecord& alma = flares.recorded[1];
    check(std::abs(gravity.t_peak_tdb - tdb_from_julian_year(2018.5558)) < 0.5 * kSecondsPerDay && gravity.radius_m == 7.5,
          "GRAVITY flare of 2018-07-22 at 7.5 r_g");
    check(std::abs(alma.t_peak_tdb - tdb_from_julian_year(2017.2754)) < 0.5 * kSecondsPerDay && alma.radius_m == 11.0,
          "ALMA hot spot of 2017-04-11 at 11 r_g");
    int with_earth_view = 0;
    for (const SceneEvent& e : scene.events) {
        for (const FlareRecord* r : {&gravity, &alma}) {
            if (e.from_earth && e.t_tdb < r->t_peak_tdb && r->t_peak_tdb - e.t_tdb < 3600.0) {
                ++with_earth_view;
            }
        }
    }
    check(with_earth_view == 2, "an event seen from the Earth before each recorded flare", with_earth_view);

    HotSpot spot;
    spot.radius_m = alma.radius_m;
    spot.sigma_m = hole.flare_diameter_m / 2.3548200450309493;
    spot.tau0 = hole.flare_tau0;
    spot.t_s_k = hole.flare_t_s_k;
    spot.strength = 1.0;
    const double base = riaf_sky_image(hole, nullptr, 72, 18.0).flux_jy;
    double added = 0.0;
    for (int k = 0; k < 8; ++k) {
        spot.phase_rad = k * kTwoPi / 8.0;
        added += (riaf_sky_image(hole, &spot, 72, 18.0).flux_jy - base) / 8.0;
    }
    std::printf("info: Sgr A* hot spot at 11 r_g adds %.3f Jy at 230 GHz (orbit average)\n", added);
    check(std::abs(added / 0.3 - 1.0) < 0.1, "ALMA hot spot flux (Wielgus 2022: <~0.3 Jy)", added);

    // Sense on the sky: the spot's light (image with the spot minus without) turns clockwise,
    // i.e. its position angle (east of north) decreases.
    auto spot_pa = [&](double phase) {
        spot.phase_rad = phase;
        const SkyImage with = riaf_sky_image(hole, &spot, 72, 18.0);
        const SkyImage without = riaf_sky_image(hole, nullptr, 72, 18.0);
        const glm::dvec2 c = (with.centroid_m * with.flux_jy - without.centroid_m * without.flux_jy) /
                             (with.flux_jy - without.flux_jy);
        return std::atan2(-c.x, c.y); // east = -west
    };
    const double turned = wrap_pi(spot_pa(0.5) - spot_pa(0.0));
    check(turned < -0.2, "hot spot moves clockwise on the sky (GRAVITY, ALMA)", turned / kDegToRad);
}

void test_sgr_a_hot_flow()
{
    Scene scene = load_scene_or_die("sgr_a.toml");
    const Body& hole = scene.bodies[0];
    check(hole.hot_flow.has_value() && hole.disk_outer_m == 0.0, "Sgr A* has a hot flow, no thin disk");
    if (!hole.hot_flow) {
        return;
    }
    // The spin axis: 150 deg from the direction to the Earth [EHT24VIII], at position angle 121 deg
    // east of north on the sky [JA25].
    const double ra = 266.41683333;
    const double dec = -29.00781611;
    const glm::dvec3 toward_earth = -unit_toward(ra, dec);
    const double inclination = std::acos(glm::dot(hole.pole, toward_earth));
    check(std::abs(inclination * kRadToDeg - 150.0) < 0.1, "Sgr A* spin 150 deg from the line of sight",
          inclination * kRadToDeg);
    const double ra_r = ra * kDegToRad;
    const double dec_r = dec * kDegToRad;
    const glm::dvec3 east(-std::sin(ra_r), std::cos(ra_r), 0.0);
    const glm::dvec3 north(-std::sin(dec_r) * std::cos(ra_r), -std::sin(dec_r) * std::sin(ra_r), std::cos(dec_r));
    const double pa = std::atan2(glm::dot(hole.pole, east), glm::dot(hole.pole, north)) * kRadToDeg;
    check(std::abs(pa - 121.0) < 0.1, "Sgr A* spin axis at position angle 121 deg", pa);

    // The orbits feel the hole's spin about its own axis (tilted KerrOrbitMotion).
    for (const char* name : {"S2", "S301", "Death Star"}) {
        const int i = scene.find(name);
        const auto* orbit = i >= 0 ? dynamic_cast<const KerrOrbitMotion*>(scene.bodies[static_cast<size_t>(i)].motion.get())
                                   : nullptr;
        check(orbit && orbit->tilted() && orbit->spin() == 0.94, "Sgr A* orbits integrated with its spin and axis");
    }

    // 230 GHz: 2.44 +- 0.23 Jy in 2017 [WI22] (n_e0_cm3 is calibrated to it by bisection on
    // this image). The ring, measured on the image blurred to the array's ~20 microarcsecond
    // resolution, comes out ~37 microarcseconds across against EHT's 51.8 +- 2.3 [EHT22]:
    // the plunging gas inside the photon ring outshines it (limitations.md). The check
    // holds the model's own value so that it does not drift unnoticed.
    const double distance_pc = 8246.7;
    const double theta_g_uas = hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS) / (distance_pc * kParsecKm) *
                               kRadToDeg * 3600.0 * 1e6;
    const RiafImage image =
        riaf_image(hole.spin, hole.gm_km3_s2, *hole.hot_flow, inclination, distance_pc, 96, 15.0, 20.0 / theta_g_uas);
    const double diameter = 2.0 * image.ring_radius_m * theta_g_uas;
    std::printf("info: Sgr A* RIAF at 230 GHz: %.3f Jy, ring %.1f microarcseconds across (EHT: 51.8 +- 2.3)\n",
                image.flux_jy, diameter);
    check(std::abs(image.flux_jy / 2.44 - 1.0) < 0.1, "Sgr A* 230 GHz flux", image.flux_jy);
    check(std::abs(diameter / 37.0 - 1.0) < 0.1, "Sgr A* 230 GHz ring diameter (microarcseconds)", diameter);
}

void test_bessel_k2()
{
    check(std::abs(bessel_k2(10.0) / 2.150981700693e-05 - 1.0) < 1e-6, "K2(10)", bessel_k2(10.0));
    check(std::abs(bessel_k2(5.0) / 5.308943712223e-03 - 1.0) < 1e-6, "K2(5)", bessel_k2(5.0));
    check(std::abs(bessel_k2(1.0) / 1.624838898635 - 1.0) < 1e-6, "K2(1)", bessel_k2(1.0));
    check(std::abs(bessel_k2(0.1) / 199.5039646419 - 1.0) < 1e-6, "K2(0.1)", bessel_k2(0.1));
}

void test_thermal_synchrotron()
{
    // Leung, Gammie & Noble 2011 eq. 72 evaluated by hand (Theta_e = 16.86, X = 55.18).
    const double j = thermal_synchrotron_emissivity(230e9, 1e6, 1e11, 30.0, kRiafSinPitch);
    check(std::abs(j / 2.29262e-16 - 1.0) < 1e-3, "Leung 2011 emissivity", j);
    check(std::abs(thermal_synchrotron_emissivity(230e9, 3e6, 1e11, 30.0, kRiafSinPitch) / j - 3.0) < 1e-12,
          "synchrotron emissivity ~ n_e");
}

// Below Theta_e ~ 0.0116 the shader's K_2(1/Theta) and the emission both underflow
// float (0/0): the reference stops emitting where the shader does, at Theta_e < 0.02.
void test_cold_electrons_do_not_emit()
{
    const double theta_015 = 0.015 * 5.9298e9; // K
    check(thermal_synchrotron_emissivity(230e9, 1e6, theta_015, 30.0, kRiafSinPitch) == 0.0,
          "no synchrotron below Theta_e = 0.02");
    check(thermal_synchrotron_emissivity(230e9, 1e6, 0.03 * 5.9298e9, 30.0, kRiafSinPitch) >= 0.0,
          "synchrotron evaluated above Theta_e = 0.02");
}

// A hot flow whose ring spans only a few pixels is point-sampled by a handful of rays
// and would flicker: it fades out between 1.5 and 4 pixels across 10 M.
void test_hot_flow_far_fade()
{
    check(riaf_far_view_fade(1.0) == 0.0, "hot flow hidden when 10 M is a pixel");
    check(riaf_far_view_fade(1.5) == 0.0 && riaf_far_view_fade(4.0) == 1.0, "hot flow fade range");
    check(riaf_far_view_fade(100.0) == 1.0, "hot flow shown up close");
    bool rising = true;
    for (double px = 1.0; px < 5.0; px += 0.1) {
        rising = rising && riaf_far_view_fade(px + 0.1) >= riaf_far_view_fade(px);
    }
    check(rising, "hot flow fades in monotonically");
}

void test_riaf_plasma()
{
    RiafParams p;
    p.n_e0_cm3 = 1e7;
    p.t_e0_k = 1e11;
    p.beta = 10.0;
    p.outer_m = 50.0;
    const RiafPlasma mid = riaf_plasma(glm::dvec3(20.0, 0.0, 0.0), 0.0, p);
    check(std::abs(mid.n_e_cm3 / 7.9433e5 - 1.0) < 1e-4, "RIAF density at 20 M", mid.n_e_cm3);
    check(std::abs(mid.t_e_k / 1.4454e10 - 1.0) < 1e-4, "RIAF electron temperature at 20 M", mid.t_e_k);
    check(std::abs(mid.b_gauss / 5.00 - 1.0) < 2e-3, "RIAF field at 20 M (beta = 10)", mid.b_gauss);
    const RiafPlasma off = riaf_plasma(glm::dvec3(std::sqrt(300.0), 0.0, 10.0), 0.0, p);
    check(std::abs(off.n_e_cm3 / (7.9433e5 * std::exp(-1.0 / 6.0)) - 1.0) < 1e-4, "RIAF vertical profile",
          off.n_e_cm3);
    check(std::abs(off.t_e_k / mid.t_e_k - 1.0) < 1e-12, "RIAF temperature depends on r only");
    check(std::abs(off.b_gauss / (mid.b_gauss * std::exp(-1.0 / 12.0)) - 1.0) < 1e-9, "RIAF field ~ sqrt(n)");
    const RiafPlasma out = riaf_plasma(glm::dvec3(60.0, 0.0, 0.0), 0.0, p);
    check(out.n_e_cm3 == 0.0 && out.t_e_k == 0.0 && out.b_gauss == 0.0, "no RIAF beyond outer_m");
}

// Kerr-Schild Cartesian point at Boyer-Lindquist (r, theta) and phi = 0.
glm::dvec3 ks_point(double r, double theta, double a)
{
    return glm::dvec3(r * std::sin(theta), a * std::sin(theta), r * std::cos(theta));
}

double ks_norm(const glm::dvec3& x, double a, const glm::dvec4& u)
{
    const KerrSchild k = kerr_schild(x, a);
    const double lu = u.w + glm::dot(k.l, glm::dvec3(u));
    return -u.w * u.w + glm::dot(glm::dvec3(u), glm::dvec3(u)) + k.f * lu * lu;
}

void test_riaf_velocity_grid()
{
    const double a = 0.94;
    const double r_isco = isco_radius(a);
    const double radii[] = {horizon_radius(a) + 0.02, 1.6, r_isco - 1e-6, r_isco + 1e-6, 6.0, 20.0, 50.0};
    const double thetas[] = {1.0, 30.0, 60.0, 90.0, 120.0, 179.0};
    bool timelike = true;
    bool inward = true;
    double worst = 0.0;
    for (double r : radii) {
        for (double th : thetas) {
            const glm::dvec3 x = ks_point(r, th * kDegToRad, a);
            const glm::dvec4 u = riaf_four_velocity(x, a);
            const double n = ks_norm(x, a, u);
            const bool finite = std::isfinite(u.x) && std::isfinite(u.y) && std::isfinite(u.z) && std::isfinite(u.w);
            timelike = timelike && finite && u.w > 0.0;
            worst = std::max(worst, std::abs(n + 1.0));
            if (r < r_isco - 1e-3) {
                // Plunging: the radius decreases along u.
                const double dr = kerr_schild_r(x + 1e-6 * glm::dvec3(u), a) - r;
                inward = inward && dr < 0.0;
            }
        }
    }
    check(timelike, "RIAF velocity finite and future directed");
    check(worst < 1e-9, "RIAF velocity normalized (u.u = -1)", worst);
    check(inward, "RIAF gas plunges inside the ISCO");

    const glm::dvec3 x6 = ks_point(6.0, 0.5 * kPi, a);
    const glm::dvec4 u6 = riaf_four_velocity(x6, a);
    const double omega = (x6.x * u6.y - x6.y * u6.x) / ((x6.x * x6.x + x6.y * x6.y) * u6.w);
    check(std::abs(omega * (std::pow(6.0, 1.5) + a) - 1.0) < 1e-9, "RIAF Keplerian outside the ISCO", omega);

    // Off the equator the plunge starts with u^r ~ sqrt(r_isco - r): continuous, not smooth.
    double jump = 0.0;
    for (double th : thetas) {
        const glm::dvec4 in = riaf_four_velocity(ks_point(r_isco - 1e-10, th * kDegToRad, a), a);
        const glm::dvec4 out = riaf_four_velocity(ks_point(r_isco + 1e-10, th * kDegToRad, a), a);
        const glm::dvec4 d = glm::abs(in - out);
        jump = std::max({jump, d.x, d.y, d.z, d.w});
    }
    check(jump < 1e-4, "RIAF velocity continuous at the ISCO", jump);
}

// Flare hot spot (riaf.cpp): Keplerian period, position and its light along a ray.
void test_hot_spot_period()
{
    const double m_km = 6.2921e6;
    check(std::abs(hot_spot_period_s(0.94, 7.0, m_km) / 60.0 - 42.77) < 0.1, "spot period at 7 r_g (GRAVITY 2018: 40 +- 8 min)",
          hot_spot_period_s(0.94, 7.0, m_km) / 60.0);
    check(std::abs(hot_spot_period_s(0.94, 9.0, m_km) / 60.0 - 61.41) < 0.1, "spot period at 9 r_g (GRAVITY 2023: 60 +- 3 min)",
          hot_spot_period_s(0.94, 9.0, m_km) / 60.0);
    check(std::abs(hot_spot_period_s(0.94, 11.0, m_km) / 60.0 - 82.25) < 0.1, "spot period at 11 r_g (ALMA: 74 +- 6 min)",
          hot_spot_period_s(0.94, 11.0, m_km) / 60.0);
}

void test_hot_spot_center()
{
    HotSpot s;
    s.radius_m = 9.0;
    s.phase_rad = 0.3;
    const double a = 0.94;
    const glm::dvec3 c0 = hot_spot_center(s, a, 0.0);
    check(std::abs(c0.z) < 1e-12 && std::abs(kerr_schild_r(c0, a) - 9.0) < 1e-9, "spot centre on the equator at R");
    const double quarter = 0.25 * kTwoPi * (std::pow(9.0, 1.5) + a); // M
    const glm::dvec3 c1 = hot_spot_center(s, a, quarter);
    const double turned = std::atan2(c0.x * c1.y - c0.y * c1.x, glm::dot(c0, c1));
    check(std::abs(turned - 0.5 * kPi) < 1e-9, "spot turns prograde by pi/2 in a quarter period", turned);
}

// A camera at 60 M in the equatorial plane looking at the hole through a spot at R = 9
// (its phase set where the light leaving it ~9 M before "now" at the hole finds it).
FlowRay spot_ray(const HotSpot* spot, const glm::dvec3& aim)
{
    const double a = 0.94;
    const glm::dvec3 cam(60.0, 0.0, 0.0);
    const ObserverFrame frame =
        static_observer_frame(cam, a, glm::dvec3(-1.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 1.0), glm::dvec3(0.0, 1.0, 0.0));
    const glm::dvec3 d = glm::normalize(aim - cam);
    const glm::dvec3 n(glm::dot(d, glm::dvec3(0.0, 1.0, 0.0)), glm::dot(d, glm::dvec3(0.0, 0.0, 1.0)), -d.x);
    return trace_flow_ray(frame_photon(cam, frame, n), a, 6.2921e11, 230e9, test_riaf(1.0), 150.0, 4000, spot);
}

void test_hot_spot_brightens_ray()
{
    HotSpot spot;
    spot.radius_m = 9.0;
    spot.t_s_k = 1e11;
    spot.strength = 1.0;
    spot.phase_rad = -9.0 / (std::pow(9.0, 1.5) + 0.94); // where it is when the light seen now leaves it
    const glm::dvec3 centre(9.0, 0.94, 0.0);
    const double base = spot_ray(nullptr, centre).tb_k;
    const double through = spot_ray(&spot, centre).tb_k - base;
    check(through > 0.03 * spot.t_s_k && through < 0.2 * spot.t_s_k, "a ray through the spot gains ~g T_s (1 - e^-0.1)",
          through / spot.t_s_k);
    const glm::dvec3 above(9.0, 0.94, 20.0);
    check(std::abs(spot_ray(&spot, above).tb_k - spot_ray(nullptr, above).tb_k) <= 1e-6 * (1.0 + base),
          "a ray 20 M from the spot does not see it");
    HotSpot off = spot;
    off.strength = 0.0;
    check(spot_ray(&off, centre).tb_k == base, "a spot of strength 0 adds nothing");
}

void test_hot_spot_ray_near_spot()
{
    HotSpot spot;
    spot.radius_m = 9.0;
    spot.t_s_k = 1e11;
    spot.strength = 1.0;
    const double a = 0.94;
    const glm::dvec3 cam(8.0, 3.0, 0.5);
    const ObserverFrame frame =
        static_observer_frame(cam, a, glm::dvec3(-1.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 1.0), glm::dvec3(0.0, 1.0, 0.0));
    bool ok = true;
    for (int i = 0; i < 64; ++i) {
        const double z = 1.0 - (i + 0.5) / 32.0;
        const double s = std::sqrt(1.0 - z * z);
        const glm::dvec3 n(s * std::cos(i * 2.39996323), s * std::sin(i * 2.39996323), z);
        const FlowRay ray = trace_flow_ray(frame_photon(cam, frame, n), a, 6.2921e11, 230e9, test_riaf(1e7), 150.0, 4000, &spot);
        ok = ok && std::isfinite(ray.tb_k) && ray.tb_k >= 0.0;
    }
    check(ok, "rays from a camera beside the hot spot are finite");
}

void test_flow_ray_from_inside()
{
    const double a = 0.94;
    const double mass_cm = 6.2921e11; // 4.261e6 M_sun
    const glm::dvec3 cam(30.0, 0.0, 0.0);
    const ObserverFrame frame =
        static_observer_frame(cam, a, glm::dvec3(-1.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 1.0), glm::dvec3(0.0, 1.0, 0.0));
    bool ok = true;
    double brightest = 0.0;
    for (int i = 0; i < 64; ++i) {
        // Fibonacci sphere of local directions.
        const double z = 1.0 - (i + 0.5) / 32.0;
        const double s = std::sqrt(1.0 - z * z);
        const double phi = i * 2.39996323;
        const glm::dvec3 n(s * std::cos(phi), s * std::sin(phi), z);
        const FlowRay ray = trace_flow_ray(frame_photon(cam, frame, n), a, mass_cm, 230e9, test_riaf(1e7), 150.0);
        ok = ok && std::isfinite(ray.tb_k) && ray.tb_k >= 0.0 && std::isfinite(ray.tau) && ray.tau >= 0.0;
        brightest = std::max(brightest, ray.tb_k);
    }
    check(ok, "RIAF rays from a camera inside the flow are finite");
    check(brightest > 0.0, "RIAF seen from inside the flow", brightest);
}

void test_flow_ray_optically_thick()
{
    // Dense enough that the ray toward the hole is opaque: the brightness approaches
    // the source function g T_e at the tau = 1 surface (Kirchhoff) and never exceeds its largest value.
    const double a = 0.94;
    const glm::dvec3 cam(40.0, 0.0, 0.0);
    const ObserverFrame frame =
        static_observer_frame(cam, a, glm::dvec3(-1.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 1.0), glm::dvec3(0.0, 1.0, 0.0));
    const FlowRay ray =
        trace_flow_ray(frame_photon(cam, frame, glm::dvec3(0.0, 0.0, 1.0)), a, 6.2921e11, 230e9, test_riaf(1e12), 150.0);
    check(ray.tau > 20.0, "dense RIAF is optically thick", ray.tau);
    check(ray.tb_k <= ray.max_source_k * (1.0 + 1e-9), "brightness below the largest source function", ray.tb_k);
    check(ray.tb_k > 0.3 * ray.tau1_source_k && ray.tb_k < 3.0 * ray.tau1_source_k,
          "thick brightness ~ source at tau = 1", ray.tb_k / ray.tau1_source_k);
}

// M87* (assets/scenes/m87.toml): the hole's orientation, its jet (ephem/jet.hpp) and hot flow.

constexpr double kM87Ra = 187.705904;
constexpr double kM87Dec = 12.391117;
constexpr double kM87DistancePc = 1.68e7;

const Body& m87_hole(const Scene& scene)
{
    const int i = scene.find("M87*");
    check(i >= 0, "M87* in its scene");
    return scene.bodies[static_cast<size_t>(std::max(i, 0))];
}

// The spin points away from the Earth [EHT19V]: the approaching jet, along the other
// pole, at [CU23]'s precession axis (position angle 288.47 deg, 17.21 deg from us).
void test_m87_orientation()
{
    const Scene scene = load_scene_or_die("m87.toml");
    const Body& hole = m87_hole(scene);
    const SkyAxes sky = sky_axes(kM87Ra, kM87Dec);
    check(std::abs(std::acos(glm::dot(hole.pole, sky.toward_observer)) * kRadToDeg - 162.79) < 0.01,
          "M87* spin 162.79 deg from the line of sight", std::acos(glm::dot(hole.pole, sky.toward_observer)) * kRadToDeg);
    const double jet_pa = wrap_two_pi(sky_position_angle(-hole.pole, sky)) * kRadToDeg;
    check(std::abs(jet_pa - 288.47) < 0.01, "M87* approaching jet at position angle 288.47 deg", jet_pa);
    check(hole.jet.has_value() && hole.jet->precesses && hole.jet->knots.size() == 12, "M87* has a precessing jet with 12 knots");
    check(hole.hot_flow.has_value(), "M87* has a hot flow");
    check_model_sky(scene, "M87* model sky (its own starlight, tools/sky/make_m87_sky.py)");
    if (hole.jet) {
        check(hole.jet->ra_deg == kM87Ra && hole.jet->dec_deg == kM87Dec, "M87* jet at the hole's place on the sky");
    }
}

// The nozzle's position angle and viewing angle follow [CU23] eqs. 3-9 with Table 1's
// parameters (written out here independently of jet_nozzle()).
void test_m87_precession()
{
    const Scene scene = load_scene_or_die("m87.toml");
    const Body& hole = m87_hole(scene);
    if (!hole.jet) {
        return;
    }
    const JetPrecession& p = hole.jet->precession;
    const SkyAxes sky = sky_axes(kM87Ra, kM87Dec);
    const glm::dvec3 axis = -hole.pole;
    const double eta_p = 288.47 * kDegToRad;
    const double theta = 17.21 * kDegToRad;
    const double phi = 1.25 * kDegToRad;
    double worst_pa = 0.0;
    double worst_view = 0.0;
    double pa_min = 1e9;
    double pa_max = -1e9;
    double view_min = 1e9;
    double view_max = -1e9;
    for (double year = 2006.0; year <= 2022.0; year += 0.05) {
        const double turn = -0.56 * (year - 2014.82);
        const double a = std::cos(theta) * std::cos(turn) * std::sin(phi) + std::sin(theta) * std::cos(phi);
        const double b = std::sin(turn) * std::sin(phi);
        const double jx = a * std::cos(eta_p) - b * std::sin(eta_p);
        const double jy = a * std::sin(eta_p) + b * std::cos(eta_p);
        const double pa = wrap_two_pi(std::atan2(jy, jx));
        const double view = std::asin(std::sqrt(jx * jx + jy * jy));

        const glm::dvec3 nozzle = jet_nozzle(p, axis, sky, year);
        const double pa_model = wrap_two_pi(sky_position_angle(nozzle, sky));
        const double view_model = std::acos(glm::dot(nozzle, sky.toward_observer));
        worst_pa = std::max(worst_pa, std::abs(wrap_pi(pa_model - pa)));
        worst_view = std::max(worst_view, std::abs(view_model - view));
        pa_min = std::min(pa_min, pa);
        pa_max = std::max(pa_max, pa);
        view_min = std::min(view_min, view);
        view_max = std::max(view_max, view);
    }
    check(worst_pa * kRadToDeg < 0.02, "M87* nozzle position angle as Cui et al. 2023 eq. 8", worst_pa * kRadToDeg);
    check(worst_view * kRadToDeg < 0.02, "M87* nozzle viewing angle as Cui et al. 2023 eq. 9", worst_view * kRadToDeg);
    // "a peak-to-peak amplitude of approximately 10 deg centred at PA 288 deg"; "the jet viewing
    // angle oscillates between 16 and 18.5 deg".
    std::printf("info: M87* nozzle: position angle %.2f-%.2f deg, viewing angle %.2f-%.2f deg\n", pa_min * kRadToDeg,
                pa_max * kRadToDeg, view_min * kRadToDeg, view_max * kRadToDeg);
    check(std::abs((pa_max - pa_min) * kRadToDeg - 8.5) < 1.0, "M87* jet position angle swings ~8-10 deg",
          (pa_max - pa_min) * kRadToDeg);
    check(std::abs(0.5 * (pa_max + pa_min) * kRadToDeg - 288.47) < 0.1, "M87* swing centred at 288.47 deg",
          0.5 * (pa_max + pa_min) * kRadToDeg);
    check(std::abs(view_min * kRadToDeg - 16.0) < 0.1 && std::abs(view_max * kRadToDeg - 18.46) < 0.1,
          "M87* viewing angle 16-18.5 deg", view_max * kRadToDeg);
    // Period 2 pi / 0.56 = 11.22 yr (Table 1: 11.24 +- 0.47).
    const double period = kTwoPi / p.omega_rad_per_year;
    const glm::dvec3 now = jet_nozzle(p, axis, sky, 2026.9);
    const glm::dvec3 later = jet_nozzle(p, axis, sky, 2026.9 + period);
    const glm::dvec3 half = jet_nozzle(p, axis, sky, 2026.9 + 0.5 * period);
    check(glm::length(later - now) < 1e-9 && std::abs(std::acos(glm::dot(now, half)) * kRadToDeg - 2.5) < 0.01,
          "M87* nozzle returns after 11.22 years, opposite (2.5 deg away) after half", period);
    check(std::abs(std::acos(glm::dot(now, axis)) * kRadToDeg - 1.25) < 1e-9, "M87* nozzle 1.25 deg from the axis");
}

// The jet's shape and speed against [AN12] and [PA19]; the knots' deprojection; the
// beaming that makes the counter-jet invisible from the Earth.
void test_m87_jet_model()
{
    const Scene scene = load_scene_or_die("m87.toml");
    const Body& hole = m87_hole(scene);
    if (!hole.jet) {
        return;
    }
    const JetParams& p = *hole.jet;
    // [AN12] Fig. 2: r = 1.43 z^0.578 (r_s); r = 0.007 z (conical); M = r_s / 2.
    const auto an12_parabola = [](double z_m) { return 2.0 * 1.43 * std::pow(0.5 * z_m, 0.578); };
    for (const double z : {3.0e3, 1.0e4, 3.0e4}) {
        check(std::abs(jet_radius(p, z) / an12_parabola(z) - 1.0) < 0.12, "M87* jet radius on the parabola [AN12]",
              jet_radius(p, z) / an12_parabola(z));
    }
    check(std::abs(jet_radius(p, 2.0e7) / (0.007 * 2.0e7) - 1.0) < 0.02, "M87* jet conical far out [AN12]",
          jet_radius(p, 2.0e7) / (0.007 * 2.0e7));
    // [PA19]: Gamma = (z / 300 R_S)^0.16.
    const double gb = jet_gamma_beta(p, 2.0e5);
    const double gamma = std::pow(1.0e5 / 300.0, 0.16);
    check(std::abs(std::sqrt(1.0 + gb * gb) / gamma - 1.0) < 1e-9, "M87* jet Gamma at 10^5 R_S [PA19]", gb);
    check(jet_gamma_beta(p, 200.0) == 0.4, "M87* jet slow near the hole [PA19]");
    check(std::abs(jet_gamma_beta(p, 1.5e7) - 1.5) < 0.15, "M87* knots at Gamma beta ~1.5 [PA19]",
          jet_gamma_beta(p, 1.5e7));

    // HST-1 starts 0.607" out in projection: 0.607" x 81.45 pc/" / sin(17.21 deg) deprojected.
    const double m_km = hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    const double expected = 0.607 / 3600.0 * kDegToRad * kM87DistancePc * kParsecKm / std::sin(17.21 * kDegToRad) / m_km;
    check(std::abs(p.knots[0].start_m / expected - 1.0) < 1e-3, "M87* HST-1 deprojected", p.knots[0].start_m);
    check(p.knots[0].line_brightness == 1.0, "M87* knots normalized to HST-1");

    // Just upstream of the first aperture, both jets have the same intrinsic ridge.
    const double theta = 17.21 * kDegToRad;
    const glm::dvec3 to_earth(std::sin(theta), 0.0, std::cos(theta)); // jet frame, the jet along +z
    const double z = 0.99 * (p.knots[0].start_m - 0.3 * (p.knots[0].end_m - p.knots[0].start_m));
    const double ridge = 0.82 * jet_radius(p, z);
    const double app = jet_emissivity(p, glm::dvec3(0.0, ridge, z), to_earth, 0.0);
    const double counter = jet_emissivity(p, glm::dvec3(0.0, ridge, -z), to_earth, 0.0);
    std::printf("info: M87* jet / counter-jet upstream of HST-1 from the Earth: %.3g\n", app / counter);
    check(app / counter > 1e3, "M87* counter-jet thousands of times fainter from the Earth", app / counter);
    // From the side the two match.
    const glm::dvec3 side(1.0, 0.0, 0.0);
    const double app_side = jet_emissivity(p, glm::dvec3(0.0, ridge, z), side, 0.0);
    const double counter_side = jet_emissivity(p, glm::dvec3(0.0, ridge, -z), side, 0.0);
    check(std::abs(app_side / counter_side - 1.0) < 1e-9, "M87* jets alike from the side", app_side / counter_side);
    // Limb brightened: the ridge outshines the axis.
    check(app > 2.0 * jet_emissivity(p, glm::dvec3(0.0, 0.0, z), to_earth, 0.0), "M87* jet limb brightened");

    // The nozzle leans: across the jet at 2000 M, the two ridges straddle z tan(1.25 deg).
    const double h = 2000.0;
    const double w = jet_radius(p, h);
    const double tan_tilt = std::tan(1.25 * kDegToRad);
    const glm::dvec3 up(0.0, 0.0, 1.0);
    double best[2] = {0.0, 0.0};
    double at[2] = {0.0, 0.0};
    for (int i = 0; i <= 4000; ++i) {
        const double x = -1.5 * w + 3.0 * w * i / 4000.0;
        const double e = jet_emissivity(p, glm::dvec3(x, 0.0, h), up, tan_tilt);
        const int k = x < h * tan_tilt ? 0 : 1;
        if (e > best[k]) {
            best[k] = e;
            at[k] = x;
        }
    }
    check(std::abs(0.5 * (at[0] + at[1]) - h * tan_tilt) < 0.01 * w, "M87* nozzle leans by 1.25 deg at 2000 M",
          0.5 * (at[0] + at[1]));
    check(jet_emissivity(p, glm::dvec3(0.0, 0.0, 3.0), up, tan_tilt) == 0.0, "M87* jet starts 5 M above the hole");
}

void test_m87_optical_epochs()
{
    const Scene scene = load_scene_or_die("m87.toml");
    const Body& hole = m87_hole(scene);
    const JetParams& p = *hole.jet;
    const JetKnot& hst = p.knots.at(0);
    const JetKnot& d = p.knots.at(1);
    const JetFeature& middle = d.features.at(1);
    const JetFeature& west = d.features.at(2);
    const double m_km = hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    const double m_per_mas = kM87DistancePc * kParsecKm * kDegToRad / 3600000.0 / m_km;
    const SkyAxes sky = sky_axes(kM87Ra, kM87Dec);
    const double pa = 290.0 * kDegToRad;
    const glm::dvec3 along = std::cos(pa) * sky.north + std::sin(pa) * sky.east;
    const glm::dvec3 across = -std::sin(pa) * sky.north + std::cos(pa) * sky.east;
    const auto position = [&](const JetFeature& f, double year) {
        const auto v = jet_feature_position(p.astrometry, f, year);
        return glm::dvec2(glm::dot(v, along), glm::dot(v, across)) / m_per_mas;
    };
    const auto motion = [&](const JetFeature& f, double year) {
        return position(f, year + 0.05) - position(f, year - 0.05);
    };
    check(std::abs(position(middle, 1995.0).x - 3610.0) < 0.01, "D-Middle projected separation in 1995 [ME13]");
    check(glm::length(10.0 * motion(middle, 2004.06) - glm::dvec2(16.2, -1.3)) < 1e-6,
          "D-Middle along and across the sky in mas/yr [ME13]");
    check(std::abs(10.0 * motion(west, 1996.0).x * 0.264 - 2.73) < 0.02 &&
          std::abs(10.0 * motion(west, 2008.0).x) < 1.0, "D-West slows from 2.73c to a near stop [ME13 Fig. 2]");
    check(std::abs(10.0 * motion(p.knots.at(4).features.at(0), 2000.0).x + 0.9) < 1e-6,
          "knot I retains its signed inward fit [ME13]");
    check(position(middle, 2050.0) == position(middle, 2007.0) && position(west, 1800.0) == position(west, 1995.0) &&
          position(west, 2050.0) == position(west, 2008.25), "astrometric fits held at their individual endpoints");
    check(std::abs(position(hst.features.at(0), 2005.0).x - 820.0) < 0.01, "HST-1 compact flare at 0.82 arcsec [PE03]");

    // Redistributing reference aperture light into compact features conserves the flux seen from Earth.
    const double theta = 17.21 * kDegToRad;
    const auto aperture_flux = [&](const JetKnot& k) {
        return (k.end_m - k.start_m) * k.line_brightness *
               std::pow(jet_doppler(p, 0.5 * (k.start_m + k.end_m), theta), 2.0 + p.spectral_index);
    };
    double d_flux = aperture_flux(d) * (1.0 - d.feature_fraction);
    for (const JetFeature& f : d.features) {
        d_flux += f.flux;
    }
    check(std::abs(d_flux / hst.features[0].flux - 117.0 / 19.8) < 1e-6 &&
          std::abs(middle.flux / west.flux - 28.6 / 23.0) < 1e-6, "D aperture total and relative subknot fluxes [PE01]");
    check(jet_line_brightness(p, 0.5 * (hst.start_m + hst.end_m), true) == 0.0,
          "HST-1 compact peak is counted only once");

    // B5-n appears in 2001 and D-Middle fades after 2007; their light returns to the
    // diffuse aperture, whose total is unchanged [ME13 sec. 3, Table 1 note b].
    const JetFeature& b5n = p.knots.at(9).features.at(1);
    check(b5n.name == "B5-n" && std::abs(position(b5n, 2001.0).x - 15180.0) < 0.01 &&
          std::abs(position(b5n, 1995.0).x - 15180.0) < 0.01, "B5-n starts at 15.18 arcsec in 2001 [ME13]");
    const auto shown = [&](const JetFeature& f, double year) { return jet_feature_visibility(p.astrometry, f, year); };
    check(shown(b5n, 1999.9) == 0.0 && shown(b5n, 2000.5) > 0.0 && shown(b5n, 2000.5) < 1.0 && shown(b5n, 2001.0) == 1.0,
          "B5-n fades in over the year before its first detection");
    check(shown(middle, 2007.0) == 1.0 && shown(middle, 2008.0) == 0.0 && shown(west, 1995.0) == 1.0 &&
          shown(west, 2008.25) == 1.0, "D-Middle fades out after 2007, the other peaks stay");
    for (const double year : {1995.0, 2007.5, 2010.0}) {
        double total = aperture_flux(d) * jet_knot_diffuse_share(p.astrometry, d, year);
        for (const JetFeature& f : d.features) {
            total += shown(f, year) * f.flux;
        }
        check(std::abs(total / hst.features[0].flux - 117.0 / 19.8) < 1e-6, "D aperture flux conserved while peaks fade");
    }

    // The approaching ridge through each aperture's peaks: knot C-1's mean offset, near the
    // centre of [PE01]'s C-1 box (rows 191-220 about the nucleus at row 200, 0.0758667"/px).
    const JetKnot& c1 = p.knots.at(10);
    const double c1_mid = 0.5 * (c1.start_m + c1.end_m);
    const double c1_ridge = jet_ridge_offset(p, c1_mid);
    check(c1.name == "C-1" && std::abs(c1_ridge - 1.4 / 3.0) < 1e-9 && std::abs(c1_ridge - 5.5 * 0.0758667) < 0.1,
          "C-1 ridge at its peaks' mean offset [ME13], near the [PE01] box centre", c1_ridge);
    check(jet_ridge_offset(p, 0.5 * (hst.start_m + hst.end_m)) == 0.0 && jet_ridge_offset(p, 0.0) == 0.0,
          "M87 ridge starts on the axis at HST-1");
    double bound_ok = 1.0;
    for (double z = p.ridge.front().x; z < p.length_m; z *= 1.01) {
        if (std::abs(jet_ridge_offset(p, z)) > jet_ridge_slope_bound(p) * (z - p.ridge.front().x) * (1.0 + 1e-12)) {
            bound_ok = 0.0;
        }
    }
    check(bound_ok == 1.0, "ridge offsets within the shader envelope's bound");
    // The diffuse jet is centred on the displaced ridge; the counter-jet is not displaced.
    const double sky_m_per_arcsec = 1000.0 * m_per_mas;
    const glm::dvec2 across_jet(0.0, sky_m_per_arcsec); // the jet frame's y on the sky
    const double ring = 0.82 * jet_radius(p, c1_mid);
    const double y_c = c1_ridge * sky_m_per_arcsec;
    const glm::dvec3 to_earth(std::sin(theta), 0.0, std::cos(theta));
    const auto at_c1 = [&](double y, double z) {
        return jet_emissivity(p, glm::dvec3(0.0, y, z), to_earth, 0.0, 0.0, {}, across_jet);
    };
    check(std::abs(at_c1(y_c + ring, c1_mid) / at_c1(y_c - ring, c1_mid) - 1.0) < 1e-9 &&
          at_c1(y_c + ring, c1_mid) > 10.0 * at_c1(-ring, c1_mid), "approaching C-1 emission about the shifted ridge");
    check(std::abs(at_c1(ring, -c1_mid) / at_c1(-ring, -c1_mid) - 1.0) < 1e-9, "counter-jet stays on the axis");

    const ObservedLightCurve& curve = scene.observed_light_curve;
    const int body = scene.find("M87*");
    const auto at = [&](const char* date, int knot = 0) {
        double t = 0.0;
        check(parse_utc(date, &t), "valid photometry test date");
        return observed_knot_factor(curve, body, knot, t);
    };
    check(std::abs(at("2005-03-27") * 19.8 - 1404.0) < 0.02 &&
          std::abs(at("2006-05-23") * 19.8 - 592.0) < 0.02 &&
          std::abs(at("2006-11-28") * 19.8 - 862.0) < 0.02, "HST-1 F606W peak, decline and rebrightening [PE11]");
    check(at("1900-01-01") == 1.0 && std::abs(at("2100-01-01") * 19.8 - 35.0) < 0.01 &&
          at("2005-03-27", 1) == 1.0, "HST-1 endpoint fluxes held, unrelated knots steady [PE01, MW21]");
}

// The color of synchrotron light: the spectrum integration reproduces the Planckian
// locus fit, and F_nu ~ nu^-0.9 is bluish white.
void test_power_law_color()
{
    const glm::dvec3 planck = spectrum_linear_srgb(
        [](double lambda, double t) {
            const double x = 1.438776877e7 / (lambda * t);
            return 1.0 / (std::pow(lambda, 5.0) * (std::exp(x) - 1.0));
        },
        5000.0);
    const glm::dvec3 locus = blackbody_linear_srgb(5000.0);
    check(glm::length(planck - locus) < 0.03, "CIE fit integration matches the Planckian locus at 5000 K",
          glm::length(planck - locus));
    const glm::dvec3 jet = power_law_linear_srgb(0.9);
    std::printf("info: synchrotron nu^-0.9 color (linear sRGB) %.3f %.3f %.3f\n", jet.r, jet.g, jet.b);
    check(jet.b == 1.0 && jet.r > 0.6, "nu^-0.9 is bluish white", jet.r);
}

// The 230 GHz flux of the hot flow from the Earth, and its ring.
void test_m87_hot_flow()
{
    const Scene scene = load_scene_or_die("m87.toml");
    const Body& hole = m87_hole(scene);
    if (!hole.hot_flow) {
        return;
    }
    const SkyAxes sky = sky_axes(kM87Ra, kM87Dec);
    const double inclination = std::acos(glm::dot(hole.pole, sky.toward_observer));
    const double theta_g_uas =
        hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS) / (kM87DistancePc * kParsecKm) * kRadToDeg * 3600.0 * 1e6;
    check(std::abs(theta_g_uas - 3.8) < 0.05, "M87* theta_g 3.8 microarcseconds [EHT19VI]", theta_g_uas);
    const RiafImage image =
        riaf_image(hole.spin, hole.gm_km3_s2, *hole.hot_flow, inclination, kM87DistancePc, 96, 15.0, 20.0 / theta_g_uas);
    const double diameter = 2.0 * image.ring_radius_m * theta_g_uas;
    std::printf("info: M87* RIAF at 230 GHz: %.3f Jy, ring %.1f microarcseconds across (EHT: 42 +- 3)\n",
                image.flux_jy, diameter);
    check(std::abs(image.flux_jy / 0.5 - 1.0) < 0.1, "M87* 230 GHz flux 0.5 Jy [EHT19V]", image.flux_jy);
    // Blurred to EHT's 20 microarcseconds the ring comes out ~22 across, against 42 +- 3 [EHT19VI]
    // (unblurred ~33): seen nearly face-on, the plunging gas inside the photon ring fills the
    // middle (limitations.md). Held at the model's own value so that it does not drift unnoticed.
    check(std::abs(diameter / 22.0 - 1.0) < 0.1, "M87* 230 GHz ring diameter (microarcseconds)", diameter);
}

} // namespace

void run_relativity_tests()
{
    test_m87_orientation();
    test_m87_precession();
    test_m87_jet_model();
    test_m87_optical_epochs();
    test_power_law_color();
    test_m87_hot_flow();
    test_kerr_conservation();
    test_kerr_tilted_conservation();
    test_tilted_orbit_matches_equatorial();
    test_kerr_schwarzschild_precession();
    test_post_keplerian();
    test_psr_b1913_scene();
    test_psr_j1141_scene();
    test_pulsar_field_lines();
    test_schwarzschild_star_light();
    test_kerr_null_geodesics();
    test_observer_frame();
    test_sgr_a_scene();
    test_death_star_precession();
    test_sgr_a_ghosts();
    test_gargantua_scene();
    test_mann_orbit();
    test_edmunds_planet();
    test_kerr_lens();
    test_gargantua_lensed_bodies();
    test_bessel_k2();
    test_thermal_synchrotron();
    test_cold_electrons_do_not_emit();
    test_hot_flow_far_fade();
    test_riaf_plasma();
    test_riaf_velocity_grid();
    test_flow_ray_from_inside();
    test_flow_ray_optically_thick();
    test_hot_spot_period();
    test_hot_spot_center();
    test_hot_spot_brightens_ray();
    test_hot_spot_ray_near_spot();
    test_near_horizon_spiral_captured();
    test_sgr_a_hot_flow();
    test_sgr_a_flares();
}
