// Tests for relativistic motion: Kerr geodesics, binary pulsars, Sgr A*.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/kerr_null.hpp"
#include "ephem/kerr_orbit.hpp"
#include "ephem/post_keplerian.hpp"
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
    // about the real orbit's normal.
    check(scene.satellite_host(companion) == pulsar && scene.satellite_host(ghost_c) == ghost_p &&
              scene.satellite_host(ghost_p) == -1,
          "PSR B1913+16 hosts (ghosts pair with ghosts)");
    check(scene.bodies[static_cast<size_t>(ghost_p)].body_gm_km3_s2 == 0.0 &&
              std::abs(scene.bodies[static_cast<size_t>(pulsar)].body_gm_km3_s2 / kSunGmKm3S2 - 1.438) < 1e-6,
          "PSR B1913+16 masses (ghosts massless)");
    check(scene.bodies[static_cast<size_t>(pulsar)].pulsar.spin_axis_orbit_of == companion,
          "PSR B1913+16 spin axis along the orbit normal");
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

void test_sgr_a_scene()
{
    Scene scene = load_scene_or_die("sgr_a.toml");
    const int s2 = scene.find("S2");
    check(scene.bodies[0].kind == BodyKind::BlackHole, "Sgr A* is a black hole");
    // Horizon of a non-spinning 4.261e6 M_sun hole: 2 GM / c^2 ~ 1.26e7 km.
    check(std::abs(scene.bodies[0].equatorial_radius_km / 1.2584e7 - 1.0) < 1e-3, "Sgr A* horizon",
          scene.bodies[0].equatorial_radius_km);

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

// The Death Star's periapsides (as marked in the app) advance by the exact
// Schwarzschild amount for its integrated turning points (~18 deg per orbit,
// vs 16.6 deg at 1PN), and the event "at periapsis" really is at one.
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
    const double expected = schwarzschild_advance(r_p, r_a);

    double worst = 0.0;
    for (size_t k = 1; k < peri.size(); ++k) {
        const glm::dvec3 a = scene.icrf_state_at(ds, peri[k - 1]).position;
        const glm::dvec3 b = scene.icrf_state_at(ds, peri[k]).position;
        const double advance = std::acos(std::clamp(glm::dot(a, b) / (glm::length(a) * glm::length(b)), -1.0, 1.0));
        worst = std::max(worst, std::abs(advance / expected - 1.0));
        // Prograde: the periapsis turns in the sense of the orbital motion.
        const glm::dvec3 h = glm::cross(a, scene.icrf_state_at(ds, peri[k - 1]).velocity);
        check(glm::dot(glm::cross(a, b), h) > 0.0, "apsidal advance is prograde");
    }
    check(worst < 0.01, "Death Star apsidal advance = exact Schwarzschild (rel. error)", worst);
    std::printf("info: Death Star r_p %.2f M, r_a %.1f M, advance %.2f deg/orbit (exact), 1PN %.2f\n", r_p, r_a,
                expected / kDegToRad, 6.0 * kPi / (0.5 * (r_p + r_a) * (1.0 - std::pow((r_a - r_p) / (r_a + r_p), 2.0))) / kDegToRad);

    // Event "Death Star at periapsis": a periapsis within the following 3 hours.
    for (const SceneEvent& e : scene.events) {
        if (e.focus == ds && e.name.find("periapsis") != std::string::npos) {
            scene.update(e.t_tdb);
            const std::vector<double> next = scene.periapsis_times(ds, e.t_tdb, e.t_tdb + 3.0 * 3600.0);
            check(next.size() == 1, "periapsis event is 0-3 h before a periapsis",
                  next.empty() ? -1.0 : (next[0] - e.t_tdb) / 3600.0);
        }
        if (e.from_orbit_normal >= 0) {
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

} // namespace

void run_relativity_tests()
{
    test_kerr_conservation();
    test_kerr_schwarzschild_precession();
    test_post_keplerian();
    test_psr_b1913_scene();
    test_psr_j1141_scene();
    test_kerr_null_geodesics();
    test_sgr_a_scene();
    test_death_star_precession();
    test_sgr_a_ghosts();
}
