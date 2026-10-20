// Tests for ephem/: Kepler, ephemeris tables, visual and transit orbits, N-body.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/coast.hpp"
#include "ephem/ephemeris.hpp"
#include "ephem/kepler.hpp"
#include "ephem/nbody.hpp"
#include "ephem/visual_orbit.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace astraxis;

namespace {

void test_kepler()
{
    const double eccentricities[] = {0.0, 0.01, 0.3, 0.7, 0.95};
    for (double e : eccentricities) {
        for (int k = -20; k <= 20; ++k) {
            const double m = k * 0.37;
            const double E = solve_kepler(m, e);
            const double residual = wrap_pi(E - e * std::sin(E) - m);
            check(std::abs(residual) < 1e-12, "Kepler equation residual", residual);
        }
    }
}

void test_kepler_propagation()
{
    // Circular orbit (mu = 1, r = 1): a quarter period turns the state by 90 deg.
    const State circle{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    const State q = propagate_kepler(circle, 1.0, kPi / 2.0);
    check(glm::length(q.position - glm::dvec3(0.0, 1.0, 0.0)) < 1e-12 &&
              glm::length(q.velocity - glm::dvec3(-1.0, 0.0, 0.0)) < 1e-12,
          "Kepler quarter circle", glm::length(q.position - glm::dvec3(0.0, 1.0, 0.0)));

    // Eccentric ellipse and hyperbola: energy and angular momentum conserved,
    // and propagating back returns the start.
    for (const double speed : {1.2, 1.6}) { // escape speed is sqrt(2) ~ 1.414
        const State s0{{1.0, 0.0, 0.0}, {0.1, speed, 0.0}};
        const State s1 = propagate_kepler(s0, 1.0, 3.7);
        const State s2 = propagate_kepler(s1, 1.0, -3.7);
        auto energy = [](const State& s) { return 0.5 * glm::dot(s.velocity, s.velocity) - 1.0 / glm::length(s.position); };
        const double de = std::abs(energy(s1) - energy(s0));
        const double dh = glm::length(glm::cross(s1.position, s1.velocity) - glm::cross(s0.position, s0.velocity));
        check(de < 1e-12 && dh < 1e-12, "Kepler propagation conserves energy and momentum", std::max(de, dh));
        check(glm::length(s2.position - s0.position) < 1e-10, "Kepler propagation round trip",
              glm::length(s2.position - s0.position));
    }

    // Eccentric orbit (ARTEMIS-P1 around the Moon in 2012, e = 0.83): plain Newton
    // cycled between two values for some steps (22,000 km off). Compare with the
    // elements propagated by the mean anomaly, over start points around the orbit.
    const double mu = 4902.800118;
    KeplerElements el;
    el.a = 10662.139;
    el.e = 0.8319;
    el.i = 0.3;
    el.arg_peri = 1.1;
    const double n = std::sqrt(mu / (el.a * el.a * el.a));
    double worst = 0.0;
    for (int k = 0; k < 64; ++k) {
        el.mean_anomaly = kTwoPi * k / 64.0;
        const State s0 = kepler_state(el, n);
        for (double dt = -60000.0; dt <= 60000.0; dt += 371.0) {
            KeplerElements later = el;
            later.mean_anomaly += n * dt;
            worst = std::max(worst, glm::length(propagate_kepler(s0, mu, dt).position - kepler_state(later, n).position));
        }
    }
    check(worst < 1e-6, "Kepler propagation on an eccentric orbit (km)", worst);
}

// Coasting from a state (mu = 1), with an optional impulse at t_burn.
class ImpulseMotion final : public MotionSource {
public:
    ImpulseMotion(const State& s0, double t_burn, const glm::dvec3& dv)
        : m_s0(s0)
        , m_t_burn(t_burn)
    {
        m_after = propagate_kepler(s0, 1.0, t_burn);
        m_after.velocity += dv;
    }
    State eval(double t) const override
    {
        return t < m_t_burn ? propagate_kepler(m_s0, 1.0, t) : propagate_kepler(m_after, 1.0, t - m_t_burn);
    }

private:
    State m_s0;
    double m_t_burn;
    State m_after;
};

void test_coast_and_delta_v()
{
    // A coast leaves nothing unexplained; an impulse adds its size, at its step
    // (to 1e-4: within the step, gravity turns the difference a little).
    const State s0{{1.0, 0.0, 0.0}, {0.1, 1.2, 0.0}};
    CoastGravity point;
    point.mu = 1.0;
    const std::vector<double> coasting =
        unexplained_delta_v(ImpulseMotion(s0, 1e9, glm::dvec3(0.0)), point, 0.0, 4.0, 40);
    check(coasting.size() == 41 && coasting.back() < 1e-10, "No delta-v while coasting", coasting.back());
    const glm::dvec3 dv(-0.2, 0.1, 0.05);
    const std::vector<double> burn = unexplained_delta_v(ImpulseMotion(s0, 2.05, dv), point, 0.0, 4.0, 40);
    check(std::abs(burn.back() - glm::length(dv)) < 1e-4 && burn[20] < 1e-10 && burn.back() - burn[21] < 1e-10,
          "The delta-v of an impulse", burn.back() - glm::length(dv));

    // The ghost coasts from the state at its epoch, both ways in time, within its window.
    const CoastMotion ghost(s0, point, 2.0, 1.0, 3.0);
    check(glm::length(ghost.eval(2.0).position - s0.position) < 1e-12 &&
              glm::length(ghost.eval(1.5).position - propagate_kepler(s0, 1.0, -0.5).position) < 1e-12,
          "Coast from a state");
    check(ghost.valid_at(1.0) && ghost.valid_at(3.0) && !ghost.valid_at(0.99) && !ghost.valid_at(3.01),
          "Coast valid within its window");

    // Integrated (a tide of nothing makes it so): the two-body path, both ways from the
    // epoch over 2 pi either way, between the knots too.
    CoastGravity integrated;
    integrated.mu = 1.0;
    integrated.tides.push_back({0.0, [](double) { return glm::dvec3(1e6, 0.0, 0.0); }});
    const CoastMotion integrated_ghost(s0, integrated, 2.0, 2.0 - 2.0 * kPi, 2.0 + 2.0 * kPi);
    double worst = 0.0;
    for (int k = 0; k <= 1000; ++k) {
        const double t = 2.0 + 2.0 * kPi * (k / 500.0 - 1.0) + 0.0017;
        if (integrated_ghost.valid_at(t)) {
            const glm::dvec3 kepler = propagate_kepler(s0, 1.0, t - 2.0).position;
            worst = std::max(worst, glm::length(integrated_ghost.eval(t).position - kepler));
        }
    }
    check(worst < 1e-9, "Integrated coast without perturbations: two-body", worst);
    const State s1 = coast(integrated, s0, 0.0, 3.0);
    check(glm::length(s1.position - propagate_kepler(s0, 1.0, 3.0).position) < 1e-10 &&
              glm::length(coast(integrated, s1, 3.0, -3.0).velocity - s0.velocity) < 1e-10,
          "Integrated coast forward and back");

    // Zonal harmonics against the closed forms: on the pole J2 adds 3 J2 mu R^2 / r^4 and
    // J4 5 J4 mu R^4 / r^6 outward; on the equator J2 pulls 3/2 J2 mu R^2 / r^4 in and J4
    // pushes 15/8 J4 mu R^4 / r^6 out.
    CoastGravity oblate;
    oblate.mu = 2.0;
    oblate.zonal = {0.0, 0.0, 0.01, 0.0, -0.002};
    oblate.radius_km = 1.0;
    oblate.pole = glm::normalize(glm::dvec3(0.3, -0.2, 1.0));
    const double r = 2.0;
    const glm::dvec3 on_pole = oblate.acceleration(r * oblate.pole, 0.0);
    const double pole_expected =
        -2.0 / (r * r) + 3.0 * 0.01 * 2.0 / std::pow(r, 4) + 5.0 * -0.002 * 2.0 / std::pow(r, 6);
    const glm::dvec3 east = glm::normalize(glm::cross(oblate.pole, glm::dvec3(1.0, 0.0, 0.0)));
    const glm::dvec3 on_equator = oblate.acceleration(r * east, 0.0);
    const double equator_expected =
        -2.0 / (r * r) - 1.5 * 0.01 * 2.0 / std::pow(r, 4) + 15.0 / 8.0 * -0.002 * 2.0 / std::pow(r, 6);
    check(glm::length(on_pole - pole_expected * oblate.pole) < 1e-15 &&
              glm::length(on_equator - equator_expected * east) < 1e-15,
          "Zonal harmonics on the pole and the equator");
    // An axisymmetric field keeps the energy and the angular momentum about the pole.
    const State tilted{{1.5, 0.2, 0.4}, {-0.1, 1.0, 0.5}};
    auto energy = [&](const State& s) {
        const double rn = glm::length(s.position);
        const double z = glm::dot(s.position, oblate.pole) / rn;
        const double p2 = 0.5 * (3.0 * z * z - 1.0);
        const double p4 = (35.0 * std::pow(z, 4) - 30.0 * z * z + 3.0) / 8.0;
        return 0.5 * glm::dot(s.velocity, s.velocity) -
               2.0 / rn * (1.0 - 0.01 * p2 / (rn * rn) - -0.002 * p4 / std::pow(rn, 4));
    };
    const State later = coast(oblate, tilted, 0.0, 30.0);
    check(std::abs(energy(later) / energy(tilted) - 1.0) < 1e-9 &&
              std::abs(glm::dot(glm::cross(later.position, later.velocity), oblate.pole) /
                           glm::dot(glm::cross(tilted.position, tilted.velocity), oblate.pole) -
                       1.0) < 1e-9,
          "Coast in an oblate field keeps the energy and the polar angular momentum",
          energy(later) / energy(tilted) - 1.0);
    // A tide pulls on the spacecraft minus on the center: none at the center, and along
    // the line to the tide-raising body twice GM d / D^3 outward to first order.
    CoastGravity tidal;
    tidal.tides.push_back({5.0, [](double) { return glm::dvec3(0.0, 100.0, 0.0); }});
    const glm::dvec3 stretch = tidal.acceleration(glm::dvec3(0.0, 0.01, 0.0), 0.0);
    check(std::abs(stretch.y / (2.0 * 5.0 * 0.01 / 1e6) - 1.0) < 1e-3, "Tide along the line", stretch.y);
}

void test_ephemeris()
{
    EphemerisTable v1;
    std::string error;
    check(v1.load(ASTRAXIS_ASSET_DIR "/ephem/voyager1.eph", &error), "voyager1.eph loads");
    check(v1.target_id() == -31 && v1.center_id() == 10, "voyager1.eph ids");

    // Horizons, Voyager 1 @ Sun, 1977-Sep-07 00:00 TDB (JD 2443393.5), ICRF km.
    const State s = v1.eval(tdb_from_jd_tdb(2443393.5));
    const glm::dvec3 horizons(1.457774659041131E+08, -3.572921412462115E+07, -1.539445114171495E+07);
    check(glm::length(s.position - horizons) < 20.0, "Voyager 1 matches Horizons (km)",
          glm::length(s.position - horizons));

    // Hermite velocity is the derivative of the interpolated position.
    const double t = tdb_from_jd_tdb(2444303.0); // 1979-03-05, near the Jupiter flyby (dense knots)
    const glm::dvec3 numeric = (v1.eval(t + 1.0).position - v1.eval(t - 1.0).position) / 2.0;
    check(glm::length(numeric - v1.eval(t).velocity) < 1e-3, "Hermite velocity consistent (km/s)",
          glm::length(numeric - v1.eval(t).velocity));

    EphemerisTable broken;
    check(!broken.load(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error), "rejects a non-ephemeris file");
}

// Akeson et al. (2021) alpha Cen AB elements reproduce their own ALMA relative
// astrometry (Table 6: 2019.6505, PA 343.7729 deg, separation 5.27869").
void test_visual_orbit_convention()
{
    VisualOrbit ab;
    const double parallax_arcsec = 0.75081;
    ab.a_km = 17.4930 / parallax_arcsec * kAuKm;
    ab.e = 0.51947;
    ab.i_deg = 79.2430;
    ab.node_deg = 205.073;
    ab.arg_peri_deg = 231.519;
    ab.t_peri_tdb = tdb_from_julian_year(1955.564);
    ab.ra_deg = 219.85892215;
    ab.dec_deg = -60.83163195;
    const double gm = (1.0788 + 0.9092) * kSunGmKm3S2;

    const State s = visual_orbit_state(ab, gm, tdb_from_julian_year(2019.6505));
    const double ra = ab.ra_deg * kDegToRad;
    const double dec = ab.dec_deg * kDegToRad;
    const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
    const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
    const double distance_km = kAuKm / parallax_arcsec * 206264.80624709636; // 1 pc / parallax
    const double n = glm::dot(s.position, north) / distance_km * 206264.80624709636;
    const double e = glm::dot(s.position, east) / distance_km * 206264.80624709636;
    const double sep = std::sqrt(n * n + e * e);
    double pa = std::atan2(e, n) * kRadToDeg;
    pa = pa < 0.0 ? pa + 360.0 : pa;
    check(std::abs(sep - 5.27869) < 0.02, "alpha Cen AB separation (arcsec)", sep);
    check(std::abs(pa - 343.7729) < 0.3, "alpha Cen AB position angle (deg)", pa);

    // Kepler's third law with these elements gives the published total mass.
    const double period_s = kTwoPi * std::sqrt(ab.a_km * ab.a_km * ab.a_km / gm);
    check(std::abs(period_s / (kDaysPerJulianYear * kSecondsPerDay) - 79.762) < 0.05, "alpha Cen AB period", period_s);
}

// TransitOrbit: whatever the fit's convention, the planet is in front of the
// star at the given mid-transit time and keeps the given period.
void test_transit_orbit_convention()
{
    TransitOrbit t;
    t.period_days = 10.0;
    t.t_transit_tdb = 1.0e6;
    t.e = 0.2;
    t.arg_peri_deg = 40.0;
    t.i_deg = 88.0;
    t.node_deg = 30.0;
    t.ra_deg = 120.0;
    t.dec_deg = 35.0;
    const double gm = kSunGmKm3S2;
    const double ra = t.ra_deg * kDegToRad;
    const double dec = t.dec_deg * kDegToRad;
    const glm::dvec3 away(std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec));
    for (double u : {90.0, 270.0}) {
        t.transit_u_deg = u;
        const VisualOrbit o = visual_orbit_from_transit(t, gm);
        const glm::dvec3 r = visual_orbit_state(o, gm, t.t_transit_tdb).position;
        const double depth = glm::dot(glm::normalize(r), away);
        check(std::abs(depth + std::sin(t.i_deg * kDegToRad)) < 1e-9, "transit orbit: in front at mid-transit", depth);
        const double period_days = kTwoPi * std::sqrt(o.a_km * o.a_km * o.a_km / gm) / kSecondsPerDay;
        check(std::abs(period_days - t.period_days) < 1e-9, "transit orbit keeps the period", period_days);
    }
    // The two conventions differ only in where omega is measured from.
    t.transit_u_deg = 270.0;
    const VisualOrbit o270 = visual_orbit_from_transit(t, gm);
    t.transit_u_deg = 90.0;
    t.arg_peri_deg -= 180.0;
    const VisualOrbit o90 = visual_orbit_from_transit(t, gm);
    const double d = glm::length(visual_orbit_state(o270, gm, 0.0).position - visual_orbit_state(o90, gm, 0.0).position);
    check(d < 1e-3, "transit orbit conventions agree (km)", d);
}

void test_nbody()
{
    // Two bodies on a Keplerian orbit: one period later they are back.
    const double gm1 = kSunGmKm3S2;
    const double gm2 = 0.5 * kSunGmKm3S2;
    const double a = 10.0 * kAuKm;
    const double e = 0.5;
    const double r_peri = a * (1.0 - e);
    const double v_peri = std::sqrt((gm1 + gm2) * (1.0 + e) / r_peri);
    std::vector<NBodyParticle> p(2);
    p[0].gm = gm1;
    p[1].gm = gm2;
    p[1].position = glm::dvec3(r_peri, 0.0, 0.0);
    p[1].velocity = glm::dvec3(0.0, v_peri, 0.0);
    // Barycentric.
    const glm::dvec3 cm_v = (p[1].velocity * gm2) / (gm1 + gm2);
    const glm::dvec3 cm_r = (p[1].position * gm2) / (gm1 + gm2);
    for (NBodyParticle& q : p) {
        q.position -= cm_r;
        q.velocity -= cm_v;
    }
    const double period = kTwoPi * std::sqrt(a * a * a / (gm1 + gm2));
    const auto knots = integrate_nbody(p, 0.0, -period, period, period / 2000.0);

    auto relative_at = [&](double t) {
        EphemerisTable t1;
        EphemerisTable t2;
        t1.set_knots(knots[0]);
        t2.set_knots(knots[1]);
        return t2.eval(t).position - t1.eval(t).position;
    };
    const double err_fwd = glm::length(relative_at(period) - glm::dvec3(r_peri, 0.0, 0.0)) / a;
    const double err_bwd = glm::length(relative_at(-period) - glm::dvec3(r_peri, 0.0, 0.0)) / a;
    check(err_fwd < 1e-6 && err_bwd < 1e-6, "N-body returns after one period", std::max(err_fwd, err_bwd));

    std::vector<NBodyParticle> end = p;
    end[0].position = knots[0].back().position;
    end[0].velocity = knots[0].back().velocity;
    end[1].position = knots[1].back().position;
    end[1].velocity = knots[1].back().velocity;
    const double de = std::abs(nbody_energy(end) / nbody_energy(p) - 1.0);
    check(de < 1e-9, "N-body energy conserved", de);

    // Thinned to 24 knots per orbit, relative to body 0: Hermite interpolation
    // between knots stays close to the full-resolution path.
    NBodyOptions output;
    output.origin = {0};
    output.stride = {1, 2000 / 24};
    const auto thin = integrate_nbody(p, 0.0, -period, period, period / 2000.0, output);
    EphemerisTable full_table;
    EphemerisTable thin_table;
    std::vector<EphemerisTable::Knot> rel(knots[1].size());
    for (size_t k = 0; k < rel.size(); ++k) {
        rel[k] = {knots[1][k].t, knots[1][k].position - knots[0][k].position,
                  knots[1][k].velocity - knots[0][k].velocity};
    }
    full_table.set_knots(rel);
    thin_table.set_knots(thin[1]);
    check(thin[1].size() < 60, "thinned N-body knots", static_cast<double>(thin[1].size()));
    double worst = 0.0;
    for (int k = 0; k < 997; ++k) {
        const double t = -period + 2.0 * period * (k + 0.5) / 997.0;
        worst = std::max(worst, glm::length(thin_table.eval(t).position - full_table.eval(t).position) / r_peri);
    }
    check(worst < 2e-3, "thinned N-body interpolation error (/ r_peri)", worst);
}

} // namespace

void run_ephem_tests()
{
    test_kepler();
    test_kepler_propagation();
    test_coast_and_delta_v();
    test_ephemeris();
    test_visual_orbit_convention();
    test_transit_orbit_convention();
    test_nbody();
}
