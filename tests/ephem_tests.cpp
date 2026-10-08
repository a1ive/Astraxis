// Tests for ephem/: Kepler, ephemeris tables, visual and transit orbits, N-body.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
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
    test_ephemeris();
    test_visual_orbit_convention();
    test_transit_orbit_convention();
    test_nbody();
}
