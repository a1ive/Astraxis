// Tests for the stellar and exoplanet systems.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/ephemeris.hpp"
#include "scene/light_curve.hpp"
#include "scene/star_catalog.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <span>
#include <utility>
#include <vector>

using namespace astraxis;

namespace {

// Mid-transit (minimum sky-projected separation, planet in front) of a planet
// around the root star, searched within +-window of t_guess; NaN if none.
double find_transit(const Scene& scene, int planet, const glm::dvec3& away, double t_guess, double window)
{
    const MotionSource& motion = *scene.bodies[static_cast<size_t>(planet)].motion;
    auto g = [&](double t) {
        const State s = motion.eval(t);
        const glm::dvec3 r = s.position - away * glm::dot(s.position, away);
        const glm::dvec3 v = s.velocity - away * glm::dot(s.velocity, away);
        return glm::dot(r, v);
    };
    const double step = window / 50.0;
    double t0 = t_guess - window;
    double g0 = g(t0);
    for (double t1 = t0 + step; t1 <= t_guess + window; t1 += step) {
        const double g1 = g(t1);
        if (g0 < 0.0 && g1 >= 0.0 && glm::dot(motion.eval(t1).position, away) < 0.0) {
            double lo = t0;
            double hi = t1;
            for (int k = 0; k < 60; ++k) {
                const double mid = 0.5 * (lo + hi);
                (g(mid) < 0.0 ? lo : hi) = mid;
            }
            return 0.5 * (lo + hi);
        }
        t0 = t1;
        g0 = g1;
    }
    return std::nan("");
}

// Mid-conjunction of `body` in front of `center` (minimum sky-projected
// separation), searched within +-window of t_guess; NaN if none.
double find_conjunction(const Scene& scene, int body, int center, const glm::dvec3& away, double t_guess,
                        double window)
{
    auto relative = [&](double t) {
        const State a = scene.icrf_state_at(body, t);
        const State b = scene.icrf_state_at(center, t);
        return State{a.position - b.position, a.velocity - b.velocity};
    };
    auto g = [&](double t) {
        const State s = relative(t);
        const glm::dvec3 r = s.position - away * glm::dot(s.position, away);
        const glm::dvec3 v = s.velocity - away * glm::dot(s.velocity, away);
        return glm::dot(r, v);
    };
    const double step = window / 50.0;
    double t0 = t_guess - window;
    double g0 = g(t0);
    for (double t1 = t0 + step; t1 <= t_guess + window; t1 += step) {
        const double g1 = g(t1);
        if (g0 < 0.0 && g1 >= 0.0 && glm::dot(relative(t1).position, away) < 0.0) {
            double lo = t0;
            double hi = t1;
            for (int k = 0; k < 60; ++k) {
                const double mid = 0.5 * (lo + hi);
                (g(mid) < 0.0 ? lo : hi) = mid;
            }
            return 0.5 * (lo + hi);
        }
        t0 = t1;
        g0 = g1;
    }
    return std::nan("");
}

// A nearby system's sky is the catalog seen from the system: from opposite the
// Sun body, without its own stars (all catalog stars at least kOwnSystemPc
// away). Returns the brightest star's V.
double check_relocated_sky(const Scene& scene, const char* what)
{
    const int sun = scene.find("Sun");
    const glm::dvec3 sun_pc = scene.icrf_state_at(sun, 0.0).position / kParsecKm;
    check(glm::length(scene.sky.viewer_pc + sun_pc) < 1e-6, what, glm::length(scene.sky.viewer_pc + sun_pc));
    double nearest = 1e9;
    double brightest = 99.0;
    for (const CatalogStar& s : scene.sky.stars) {
        nearest = s.distance_pc > 0.0 ? std::min(nearest, s.distance_pc) : nearest;
        brightest = std::min(brightest, s.vmag);
    }
    std::printf("info: %s: %zu stars, nearest %.2f pc, brightest V = %.2f\n", what, scene.sky.stars.size(), nearest,
                brightest);
    check(scene.sky.stars.size() > 8000 && nearest >= kOwnSystemPc, what, nearest);
    return brightest;
}

void test_alpha_centauri_scene()
{
    Scene scene = load_scene_or_die("alpha_centauri.toml");
    const int a = scene.find("Alpha Centauri A");
    const int b = scene.find("Alpha Centauri B");
    const int prox = scene.find("Proxima Centauri");
    const int sun = scene.find("Sun");

    // The integrated (three-body) AB orbit still reproduces the ALMA astrometry.
    const double t_alma = tdb_from_julian_year(2019.6505);
    const glm::dvec3 rel = scene.icrf_state_at(b, t_alma).position - scene.icrf_state_at(a, t_alma).position;
    const double ra = 219.85892215 * kDegToRad;
    const double dec = -60.83163195 * kDegToRad;
    const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
    const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
    const double arcsec_per_km = 0.75081 / kAuKm; // parallax / au
    const double n = glm::dot(rel, north) * arcsec_per_km;
    const double e = glm::dot(rel, east) * arcsec_per_km;
    check(std::abs(std::sqrt(n * n + e * e) - 5.27869) < 0.02, "alpha Cen N-body separation", std::sqrt(n * n + e * e));

    // Proxima ~12,947 au from the AB barycentre at the Hipparcos epoch [Kervella 2017].
    const double d = glm::length(scene.icrf_state_at(prox, tdb_from_julian_year(1991.25)).position) / kAuKm;
    check(std::abs(d - 12947.0) < 50.0, "Proxima distance (au)", d);

    // The Sun is ~1.33 pc away, in Cassiopeia.
    const double d_sun = glm::length(scene.icrf_state_at(sun, 0.0).position) / kParsecKm;
    check(std::abs(d_sun - 1.3319) < 1e-3, "Sun distance from alpha Cen (pc)", d_sun);

    // The sky from here: A, B and Proxima are bodies, not catalog stars; Sirius
    // (2.64 pc from the Sun, 2.92 pc from here) is still the brightest, dimmed.
    const double brightest = check_relocated_sky(scene, "alpha Cen sky");
    check(std::abs(brightest + 1.22) < 0.01, "Sirius from alpha Cen", brightest);

    // Proxima's planets [Suarez Mascareno et al. 2025, Table 3; Damasso et al. 2020,
    // Table 1]: circular orbits whose time of inferior conjunction T0 puts the planet
    // in front of the star, so Proxima's reflex RV is -K sin(2 pi (t - T0) / P)
    // (positive = receding) and the planet's own line-of-sight velocity relative
    // to Proxima is +v sin i sin(2 pi (t - T0) / P).
    struct ProximaPlanet {
        const char* name;
        double t0_bjd;
        double period_days;
        double a_au;
    };
    const ProximaPlanet planets[] = {
        {"Proxima Centauri d", 2460557.55, 5.12338, 0.02881},
        {"Proxima Centauri b", 2460548.59, 11.18465, 0.04848},
        {"Proxima Centauri c (candidate)", 2455892.0, 1900.0, 1.48},
    };
    const glm::dvec3 away = unit_from_ra_dec(217.42894222 * kDegToRad, -62.67949019 * kDegToRad);
    const double sin_i = std::sin(47.0 * kDegToRad);
    scene.update(tdb_from_julian_year(2026.0)); // lighting goes by where the stars are
    for (const ProximaPlanet& p : planets) {
        const int index = scene.find(p.name);
        check(index >= 0 && scene.bodies[static_cast<size_t>(index)].parent == prox, p.name);
        check(scene.lighting_star(index) == prox, "lit by Proxima");
        const MotionSource& orbit = *scene.bodies[static_cast<size_t>(index)].motion;
        const double period = p.period_days * kSecondsPerDay;
        const double t0 = (p.t0_bjd - kJ2000Jd) * kSecondsPerDay;
        const double a_km = p.a_au * kAuKm;
        const double amplitude = kTwoPi / period * a_km * sin_i;
        double worst = 0.0;
        for (int k = 0; k < 12; ++k) {
            const double phase = k / 12.0;
            const double t = t0 + (7.0 + phase) * period;
            const double planet_los = glm::dot(orbit.eval(t).velocity, away);
            worst = std::max(worst, std::abs(planet_los - amplitude * std::sin(kTwoPi * phase)) / amplitude);
        }
        check(worst < 1e-6, "Proxima planet RV curve (relative error)", worst);
        // In front of the star at conjunction.
        const double depth = glm::dot(orbit.eval(t0).position, away);
        check(std::abs(depth / (a_km * sin_i) + 1.0) < 1e-6, "Proxima planet in front at conjunction", depth);
    }

    // A's candidate planet [Beichman, Sanghi et al. 2025]: the chosen orbit from
    // their Zenodo set puts it where Table 2 has S1 (JWST, 2024-08-10) and C1
    // (VLT/NEAR, 2019-06-01), within the measurement errors, in arcsec east and
    // north of A. Its inclination to the AB orbit is the file's imut (70.56 deg,
    // Xuan & Wyatt 2020 eq. 19 with Akeson et al. 2021's elements).
    const int ab = scene.find("Alpha Centauri Ab (candidate)");
    check(ab >= 0 && scene.bodies[static_cast<size_t>(ab)].parent == a, "alpha Cen Ab orbits A");
    check(scene.lighting_star(ab) == a, "alpha Cen Ab lit by A");
    struct Sighting {
        const char* what;
        double jd;
        double east;
        double north;
        double sigma;
    };
    const Sighting sightings[] = {
        {"alpha Cen Ab at C1 (2019)", 2458635.5, -0.64, -0.56, 0.05},
        {"alpha Cen Ab at S1 (2024)", 2460532.5, 1.50, 0.17, 0.13},
    };
    for (const Sighting& s : sightings) {
        const double t = (s.jd - kJ2000Jd) * kSecondsPerDay;
        const glm::dvec3 offset = scene.icrf_state_at(ab, t).position - scene.icrf_state_at(a, t).position;
        const double de = glm::dot(offset, east) * arcsec_per_km - s.east;
        const double dn = glm::dot(offset, north) * arcsec_per_km - s.north;
        check(std::hypot(de, dn) < s.sigma, s.what, std::hypot(de, dn));
    }
    const double t_ab = tdb_from_julian_year(2024.6);
    const State planet = scene.icrf_state_at(ab, t_ab);
    const State sa = scene.icrf_state_at(a, t_ab);
    const State sb = scene.icrf_state_at(b, t_ab);
    const glm::dvec3 n_planet = glm::normalize(glm::cross(planet.position - sa.position, planet.velocity - sa.velocity));
    const glm::dvec3 n_binary = glm::normalize(glm::cross(sb.position - sa.position, sb.velocity - sa.velocity));
    const double mutual = std::acos(glm::dot(n_planet, n_binary)) * kRadToDeg;
    check(std::abs(mutual - 70.56) < 0.1, "alpha Cen Ab mutual inclination (deg)", mutual);

    // B's false positive [Dumusque et al. 2012; Rajpaul et al. 2016]: drawn only
    // during its own event, in front of B at the conjunction from the RVs
    // [Demory et al. 2015], and (coplanar with AB) passing beside the star.
    const int bb = scene.find("Alpha Centauri Bb (false positive)");
    check(bb >= 0 && scene.bodies[static_cast<size_t>(bb)].parent == b, "alpha Cen Bb orbits B");
    int bb_event = -1;
    for (size_t k = 0; k < scene.events.size(); ++k) {
        bb_event = scene.events[k].focus == b ? static_cast<int>(k) : bb_event;
    }
    check(bb_event >= 0 && scene.bodies[static_cast<size_t>(bb)].event == bb_event, "alpha Cen Bb shown by its event");
    scene.set_active_event(-1);
    scene.update(tdb_from_julian_year(2010.2));
    check(!scene.bodies[static_cast<size_t>(bb)].visible, "alpha Cen Bb hidden outside its event");
    scene.set_active_event(bb_event);
    scene.update(tdb_from_julian_year(2010.2));
    check(scene.bodies[static_cast<size_t>(bb)].visible, "alpha Cen Bb shown during its event");
    scene.set_active_event(-1);
    const glm::dvec3 away_ab = unit_from_ra_dec(ra, dec);
    const double t_conj = (2455281.24 - kJ2000Jd) * kSecondsPerDay;
    const double t_found = find_conjunction(scene, bb, b, away_ab, t_conj, 0.5 * kSecondsPerDay);
    check(std::abs(t_found - t_conj) < 60.0, "alpha Cen Bb conjunction (s)", t_found - t_conj);
    const glm::dvec3 rel_bb = scene.icrf_state_at(bb, t_conj).position - scene.icrf_state_at(b, t_conj).position;
    const double miss = glm::length(rel_bb - away_ab * glm::dot(rel_bb, away_ab)) /
                        scene.bodies[static_cast<size_t>(b)].equatorial_radius_km;
    check(glm::dot(rel_bb, away_ab) < 0.0 && std::abs(miss - 1.94) < 0.02, "alpha Cen Bb passes B (stellar radii)",
          miss);
}

// TRAPPIST-1 [Agol et al. 2021]: integrated from the Table 2 elements, the
// planets transit when the paper's posterior says they do (times_obs_and_posterior.txt
// in github.com/ericagol/TRAPPIST1_Spitzer: first, middle and last fitted transit,
// BJD_TDB - 2450000), with the Table 5 impact parameters.
void test_trappist1_scene()
{
    Scene scene = load_scene_or_die("trappist1.toml");
    const glm::dvec3 away = unit_toward(346.622368729, -5.041399251);
    const double r_star = scene.bodies[0].equatorial_radius_km;
    struct Planet {
        const char* name;
        double times[3];
        double impact; // b0 [Table 5]
    };
    const Planet planets[] = {
        {"TRAPPIST-1 b", {7322.517901, 7804.486285, 8769.942058}, 0.095},
        {"TRAPPIST-1 c", {7282.805834, 7796.232104, 8777.058537}, 0.109},
        {"TRAPPIST-1 d", {7560.801874, 7815.940948, 8374.794749}, 0.063},
        {"TRAPPIST-1 e", {7312.713901, 7831.152074, 8770.478717}, 0.191},
        {"TRAPPIST-1 f", {7321.522392, 7827.886782, 8785.389199}, 0.312},
        {"TRAPPIST-1 g", {7294.772224, 7961.825511, 8777.174343}, 0.379},
        {"TRAPPIST-1 h", {7662.550861, 7981.632640, 8769.837235}, 0.378},
    };
    double worst_min = 0.0;
    for (const Planet& p : planets) {
        const int index = scene.find(p.name);
        check(index > 0, p.name);
        if (index <= 0) {
            continue;
        }
        for (double t_post : p.times) {
            const double t_ref = tdb_from_jd_tdb(2450000.0 + t_post);
            const double t = find_transit(scene, index, away, t_ref, 0.1 * kSecondsPerDay);
            const double dt_min = std::isnan(t) ? 1e9 : (t - t_ref) / 60.0;
            worst_min = std::max(worst_min, std::abs(dt_min));
            if (std::abs(dt_min) > 3.0) {
                std::printf("  %s transit at %.6f: off by %.2f min\n", p.name, t_post, dt_min);
            }
        }
        const double t = find_transit(scene, index, away, tdb_from_jd_tdb(2450000.0 + p.times[0]), 0.1 * kSecondsPerDay);
        const glm::dvec3 r = scene.bodies[static_cast<size_t>(index)].motion->eval(t).position;
        const double impact = glm::length(r - away * glm::dot(r, away)) / r_star;
        check(std::abs(impact - p.impact) < 0.02, "TRAPPIST-1 impact parameter", impact);
    }
    std::printf("info: TRAPPIST-1 transit times vs Agol 2021 posterior: worst %.2f min\n", worst_min);
    check(worst_min < 3.0, "TRAPPIST-1 transit times match Agol 2021 (min)", worst_min);

    // The output is thinned to ~24 knots per orbit (span 2000-2050).
    const auto& b = static_cast<const EphemerisMotion&>(*scene.bodies[static_cast<size_t>(scene.find("TRAPPIST-1 b"))].motion);
    const double per_orbit = 1.510826 * static_cast<double>(b.table().knots().size()) / (50.0 * kDaysPerJulianYear);
    check(per_orbit > 20.0 && per_orbit < 30.0, "TRAPPIST-1 b knots per orbit", per_orbit);

    // The sky from 12.5 pc: Sirius is now far away; Canopus (95 pc) is the brightest.
    const double brightest = check_relocated_sky(scene, "TRAPPIST-1 sky");
    check(std::abs(brightest + 0.57) < 0.01, "Canopus from TRAPPIST-1", brightest);
}

// Kepler-223 [Mills et al. 2016]: integrated from the best-fit initial
// conditions (Extended Data Table 3), the planets reproduce the transit-timing
// variations measured in each Kepler quarter (Extended Data Table 1: mean time
// of the quarter's transits, then the TTV relative to the given linear
// ephemeris with its -3, -1, +1, +3 sigma offsets; days, BJD - 2454900).
void test_kepler223_scene()
{
    Scene scene = load_scene_or_die("kepler223.toml");
    check_model_sky(scene, "Kepler-223 model sky");
    const glm::dvec3 away = unit_toward(298.318417682, 47.279530121);
    struct Quarter {
        double t_mean, m3, m1, ttv, p1, p3;
    };
    struct Planet {
        const char* name;
        double period_days;
        double t0;
        Quarter quarters[11];
    };
    const Planet planets[] = {
        {"Kepler-223 b", 7.3840154, 70.49489,
         {{123.32662, -0.0354, -0.0058, -0.0006, 0.0059, 0.0316},
          {239.02516, -0.0517, -0.0103, 0.0137, 0.0101, 0.0423},
          {416.51724, -0.0200, -0.0061, -0.0010, 0.0062, 0.0210},
          {521.69775, -0.0628, -0.0123, 0.0068, 0.0113, 0.0342},
          {699.18988, -0.0470, -0.0088, -0.0010, 0.0084, 0.0370},
          {797.79657, -0.0417, -0.0097, -0.0123, 0.0088, 0.0243},
          {886.54260, -0.0343, -0.0072, -0.0187, 0.0074, 0.0617},
          {1073.89526, -0.0500, -0.0118, 0.0730, 0.0140, 0.1150},
          {1162.64136, -0.0542, -0.0071, 0.0692, 0.0071, 0.0708},
          {1251.38745, -0.0217, -0.0062, 0.0507, 0.0065, 0.0333},
          {1458.46155, -0.0379, -0.0129, 0.0659, 0.0090, 0.0241}}},
        {"Kepler-223 c", 9.8487130, 71.37624,
         {{116.10564, -0.0362, -0.0103, -0.0168, 0.0133, 0.0518},
          {242.86336, -0.0683, -0.0405, 0.0023, 0.0077, 0.0747},
          {420.32413, -0.0254, -0.0077, 0.0264, 0.0076, 0.0476},
          {509.05453, -0.0266, -0.0077, 0.0166, 0.0090, 0.0304},
          {701.30371, -0.0585, -0.0121, -0.0155, 0.0126, 0.0535},
          {790.03418, -0.0302, -0.0075, -0.0318, 0.0084, 0.0268},
          {886.15869, -0.0173, -0.0046, -0.0737, 0.0048, 0.0537},
          {1071.01367, -0.1404, -0.0078, -0.0766, 0.0066, 0.0226},
          {1148.65283, -0.0411, -0.0067, -0.0959, 0.0064, 0.0349},
          {1252.17163, -0.0392, -0.0067, -0.0418, 0.0068, 0.0548},
          {1470.30054, -0.0361, -0.0056, -0.0449, 0.0051, 0.0179}}},
        {"Kepler-223 d", 14.7883997, 109.76775,
         {{132.10997, -0.0416, -0.0058, 0.0376, 0.0054, 0.0134},
          {248.65308, -0.0221, -0.0062, -0.0169, 0.0063, 0.0229},
          {427.57138, -0.0351, -0.0084, -0.0099, 0.0070, 0.0169},
          {519.49268, -0.0285, -0.0070, 0.0035, 0.0066, 0.0245},
          {711.54260, -0.0260, -0.0086, 0.0240, 0.0094, 0.0420},
          {800.18097, -0.0226, -0.0060, 0.0256, 0.0057, 0.0194},
          {898.66815, -0.0192, -0.0057, 0.0212, 0.0055, 0.0238},
          {1077.35193, -0.0530, -0.0080, 0.0020, 0.0077, 0.0210},
          {1169.50781, -0.0354, -0.0085, -0.0236, 0.0093, 0.0286},
          {1271.27771, -0.0272, -0.0131, -0.0578, 0.0132, 0.0328},
          {1483.43542, -0.0298, -0.0061, -0.1612, 0.0057, 0.0302}}},
        {"Kepler-223 e", 19.7213435, 68.10686,
         {{135.47421, -0.0303, -0.0060, -0.0067, 0.0053, 0.0187},
          {238.21753, -0.0232, -0.0067, 0.0022, 0.0072, 0.0458},
          {433.78842, -0.0542, -0.0095, 0.0302, 0.0084, 0.0298},
          {524.27625, -0.0244, -0.0061, -0.0106, 0.0063, 0.0296},
          {709.21222, -0.0432, -0.0071, 0.0022, 0.0065, 0.0208},
          {797.82037, -0.0240, -0.0060, 0.0090, 0.0061, 0.0240},
          {893.81256, -0.0357, -0.0216, 0.0297, 0.0242, 0.0513},
          {1079.88989, -0.0662, -0.0083, 0.0602, 0.0078, 0.0308},
          {1170.71301, -0.1067, -0.0118, 0.1167, 0.0110, 0.0453},
          {1263.01343, -0.0252, -0.0049, 0.1352, 0.0049, 0.0188},
          {1469.48169, -0.0393, -0.0097, 0.2283, 0.0100, 0.0467}}},
    };
    double residual_sq = 0.0;
    double signal_sq = 0.0;
    int within_3sigma = 0;
    int count = 0;
    for (const Planet& p : planets) {
        const int index = scene.find(p.name);
        check(index > 0, p.name);
        if (index <= 0) {
            continue;
        }
        for (const Quarter& q : p.quarters) {
            // Mean modelled TTV of the transits within half a quarter of the mean time.
            double sum = 0.0;
            int n_transits = 0;
            const int first = static_cast<int>(std::ceil((q.t_mean - 45.0 - p.t0) / p.period_days));
            const int last = static_cast<int>(std::floor((q.t_mean + 45.0 - p.t0) / p.period_days));
            for (int n = first; n <= last; ++n) {
                const double linear = 2454900.0 + p.t0 + n * p.period_days;
                const double t = find_transit(scene, index, away, tdb_from_jd_tdb(linear), 0.5 * kSecondsPerDay);
                if (!std::isnan(t)) {
                    sum += t / kSecondsPerDay + kJ2000Jd - linear;
                    ++n_transits;
                }
            }
            check(n_transits > 0, "Kepler-223 transits found in a quarter", q.t_mean);
            if (n_transits == 0) {
                continue;
            }
            const double model = sum / n_transits;
            ++count;
            residual_sq += (model - q.ttv) * (model - q.ttv);
            signal_sq += q.ttv * q.ttv;
            if (model >= q.ttv + q.m3 && model <= q.ttv + q.p3) {
                ++within_3sigma;
            } else {
                std::printf("  %s quarter at %.1f: model TTV %+.4f d, measured %+.4f d\n", p.name, q.t_mean, model,
                            q.ttv);
            }
        }
    }
    // The binned TTVs scatter around the photodynamical model by 1-3 sigma
    // ([MI16] Fig. 1), but the model must explain the signal (up to 0.23 d).
    const double rms_residual = std::sqrt(residual_sq / std::max(count, 1));
    const double rms_signal = std::sqrt(signal_sq / std::max(count, 1));
    std::printf("info: Kepler-223 quarterly TTVs: %d of %d within 3 sigma, rms %.4f d (signal %.4f d)\n",
                within_3sigma, count, rms_residual, rms_signal);
    check(count == 44 && within_3sigma == count, "Kepler-223 TTVs within 3 sigma", within_3sigma);
    check(rms_residual < 0.25 * rms_signal, "Kepler-223 TTV residuals small against the signal", rms_residual);
}

// Kepler-47 [Orosz et al. 2019]: integrated from the Table 7 barycentric
// initial conditions (BJD - 2455000 below). In the Kepler years the planets
// transit the primary when the paper's model does (Table 2, "model time")
// to within the light-travel time across the orbits (2-8 min; ELC corrects
// for it, we do not). The model's predictions for 2019-2026 (Table 10) fall
// within a quarter of the predicted transit duration: c drifts ahead of ELC by
// about 0.2 primary radii on the sky (cause unknown; 1-sigma 0.04-0.26 d), b
// and d stay within minutes. The binary eclipses drift 0.6 s per orbit from
// the linear ephemeris of Eq. (1), most likely the GR terms ELC includes
// (Sect. 5.1) and we leave out.
void test_kepler47_scene()
{
    Scene scene = load_scene_or_die("kepler47.toml");
    check_model_sky(scene, "Kepler-47 model sky");
    const glm::dvec3 away = unit_toward(295.297909668, 46.920474260);
    const int a = scene.find("Kepler-47 A");
    const int b = scene.find("Kepler-47 B");
    check(a > 0 && b > 0, "Kepler-47 stars");
    if (a <= 0 || b <= 0) {
        return;
    }
    struct Transit {
        const char* name;
        double time;
        double duration_h; // 0: Kepler-era model time
    };
    const Transit transits[] = {
        {"Kepler-47 b", -30.80876, 0.0},  {"Kepler-47 b", 352.25805, 0.0},  {"Kepler-47 b", 1408.96429, 0.0},
        {"Kepler-47 d", 604.44003, 0.0},  {"Kepler-47 d", 1350.36740, 0.0}, {"Kepler-47 c", 246.65033, 0.0},
        {"Kepler-47 c", 550.47748, 0.0},  {"Kepler-47 c", 850.98582, 0.0},  {"Kepler-47 c", 1154.78852, 0.0},
        {"Kepler-47 b", 3522.04361, 4.86}, {"Kepler-47 b", 4961.81627, 3.34}, {"Kepler-47 b", 6402.38804, 10.99},
        {"Kepler-47 d", 3590.46599, 5.02}, {"Kepler-47 d", 5456.14904, 34.28}, {"Kepler-47 d", 6390.01925, 6.46},
        {"Kepler-47 c", 3575.38963, 5.78}, {"Kepler-47 c", 4483.55306, 6.61}, {"Kepler-47 c", 5996.00286, 6.10},
    };
    double worst_kepler_min = 0.0;
    double worst_fraction = 0.0;
    for (const Transit& p : transits) {
        const int index = scene.find(p.name);
        check(index > 0, p.name);
        if (index <= 0) {
            continue;
        }
        const double t_ref = tdb_from_jd_tdb(2455000.0 + p.time);
        const double t = find_conjunction(scene, index, a, away, t_ref, 0.5 * kSecondsPerDay);
        const double dt_min = std::isnan(t) ? 1e9 : (t - t_ref) / 60.0;
        if (p.duration_h == 0.0) {
            worst_kepler_min = std::max(worst_kepler_min, std::abs(dt_min));
        } else {
            worst_fraction = std::max(worst_fraction, std::abs(dt_min) / (60.0 * p.duration_h));
        }
    }
    std::printf("info: Kepler-47 transits vs Orosz 2019: Kepler era worst %.2f min, 2019-2026 worst %.2f durations\n",
                worst_kepler_min, worst_fraction);
    check(worst_kepler_min < 8.0, "Kepler-47 transits match Orosz 2019 Table 2 (min)", worst_kepler_min);
    check(worst_fraction < 0.25, "Kepler-47 transits within the Table 10 predictions (durations)", worst_fraction);

    // Primary eclipses: B in front of A, P = 7.44837568 d, T0 = BJD 2454963.246137.
    double worst_eclipse_min = 0.0;
    for (int n : {1, 190, 400, 816}) {
        const double t_ref = tdb_from_jd_tdb(2454963.246137 + 7.44837568 * n);
        const double t = find_conjunction(scene, b, a, away, t_ref, 0.1 * kSecondsPerDay);
        const double dt_min = std::isnan(t) ? 1e9 : (t - t_ref) / 60.0;
        worst_eclipse_min = std::max(worst_eclipse_min, std::abs(dt_min));
    }
    std::printf("info: Kepler-47 primary eclipses vs linear ephemeris 2009-2025: worst %.2f min\n",
                worst_eclipse_min);
    check(worst_eclipse_min < 10.0, "Kepler-47 primary eclipses near the linear ephemeris (min)", worst_eclipse_min);

    // Both stars light the planets, the red dwarf with ~2% of the light; the Sun not at all.
    scene.update(tdb_from_jd_tdb(2461041.5));
    StarLight lights[kMaxStarLights];
    const int count = scene.lighting_stars(scene.find("Kepler-47 c"), lights, kMaxStarLights);
    check(count == 2 && lights[0].star == a && lights[1].star == b, "Kepler-47 c lit by A and B", count);
    if (count == 2) {
        check(lights[1].relative_flux > 0.01 && lights[1].relative_flux < 0.03, "Kepler-47 B adds ~2% light",
              lights[1].relative_flux);
    }
}

// TIC 168789840 [Powell et al. 2021]: nested Keplerian orbits. Each binary's
// secondary passes in front of its primary at the Table 2 primary eclipses,
// close enough on the sky to eclipse; B stands at the Table 4 speckle position
// relative to AC; A and C keep the Table 8 period and eccentricity.
void test_tic168789840_scene()
{
    Scene scene = load_scene_or_die("tic168789840.toml");
    check_model_sky(scene, "TIC 168789840 model sky");
    const glm::dvec3 away = unit_toward(63.5201633, -31.9228995);
    struct Binary {
        const char* primary;
        const char* secondary;
        double period_days;
        double t0_bjd;
    };
    const Binary binaries[] = {
        {"TIC 168789840 A1", "TIC 168789840 A2", 1.570013, 2458412.3855},
        {"TIC 168789840 B1", "TIC 168789840 B2", 8.217111, 2458411.9008},
        {"TIC 168789840 C1", "TIC 168789840 C2", 1.305883, 2458413.6822},
    };
    for (const Binary& bin : binaries) {
        const int p = scene.find(bin.primary);
        const int s = scene.find(bin.secondary);
        check(p > 0 && s > 0, bin.secondary);
        if (p <= 0 || s <= 0) {
            continue;
        }
        for (int n : {0, 1000}) {
            const double t_ref = tdb_from_jd_tdb(bin.t0_bjd + bin.period_days * n);
            const double t = find_conjunction(scene, s, p, away, t_ref, 0.2 * kSecondsPerDay);
            const double dt_min = std::isnan(t) ? 1e9 : (t - t_ref) / 60.0;
            check(std::abs(dt_min) < 0.1, "TIC 168789840 primary eclipse on the ephemeris (min)", dt_min);
            if (!std::isnan(t)) {
                const glm::dvec3 r = scene.icrf_state_at(s, t).position - scene.icrf_state_at(p, t).position;
                const double sky = glm::length(r - away * glm::dot(r, away));
                const double radii = scene.bodies[static_cast<size_t>(p)].equatorial_radius_km +
                                     scene.bodies[static_cast<size_t>(s)].equatorial_radius_km;
                check(sky < radii, "TIC 168789840 binaries eclipse", sky / radii);
            }
        }
    }

    // B relative to AC at the speckle epoch: PA 257.74 deg, 0.4230" (at 584 pc).
    const int ac = scene.find("TIC 168789840 AC");
    const int b = scene.find("TIC 168789840 B");
    check(ac > 0 && b > 0, "TIC 168789840 AC and B");
    if (ac > 0 && b > 0) {
        const double t = tdb_from_julian_year(2020.8236);
        const glm::dvec3 r = scene.icrf_state_at(b, t).position - scene.icrf_state_at(ac, t).position;
        const double ra = 63.5201633 * kDegToRad;
        const double dec = -31.9228995 * kDegToRad;
        const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
        const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
        const double pa = wrap_two_pi(std::atan2(glm::dot(r, east), glm::dot(r, north))) * kRadToDeg;
        const double sep = std::hypot(glm::dot(r, east), glm::dot(r, north)) / (584.0 * kAuKm) * 1000.0;
        check(std::abs(pa - 257.74) < 0.05, "TIC 168789840 B position angle (deg)", pa);
        check(std::abs(sep - 423.0) < 0.5, "TIC 168789840 B separation (mas)", sep);
    }

    // A and C: closest at periastron (BJD 2457662), a (1 - e) apart, every 3.7 years.
    const int a = scene.find("TIC 168789840 A");
    const int c = scene.find("TIC 168789840 C");
    check(a > 0 && c > 0, "TIC 168789840 A and C");
    if (a > 0 && c > 0) {
        auto separation = [&](double jd) {
            const double t = tdb_from_jd_tdb(jd);
            return glm::length(scene.icrf_state_at(c, t).position - scene.icrf_state_at(a, t).position) / kAuKm;
        };
        const double a_au = 3.67647;
        for (double jd : {2457662.0, 2457662.0 + 1351.425 * 3}) {
            check(std::abs(separation(jd) - a_au * (1.0 - 0.28)) < 1e-3, "TIC 168789840 AC periastron (au)",
                  separation(jd));
        }
        check(std::abs(separation(2457662.0 + 0.5 * 1351.425) - a_au * 1.28) < 1e-3, "TIC 168789840 AC apastron (au)",
              separation(2457662.0 + 0.5 * 1351.425));
    }
}

// The light of overlapping disks (scene/light_curve.hpp) against closed forms:
// a lone limb-darkened disk gives pi R^2 (1 - u / 3); for uniform disks the
// hidden part is the lens where they overlap; a disk wholly behind a larger
// one gives nothing, and one in front of it hides its own area.
void test_light_curve_model()
{
    const double planck = planck_shape(5772.0, 550.0);
    const StarDisk one[] = {{glm::dvec3(0.0), 1000.0, 5772.0}};
    const double lone = disks_full_flux(one, 550.0, 0.6);
    check(std::abs(lone / (kPi * 1000.0 * 1000.0 * (1.0 - 0.6 / 3.0) * planck) - 1.0) < 1e-4,
          "light curve: limb-darkened disk", lone);

    // Equal uniform disks of radius 1 at distance d: the lens is 2 acos(d/2) - (d/2) sqrt(4 - d^2).
    const glm::dvec3 toward(0.0, 0.0, 1.0);
    for (double d : {0.3, 1.0, 1.7}) {
        const StarDisk pair[] = {{glm::dvec3(0.0), 1.0, 6000.0}, {glm::dvec3(d, 0.0, 5.0), 1.0, 6000.0}};
        const double lens = 2.0 * std::acos(0.5 * d) - 0.5 * d * std::sqrt(4.0 - d * d);
        // Of the back disk (the pair's full light is twice one disk's).
        const double hidden = 2.0 - 2.0 * disks_flux(pair, toward, 550.0, 0.0) / disks_full_flux(pair, 550.0, 0.0);
        check(std::abs(hidden - lens / kPi) < 2e-3, "light curve: lens of two uniform disks", hidden - lens / kPi);
    }
    const StarDisk behind[] = {{glm::dvec3(0.0), 3.0, 4000.0}, {glm::dvec3(0.5, 0.2, -10.0), 1.0, 9000.0}};
    const double big = disks_full_flux(std::span(behind, 1), 550.0, 0.6);
    check(std::abs(disks_flux(behind, toward, 550.0, 0.6) / big - 1.0) < 1e-9, "light curve: disk wholly behind");
}

// A sphere of radius R as a RocheSurface (the same latitude-longitude mesh
// about the x axis), uniformly bright.
RocheSurface faceted_sphere(double radius_km, double temperature_k)
{
    constexpr int kRingCount = 64;
    constexpr int kSegments = 128;
    RocheSurface s;
    auto vertex = [&](int k, int m) {
        const double theta = kPi * k / kRingCount;
        const double phi = kTwoPi * m / kSegments;
        return radius_km *
               glm::dvec3(std::cos(theta), std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi));
    };
    auto add = [&](const glm::dvec3& p0, const glm::dvec3& p1, const glm::dvec3& p2) {
        const glm::dvec3 c = glm::cross(p1 - p0, p2 - p0);
        if (glm::length(c) <= 0.0) {
            return;
        }
        RocheSurface::Facet f;
        f.center = (p0 + p1 + p2) / 3.0;
        f.normal = glm::normalize(c);
        f.area_km2 = 0.5 * glm::length(c);
        f.temperature_k = temperature_k;
        f.radiance = planck_shape(temperature_k, kRocheRadianceWavelengthNm);
        s.facets.push_back(f);
    };
    for (int k = 0; k < kRingCount; ++k) {
        for (int m = 0; m < kSegments; ++m) {
            add(vertex(k, m), vertex(k + 1, m), vertex(k + 1, m + 1));
            add(vertex(k, m), vertex(k + 1, m + 1), vertex(k, m + 1));
            s.mesh.positions.push_back(glm::vec3(vertex(k, m)));
        }
    }
    s.mesh.positions.push_back(glm::vec3(vertex(kRingCount, 0)));
    s.l1_km = s.back_km = s.side_km = s.pole_km = s.volume_radius_km = radius_km;
    return s;
}

// A Roche-lobe-filling star (scene/roche.hpp): the lobe of Kim's (1989)
// Wilson-Devinney solution of Algol (Table 3: q = 0.227, B's radii at the pole,
// side and back 0.241, 0.251, 0.284 of the separation); its volume against
// Eggleton's (1983) fit, good to 1%. The light of such a surface (the facets
// facing the observer) and the outline it hides, on a faceted sphere: as the
// limb-darkened disk, alone, behind a disk and in front of one.
void test_roche_surface()
{
    const double q = 1.0 / 0.227;
    const double pole = roche_radius(q, glm::dvec3(0.0, 0.0, 1.0));
    const double side = roche_radius(q, glm::dvec3(0.0, 1.0, 0.0));
    const double back = roche_radius(q, glm::dvec3(-1.0, 0.0, 0.0));
    std::printf("info: Roche lobe for q = 0.227: pole %.4f, side %.4f, back %.4f, L1 %.4f\n", pole, side, back,
                roche_l1(q));
    check(std::abs(pole - 0.241) < 0.0015 && std::abs(side - 0.251) < 0.0015 && std::abs(back - 0.284) < 0.0015,
          "Roche lobe radii of Kim's Algol B", pole);

    RocheStar star;
    star.separation_km = 1.0e7;
    star.temperature_k = 5000.0;
    for (double mass_ratio : {1.0 / 0.221, 1.0, 0.3}) {
        star.mass_ratio = mass_ratio;
        const RocheSurface s = make_roche_surface(star);
        const double r = 1.0 / mass_ratio;
        const double eggleton = 0.49 * std::cbrt(r * r) / (0.6 * std::cbrt(r * r) + std::log(1.0 + std::cbrt(r)));
        check(std::abs(s.volume_radius_km / star.separation_km / eggleton - 1.0) < 0.01,
              "Roche lobe volume against Eggleton's fit", s.volume_radius_km / star.separation_km - eggleton);
        check(std::abs(s.l1_km / star.separation_km - roche_l1(mass_ratio)) < 1e-12 && s.l1_km > s.back_km &&
                  s.back_km > s.side_km && s.side_km > s.pole_km,
              "Roche lobe longest toward L1, flattened at the poles");
    }

    // The facet sum and the outline, on a sphere made of facets like a lobe's.
    const RocheSurface round = faceted_sphere(1.0e6, 5000.0);
    StarDisk lobe;
    lobe.position = glm::dvec3(0.0, 0.0, 10.0 * round.l1_km);
    lobe.radius_km = round.l1_km;
    lobe.temperature_k = star.temperature_k;
    lobe.roche = &round;
    const StarDisk sphere{lobe.position, round.volume_radius_km, star.temperature_k};
    const glm::dvec3 toward(0.0, 0.0, 1.0);
    const double alone = disks_flux(std::span(&lobe, 1), toward, 550.0, 0.6);
    const double disk = disks_flux(std::span(&sphere, 1), toward, 550.0, 0.6);
    check(std::abs(alone / disk - 1.0) < 3e-3, "Roche surface's light (faceted sphere)", alone / disk - 1.0);
    const double behind_r = 1.5 * round.volume_radius_km;
    for (double offset : {0.5, 1.0, 1.5}) {
        const StarDisk back_disk{glm::dvec3(offset * behind_r, 0.0, 0.0), behind_r, 8000.0};
        const double back_full = disks_full_flux(std::span(&back_disk, 1), 550.0, 0.6);
        const StarDisk with_lobe[] = {back_disk, lobe};
        const StarDisk with_sphere[] = {back_disk, sphere};
        const double hidden_by_lobe = back_full - (disks_flux(with_lobe, toward, 550.0, 0.6) - alone);
        const double hidden_by_sphere = back_full - (disks_flux(with_sphere, toward, 550.0, 0.6) - disk);
        check(std::abs(hidden_by_lobe / hidden_by_sphere - 1.0) < 3e-3, "Roche outline in front of a disk",
              hidden_by_lobe / hidden_by_sphere - 1.0);
        const StarDisk front_disk{glm::dvec3(offset * behind_r, 0.0, 20.0 * round.l1_km), behind_r, 8000.0};
        const StarDisk lobe_behind[] = {front_disk, lobe};
        const StarDisk sphere_behind[] = {front_disk, sphere};
        const double lobe_left = disks_flux(lobe_behind, toward, 550.0, 0.6);
        const double sphere_left = disks_flux(sphere_behind, toward, 550.0, 0.6);
        check(std::abs(lobe_left / sphere_left - 1.0) < 3e-3, "Roche surface behind a disk",
              lobe_left / sphere_left - 1.0);
    }
}

// Algol [Baron et al. 2012, Kim 1989, Kolbas et al. 2015]: B and C stand where
// CHARA saw them relative to A (Table 3, 2006-2010); B passes in front of A at
// the primary minima TESS measured, once the light-time effect of the outer
// orbit is added; B fills its Roche lobe, as elongated as CHARA saw it; the
// eclipse is as deep and long as the GCVS has it (V = 2.12 .. 3.39, 14% of the
// period); A gives most of the visible light [KO15 Sect. 5.1].
void test_algol_scene()
{
    test_light_curve_model();
    test_roche_surface();
    Scene scene = load_scene_or_die("algol.toml");
    const double brightest = check_relocated_sky(scene, "Algol sky from the system");
    check(brightest > -2.0, "Algol sky without the catalog's Algol (V)", brightest); // its copy would be -4.6
    const int root = scene.find("Algol");
    const int ab = scene.find("Algol AB");
    const int a = scene.find("Algol A");
    const int b = scene.find("Algol B");
    const int c = scene.find("Algol C");
    check(root == 0 && ab > 0 && a > 0 && b > 0 && c > 0, "Algol bodies");
    check(scene.light_curve.enabled() && scene.light_curve.stars.size() == 3, "Algol light curve of A, B and C");
    if (ab <= 0 || a <= 0 || b <= 0 || c <= 0) {
        return;
    }
    const double distance_km = 28.8184 * kParsecKm;
    const double ra = 47.0422185 * kDegToRad;
    const double dec = 40.9556467 * kDegToRad;
    const glm::dvec3 away = unit_toward(47.0422185, 40.9556467);
    const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
    const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
    auto offset_mas = [&](int body, double t) { // (east, north) from A on the sky
        const glm::dvec3 r = scene.icrf_state_at(body, t).position - scene.icrf_state_at(a, t).position;
        return glm::dvec2(glm::dot(r, east), glm::dot(r, north)) / distance_km * (kRadToDeg * 3600.0e3);
    };

    // Kepler's third law for the masses gives [BA12]'s semi-major axes.
    const double a_inner =
        glm::length(scene.icrf_state_at(b, 0.0).position - scene.icrf_state_at(a, 0.0).position) / distance_km;
    check(std::abs(a_inner * kRadToDeg * 3600.0e3 - 2.15) < 0.01, "Algol inner a (mas)",
          a_inner * kRadToDeg * 3600.0e3);
    {
        OrbitElements el;
        const State ra_c = scene.icrf_state_at(c, 0.0);
        const State ra_ab = scene.icrf_state_at(ab, 0.0);
        osculating_elements(State{ra_c.position - ra_ab.position, ra_c.velocity - ra_ab.velocity}, 7.471710374e11,
                            &el);
        const double a_mas = 0.5 * (el.periapsis_km + el.apoapsis_km) / distance_km * kRadToDeg * 3600.0e3;
        check(std::abs(a_mas - 93.43) < 0.2, "Algol outer a (mas)", a_mas);
    }

    // [BA12] Table 3: MJD, rho (mas), theta (deg east of north) and the major
    // semi-axis of the 1-sigma error ellipse (mas). C moves slowly: within 2
    // sigma (2007-11-23 is off by 1.7 sigma, the rest by less than 1). B's
    // phase comes from the TESS ephemeris carried back 10-13 years, across
    // Algol's wandering O-C [JE21]: 0.1 d would move it 0.5 mas along the orbit.
    struct Position {
        double mjd, rho, theta, sigma;
    };
    const Position b_seen[] = {
        {54017.535773, 1.630, 40.74, 0.214},   {54427.234904, 2.200, 44.83, 0.330},
        {54698.455342, 1.829, -144.44, 0.318}, {55056.508177, 2.115, -133.91, 0.279},
        {55063.416862, 1.516, 51.64, 0.174},   {55416.495590, 2.130, 43.59, 0.068},
    };
    const Position c_seen[] = {
        {54017.535773, 44.54, -61.91, 0.217}, {54377.358738, 18.86, 108.01, 0.374},
        {54427.234904, 61.06, 128.69, 0.508}, {55055.480222, 17.18, 103.62, 0.124},
        {55416.495590, 67.73, -55.74, 0.111},
    };
    auto miss_mas = [&](int body, const Position& p) {
        const double th = p.theta * kDegToRad;
        const glm::dvec2 measured(p.rho * std::sin(th), p.rho * std::cos(th));
        return glm::length(offset_mas(body, tdb_from_jd_tdb(p.mjd + 2400000.5)) - measured);
    };
    double worst_b = 0.0;
    for (const Position& p : b_seen) {
        worst_b = std::max(worst_b, miss_mas(b, p));
    }
    check(worst_b < 0.4, "Algol B on the CHARA positions (mas)", worst_b);
    double worst_c = 0.0; // in sigma
    for (const Position& p : c_seen) {
        worst_c = std::max(worst_c, miss_mas(c, p) / p.sigma);
    }
    std::printf("info: Algol on the CHARA positions: B within %.2f mas, C within %.2f sigma\n", worst_b, worst_c);
    check(worst_c < 2.0, "Algol C on the CHARA positions (sigma)", worst_c);

    // TESS minima (tools/algol/algol_ephemeris.py): the geometric conjunction of
    // B in front of A, delayed by AB's depth around the barycenter.
    double worst_min = 0.0;
    for (double bjd : {2458793.97822, 2458811.18274, 2459883.57643, 2459909.38099, 2460614.74683, 2460634.81733}) {
        const double t_seen = tdb_from_jd_tdb(bjd);
        const double t = find_conjunction(scene, b, a, away, t_seen, 0.1 * kSecondsPerDay);
        const double delay = glm::dot(scene.icrf_state_at(ab, t).position, away) / kSpeedOfLightKmS;
        const double dt_min = std::isnan(t) ? 1e9 : (t + delay - t_seen) / 60.0;
        worst_min = std::max(worst_min, std::abs(dt_min));
    }
    check(worst_min < 1.5, "Algol primary minima at the TESS times (min)", worst_min);

    // B is its Roche lobe [KI89] for [BA12]'s masses, turning with the orbit
    // (its x axis on A): as large as a sphere of 3.44 R_sun.
    const Body& star_b = scene.bodies[static_cast<size_t>(b)];
    check(star_b.roche != nullptr, "Algol B fills its Roche lobe");
    if (!star_b.roche) {
        return;
    }
    const double volume_r = star_b.roche->volume_radius_km / kSunRadiusKm;
    check(std::abs(volume_r - 3.44) < 0.01, "Algol B's volume radius (R_sun)", volume_r);
    {
        const double t = tdb_from_jd_tdb(2459743.3);
        const glm::dvec3 to_a = scene.icrf_state_at(a, t).position - scene.icrf_state_at(b, t).position;
        const double off = std::acos(glm::dot(scene.body_to_icrf(star_b, t)[0], glm::normalize(to_a)));
        check(off < 1e-6, "Algol B's long axis on A (rad)", off);
    }

    // Light of the three stars seen from Earth at t, with and without eclipses.
    const glm::dvec3 toward_earth = -away;
    auto flux = [&](double t, bool hide) {
        std::vector<StarDisk> disks;
        for (int s : {a, b, c}) {
            disks.push_back(star_disk_at(scene, s, t));
        }
        if (hide) {
            return disks_flux(disks, toward_earth, 550.0, 0.6);
        }
        double sum = 0.0;
        for (const StarDisk& d : disks) {
            sum += disks_flux(std::span(&d, 1), toward_earth, 550.0, 0.6);
        }
        return sum;
    };

    // Depth of the eclipses in V below the peak light, and the length of
    // the primary eclipse (while anything is hidden).
    const double period_s = 2.8673453 * kSecondsPerDay;
    const double t_min = tdb_from_jd_tdb(2459743.07502 + 555.0 * 2.8673453);
    double peak = 1e9;
    double faintest_out = 0.0; // out of eclipse
    double first = 0.0;
    double last = 0.0;
    for (double dt = -0.5; dt < 0.5; dt += 0.001) { // fraction of the period, around the primary
        const double t = t_min + dt * period_s;
        const double m = light_curve_magnitude(scene, t);
        peak = std::min(peak, m);
        if (flux(t, true) < flux(t, false) * (1.0 - 1e-6)) {
            if (dt < 0.25 && dt > -0.25) {
                first = first == 0.0 ? dt : first;
                last = dt;
            }
        } else {
            faintest_out = std::max(faintest_out, m);
        }
    }
    const double primary = light_curve_magnitude(scene, t_min) - peak;
    const double secondary = light_curve_magnitude(scene, t_min + 0.5 * period_s) - peak;
    std::printf("info: Algol V eclipse depths: primary %.3f, secondary %.3f mag; out of eclipse within %.3f mag;"
                " primary eclipse %.3f of the period\n",
                primary, secondary, faintest_out - peak, last - first);
    check(std::abs(primary - (3.39 - 2.12)) < 0.1, "Algol primary eclipse depth (mag)", primary);
    check(secondary > 0.02 && secondary < 0.1, "Algol secondary eclipse depth (mag)", secondary);
    check(std::abs((last - first) - 0.14) < 0.01, "Algol eclipse duration (fraction of the period)", last - first);
    // Out of eclipse, B's tidal shape and the side A warms vary the light by a
    // few hundredths of a magnitude.
    check(faintest_out - peak > 0.02 && faintest_out - peak < 0.06,
          "Algol out-of-eclipse variation (mag)", faintest_out - peak);

    // A's share of the light at quadrature: between the light-curve solution
    // (0.814 in V, [KI89]) and the spectra (0.915 at 5015 A [KO15] Table 2).
    {
        const double t = t_min + 0.25 * period_s;
        const StarDisk star_a = star_disk_at(scene, a, t);
        const double share = disks_flux(std::span(&star_a, 1), toward_earth, 550.0, 0.6) / flux(t, true);
        check(share > 0.814 && share < 0.915, "Algol A's share of the V light", share);
    }

    // B's disk as CHARA saw it [BA12 Sect. 3]: fitted as a uniform ellipse,
    // 1.04-1.22 times longer than wide over their dates. Here from the second
    // moments of the uniformly bright outline (the facets facing Earth, each
    // weighted by its area on the sky): 1.11-1.27. Moments weigh the pointed
    // end toward A more than a fit to the visibilities would (with B's light in
    // the H band, warmed toward A, they reach 1.34), so only the ranges are
    // compared.
    double least = 1e9;
    double most = 0.0;
    for (const Position& p : b_seen) {
        const double t = tdb_from_jd_tdb(p.mjd + 2400000.5);
        const glm::dmat3 axes = scene.body_to_icrf(star_b, t);
        double sxx = 0.0, syy = 0.0, sxy = 0.0, sx = 0.0, sy = 0.0, sw = 0.0;
        for (const RocheSurface::Facet& f : star_b.roche->facets) {
            const double mu = glm::dot(axes * f.normal, toward_earth);
            if (mu <= 0.0) {
                continue;
            }
            const glm::dvec3 p_sky = axes * f.center;
            const double x = glm::dot(p_sky, east);
            const double y = glm::dot(p_sky, north);
            const double w = f.area_km2 * mu;
            sw += w;
            sx += w * x;
            sy += w * y;
            sxx += w * x * x;
            syy += w * y * y;
            sxy += w * x * y;
        }
        const double mx = sx / sw;
        const double my = sy / sw;
        const double vxx = sxx / sw - mx * mx;
        const double vyy = syy / sw - my * my;
        const double vxy = sxy / sw - mx * my;
        const double mean = 0.5 * (vxx + vyy);
        const double spread = std::sqrt(0.25 * (vxx - vyy) * (vxx - vyy) + vxy * vxy);
        const double aspect = std::sqrt((mean + spread) / (mean - spread));
        least = std::min(least, aspect);
        most = std::max(most, aspect);
    }
    std::printf("info: Algol B's aspect ratio on the CHARA dates: %.3f .. %.3f\n", least, most);
    check(least > 1.0 && least < 1.22 && most > 1.04 && most < 1.3, "Algol B's elongation on the sky", most);
}

// Time of the light curve's maximum (least magnitude) within +-window of t_guess.
double find_light_maximum(const Scene& scene, double t_guess, double window)
{
    double lo = t_guess - window;
    double hi = t_guess + window;
    // Coarse scan, then golden-section search around the best sample.
    double best = lo;
    double best_mag = 1e9;
    for (double t = lo; t <= hi; t += window / 100.0) {
        const double m = light_curve_magnitude(scene, t);
        if (m < best_mag) {
            best_mag = m;
            best = t;
        }
    }
    lo = best - window / 100.0;
    hi = best + window / 100.0;
    for (int k = 0; k < 60; ++k) {
        const double a = hi - (hi - lo) * 0.618033988749895;
        const double b = lo + (hi - lo) * 0.618033988749895;
        if (light_curve_magnitude(scene, a) < light_curve_magnitude(scene, b)) {
            hi = b;
        } else {
            lo = a;
        }
    }
    return 0.5 * (lo + hi);
}

// Delta Cephei [Merand et al. 2015, Engle et al. 2014, Anderson et al. 2015]:
// the pulsation follows the SPIPS model (radii, mean radius) on Engle's
// ephemeris with its period change; the light curve peaks on Engle's maxima
// and on SPIPS's own ephemeris, with the GCVS amplitude; the Cepheid's radial
// velocity is the spectroscopic orbit's.
void test_delta_cephei_scene()
{
    Scene scene = load_scene_or_die("delta_cephei.toml");
    check_model_sky(scene, "Delta Cephei model sky");
    const int a = scene.find("Delta Cephei A");
    const int b = scene.find("Delta Cephei B");
    check(a > 0 && b > 0 && scene.bodies[static_cast<size_t>(a)].pulsation != nullptr, "Delta Cephei A pulsates");
    check(scene.light_curve.enabled() && scene.light_curve.stars.size() == 2, "Delta Cephei light curve of A and B");
    if (a <= 0 || b <= 0 || !scene.bodies[static_cast<size_t>(a)].pulsation) {
        return;
    }
    const Pulsation& pulsation = *scene.bodies[static_cast<size_t>(a)].pulsation;

    // [ME15]: 40.1 .. 44.7 R_sun (diameters 1.362 .. 1.518 mas at 274 pc), on average 43.0.
    double r_min = 1e30, r_max = 0.0, r_sum = 0.0, t_min = 1e9, t_max = 0.0;
    for (const Pulsation::Sample& s : pulsation.samples) {
        r_min = std::min(r_min, s.radius_km);
        r_max = std::max(r_max, s.radius_km);
        r_sum += s.radius_km;
        t_min = std::min(t_min, s.temperature_k);
        t_max = std::max(t_max, s.temperature_k);
    }
    const double r_mean = r_sum / pulsation.samples.size() / kSunRadiusKm;
    check(std::abs(r_mean - 43.0) < 0.15, "Delta Cephei mean radius (R_sun)", r_mean);
    check(std::abs(r_min / kSunRadiusKm - 40.14) < 0.05 && std::abs(r_max / kSunRadiusKm - 44.72) < 0.05,
          "Delta Cephei radius range", r_max / r_min);
    check(t_min > 5300.0 && t_min < 5330.0 && t_max > 6330.0 && t_max < 6360.0, "Delta Cephei T_eff range",
          t_max - t_min);

    // [EN14]'s maxima: t(E) = T0 + P0 E + P0 Pdot E^2 / 2, Pdot = -0.1006 s/yr.
    const double t0 = tdb_from_jd_tdb(2455479.905);
    const double p0 = 5.366208 * kSecondsPerDay;
    const double pdot = -0.1006 / (kDaysPerJulianYear * kSecondsPerDay);
    double worst_cycles = 0.0;
    for (double e : {-25000.0, -5000.0, 1091.0, 5000.0}) {
        const double t = t0 + p0 * e + 0.5 * p0 * pdot * e * e;
        worst_cycles = std::max(worst_cycles, std::abs(pulsation.cycles(t) - e));
    }
    check(worst_cycles < 1e-4, "Delta Cephei cycles on the quadratic ephemeris", worst_cycles);

    // The light curve (a black body at 550 nm) peaks on [EN14]'s maxima, 2010 and 2026.
    for (double e : {0.0, 1091.0}) {
        const double t = t0 + p0 * e + 0.5 * p0 * pdot * e * e;
        const double dt_d = (find_light_maximum(scene, t, 0.5 * kSecondsPerDay) - t) / kSecondsPerDay;
        check(std::abs(dt_d) < 0.03, "Delta Cephei light maximum on the ephemeris (d)", dt_d);
    }
    // [ME15]'s own ephemeris (bolometric maximum, model phase 0, at MJD 48304.7362421 with
    // P = 5.3662906 d and -0.069 s/yr) against the scene's model phase 0 near 2010.
    {
        const double me_t0 = tdb_from_jd_tdb(48304.7362421 + 2400000.5);
        const double me_p = 5.3662906 * kSecondsPerDay;
        const double me_pdot = -0.069 / (kDaysPerJulianYear * kSecondsPerDay);
        const double e = std::round((t0 - me_t0) / me_p);
        const double me_max = me_t0 + me_p * e + 0.5 * me_p * me_pdot * e * e;
        // Where the scene is at model phase 0: cycles + epoch_phase is a whole number.
        const double phase = pulsation.epoch_phase + pulsation.cycles(me_max);
        const double dt_d = (phase - std::round(phase)) * 5.366208;
        check(std::abs(dt_d) < 0.05, "Delta Cephei SPIPS phase 0 against the scene (d)", dt_d);
    }

    // V amplitude: 3.48 .. 4.37 [GCVS] (the model's black body: 0.80).
    double faintest = 0.0;
    double brightest = 99.0;
    for (double t = t0; t < t0 + p0; t += p0 / 400.0) {
        const double m = light_curve_magnitude(scene, t);
        faintest = std::max(faintest, m);
        brightest = std::min(brightest, m);
    }
    std::printf("info: Delta Cephei V amplitude %.3f mag (brightest %.3f)\n", faintest - brightest, brightest);
    check(std::abs((faintest - brightest) - (4.37 - 3.48)) < 0.12, "Delta Cephei V amplitude (mag)",
          faintest - brightest);
    check(brightest < 0.01, "Delta Cephei light curve reaches its reference", brightest);

    // [AN15] Table 3: the Cepheid's radial velocity K (cos(nu + omega) + e cos omega).
    const glm::dvec3 away = unit_toward(337.2927611, 58.4151938);
    const double t_peri = tdb_from_jd_tdb(55649.68 + 2400000.5);
    const double period = 2201.87 * kSecondsPerDay;
    const double e = 0.674;
    const double omega = 246.77 * kDegToRad;
    double worst_rv = 0.0;
    for (double f = 0.0; f < 1.0; f += 0.01) {
        const double t = t_peri + f * period;
        const double m = kTwoPi * f;
        double ea = m;
        for (int k = 0; k < 50; ++k) {
            ea -= (ea - e * std::sin(ea) - m) / (1.0 - e * std::cos(ea));
        }
        const double nu = 2.0 * std::atan2(std::sqrt(1.0 + e) * std::sin(0.5 * ea), std::sqrt(1.0 - e) * std::cos(0.5 * ea));
        const double rv = 1.509 * (std::cos(nu + omega) + e * std::cos(omega));
        worst_rv = std::max(worst_rv, std::abs(glm::dot(scene.icrf_state_at(a, t).velocity, away) - rv));
    }
    check(worst_rv < 0.02, "Delta Cephei radial velocity on the orbit (km/s)", worst_rv);
}

// Kepler-64 [Schwamb et al. 2013]: integrated from the Table 7 joint solution,
// the binary eclipses and the planet transits Aa when they were measured
// (Table 3, 1-sigma 0.9 and 7 min; Table 2, 6-9 min; BJD - 2455000), and Ba +
// Bb stand where Gaia sees them.
void test_kepler64_scene()
{
    Scene scene = load_scene_or_die("kepler64.toml");
    check_model_sky(scene, "Kepler-64 model sky");
    const glm::dvec3 away = unit_toward(298.215070011, 39.955103210);
    const int aa = scene.find("Kepler-64 Aa");
    const int ab = scene.find("Kepler-64 Ab");
    const int b = scene.find("Kepler-64 b");
    check(aa > 0 && ab > 0 && b > 0, "Kepler-64 bodies");
    if (aa <= 0 || ab <= 0 || b <= 0) {
        return;
    }
    struct Event {
        int body;
        int center;
        double time;
    };
    const Event primary[] = {{ab, aa, -32.18064}, {ab, aa, 447.82520}, {ab, aa, 927.83150}};
    const Event secondary[] = {{aa, ab, -24.32048}, {aa, ab, 455.67885}, {aa, ab, 915.68830}};
    const Event transits[] = {{b, aa, 70.80674}, {b, aa, 344.11218}, {b, aa, 613.17869}, {b, aa, 885.91042}};
    auto worst = [&](const auto& events) {
        double w = 0.0;
        for (const Event& e : events) {
            const double t_ref = tdb_from_jd_tdb(2455000.0 + e.time);
            const double t = find_conjunction(scene, e.body, e.center, away, t_ref, 0.5 * kSecondsPerDay);
            w = std::max(w, std::isnan(t) ? 1e9 : std::abs(t - t_ref) / 60.0);
        }
        return w;
    };
    const double w_primary = worst(primary);
    const double w_secondary = worst(secondary);
    const double w_transit = worst(transits);
    std::printf("info: Kepler-64 vs Schwamb 2013: primary eclipses %.2f, secondary %.2f, transits %.2f min\n",
                w_primary, w_secondary, w_transit);
    check(w_primary < 3.0, "Kepler-64 primary eclipses (min)", w_primary);
    check(w_secondary < 20.0, "Kepler-64 secondary eclipses (min)", w_secondary);
    check(w_transit < 20.0, "Kepler-64 planet transits (min)", w_transit);

    // Ba + Bb: 0.7043" at position angle 123.29 deg from the A system (J2016.0, 1906.58 pc).
    const int bb = scene.find("Kepler-64 B");
    check(bb > 0, "Kepler-64 B");
    if (bb > 0) {
        const glm::dvec3 r = scene.icrf_state_at(bb, tdb_from_julian_year(2016.0)).position;
        const double ra = 298.215070011 * kDegToRad;
        const double dec = 39.955103210 * kDegToRad;
        const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
        const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
        const double pa = wrap_two_pi(std::atan2(glm::dot(r, east), glm::dot(r, north))) * kRadToDeg;
        const double sep = std::hypot(glm::dot(r, east), glm::dot(r, north)) / (1906.58 * kAuKm);
        check(std::abs(pa - 123.29) < 0.05, "Kepler-64 B position angle (deg)", pa);
        check(std::abs(sep - 0.7043) < 0.001, "Kepler-64 B separation (arcsec)", sep);
    }
}

// PSR B1620-26 [Thorsett et al. 1999]: the pulsar's line-of-sight motion
// reproduces the Table 1 timing orbit (Roemer delay x (1 - e) sin omega at
// periastron, semi-amplitude K), the inner binary moves 6.4 lt-s around the
// triple's barycentre and recedes through its node at JD 2449104 (Table 2),
// the two orbits are 40 deg apart [Sigurdsson & Thorsett 2005], and the sky is
// M4 seen from inside.
void test_psr_b1620_scene()
{
    Scene scene = load_scene_or_die("psr_b1620.toml");
    const glm::dvec3 away = unit_toward(245.9092575, -26.5316025);
    const int pulsar = scene.find("PSR B1620-26 A");
    const int wd = scene.find("PSR B1620-26 B");
    const int ab = scene.find("PSR B1620-26 AB");
    const int planet = scene.find("PSR B1620-26 b");
    check(pulsar > 0 && wd > 0 && ab > 0 && planet > 0, "PSR B1620-26 bodies");
    if (pulsar <= 0 || wd <= 0 || ab <= 0 || planet <= 0) {
        return;
    }
    const double lt_s = kSpeedOfLightKmS; // km per light-second
    auto pulsar_depth = [&](double t) {
        const State p = scene.icrf_state_at(pulsar, t);
        const State c = scene.icrf_state_at(ab, t);
        return std::make_pair(glm::dot(p.position - c.position, away), glm::dot(p.velocity - c.velocity, away));
    };
    const double x = 64.809460;
    const double e = 0.02531545;
    const double w = 117.1291 * kDegToRad;
    const double t_peri = tdb_from_jd_tdb(2448728.76242);
    const double depth_peri = pulsar_depth(t_peri).first / lt_s;
    check(std::abs(depth_peri - x * (1.0 - e) * std::sin(w)) < 1e-4, "PSR B1620-26 Roemer delay at periastron (s)",
          depth_peri);
    double k_max = 0.0;
    for (int i = 0; i < 2000; ++i) {
        k_max = std::max(k_max, pulsar_depth(t_peri + 191.44281 * kSecondsPerDay * i / 2000.0).second);
    }
    // Maximum radial velocity K (1 + e cos omega), K = 2 pi x c / (P sqrt(1 - e^2)).
    const double k = kTwoPi * x * lt_s / (191.44281 * kSecondsPerDay * std::sqrt(1.0 - e * e));
    check(std::abs(k_max / (k * (1.0 + e * std::cos(w))) - 1.0) < 1e-4, "PSR B1620-26 pulsar velocity amplitude",
          k_max);

    // The binary around the triple's barycentre (the root): node at JD 2449104, 6.4 lt-s.
    const double t_node = tdb_from_jd_tdb(2449104.0);
    const double quarter = 61.8 * kDaysPerJulianYear * kSecondsPerDay / 4.0;
    const State at_node = scene.icrf_state_at(ab, t_node);
    const double depth_q = glm::dot(scene.icrf_state_at(ab, t_node + quarter).position, away) / lt_s;
    check(std::abs(glm::dot(at_node.position, away)) / lt_s < 1e-6 && glm::dot(at_node.velocity, away) > 0.0,
          "PSR B1620-26 binary recedes through its node in 1993");
    check(std::abs(depth_q - 6.4) < 0.01, "PSR B1620-26 binary outer amplitude (lt-s)", depth_q);

    scene.update(tdb_from_julian_year(2026.0));
    const double tilt = std::acos(glm::dot(scene.orbit_normal(wd), scene.orbit_normal(planet))) * kRadToDeg;
    check(std::abs(tilt - 40.0) < 0.05, "PSR B1620-26 relative inclination (deg)", tilt);

    // The sky: ~13,500 members, many brighter than the naked-eye limit, the
    // brightest far brighter than Sirius.
    int naked_eye = 0;
    double brightest = 99.0;
    for (const CatalogStar& s : scene.sky.stars) {
        naked_eye += s.vmag < 6.5 ? 1 : 0;
        brightest = std::min(brightest, s.vmag);
    }
    std::printf("info: M4 from PSR B1620-26: %zu stars, %d brighter than V = 6.5, brightest V = %.2f\n",
                scene.sky.stars.size(), naked_eye, brightest);
    check(scene.sky.stars.size() > 13000, "M4 sky star count", static_cast<double>(scene.sky.stars.size()));
    check(naked_eye > 5000 && brightest < -4.0, "M4 sky is bright", brightest);
    // Behind the members, the model's field stars and Milky Way from M4.
    check_model_sky(scene, "M4 field sky");

    // The pulsar's spin axis is the inner orbit's normal.
    check(scene.bodies[static_cast<size_t>(pulsar)].pulsar.enabled &&
              scene.bodies[static_cast<size_t>(pulsar)].pulsar.spin_axis_orbit_of == wd,
          "PSR B1620-26 spin axis");
}

// WR 104's dust as seen from Earth at t_tdb: sky offsets (east, north; mas)
// from the pair's barycenter and the light of each grain (pinwheel_grain_brightness).
struct SkyGrain {
    double east = 0.0;
    double north = 0.0;
    double brightness = 0.0;
};

std::vector<SkyGrain> wr104_sky_grains(const Pinwheel& pinwheel, double t_tdb, double distance_pc)
{
    const double ra = 270.5171852 * kDegToRad;
    const double dec = -23.6283801 * kDegToRad;
    const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
    const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
    const double mas_per_km = 1000.0 / (distance_pc * kAuKm);
    std::vector<SkyGrain> out;
    out.reserve(pinwheel.grains.size());
    for (const PinwheelGrain& g : pinwheel.grains) {
        const glm::dvec3 r = pinwheel_grain_offset(pinwheel, g, t_tdb);
        out.push_back({glm::dot(r, east) * mas_per_km, glm::dot(r, north) * mas_per_km,
                       pinwheel_grain_brightness(pinwheel, glm::length(r))});
    }
    return out;
}

// Position angle (deg, north through east) at the origin of the Archimedean
// spiral that best follows the grains within max_mas of the center, for a
// spiral that winds outward toward larger position angles (sense = +1, as
// WR 104's [TU08 Fig. 8]) or smaller ones (-1): the brightness-weighted
// circular mean of PA - sense * 360 rho / step. `coherence` (0..1) is the
// length of that mean: near 1 when the light lies along such a spiral.
double fit_spiral_theta0(const std::vector<SkyGrain>& grains, double step_mas, double max_mas, double sense,
                         double* coherence)
{
    double c = 0.0;
    double s = 0.0;
    double w = 0.0;
    for (const SkyGrain& g : grains) {
        const double rho = std::hypot(g.east, g.north);
        if (rho > max_mas || g.brightness <= 0.0) {
            continue;
        }
        const double a = std::atan2(g.east, g.north) - sense * kTwoPi * rho / step_mas;
        c += g.brightness * std::cos(a);
        s += g.brightness * std::sin(a);
        w += g.brightness;
    }
    *coherence = std::hypot(c, s) / w;
    return wrap_two_pi(std::atan2(s, c)) * kRadToDeg;
}

// Axis ratio (minor / major, from the second moments) of the light of the
// grains between rho0 and rho1 mas on the sky.
double sky_axis_ratio(const std::vector<SkyGrain>& grains, double rho0, double rho1)
{
    double xx = 0.0;
    double yy = 0.0;
    double xy = 0.0;
    for (const SkyGrain& g : grains) {
        const double rho = std::hypot(g.east, g.north);
        if (rho < rho0 || rho > rho1) {
            continue;
        }
        xx += g.brightness * g.east * g.east;
        yy += g.brightness * g.north * g.north;
        xy += g.brightness * g.east * g.north;
    }
    const double mean = 0.5 * (xx + yy);
    const double diff = std::sqrt(0.25 * (xx - yy) * (xx - yy) + xy * xy);
    return std::sqrt((mean - diff) / (mean + diff));
}

// WR 104 [HI24] [RI26] [TU08] [SO18]: the binary's spectroscopic orbit and
// masses, and the dust spiral (scene/pinwheel.hpp) as the images show it.
void test_wr104_scene()
{
    Scene scene = load_scene_or_die("wr104.toml");
    check_model_sky(scene, "WR 104 model sky");
    const glm::dvec3 away = unit_toward(270.5171852, -23.6283801);
    const int wr = scene.find("WR 104 WC9d");
    const int ob = scene.find("WR 104 OB");
    const int pair = scene.find("WR 104 WC+OB");
    const int b = scene.find("WR 104 B");
    check(wr > 0 && ob > 0 && pair > 0 && b > 0, "WR 104 bodies");
    if (wr <= 0 || ob <= 0 || pair <= 0 || b <= 0) {
        return;
    }
    const double period_days = 241.54;
    const double e0_jd = 2449439.7;

    // [HI24]: the WR star in front at E0 (inferior conjunction), still after
    // 60 orbits; radial-velocity amplitudes 63.8 and 23.7 km/s.
    for (int n : {0, 60}) {
        const double t_ref = tdb_from_jd_tdb(e0_jd + period_days * n);
        const double t = find_conjunction(scene, wr, ob, away, t_ref, 2.0 * kSecondsPerDay);
        const double dt_days = std::isnan(t) ? 1e9 : (t - t_ref) / kSecondsPerDay;
        check(std::abs(dt_days) < 0.01, "WR 104 WR star in front on [HI24]'s E0 (days)", dt_days);
    }
    double k_wr = 0.0;
    double k_ob = 0.0;
    for (int k = 0; k < 400; ++k) {
        const double t = tdb_from_jd_tdb(e0_jd + period_days * k / 400.0);
        const glm::dvec3 v0 = scene.icrf_state_at(pair, t).velocity;
        k_wr = std::max(k_wr, glm::dot(scene.icrf_state_at(wr, t).velocity - v0, away));
        k_ob = std::max(k_ob, glm::dot(scene.icrf_state_at(ob, t).velocity - v0, away));
    }
    check(std::abs(k_wr - 63.8) < 0.05, "WR 104 K of the WR star (km/s)", k_wr);
    check(std::abs(k_ob - 23.7) < 0.05, "WR 104 K of the OB star (km/s)", k_ob);
    // [RI26] Sect. 5.1: 14.9 and 40.2 M_sun at i = 41.88 deg (their numbers are
    // 2.4% below what their formula gives for it: 15.26, 41.09).
    const double m_wr = scene.bodies[static_cast<size_t>(wr)].body_gm_km3_s2 / kSunGmKm3S2;
    const double m_ob = scene.bodies[static_cast<size_t>(ob)].body_gm_km3_s2 / kSunGmKm3S2;
    std::printf("info: WR 104 masses %.2f + %.2f M_sun\n", m_wr, m_ob);
    check(std::abs(m_wr / 14.9 - 1.0) < 0.03, "WR 104 WR mass against [RI26] (M_sun)", m_wr);
    check(std::abs(m_ob / 40.2 - 1.0) < 0.03, "WR 104 OB mass against [RI26] (M_sun)", m_ob);

    // B relative to the pair at the [SO18] epoch: PA 74.8 deg, 975.5 mas.
    {
        const double t = tdb_from_jd_tdb(2457590.6);
        const glm::dvec3 r = scene.icrf_state_at(b, t).position - scene.icrf_state_at(pair, t).position;
        const double ra = 270.5171852 * kDegToRad;
        const double dec = -23.6283801 * kDegToRad;
        const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
        const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
        const double pa = wrap_two_pi(std::atan2(glm::dot(r, east), glm::dot(r, north))) * kRadToDeg;
        const double sep = std::hypot(glm::dot(r, east), glm::dot(r, north)) / (2580.0 * kAuKm) * 1000.0;
        check(std::abs(pa - 74.8) < 0.05, "WR 104 B position angle (deg)", pa);
        check(std::abs(sep - 975.5) < 0.5, "WR 104 B separation (mas)", sep);
    }

    const Pinwheel& pw = scene.pinwheel;
    check(pw.enabled() && pw.center == pair && pw.primary == wr && pw.secondary == ob, "WR 104 pinwheel");
    if (!pw.enabled()) {
        return;
    }
    // [SO18] Table 7: 66 +- 3 mas per turn and r_nuc = 12 +- 2 mas at 2.58 kpc;
    // [TU08]: 0.28 +- 0.02 mas/day.
    const double mas_per_km = 1000.0 / (2580.0 * kAuKm);
    const double step_mas = pw.step_km() * mas_per_km;
    check(std::abs(step_mas - 66.0) < 3.0, "WR 104 spiral step (mas)", step_mas);
    check(std::abs(step_mas / period_days - 0.28) < 0.02, "WR 104 dust motion (mas/day)", step_mas / period_days);
    check(std::abs(pw.nucleation_radius_km * mas_per_km - 12.0) < 0.1, "WR 104 nucleation radius (mas)",
          pw.nucleation_radius_km * mas_per_km);
    check(std::abs(pinwheel_dust_temperature(pw, 2.0 * pw.nucleation_radius_km) - 2000.0 * std::pow(2.0, -0.45)) <
              1e-9,
          "WR 104 dust temperature law");

    // Grains coast at the wind speed, and the youngest ring circles the
    // current WR -> OB line at half the opening angle.
    const double t98 = tdb_from_jd_tdb(2450918.0); // [TU08]'s first epoch, 1998 Apr 14
    {
        const glm::dvec3 axis = glm::normalize(scene.icrf_state_at(ob, t98).position -
                                               scene.icrf_state_at(wr, t98).position);
        double youngest = 1e30;
        for (const PinwheelGrain& g : pw.grains) {
            youngest = std::min(youngest, pinwheel_grain_age(pw, g, t98));
        }
        glm::dvec3 ring_axis(0.0);
        double worst_speed = 0.0;
        for (const PinwheelGrain& g : pw.grains) {
            const double age = pinwheel_grain_age(pw, g, t98);
            if (age < youngest + 1.0) {
                ring_axis += glm::dvec3(g.direction);
            }
            if (age > kSecondsPerDay && age < pw.turns * pw.period_s - kSecondsPerDay) {
                const glm::dvec3 d = pinwheel_grain_offset(pw, g, t98 + kSecondsPerDay) - pinwheel_grain_offset(pw, g, t98);
                worst_speed = std::max(worst_speed, std::abs(glm::length(d) / kSecondsPerDay - pw.speed_km_s));
            }
        }
        const double lag_deg = std::acos(std::min(1.0, glm::dot(glm::normalize(ring_axis), axis))) * kRadToDeg;
        const double ring_deg = 360.0 * youngest / pw.period_s;
        check(std::abs(lag_deg - ring_deg) < 0.01, "WR 104 youngest ring around the WR -> OB line (deg)",
              lag_deg - ring_deg);
        // (directions are floats: unit to ~1e-7)
        check(worst_speed < 1e-6 * pw.speed_km_s, "WR 104 grains coast at the wind speed (km/s)", worst_speed);
    }

    // As seen from Earth the pair turns clockwise ([TU08] Fig. 8: the spiral
    // winds counterclockwise outward): the WR -> OB line's position angle falls.
    {
        const double ra = 270.5171852 * kDegToRad;
        const double dec = -23.6283801 * kDegToRad;
        const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
        const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
        auto pa = [&](double t) {
            const glm::dvec3 r = scene.icrf_state_at(ob, t).position - scene.icrf_state_at(wr, t).position;
            return std::atan2(glm::dot(r, east), glm::dot(r, north));
        };
        const double turn = wrap_pi(pa(t98 + kSecondsPerDay) - pa(t98)) * kRadToDeg;
        check(turn < -0.5, "WR 104 turns clockwise on the sky (deg/day)", turn);
    }
    // At [TU08]'s first epoch the light lies along such a spiral, starting at
    // PA theta_0 = 269 +- 1 deg ([TU08] Sect. 2.2; the node of the orbit is
    // fitted to it, within [HI24]'s phase uncertainty of ~8 deg then).
    const std::vector<SkyGrain> sky98 = wr104_sky_grains(pw, t98, 2580.0);
    double coherence = 0.0;
    double coherence_reversed = 0.0;
    const double theta0 = fit_spiral_theta0(sky98, step_mas, 100.0, 1.0, &coherence);
    fit_spiral_theta0(sky98, step_mas, 100.0, -1.0, &coherence_reversed);
    std::printf("info: WR 104 spiral at 1998 Apr 14: theta_0 = %.1f deg, coherence %.2f (reversed %.2f)\n", theta0,
                coherence, coherence_reversed);
    check(coherence > 0.4 && coherence > 1.5 * coherence_reversed,
          "WR 104 light along a spiral winding counterclockwise outward", coherence);
    check(std::abs(wrap_pi((theta0 - 269.0) * kDegToRad)) * kRadToDeg < 1.0, "WR 104 spiral theta_0 (deg)", theta0);

    // The wide cone keeps the projected coils as round as [TU08]'s flat spiral
    // fits allow (i <= 16 deg: cos 16 deg = 0.961) despite the 42 deg tilt; a
    // 60 deg cone ([RI26] Fig. 6) would not.
    const double round = sky_axis_ratio(sky98, step_mas, 3.0 * step_mas);
    Pinwheel narrow = pw;
    narrow.opening_angle_deg = 60.0;
    build_pinwheel_grains(scene, narrow, 200, 24);
    const double round60 = sky_axis_ratio(wr104_sky_grains(narrow, t98, 2580.0), step_mas, 3.0 * step_mas);
    std::printf("info: WR 104 coils 1-3 on the sky, axis ratio %.3f (%.0f deg cone), %.3f (60 deg)\n", round,
                pw.opening_angle_deg, round60);
    check(round > std::cos(16.0 * kDegToRad), "WR 104 coils round on the sky", round);

    // In visible light ([HA04]: 0.01 um grains scatter a negligible share): at
    // the nucleation radius the dust's own glow outshines the scattered
    // starlight; a turn out, cooled below 1,000 K, it is the other way round.
    {
        const double share = pinwheel_star_visible_share(scene, pw);
        const PinwheelVisibleLight inner = pinwheel_visible_light(pw, pw.nucleation_radius_km, share);
        const PinwheelVisibleLight outer = pinwheel_visible_light(pw, pw.step_km(), share);
        std::printf("info: WR 104 visible light: stars %.3f in 400-700 nm; dust at r_nuc %.2e glow, %.2e "
                    "scattered; at a turn %.2e, %.2e\n",
                    share, inner.thermal, inner.scattered, outer.thermal, outer.scattered);
        check(share > 0.01 && share < 0.1, "WR 104 stars' visible share", share);
        check(inner.thermal > 10.0 * inner.scattered, "WR 104 dust glows in visible light where it forms",
              inner.thermal / inner.scattered);
        check(outer.scattered > 10.0 * outer.thermal, "WR 104 dust a turn out only scatters visible light",
              outer.scattered / outer.thermal);
    }

    // [SO18] Table 3: the azimuthally averaged K-band profile falls as rho^-3.0
    // +- 0.2 between their resolution limit and the noise floor.
    {
        constexpr int kBins = 12;
        const double rho0 = 40.0;
        const double rho1 = 400.0;
        double light[kBins] = {};
        for (const SkyGrain& g : sky98) {
            const double rho = std::hypot(g.east, g.north);
            const int bin = static_cast<int>(std::floor(kBins * std::log(rho / rho0) / std::log(rho1 / rho0)));
            if (bin >= 0 && bin < kBins) {
                light[bin] += g.brightness;
            }
        }
        double sx = 0.0;
        double sy = 0.0;
        double sxx = 0.0;
        double sxy = 0.0;
        for (int k = 0; k < kBins; ++k) {
            const double lo = rho0 * std::pow(rho1 / rho0, double(k) / kBins);
            const double hi = rho0 * std::pow(rho1 / rho0, double(k + 1) / kBins);
            const double x = std::log(std::sqrt(lo * hi));
            const double y = std::log(light[k] / (kPi * (hi * hi - lo * lo)));
            sx += x;
            sy += y;
            sxx += x * x;
            sxy += x * y;
        }
        const double slope = (kBins * sxy - sx * sy) / (kBins * sxx - sx * sx);
        check(std::abs(slope + 3.0) < 0.2, "WR 104 radial profile slope", slope);
    }
}

} // namespace

void run_stellar_system_tests()
{
    test_alpha_centauri_scene();
    test_trappist1_scene();
    test_kepler223_scene();
    test_kepler47_scene();
    test_tic168789840_scene();
    test_algol_scene();
    test_delta_cephei_scene();
    test_kepler64_scene();
    test_psr_b1620_scene();
    test_wr104_scene();
}
