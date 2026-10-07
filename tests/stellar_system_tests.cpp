// Tests for the stellar and exoplanet systems.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/ephemeris.hpp"
#include "scene/star_catalog.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>

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
        {"Proxima Centauri c (unconfirmed)", 2455892.0, 1900.0, 1.48},
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

// Kepler-64 [Schwamb et al. 2013]: integrated from the Table 7 joint solution,
// the binary eclipses and the planet transits Aa when they were measured
// (Table 3, 1-sigma 0.9 and 7 min; Table 2, 6-9 min; BJD - 2455000), and Ba +
// Bb stand where Gaia sees them.
void test_kepler64_scene()
{
    Scene scene = load_scene_or_die("kepler64.toml");
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
    check(scene.sky.milky_way.empty(), "M4 sky without the Milky Way seen from Earth");

    // The pulsar's spin axis is the inner orbit's normal.
    check(scene.bodies[static_cast<size_t>(pulsar)].pulsar.enabled &&
              scene.bodies[static_cast<size_t>(pulsar)].pulsar.spin_axis_orbit_of == wd,
          "PSR B1620-26 spin axis");
}

} // namespace

void run_stellar_system_tests()
{
    test_alpha_centauri_scene();
    test_trappist1_scene();
    test_kepler223_scene();
    test_kepler47_scene();
    test_tic168789840_scene();
    test_kepler64_scene();
    test_psr_b1620_scene();
}
