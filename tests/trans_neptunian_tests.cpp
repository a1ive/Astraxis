// Tests for the trans-Neptunian scene: its N-body integration against JPL Horizons, and
// the resonances, closest approaches and orbits its captions describe.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

using namespace astraxis;

namespace {

// [GM] in the scene: the Sun with the inner planets, and Neptune's system.
constexpr double kGmSun = 132713233263.351411;
constexpr double kGmNeptune = 6836527.10058;
constexpr double kYearS = kDaysPerJulianYear * kSecondsPerDay;

// Heliocentric osculating elements in the J2000 ecliptic (radians, km).
struct Elements {
    double a = 0.0;
    double e = 0.0;
    double inc = 0.0;
    double lambda = 0.0; // mean longitude
    double varpi = 0.0;  // longitude of perihelion
    double node = 0.0;   // longitude of the ascending node
};

Elements elements(const State& s, double mu)
{
    const double eps = 84381.412 / 3600.0 * kDegToRad;
    auto to_ecliptic = [&](const glm::dvec3& v) {
        return glm::dvec3(v.x, std::cos(eps) * v.y + std::sin(eps) * v.z, -std::sin(eps) * v.y + std::cos(eps) * v.z);
    };
    const glm::dvec3 r = to_ecliptic(s.position);
    const glm::dvec3 v = to_ecliptic(s.velocity);
    const glm::dvec3 h = glm::cross(r, v);
    const glm::dvec3 hn = glm::normalize(h);
    const glm::dvec3 ev = glm::cross(v, h) / mu - r / glm::length(r);
    Elements el;
    el.e = glm::length(ev);
    el.a = 1.0 / (2.0 / glm::length(r) - glm::dot(v, v) / mu);
    el.inc = std::acos(hn.z);
    const double node = std::atan2(h.x, -h.y);
    const glm::dvec3 n(std::cos(node), std::sin(node), 0.0);
    const glm::dvec3 en = ev / el.e;
    const glm::dvec3 rn = glm::normalize(r);
    const double omega = std::atan2(glm::dot(glm::cross(n, en), hn), glm::dot(n, en));
    const double f = std::atan2(glm::dot(glm::cross(en, rn), hn), glm::dot(en, rn));
    const double big_e = 2.0 * std::atan(std::sqrt((1.0 - el.e) / (1.0 + el.e)) * std::tan(f / 2.0));
    const double mean = big_e - el.e * std::sin(big_e);
    el.node = node;
    el.varpi = node + omega;
    el.lambda = el.varpi + mean;
    return el;
}

double wrap_deg(double deg)
{
    deg = std::fmod(deg, 360.0);
    return deg < 0.0 ? deg + 360.0 : deg;
}

struct Libration {
    double min_deg = 360.0; // range of the smoothed angle, measured from `cut`
    double max_deg = 0.0;
    double turns = 0.0;     // net turns of the smoothed angle (circulation shows as |turns| >= 1)
    double period_years = 0.0;
};

// The resonant angle p lambda - q lambda_N - (p - q - nodes) varpi - nodes Omega, sampled
// every 25 years over the integration and smoothed over 500 years (the short-period terms
// average out). Angles are reported in [cut, cut + 360).
Libration libration(const Scene& scene, int body, int neptune, int p, int q, double cut_deg, int nodes = 0)
{
    const double t0 = tdb_from_julian_year(-7900.0);
    const double t1 = tdb_from_julian_year(51900.0);
    const double dt = 25.0 * kYearS;
    std::vector<std::complex<double>> z;
    std::vector<double> t;
    for (double tk = t0; tk <= t1; tk += dt) {
        const Elements o = elements(scene.icrf_state_at(body, tk), kGmSun);
        const Elements n = elements(scene.icrf_state_at(neptune, tk), kGmSun + kGmNeptune);
        z.push_back(std::polar(1.0, p * o.lambda - q * n.lambda - (p - q - nodes) * o.varpi - nodes * o.node));
        t.push_back(tk);
    }
    constexpr int kHalf = 10;
    Libration out;
    double prev = 0.0;
    double unwrapped = 0.0;
    std::vector<double> rise; // upward crossings of the mean
    std::vector<double> smooth;
    for (size_t k = kHalf; k + kHalf < z.size(); ++k) {
        std::complex<double> sum = 0.0;
        for (size_t j = k - kHalf; j <= k + kHalf; ++j) {
            sum += z[j];
        }
        const double deg = std::arg(sum) * kRadToDeg;
        if (!smooth.empty()) {
            unwrapped += wrap_deg(deg - prev + 180.0) - 180.0;
        }
        prev = deg;
        const double shifted = wrap_deg(deg - cut_deg) + cut_deg;
        smooth.push_back(shifted);
        out.min_deg = std::min(out.min_deg, shifted);
        out.max_deg = std::max(out.max_deg, shifted);
    }
    out.turns = unwrapped / 360.0;
    const double mid = 0.5 * (out.min_deg + out.max_deg);
    for (size_t k = 1; k < smooth.size(); ++k) {
        if (smooth[k - 1] < mid && smooth[k] >= mid) {
            rise.push_back(t[k + kHalf]);
        }
    }
    if (rise.size() >= 2) {
        out.period_years = (rise.back() - rise.front()) / (rise.size() - 1) / kYearS;
    }
    return out;
}

// The integration takes a while, so the tests share one copy.
Scene& tno_scene()
{
    static Scene scene = [] {
        const auto start = std::chrono::steady_clock::now();
        Scene loaded = load_scene_or_die("trans_neptunian.toml");
        const std::chrono::duration<double> took = std::chrono::steady_clock::now() - start;
        std::printf("info: trans_neptunian.toml loads in %.2f s (60,000 years of N-body)\n", took.count());
        return loaded;
    }();
    return scene;
}

void check_position(const Scene& scene, const char* name, double jd, const glm::dvec3& horizons_km, double tol_au)
{
    const int b = scene.find(name);
    const glm::dvec3 p = scene.icrf_state_at(b, tdb_from_jd_tdb(jd)).position;
    const double off = glm::length(p - horizons_km) / kAuKm;
    const std::string what = std::string(name) + " near Horizons (au)";
    check(off < tol_au, what.c_str(), off);
}

void test_horizons()
{
    const Scene& scene = tno_scene();
    // JPL Horizons, heliocentric ICRF (DE441 barycenters; Sedna JPL#51, 2001 QR322 JPL#18),
    // 1950-01-01 and 2100-01-01 00:00 TDB: three quarters of a century either side of the
    // epoch. The scene leaves the inner planets out and the small bodies massless.
    const double jd1950 = 2433282.5;
    const double jd2100 = 2488069.5;
    check_position(scene, "Neptune", jd1950, {-4.352194623448150e9, -1.204353882831413e9, -3.846359404596304e8}, 0.01);
    check_position(scene, "Neptune", jd2100, {-4.347198940227345e9, 1.100540332667645e9, 5.587006007437096e8}, 0.01);
    check_position(scene, "Pluto system", jd1950, {-3.969438644326406e9, 3.031115209821407e9, 2.141548453110631e9}, 0.01);
    check_position(scene, "Pluto system", jd2100, {5.935341284192660e9, 4.259877013287137e9, -4.592031002406538e8}, 0.01);
    check_position(scene, "Sedna", jd1950, {1.484315903508949e10, 7.264381199548445e9, -8.174714073940578e7}, 0.01);
    check_position(scene, "Sedna", jd2100, {-4.631554446364971e9, 1.022975817730068e10, 3.069951264715179e9}, 0.01);
    check_position(scene, "2001 QR322", jd1950, {-2.256167498071764e9, -3.588928366241788e9, -1.443338114156276e9}, 0.01);
    check_position(scene, "2001 QR322", jd2100, {-4.135863474910993e9, -1.780483764158350e9, -6.802158952573375e8}, 0.01);
}

void test_resonances()
{
    const Scene& scene = tno_scene();
    const int neptune = scene.find("Neptune");
    struct Case {
        const char* name;
        int p, q;
        double cut, lo, hi; // the smoothed angle stays in [lo, hi] (degrees, from cut)
    };
    // [CH03], [DES]; ranges from this integration with a margin.
    const Case cases[] = {
        {"2001 QR322", 1, 1, -180.0, 20.0, 110.0},  // L4: about +60
        {"2008 LC18", 1, 1, 0.0, 250.0, 340.0},     // L5: about -60
        {"1995 DA2", 4, 3, 0.0, 100.0, 270.0},
        {"Pluto system", 3, 2, 0.0, 90.0, 270.0},
        {"Orcus", 3, 2, 0.0, 100.0, 265.0},
        {"1996 TR66", 2, 1, -180.0, 40.0, 95.0},    // asymmetric, ahead of Neptune
        {"2000 JG81", 2, 1, 0.0, 10.0, 350.0},      // symmetric: about 180, wide
        {"1999 DE9", 5, 2, 0.0, 45.0, 320.0},
        {"Gonggong", 10, 3, 0.0, 125.0, 235.0},
    };
    for (const Case& c : cases) {
        const Libration l = libration(scene, scene.find(c.name), neptune, c.p, c.q, c.cut);
        std::printf("info: %s %d:%d angle %.1f .. %.1f deg, net %.2f turns, period ~%.0f yr\n", c.name, c.p, c.q,
                    l.min_deg, l.max_deg, l.turns, l.period_years);
        const std::string what = std::string(c.name) + " librates in its resonance";
        check(l.min_deg >= c.lo && l.max_deg <= c.hi && std::abs(l.turns) < 0.5, what.c_str(), l.max_deg - l.min_deg);
    }

    // 2026 ED138 at the edge of the 12:7 [HO26]: 12 lambda - 7 lambda_N - 3 varpi - 2 Omega
    // librates about 180 deg here ("swings within it here").
    const Libration ed138 = libration(scene, scene.find("2026 ED138"), neptune, 12, 7, 0.0, 2);
    std::printf("info: 2026 ED138 12:7 angle %.1f .. %.1f deg, net %.2f turns, period ~%.0f yr\n", ed138.min_deg,
                ed138.max_deg, ed138.turns, ed138.period_years);
    check(ed138.min_deg > 80.0 && ed138.max_deg < 280.0 && std::abs(ed138.turns) < 0.5,
          "2026 ED138 librates in the 12:7", ed138.max_deg - ed138.min_deg);
    // Haumea and 2026 ED138 near the resonance's semi-major axis, 30.07 (12/7)^(2/3) = 43.06 au.
    for (const char* name : {"Haumea", "2026 ED138"}) {
        double sum = 0.0;
        int n = 0;
        for (double y = -7900.0; y <= 51900.0; y += 50.0, ++n) {
            sum += elements(scene.icrf_state_at(scene.find(name), tdb_from_julian_year(y)), kGmSun).a;
        }
        const double a = sum / n / kAuKm;
        const std::string what = std::string(name) + " mean a near the 12:7 (au)";
        check(std::abs(a - 43.06) < 0.2, what.c_str(), a);
    }

    // The Trojans swing about their points over about 9,500 years (the caption and trails).
    for (const char* name : {"2001 QR322", "2008 LC18"}) {
        const Libration l = libration(scene, scene.find(name), neptune, 1, 1, name[3] == '1' ? -180.0 : 0.0);
        const std::string what = std::string(name) + " libration period about 9,500 yr";
        check(l.period_years > 8500.0 && l.period_years < 10500.0, what.c_str(), l.period_years);
    }

    // Pluto: about 180 deg, +-80 deg, over about 20,000 years [CO65] (the captions).
    const Libration pluto = libration(scene, scene.find("Pluto system"), neptune, 3, 2, 0.0);
    const double amplitude = 0.5 * (pluto.max_deg - pluto.min_deg);
    check(std::abs(0.5 * (pluto.max_deg + pluto.min_deg) - 180.0) < 5.0, "Pluto librates about 180 deg",
          0.5 * (pluto.max_deg + pluto.min_deg));
    check(amplitude > 70.0 && amplitude < 90.0, "Pluto's libration amplitude about 80 deg", amplitude);
    check(pluto.period_years > 18000.0 && pluto.period_years < 22000.0, "Pluto's libration period about 20,000 yr",
          pluto.period_years);
}

// Closest approach to Neptune over the integration, sampled every 2 years (au).
double closest_to_neptune(const Scene& scene, const char* name)
{
    const int b = scene.find(name);
    const int n = scene.find("Neptune");
    double best = 1e300;
    for (double y = -7900.0; y <= 51900.0; y += 2.0) {
        const double t = tdb_from_julian_year(y);
        best = std::min(best, glm::length(scene.icrf_state_at(b, t).position - scene.icrf_state_at(n, t).position));
    }
    return best / kAuKm;
}

void test_orbits()
{
    const Scene& scene = tno_scene();
    // "the two never come closer than 17 au" [INT]
    const double pluto = closest_to_neptune(scene, "Pluto system");
    std::printf("info: Pluto comes within %.2f au of Neptune\n", pluto);
    check(pluto > 17.0, "Pluto stays more than 17 au from Neptune", pluto);
    // Pluto's perihelion is inside Neptune's orbit.
    const Elements p = elements(scene.icrf_state_at(scene.find("Pluto system"), 0.0), kGmSun);
    check(p.a * (1.0 - p.e) / kAuKm < 30.0, "Pluto's perihelion inside 30 au", p.a * (1.0 - p.e) / kAuKm);

    // "never come within 75 au of the Sun" [BR04], [TS14]: the closest distance over the
    // span (osculating heliocentric perihelia of such wide orbits swing with the Sun's
    // own motion around the barycenter, so the distance itself is sampled).
    for (const char* name : {"Sedna", "2012 VP113"}) {
        double r_min = 1e300;
        for (double y = -7900.0; y <= 51900.0; y += 2.0) {
            const double r = glm::length(scene.icrf_state_at(scene.find(name), tdb_from_julian_year(y)).position);
            r_min = std::min(r_min, r / kAuKm);
        }
        std::printf("info: %s comes within %.1f au of the Sun\n", name, r_min);
        const std::string what = std::string(name) + " stays beyond 75 au";
        check(r_min > 75.0, what.c_str(), r_min);
    }

    // The classical belt: Arrokoth cold (near-circular, in the plane), Makemake hot ("up to
    // 30 degrees"), both between the 3:2 (39.4 au) and 2:1 (47.8 au).
    const Elements arrokoth = elements(scene.icrf_state_at(scene.find("Arrokoth"), 0.0), kGmSun);
    const Elements makemake = elements(scene.icrf_state_at(scene.find("Makemake"), 0.0), kGmSun);
    check(arrokoth.inc * kRadToDeg < 5.0 && arrokoth.e < 0.1, "Arrokoth cold classical", arrokoth.inc * kRadToDeg);
    check(makemake.inc * kRadToDeg > 25.0 && makemake.inc * kRadToDeg < 30.0, "Makemake inclined ~29 deg",
          makemake.inc * kRadToDeg);
    for (const Elements& el : {arrokoth, makemake}) {
        check(el.a / kAuKm > 39.5 && el.a / kAuKm < 47.8, "classical between 3:2 and 2:1", el.a / kAuKm);
    }
}

// The moons: [ELEM] and [PR24] put each orbit where the fits say.
void test_moons()
{
    Scene& scene = tno_scene();
    const double t = tdb_from_julian_year(2026.0);
    scene.update(t);
    auto normal = [&](const char* moon, const char* host) {
        const State m = scene.icrf_state_at(scene.find(moon), t);
        const State h = scene.icrf_state_at(scene.find(host), t);
        return glm::normalize(glm::cross(m.position - h.position, m.velocity - h.velocity));
    };
    auto deg = [](const glm::dvec3& a, const glm::dvec3& b) {
        return std::acos(std::clamp(glm::dot(a, b), -1.0, 1.0)) * kRadToDeg;
    };
    // [PR24] Table 2: the orbit poles of Namaka (RA 292.1, Dec -0.6) and Hi'iaka (283.00,
    // -10.24) from the ecliptic elements, and Hi'iaka close to the occultation pole of
    // Haumea's template (285.1, -10.6): 1.0 deg from the fitted equator, which is itself
    // 2 deg from that pole.
    const double namaka = deg(normal("Namaka", "Haumea"), unit_toward(292.1, -0.6));
    const double hiiaka = deg(normal("Hi'iaka", "Haumea"), unit_toward(283.00, -10.24));
    const double equator = deg(normal("Hi'iaka", "Haumea"), unit_toward(285.1, -10.6));
    check(namaka < 0.3, "Namaka's orbit pole as fitted (deg)", namaka);
    check(hiiaka < 0.3, "Hi'iaka's orbit pole as fitted (deg)", hiiaka);
    check(equator < 3.0, "Hi'iaka orbits near Haumea's equator (deg)", equator);

    // Pluto and Charon 19,600 km apart, on either side of the barycenter; Hydra in
    // Pluto's equator.
    const glm::dvec3 pluto = scene.icrf_state_at(scene.find("Pluto"), t).position;
    const glm::dvec3 charon = scene.icrf_state_at(scene.find("Charon"), t).position;
    const glm::dvec3 bary = scene.icrf_state_at(scene.find("Pluto system"), t).position;
    check(std::abs(glm::length(charon - pluto) - 19600.0) < 1.0, "Pluto-Charon distance", glm::length(charon - pluto));
    check(glm::length(pluto * 0.8912352640 + charon * 0.1087647360 - bary) < 1.0, "barycenter between them");
    check(deg(normal("Hydra", "Pluto system"), unit_toward(132.993, -6.163)) < 0.5, "Hydra in Pluto's equator");
}

// Ecliptic longitude (deg, 0..360) of an ICRF direction.
double ecliptic_longitude(const glm::dvec3& icrf)
{
    const double eps = 84381.448 / 3600.0 * kDegToRad;
    const double y = std::cos(eps) * icrf.y + std::sin(eps) * icrf.z;
    return wrap_deg(std::atan2(y, icrf.x) * kRadToDeg);
}

int belt_index(const Scene& scene, const char* name)
{
    for (size_t k = 0; k < scene.belts.size(); ++k) {
        if (scene.belts[k].name == name) {
            return static_cast<int>(k);
        }
    }
    check(false, (std::string("belt ") + name + " exists").c_str());
    return 0;
}

// The known objects ([BELT]): what the captions say about them in 2026, and which events
// show which family.
void test_belts()
{
    Scene& scene = tno_scene();
    const double t = tdb_from_julian_year(2026.0);
    const double neptune = ecliptic_longitude(scene.icrf_state_at(scene.find("Neptune"), t).position);
    auto from_neptune = [&](const glm::dvec3& p) { return wrap_deg(ecliptic_longitude(p) - neptune + 180.0) - 180.0; };

    // "The 443 objects known to share Pluto's 3:2 resonance ... The 272 now inside 36 au,
    // near perihelion, all lie more than 20 degrees from Neptune's direction ... None is
    // within 12 au of Neptune."
    const SceneBelt& plutinos = scene.belts[static_cast<size_t>(belt_index(scene, "Plutinos"))];
    const glm::dvec3 neptune_at = scene.icrf_state_at(scene.find("Neptune"), t).position;
    double closest = 360.0;
    double nearest_au = 1e9;
    int inner = 0;
    for (size_t k = 0; k < plutinos.size(); ++k) {
        const glm::dvec3 p = belt_position(plutinos, k, t);
        nearest_au = std::min(nearest_au, glm::length(p - neptune_at) / kAuKm);
        if (glm::length(p) < 36.0 * kAuKm) {
            ++inner;
            closest = std::min(closest, std::abs(from_neptune(p)));
        }
    }
    check(plutinos.size() == 443 && inner == 272, "443 plutinos, 272 inside 36 au (the caption)", inner);
    check(closest > 20.0, "no plutino near perihelion within 20 deg of Neptune", closest);
    check(nearest_au > 12.0, "no plutino within 12 au of Neptune", nearest_au);

    // "Of the 33 Trojans known, 28 lead Neptune."
    const SceneBelt& trojans = scene.belts[static_cast<size_t>(belt_index(scene, "Neptune Trojans"))];
    int leading = 0;
    for (size_t k = 0; k < trojans.size(); ++k) {
        leading += from_neptune(belt_position(trojans, k, t)) > 0.0 ? 1 : 0;
    }
    check(trojans.size() == 33 && leading == 28, "33 Trojans, 28 leading (the caption)", leading);
    // "In the 2000s and 2010s the trailing cloud lay in front of the crowded center of the
    // Milky Way": L5, 60 deg behind Neptune, within 25 deg of where the Galactic plane
    // crosses the ecliptic toward the center (ecliptic longitude ~266 deg).
    for (double year = 2000.0; year < 2020.0; year += 1.0) {
        const double l5 = ecliptic_longitude(scene.icrf_state_at(scene.find("Neptune"), tdb_from_julian_year(year)).position) - 60.0;
        const double off = std::abs(wrap_deg(l5 - 266.0 + 180.0) - 180.0);
        check(off < 25.0, "L5 toward the Galactic center in 2000-2019", off);
    }

    // "Few are known toward Taurus and Gemini" (ecliptic longitude 90-120 deg, where the
    // Galactic anticenter crosses it): fewer than 40% of the average per 30 deg.
    int total = 0;
    int gap = 0;
    int scattered = 0;
    int inside = 0;
    for (size_t b = 0; b < scene.belts.size(); ++b) {
        const SceneBelt& belt = scene.belts[b];
        const bool far_family = belt.name == "Scattered disk" || belt.name == "Detached objects";
        for (size_t k = 0; k < belt.size(); ++k) {
            const glm::dvec3 p = belt_position(belt, k, t);
            const double lon = ecliptic_longitude(p);
            ++total;
            gap += lon >= 90.0 && lon < 120.0 ? 1 : 0;
            if (far_family) {
                ++scattered;
                inside += glm::length(p) < belt.elements[k * SceneBelt::kStride] * kAuKm ? 1 : 0;
            }
        }
    }
    const double gap_share = gap / (total / 12.0);
    std::printf("info: %d known objects; toward Taurus-Gemini %.2f of the average; %d of %d scattered or detached "
                "inside their semi-major axis\n",
                total, gap_share, inside, scattered);
    check(gap_share < 0.4, "few known toward Taurus and Gemini", gap_share);
    check(inside > 0.9 * scattered, "nearly every scattered or detached object inside a", inside / double(scattered));

    // Events choose their families; without an event, all are drawn.
    int event = -1;
    for (size_t k = 0; k < scene.events.size(); ++k) {
        event = scene.events[k].name == "Every known plutino" ? static_cast<int>(k) : event;
    }
    scene.set_active_event(event);
    const bool only_plutinos = scene.belt_shown(static_cast<size_t>(belt_index(scene, "Plutinos"))) &&
                               !scene.belt_shown(static_cast<size_t>(belt_index(scene, "Classical belt")));
    scene.set_active_event(-1);
    check(event >= 0 && only_plutinos, "the plutino event shows the plutinos alone");
    check(scene.belt_shown(static_cast<size_t>(belt_index(scene, "Classical belt"))), "all belts without an event");
}

} // namespace

void run_trans_neptunian_tests()
{
    test_horizons();
    test_resonances();
    test_orbits();
    test_moons();
    test_belts();
}
