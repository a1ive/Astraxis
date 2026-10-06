// Tests for body appearance: comet tails, atmospheres, plumes.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "scene/comet.hpp"
#include "scene/scene_loader.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace astraxis;

namespace {

// Comets in the ISEE-3 scene ([SBDB], [SLV], [ESO] in isee3.toml). The ion tail points
// along the solar wind as the comet sees it: a few degrees off the anti-sun direction,
// toward where the comet came from. ICE flew through Giacobini-Zinner's tail 7,800 km
// behind the nucleus (Bame et al. 1986, Science 232, 356): at closest approach it should
// be downstream and within the tail's radius (1.2e4 km across, [SLV]) of the axis.
void test_comet_tails()
{
    Scene scene = load_scene_or_die("isee3.toml");
    const int gz = scene.find("Giacobini-Zinner");
    const int halley = scene.find("Halley");
    const int craft = scene.find("ISEE-3");
    const Body::Comet& gz_comet = scene.bodies[static_cast<size_t>(gz)].comet;
    const Body::Comet& halley_comet = scene.bodies[static_cast<size_t>(halley)].comet;
    check(gz_comet.enabled && halley_comet.enabled, "isee3.toml comets");
    check(gz_comet.m1 == 13.2 && gz_comet.k1 == 8.25 && halley_comet.m1 == 5.5 && halley_comet.k1 == 8.0,
          "comet magnitude laws [SBDB]");
    check(std::abs(comet_total_magnitude(5.5, 8.0, 0.1) - (5.5 - 8.0)) < 1e-12, "total magnitude law");

    auto utc = [](const CalendarDateTime& c) { return tdb_from_utc_jd(jd_from_calendar(c)); };
    const double t_ca = utc({1985, 9, 11, 11, 2, 43});
    const State sun = scene.sun_icrf_state_at(t_ca);
    const State comet = scene.icrf_state_at(gz, t_ca);
    const glm::dvec3 r = comet.position - sun.position;
    const glm::dvec3 v = comet.velocity - sun.velocity;
    const glm::dvec3 axis = ion_tail_direction(r, v, gz_comet.solar_wind_km_s);

    // Aberration: atan(v_perp / v_wind) away from the anti-sun direction, against the motion.
    const glm::dvec3 out = glm::normalize(r);
    const glm::dvec3 v_perp = v - glm::dot(v, out) * out;
    const double expected = std::atan2(glm::length(v_perp), gz_comet.solar_wind_km_s - glm::dot(v, out));
    const double aberration = std::acos(std::clamp(glm::dot(axis, out), -1.0, 1.0));
    check(std::abs(aberration - expected) < 1e-9, "ion tail aberration angle", aberration / kDegToRad);
    check(glm::dot(axis, v_perp) < 0.0, "ion tail lags behind the motion");

    const glm::dvec3 rel = scene.icrf_state_at(craft, t_ca).position - comet.position;
    const double downstream = glm::dot(rel, axis);
    const double off_axis = glm::length(rel - downstream * axis);
    check(downstream > 0.0, "ICE passed on the tail side (km)", downstream);
    check(off_axis < 0.5 * gz_comet.ion_tail_diameter_km, "ICE within the ion tail's radius (km)", off_axis);
    const double off_anti_sun = glm::length(rel - glm::dot(rel, out) * out);
    std::printf("info: G-Z tail %.2f deg from anti-sun; ICE %.0f km downstream, %.0f km off its axis (%.0f km off the "
                "anti-sun line)\n",
                aberration / kDegToRad, downstream, off_axis, off_anti_sun);
}

// Finson-Probstein dust (comet.hpp). A grain of radiation pressure ratio beta, released
// without ejection speed, drifts from the beta = 0 one by the radiation pressure
// acceleration: beta GM / r^2 outward, 1/2 a t^2 after a short time. Halley's tail on
// 1986-03-10: grains on the anti-sun side, curving back along the orbit (lagging the
// ion tail, which only aberrates by v / 400 km/s).
void test_dust_tail()
{
    const double gm = 1.32712440041279419e11; // [GM] in isee3.toml
    const State comet{glm::dvec3(0.6 * kAuKm, 0.0, 0.0), glm::dvec3(0.0, 50.0, 5.0)};
    const double dt = 6.0 * 3600.0;
    const glm::dvec3 drift = dust_grain_state(comet, glm::dvec3(0.0), 0.5, gm, dt).position -
                             dust_grain_state(comet, glm::dvec3(0.0), 0.0, gm, dt).position;
    const double expected = 0.5 * 0.5 * gm / (0.6 * kAuKm * 0.6 * kAuKm) * dt * dt;
    check(std::abs(glm::length(drift) / expected - 1.0) < 1e-3, "dust drift by radiation pressure (km)",
          glm::length(drift));
    check(glm::dot(glm::normalize(drift), glm::dvec3(1.0, 0.0, 0.0)) > 0.9999, "dust drift away from the sun");

    Scene scene = load_scene_or_die("isee3.toml");
    const int halley = scene.find("Halley");
    const double t = tdb_from_utc_jd(jd_from_calendar({1986, 3, 10, 0, 0, 0}));
    DustTail tail;
    tail.update(scene, halley, t);
    const auto& grains = tail.grains();
    const State sun = scene.sun_icrf_state_at(t);
    const State c = scene.icrf_state_at(halley, t);
    const glm::dvec3 out = glm::normalize(c.position - sun.position);
    const glm::dvec3 v = c.velocity - sun.velocity;
    const glm::dvec3 v_perp = glm::normalize(v - glm::dot(v, out) * out);
    glm::dvec3 mean(0.0);
    double total = 0.0;
    double far = 0.0;
    for (const DustGrain& g : grains) {
        mean += g.weight * (g.position - c.position);
        total += g.weight;
        far = std::max(far, glm::length(g.position - c.position));
    }
    mean /= total;
    const double anti_sun = glm::dot(mean, out);
    const double lag = glm::dot(mean, v_perp);
    const glm::dvec3 ion = ion_tail_direction(c.position - sun.position, v, 400.0);
    check(grains.size() > DustTail::kSlots * DustTail::kBetas * DustTail::kDirections * 9 / 10, "Halley dust grains",
          static_cast<double>(grains.size()));
    check(anti_sun > 0.0, "Halley dust tail on the anti-sun side (km)", anti_sun);
    check(lag / anti_sun < glm::dot(ion, v_perp) / glm::dot(ion, out), "Halley dust tail lags the ion tail");
    // Cost of a frame: the next update a minute later (cached slots, warm-started anomalies).
    const auto clock_start = std::chrono::steady_clock::now();
    for (int k = 1; k <= 20; ++k) {
        tail.update(scene, halley, t + 60.0 * k);
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - clock_start).count();
    std::printf("info: dust tail update %.2f ms per frame\n", ms / 20.0);
    // Another output in the same frame: the grains are kept, not propagated again.
    const size_t kept = tail.grains().size();
    const auto repeat_start = std::chrono::steady_clock::now();
    tail.update(scene, halley, t + 60.0 * 20);
    const double repeat_ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - repeat_start).count();
    check(tail.grains().size() == kept, "dust tail kept at the same time");
    std::printf("info: dust tail update at the same time %.4f ms\n", repeat_ms);
    // Giacobini-Zinner at ICE's pass: also a tail on the anti-sun side.
    {
        const int gz = scene.find("Giacobini-Zinner");
        const double t_gz = tdb_from_utc_jd(jd_from_calendar({1985, 9, 11, 11, 0, 0}));
        DustTail gz_tail;
        gz_tail.update(scene, gz, t_gz);
        const State gs = scene.sun_icrf_state_at(t_gz);
        const State gc = scene.icrf_state_at(gz, t_gz);
        glm::dvec3 gz_mean(0.0);
        double gz_total = 0.0;
        for (const DustGrain& g : gz_tail.grains()) {
            gz_mean += g.weight * (g.position - gc.position);
            gz_total += g.weight;
        }
        gz_mean /= gz_total;
        const double gz_anti_sun = glm::dot(gz_mean, glm::normalize(gc.position - gs.position));
        check(gz_anti_sun > 0.8 * glm::length(gz_mean), "Giacobini-Zinner dust tail on the anti-sun side (km)",
              gz_anti_sun);
    }
    std::printf("info: Halley dust 1986-03-10: %zu grains, weighted center %.3g km anti-sun, %.3g km along -v "
                "(%.1f deg behind the anti-sun line), farthest %.3g km\n",
                grains.size(), anti_sun, -lag, std::atan2(-lag, anti_sun) / kDegToRad, far);
}

// Masses shown in the info panel (body_gm_km3_s2).
// Atmospheres: per-channel optical depths from the 550 nm values, Earth's
// against Bruneton's Rayleigh formula, Pluto's haze phase function against
// the New Horizons phase ratios, the same parameters in every scene, and the
// loader's checks.
void test_atmospheres()
{
    const Scene scene = load_scene_or_die("solar_system.toml");
    auto atmosphere = [&](const char* name) -> const Atmosphere& {
        return scene.bodies[static_cast<size_t>(scene.find(name))].atmosphere;
    };
    for (const char* name : {"Earth", "Venus", "Mars", "Titan", "Pluto", "Triton"}) {
        check(atmosphere(name).enabled, name);
    }
    check(!atmosphere("Moon").enabled && !atmosphere("Jupiter").enabled, "no atmosphere on the Moon and Jupiter");

    // Bruneton 2017 (demo.cc): beta_R = 1.24062e-6 lambda^-4 m^-1 (lambda in um), H_R = 8 km.
    const Atmosphere& earth = atmosphere("Earth");
    double worst = 0.0;
    for (int c = 0; c < 3; ++c) {
        const double um = kChannelWavelengthsNm[c] / 1000.0;
        const double expected = 1.24062e-6 * std::pow(um, -4.0) * 8000.0;
        worst = std::max(worst, std::abs(earth.rayleigh_depth[c] / expected - 1.0));
    }
    check(worst < 1e-3, "Earth Rayleigh depths match Bruneton's formula", worst);
    check(std::abs(earth.haze_depth.r - earth.haze_depth.b) < 1e-12, "Earth aerosols are gray (alpha = 0)");

    // Titan: opacity ~ lambda^-2.34 above 80 km (Tomasko et al. 2008).
    const Atmosphere& titan = atmosphere("Titan");
    const double ratio = titan.haze_depth.b / titan.haze_depth.r;
    check(std::abs(ratio - std::pow(440.0 / 680.0, -2.34)) < 1e-9, "Titan haze Angstrom exponent", ratio);
    // Triton: scattering optical depth ~ lambda^-2 (Ohno et al. 2021), haze to ~30 km inside the shell.
    const Atmosphere& triton = atmosphere("Triton");
    check(std::abs(triton.haze_depth.b / triton.haze_depth.r - std::pow(440.0 / 680.0, -2.0)) < 1e-9 &&
              triton.height_km >= 30.0,
          "Triton haze");
    check(titan.deck_altitude_km == 80.0 && atmosphere("Venus").deck_altitude_km == 74.0, "deck altitudes");

    // Pluto: Henyey-Greenstein ratio between 167 and 20 deg phase (scattering angles
    // 13 and 160 deg) within the haze I/F ratios of Cheng et al. 2017, Table 4
    // (peak: 0.3 / 0.02 = 15; at 45 km: 0.15 / 0.004 = 37).
    const double g = atmosphere("Pluto").haze_g;
    auto hg = [g](double theta_deg) {
        const double d = 1.0 + g * g - 2.0 * g * std::cos(theta_deg * kDegToRad);
        return (1.0 - g * g) / (d * std::sqrt(d));
    };
    const double pluto_ratio = hg(13.0) / hg(160.0);
    check(pluto_ratio > 15.0 && pluto_ratio < 37.0, "Pluto haze forward/back ratio", pluto_ratio);

    // Scenes that repeat a body carry the same atmosphere.
    for (const char* file : {"earth_moon.toml", "jwst.toml", "mercury.toml", "parker.toml", "saturn.toml"}) {
        const Scene other = load_scene_or_die(file);
        for (const Body& b : other.bodies) {
            const int k = scene.find(b.name);
            if (k < 0 || !scene.bodies[static_cast<size_t>(k)].atmosphere.enabled) {
                continue;
            }
            const Atmosphere& a = b.atmosphere;
            const Atmosphere& ref = scene.bodies[static_cast<size_t>(k)].atmosphere;
            check(a.enabled && a.height_km == ref.height_km && a.deck_altitude_km == ref.deck_altitude_km &&
                      a.rayleigh_depth == ref.rayleigh_depth && a.haze_depth == ref.haze_depth &&
                      a.haze_albedo == ref.haze_albedo && a.haze_g == ref.haze_g &&
                      a.haze_scale_height_km == ref.haze_scale_height_km &&
                      b.texture == scene.bodies[static_cast<size_t>(k)].texture,
                  (std::string(file) + ": " + b.name + " atmosphere as in solar_system.toml").c_str());
        }
    }

    std::string error;
    Scene bad;
    const char* head = R"(name = 'x'
[[bodies]]
)";
    auto bad_scene = [&](const char* body) {
        return !load_scene_string(std::string(head) + body, "bad_atmosphere", bad, &error);
    };
    check(bad_scene(R"(name = 'S'
kind = 'star'
radii_km = [1.0]
[bodies.atmosphere]
height_km = 1.0
)"),
          "rejects an atmosphere on a star");
    check(bad_scene(R"(name = 'P'
radii_km = [1.0]
[bodies.atmosphere]
height_km = 1.0
haze_depth = 0.1
)") && error.find("haze_scale_height_km") != std::string::npos,
          "requires a haze scale height");
    check(bad_scene(R"(name = 'P'
radii_km = [1.0]
[bodies.atmosphere]
height_km = 1.0
haze_albedo = [0.5, 1.5, 0.5]
)"),
          "rejects a haze albedo above 1");
}

// Plumes: vents from the IAU Gazetteer (west longitudes converted to east),
// the tiger stripes 30-36 km apart and perpendicular to the line through their
// centers, surface points on the ellipsoid, the same plumes in every scene,
// and the loader's checks.
void test_plumes()
{
    const Scene scene = load_scene_or_die("solar_system.toml");
    auto body = [&](const char* name) -> const Body& { return scene.bodies[static_cast<size_t>(scene.find(name))]; };
    auto plume = [&](const Body& b, const char* name) -> const Plume* {
        for (const Plume& p : b.plumes) {
            if (p.name == name) {
                return &p;
            }
        }
        return nullptr;
    };

    const Plume* pele = plume(body("Io"), "Pele");
    const Plume* tvashtar = plume(body("Io"), "Tvashtar");
    check(pele && pele->type == Plume::Type::Umbrella && std::abs(pele->lat_lon_deg.y - (360.0 - 255.28)) < 1e-9 &&
              pele->height_km == 300.0 && pele->radius_km == 600.0,
          "Pele: 255.28 W, 300 km high, 1,200 km across");
    check(tvashtar && std::abs(tvashtar->lat_lon_deg.y - (360.0 - 123.53)) < 1e-9 &&
              tvashtar->radius_km == 2.0 * tvashtar->height_km,
          "Tvashtar: 123.53 W, ballistic reach 2 H");

    // Tiger stripes: midpoints near the Gazetteer centers, lengths, spacing and orientation.
    const Body& enceladus = body("Enceladus");
    const char* stripes[4] = {"Alexandria Sulcus", "Cairo Sulcus", "Baghdad Sulcus", "Damascus Sulcus"};
    const double centers[4][2] = {{-75.63, 360.0 - 137.56}, {-81.62, 360.0 - 154.48},
                                  {-86.91, 360.0 - 230.54}, {-80.59, 360.0 - 285.87}};
    const double lengths[4] = {111.0, 165.0, 176.0, 125.0};
    glm::dvec3 mid[4];
    glm::dvec3 along[4];
    bool ok = true;
    for (int k = 0; k < 4; ++k) {
        const Plume* p = plume(enceladus, stripes[k]);
        if (!p || p->type != Plume::Type::Jets) {
            ok = false;
            continue;
        }
        const glm::dvec3 a = body_surface_point(enceladus, p->lat_lon_deg.x, p->lat_lon_deg.y);
        const glm::dvec3 b = body_surface_point(enceladus, p->end_lat_lon_deg.x, p->end_lat_lon_deg.y);
        // The chord's midpoint, raised to the surface (the sagitta is up to 15 km).
        mid[k] = glm::normalize(a + b) * glm::length(body_surface_point(enceladus, centers[k][0], centers[k][1]));
        along[k] = glm::normalize(b - a);
        const glm::dvec3 c = body_surface_point(enceladus, centers[k][0], centers[k][1]);
        ok = ok && glm::length(mid[k] - c) < 3.0 && std::abs(glm::length(b - a) / lengths[k] - 1.0) < 0.1;
    }
    check(ok, "tiger stripes centered on the Gazetteer centers, with their lengths");
    double worst_spacing = 0.0;
    double worst_angle = 0.0;
    for (int k = 0; k + 1 < 4; ++k) {
        const glm::dvec3 step = mid[k + 1] - mid[k];
        worst_spacing = std::max(worst_spacing, std::abs(glm::length(step) - 33.0));
        worst_angle = std::max(worst_angle, std::abs(glm::dot(glm::normalize(step), along[k])));
    }
    check(worst_spacing < 5.0, "tiger stripes ~35 km apart (km off)", worst_spacing);
    check(worst_angle < 0.2, "tiger stripes perpendicular to the line of centers (cos)", worst_angle);

    const Body& triton = body("Triton");
    const Plume* hili = plume(triton, "Hili");
    const Plume* mahilani = plume(triton, "Mahilani");
    check(hili && mahilani && hili->type == Plume::Type::Geyser && hili->height_km == 8.0 &&
              hili->tail_azimuth_deg == 270.0 && mahilani->lat_lon_deg == glm::dvec2(-50.5, 359.5),
          "Triton's geysers: 8 km, tails to the west");

    // Surface points lie on the ellipsoid; normals are unit and outward.
    double worst_surface = 0.0;
    for (const Body* b : {&body("Io"), &enceladus, &triton}) {
        for (double lat = -80.0; lat <= 80.0; lat += 40.0) {
            for (double lon = 0.0; lon < 360.0; lon += 60.0) {
                const glm::dvec3 p = body_surface_point(*b, lat, lon);
                const glm::dvec3 q = p / glm::dvec3(b->equatorial_radius_km, b->equatorial_radius_b_km, b->polar_radius_km);
                const glm::dvec3 n = body_surface_normal(*b, p);
                worst_surface = std::max({worst_surface, std::abs(glm::dot(q, q) - 1.0),
                                          std::abs(glm::length(n) - 1.0), glm::dot(n, glm::normalize(p)) > 0.99 ? 0.0 : 1.0});
            }
        }
    }
    check(worst_surface < 1e-12, "surface points and normals", worst_surface);

    for (const char* file : {"jupiter.toml", "saturn.toml"}) {
        const Scene other = load_scene_or_die(file);
        for (const Body& b : other.bodies) {
            const int k = scene.find(b.name);
            if (k < 0 || scene.bodies[static_cast<size_t>(k)].plumes.empty()) {
                continue;
            }
            const std::vector<Plume>& ref = scene.bodies[static_cast<size_t>(k)].plumes;
            bool same = b.plumes.size() == ref.size();
            for (size_t i = 0; same && i < ref.size(); ++i) {
                const Plume& x = b.plumes[i];
                const Plume& y = ref[i];
                same = x.name == y.name && x.type == y.type && x.lat_lon_deg == y.lat_lon_deg &&
                       x.end_lat_lon_deg == y.end_lat_lon_deg && x.count == y.count && x.height_km == y.height_km &&
                       x.radius_km == y.radius_km && x.spread_deg == y.spread_deg &&
                       x.optical_depth == y.optical_depth && x.albedo == y.albedo && x.g == y.g;
            }
            check(same, (std::string(file) + ": " + b.name + " plumes as in solar_system.toml").c_str());
        }
    }

    std::string error;
    Scene bad;
    const char* head = R"(name = 'x'
[[bodies]]
name = 'P'
radii_km = [100.0]
[[bodies.plumes]]
name = 'v'
lat_lon_deg = [0.0, 0.0]
height_km = 10.0
optical_depth = 0.1
)";
    auto bad_scene = [&](const char* rest) {
        return !load_scene_string(std::string(head) + rest, "bad_plume", bad, &error);
    };
    check(bad_scene("type = 'fountain'\n"), "rejects an unknown plume type");
    check(bad_scene("type = 'jets'\nend_lat_lon_deg = [1.0, 1.0]\ncount = 13\nspread_deg = 5.0\n"),
          "rejects too many jets");
    check(bad_scene("type = 'geyser'\nradius_km = 0.5\ntail_km = 10.0\n") &&
              error.find("tail_azimuth_deg") != std::string::npos,
          "requires a geyser's tail azimuth");
    check(!bad_scene("type = 'umbrella'\n") && bad.bodies[0].plumes[0].radius_km == 20.0,
          "umbrella radius defaults to 2 H");
}

} // namespace

void run_appearance_tests()
{
    test_comet_tails();
    test_dust_tail();
    test_atmospheres();
    test_plumes();
}
