// Tests for body appearance: comet tails, atmospheres, plumes, night lights and aurora.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "scene/aurora.hpp"
#include "scene/aurora_model.hpp"
#include "scene/jovian_plasma_sheet.hpp"
#include "scene/satellite_aurora.hpp"
#include "scene/saturn_aurora.hpp"
#include "scene/mars_aurora.hpp"
#include "scene/magnetic_field.hpp"
#include "scene/magnetosphere.hpp"
#include "scene/comet.hpp"
#include "scene/corona.hpp"
#include "scene/relief.hpp"
#include "scene/scene_loader.hpp"
#include "scene/simulation.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <string>
#include <tuple>
#include <utility>
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
    // Siding Spring in the Mars scene: a child of Mars (the root), whose grains still
    // move under the Sun's gravity, not Mars'.
    {
        Scene mars = load_scene_or_die("mars.toml");
        const int ss = mars.find("Siding Spring");
        const double t_ss = tdb_from_utc_jd(jd_from_calendar({2014, 10, 19, 18, 0, 0}));
        DustTail ss_tail;
        ss_tail.update(mars, ss, t_ss);
        const State ss_sun = mars.sun_icrf_state_at(t_ss);
        const State mc = mars.icrf_state_at(ss, t_ss);
        glm::dvec3 ss_mean(0.0);
        double ss_total = 0.0;
        for (const DustGrain& g : ss_tail.grains()) {
            ss_mean += g.weight * (g.position - mc.position);
            ss_total += g.weight;
        }
        ss_mean /= ss_total;
        const double ss_anti_sun = glm::dot(ss_mean, glm::normalize(mc.position - ss_sun.position));
        check(!ss_tail.grains().empty() && ss_anti_sun > 0.8 * glm::length(ss_mean),
              "Siding Spring dust tail on the anti-sun side, in a Mars-centered scene (km)", ss_anti_sun);
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
    for (const char* file :
         {"didymos.toml", "earth_moon.toml", "mercury.toml", "parker.toml", "saturn.toml", "space_telescopes.toml",
          "venus.toml"}) {
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
                      a.haze_scale_height_km == ref.haze_scale_height_km && a.refractivity == ref.refractivity &&
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


// Earth's night lights and aurora: the same in every scene with the Earth;
// the Kp record (GFZ) at the great storms; Starkov's ovals (Sigernes et al.
// 2011); the hemispheric power, the IGRF dipole, the line yields of Mathieu
// et al. (2025) and the emission they make; the simulation's Kp window.
void test_night_and_aurora()
{
    const Scene scene = load_scene_or_die("solar_system.toml");
    const Body& earth = scene.bodies[static_cast<size_t>(scene.find("Earth"))];
    check(earth.aurora.kp != nullptr && earth.night_lights.gain > 0.0 && !earth.night_lights.texture.empty(),
          "Earth has night lights and an aurora");
    for (const char* file : {"didymos.toml", "earth_moon.toml", "isee3.toml", "mercury.toml", "parker.toml",
                             "rosetta.toml", "space_telescopes.toml", "venus.toml"}) {
        const Scene other = load_scene_or_die(file);
        const Body& e = other.bodies[static_cast<size_t>(other.find("Earth"))];
        check(e.aurora.kp != nullptr && e.night_lights.texture == earth.night_lights.texture &&
                  e.night_lights.color == earth.night_lights.color && e.night_lights.gain == earth.night_lights.gain,
              (std::string(file) + ": Earth's night lights and aurora as in solar_system.toml").c_str());
    }
    for (const char* name : {"Mars", "Venus", "Moon"}) {
        const Body& b = scene.bodies[static_cast<size_t>(scene.find(name))];
        check(!b.aurora.kp && b.night_lights.texture.empty(), name);
    }

    // The Kp record: 1932 on, 3-hour UT intervals.
    const KpRecord& kp = *earth.aurora.kp;
    auto utc = [](int y, int mo, int d, int h) { return tdb_from_utc_jd(jd_from_calendar({y, mo, d, h, 0, 0})); };
    check(std::abs(kp.start_tdb() - utc(1932, 1, 1, 0)) < 1e-6, "Kp record starts 1932-01-01 00 UT");
    check(kp.end_tdb() > utc(2026, 9, 1, 0), "Kp record reaches 2026-09");
    // Kp 9o: 1989-03-13 21-24 UT (the Quebec blackout), 2003-10-29 06-09 UT,
    // 2024-05-11 00-03 UT (the Gannon storm).
    check(kp.at(utc(1989, 3, 13, 22)) == 9.0 && kp.at(utc(2003, 10, 29, 7)) == 9.0 &&
              kp.at(utc(2024, 5, 11, 1)) == 9.0,
          "Kp 9 in the great storms");
    check(kp.at(utc(2024, 5, 10, 22)) == 26.0 / 3.0, "Kp 9- on 2024-05-10 21-24 UT", kp.at(utc(2024, 5, 10, 22)));
    check(kp.median() == 2.0, "median Kp 2o", kp.median());
    check(kp.at(utc(1900, 1, 1, 0)) == kp.median() && kp.at(utc(2100, 1, 1, 0)) == kp.median(),
          "the median outside the record");
    const double t0 = utc(2024, 5, 10, 21);
    check(std::abs(kp.mean(t0 + 600.0, t0 + 3600.0) - kp.at(t0 + 60.0)) < 1e-9, "mean within one interval");
    const double two = 0.5 * (kp.at(t0 + 60.0) + kp.at(t0 + 3.0 * 3600.0 + 60.0));
    check(std::abs(kp.mean(t0, t0 + 6.0 * 3600.0) - two) < 1e-6, "mean over two intervals",
          kp.mean(t0, t0 + 21600.0));
    check(std::abs(kp.mean(utc(2101, 1, 1, 0), utc(2102, 1, 1, 0)) - kp.median()) < 1e-12, "mean beyond the record");

    // Starkov's ovals: |AL| from Kp (Sigernes eq. 1); poleward < equatorward
    // <= diffuse everywhere; they grow with Kp, mostly on the night side
    // (antisunward), and the oval is broader at midnight than at noon.
    check(std::abs(starkov_al_nt(3.0) - (18.0 - 36.9 + 244.8 - 54.0)) < 1e-9, "|AL| from Kp");
    bool ordered = true;
    for (double k = 0.0; k <= 9.0; k += 0.5) {
        for (double t = 0.0; t < 24.0; t += 0.5) {
            const AuroralBoundaries b = starkov_boundaries(k, t);
            ordered = ordered && b.poleward_deg < b.equatorward_deg && b.equatorward_deg < b.diffuse_deg + 1e-9;
        }
    }
    check(ordered, "poleward < equatorward <= diffuse boundary");
    const AuroralBoundaries quiet = starkov_boundaries(1.0, 0.0);
    const AuroralBoundaries storm = starkov_boundaries(8.0, 0.0);
    check(storm.equatorward_deg > quiet.equatorward_deg + 5.0 && storm.diffuse_deg > quiet.diffuse_deg + 5.0,
          "the night-side oval moves equatorward with Kp", storm.equatorward_deg - quiet.equatorward_deg);
    const AuroralBoundaries midnight = starkov_boundaries(3.0, 0.0);
    const AuroralBoundaries noon = starkov_boundaries(3.0, 12.0);
    check(midnight.equatorward_deg > noon.equatorward_deg + 5.0, "the oval is displaced antisunward");
    check(midnight.equatorward_deg - midnight.poleward_deg > noon.equatorward_deg - noon.poleward_deg,
          "Kp 3 oval broader at midnight than at noon [Sigernes Fig. 1]");

    // IGRF-14 dipole: the north geomagnetic pole from (g10, g11, h11) at the
    // epochs, and the 2025-30 secular variation after 2025 (held after 2030).
    auto pole_at = [](double year) { return geomagnetic_north_pole(tdb_from_julian_year(year)); };
    auto lat_lon = [](const glm::dvec3& p) {
        return glm::dvec2(std::asin(p.z) * kRadToDeg, std::atan2(p.y, p.x) * kRadToDeg);
    };
    const glm::dvec3 p2020 = glm::normalize(glm::dvec3(1451.37, -4653.35, 29403.41));
    check(glm::length(pole_at(2020.0) - p2020) < 1e-9, "IGRF 2020 dipole pole");
    const glm::dvec2 ll2025 = lat_lon(pole_at(2025.0));
    check(std::abs(ll2025.x - 80.79) < 0.01 && std::abs(ll2025.y + 72.76) < 0.01, "geomagnetic pole 2025.0 (deg)",
          ll2025.x);
    const glm::dvec3 p2030 = glm::normalize(glm::dvec3(1410.3 - 50.0, -4545.5 + 107.5, 29350.0 - 63.0));
    check(glm::length(pole_at(2030.0) - p2030) < 1e-9 && glm::length(pole_at(2040.0) - p2030) < 1e-12,
          "secular variation to 2030, then held");
    check(lat_lon(pole_at(1950.0)).x < lat_lon(pole_at(2025.0)).x - 2.0, "the pole drifted north since 1950");

    // Longyearbyen (78.2 N, 16.0 E) lies under the Kp-3 oval near magnetic
    // noon at 08:50 UT on 2009-12-24 [Sigernes Fig. 1]. Their corrected
    // geomagnetic pole puts noon there then; the centered dipole about an
    // hour later.
    {
        Scene s = load_scene_or_die("solar_system.toml");
        const double t = utc(2009, 12, 24, 8) + 50.0 * 60.0;
        s.update(t);
        const Body& e = s.bodies[static_cast<size_t>(s.find("Earth"))];
        const glm::dmat3 to_body = glm::transpose(e.orientation);
        const glm::dvec3 sun = glm::normalize(to_body * (s.sun_position() - e.world_position));
        const glm::dvec3 pole = geomagnetic_north_pole(t);
        const double lat = 78.2 * kDegToRad;
        const double lon = 16.0 * kDegToRad;
        const glm::dvec3 site(std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat));
        const double mlt = magnetic_local_time(site, pole, sun);
        check(std::abs(mlt - 12.0) < 1.0, "Longyearbyen near magnetic noon (h)", mlt);
        const double colat = std::acos(glm::dot(site, pole)) * kRadToDeg;
        const AuroralBoundaries b = starkov_boundaries(3.0, mlt);
        check(colat > b.poleward_deg && colat < b.equatorward_deg, "Longyearbyen under the Kp-3 oval", colat);
        const glm::dvec3 dusk_site = glm::normalize(glm::cross(pole, sun));
        check(std::abs(magnetic_local_time(dusk_site, pole, sun) - 18.0) < 1e-9, "18 h on the dusk side");
    }

    // Electron hemispheric power (Emery et al. 2008).
    check(std::abs(electron_hemispheric_power_gw(0.0) - 2.5103) < 1e-12 &&
              std::abs(electron_hemispheric_power_gw(3.0) - (2.5103 + 29.7621 - 19.6425 + 8.1891)) < 1e-9,
          "hemispheric power fit");
    bool rising = true;
    for (double k = 0.0; k < 9.0; k += 0.25) {
        rising = rising && electron_hemispheric_power_gw(k + 0.25) > electron_hemispheric_power_gw(k);
    }
    check(rising, "hemispheric power grows with Kp");

    // Line yields: Table 2 at 1 and 10 keV, in between in log energy.
    const AuroraLineYields y1 = aurora_line_yields(1.0);
    const AuroraLineYields y10 = aurora_line_yields(10.0);
    check(std::abs(y1.green - 0.325) < 1e-12 && std::abs(y1.red - 0.525) < 1e-12 &&
              std::abs(y1.blue - 0.072) < 1e-12 && std::abs(y10.green - 0.445) < 1e-12 &&
              std::abs(y10.red - 0.082) < 1e-12 && std::abs(y10.violet - 0.066) < 1e-12,
          "line yields [Mathieu Table 2]");
    const AuroraLineYields y3 = aurora_line_yields(std::sqrt(10.0));
    check(std::abs(y3.green - std::sqrt(0.325 * 0.445)) < 1e-12, "geometric mean halfway in log energy");

    // The emission: the hemispheric power over the oval (and the diffuse
    // aurora) gives a few erg cm^-2 s^-1 and about a kR of green light at
    // Kp 3, growing with Kp; green outshines red (in luminance).
    const AuroraEmission e0 = aurora_emission(0.0);
    const AuroraEmission e3 = aurora_emission(3.0);
    const AuroraEmission e9 = aurora_emission(9.0);
    check(e3.oval_flux > 1.0 && e3.oval_flux < 5.0, "Kp 3 oval energy flux (erg cm^-2 s^-1)", e3.oval_flux);
    check(e0.oval_flux < e3.oval_flux && e3.oval_flux < e9.oval_flux, "oval flux grows with Kp");
    const double green_kr = kGreenOverBlue * e3.blue_kr;
    check(green_kr > 0.3 && green_kr < 3.0, "Kp 3 green column (kR)", green_kr);
    auto luma = [](const glm::dvec3& c) { return 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b; };
    const double red_kr = e3.blue_kr * red_over_blue(e3.blue_kr);
    check(e3.low_rgb.g > e3.low_rgb.r && luma(e3.low_rgb) > 3.0 * luma(red_kr * kRedPerKr), "green dominates");
    check(std::abs(red_over_blue(0.1) - 2.0) < 1e-12 && std::abs(red_over_blue(3.0) - 0.4) < 1e-12 &&
              red_over_blue(0.01) == 2.0 && red_over_blue(10.0) == red_over_blue(3.0),
          "I630/I427.8 from 2 to 0.4 [Dashkevich 2006]");
    check(std::abs(e3.oval[0].y - starkov_boundaries(3.0, 0.0).equatorward_deg * kDegToRad) < 1e-6 &&
              std::abs(e3.oval[24].z - starkov_boundaries(3.0, 12.0).poleward_deg * kDegToRad) < 1e-6,
          "oval samples at midnight and noon");

    // The simulation averages Kp over a window that grows with the warp.
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/earth_moon.toml", &error), "earth_moon loads");
    sim.stop_tour();
    const auto storm_event = std::find_if(sim.scene().events.begin(), sim.scene().events.end(),
        [](const SceneEvent& e) { return e.name == "Geomagnetic storm: Kp 9 (2024-05-11)"; });
    check(storm_event != sim.scene().events.end() && storm_event->aurora_visible && !storm_event->aurora_ultraviolet,
          "geomagnetic storm tour enables visible aurora");
    check(!sim.scene().events.front().aurora_visible && !sim.scene().events.front().aurora_ultraviolet,
          "ordinary Earth-Moon event leaves aurora off");
    sim.clock().t_tdb = utc(2024, 5, 11, 1) + 30.0 * 60.0;
    sim.clock().warp = 60.0;
    sim.update(0.0);
    check(sim.aurora_kp() > 8.5, "Kp near 9 during the Gannon storm at low warp", sim.aurora_kp());
    sim.clock().warp = 86400.0 * 30.0;
    sim.update(0.0);
    check(sim.aurora_kp() < 4.0, "averaged down at a month per second", sim.aurora_kp());

    std::string err;
    Scene bad;
    check(!load_scene_string("name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [1.0]\n"
                             "[bodies.night_lights]\ntexture = 'a.jpg'\ngain = 1.0\n",
                             "bad_lights", bad, &err),
          "night lights need the body's texture");
    check(!load_scene_string("name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [1.0]\n"
                             "[bodies.aurora]\nkp = 'geomag/missing.bin'\n",
                             "bad_aurora", bad, &err),
          "rejects a missing Kp record");
}

// Relief maps and the sea's glint (scene/relief.hpp, tools/textures/relief_map.py).
void test_relief()
{
    // The encoding: zero at mid-scale, odd, the square of the offset.
    check(relief_decode(0.5) == 0.0 && relief_decode(1.0) == kReliefMaxSlope && relief_decode(0.0) == -kReliefMaxSlope,
          "relief encoding: its ends and middle");
    check(std::abs(relief_decode(0.75) - 0.25 * kReliefMaxSlope) < 1e-12 &&
              std::abs(relief_decode(0.25) + 0.25 * kReliefMaxSlope) < 1e-12,
          "relief encoding: the square of the offset");

    // Ground that rises toward the east faces west, by the angle of the slope.
    const glm::dvec3 pole(0.0, 0.0, 1.0);
    const glm::dvec3 n = glm::normalize(glm::dvec3(1.0, 0.3, 0.4));
    const glm::dvec3 east = glm::normalize(glm::cross(pole, n));
    const glm::dvec3 north = glm::cross(n, east);
    const glm::dvec3 tilted = relief_normal(n, pole, 0.2, 0.0);
    check(std::abs(glm::dot(tilted, east) + std::sin(std::atan(0.2))) < 1e-12 && std::abs(glm::dot(tilted, north)) < 1e-12,
          "a slope up to the east tilts the normal west", glm::dot(tilted, east));
    const glm::dvec3 tilted_north = relief_normal(n, pole, 0.0, -0.5);
    check(std::abs(glm::dot(tilted_north, north) - std::sin(std::atan(0.5))) < 1e-12,
          "a slope down to the north tilts it north", glm::dot(tilted_north, north));
    check(glm::length(relief_normal(n, pole, 0.0, 0.0) - n) < 1e-15 && relief_normal(pole, pole, 0.3, 0.3) == pole,
          "level ground and the pole keep the smooth normal");
    // Near the terminator a slope of a few degrees is most of the light: the sun 3 deg up in
    // the east, a wall of 10 deg facing it and one facing away.
    const glm::dvec3 sun = glm::normalize(std::sin(glm::radians(3.0)) * n + std::cos(glm::radians(3.0)) * east);
    const double flat = glm::dot(n, sun);
    const double toward = glm::dot(relief_normal(n, pole, -std::tan(glm::radians(10.0)), 0.0), sun);
    const double away = glm::dot(relief_normal(n, pole, std::tan(glm::radians(10.0)), 0.0), sun);
    check(std::abs(toward / flat - std::sin(glm::radians(13.0)) / std::sin(glm::radians(3.0))) < 1e-9 && away < 0.0,
          "a 10 deg wall facing the low sun is 4 times brighter, the far one is dark", toward / flat);

    // The sea. Cox & Munk: 0.003 + 0.00512 U; water reflects 2% head-on.
    check(std::abs(sea_slope_variance(10.0) - 0.0542) < 1e-12 && sea_slope_variance(0.0) == 0.003,
          "Cox-Munk mean square slope");
    const double f0 = fresnel_normal_reflectance(1.333);
    check(std::abs(f0 - 0.02037) < 2e-5, "water's reflectance at normal incidence", f0);
    const double variance = sea_slope_variance(7.0);
    const glm::dvec3 up(0.0, 0.0, 1.0);
    // Sun overhead, seen from overhead: every facet that is level mirrors it, F / (4 variance).
    const double peak = sea_glint(up, up, up, variance, f0);
    check(std::abs(peak - f0 / (4.0 * variance)) < 1e-9 * peak, "the glint's peak with the sun overhead", peak);
    check(peak > 0.1 && peak < 0.2, "comparable to the sea's own albedo", peak);
    // What the sea reflects in all, with the sun 30 deg from the zenith: sum of
    // f mu_view over the sky = glint / (pi mu_sun) mu_view. It is the Fresnel reflectance
    // there, a little over the 2% of normal incidence: no light made or lost.
    const glm::dvec3 sun30(std::sin(glm::radians(30.0)), 0.0, std::cos(glm::radians(30.0)));
    double reflected = 0.0;
    double brightest = 0.0;
    double brightest_zenith_deg = 0.0;
    const int steps = 400;
    for (int i = 0; i < steps; ++i) {
        const double mu = (i + 0.5) / steps;
        const double sine = std::sqrt(1.0 - mu * mu);
        for (int j = 0; j < steps; ++j) {
            const double phi = kTwoPi * (j + 0.5) / steps;
            const glm::dvec3 v(sine * std::cos(phi), sine * std::sin(phi), mu);
            const double g = sea_glint(up, sun30, v, variance, f0);
            reflected += g / (kPi * sun30.z) * mu * (kTwoPi / steps / steps);
            if (g > brightest) {
                brightest = g;
                brightest_zenith_deg = glm::degrees(std::atan2(v.x, v.z));
            }
        }
    }
    check(reflected > f0 && reflected < 1.15 * f0, "the sea reflects its Fresnel share of the sunlight", reflected / f0);
    check(std::abs(brightest_zenith_deg + 30.0) < 4.0, "brightest near the mirror direction", brightest_zenith_deg);
    // Reciprocity of the reflectance (the glint is pi f mu_sun), and nothing from below.
    const glm::dvec3 v60(-std::sin(glm::radians(60.0)) * 0.8, std::sin(glm::radians(60.0)) * 0.6, std::cos(glm::radians(60.0)));
    check(std::abs(sea_glint(up, sun30, v60, variance, f0) / sun30.z - sea_glint(up, v60, sun30, variance, f0) / v60.z) <
              1e-12,
          "the sea's reflectance is reciprocal");
    check(sea_glint(up, -sun30, v60, variance, f0) == 0.0 && sea_glint(up, sun30, -v60, variance, f0) == 0.0,
          "no glint with the sun or the viewer below the horizon");
    // A calmer sea: a smaller, brighter glint.
    check(sea_glint(up, up, up, sea_slope_variance(2.0), f0) > 2.5 * peak, "a calm sea's glint is brighter at its center");

    // Scenes: the Earth's relief and sea, and the airless or exposed terrain maps.
    const Scene solar = load_scene_or_die("solar_system.toml");
    const Body& earth = solar.bodies[static_cast<size_t>(solar.find("Earth"))];
    check(!earth.relief.texture.empty() && earth.relief.exaggeration == 1.0 && earth.sea.wind_m_s == 7.0 &&
              earth.sea.refractive_index == 1.333,
          "the Earth has relief and a sea");
    for (const char* name : {"Moon", "Mars"}) {
        const Body& b = solar.bodies[static_cast<size_t>(solar.find(name))];
        check(!b.relief.texture.empty() && b.relief.exaggeration == 1.0 && b.sea.wind_m_s < 0.0, name);
    }
    auto check_dem = [&](const char* name, const Scene& system, double left_lon) {
        const Body& b = solar.bodies[static_cast<size_t>(solar.find(name))];
        const Body& other = system.bodies[static_cast<size_t>(system.find(name))];
        check(!b.relief.texture.empty() && b.relief.texture == other.relief.texture &&
                  b.relief.exaggeration == 1.0 && b.sea.wind_m_s < 0.0 &&
                  b.texture_left_lon_deg == left_lon && !b.texture_west_positive,
              (b.name + ": relief follows its template and albedo coordinates").c_str());
        const TestImage relief = load_test_image(b.relief.texture.c_str());
        const TestImage albedo = load_test_image(b.texture.c_str());
        check(relief.width > 0 && relief.width == albedo.width && relief.height == albedo.height,
              (b.name + ": the relief map matches the albedo map").c_str());
        return relief;
    };
    const TestImage mercury = check_dem("Mercury", load_scene_or_die("mercury.toml"), -180.0);
    const TestImage venus = check_dem("Venus", load_scene_or_die("venus.toml"), -180.0);
    const Scene saturn = load_scene_or_die("saturn.toml");
    const TestImage enceladus = check_dem("Enceladus", saturn, 0.0);
    const TestImage mimas = check_dem("Mimas", saturn, -180.0);
    const TestImage tethys = check_dem("Tethys", saturn, -180.0);
    const TestImage dione = check_dem("Dione", saturn, -180.0);
    const TestImage rhea = check_dem("Rhea", saturn, -180.0);
    const Scene earth_moon = load_scene_or_die("earth_moon.toml");
    check(earth_moon.bodies[static_cast<size_t>(earth_moon.find("Moon"))].relief.texture ==
                  solar.bodies[static_cast<size_t>(solar.find("Moon"))].relief.texture &&
              earth_moon.bodies[static_cast<size_t>(earth_moon.find("Earth"))].sea.wind_m_s == earth.sea.wind_m_s,
          "the Earth-Moon scene has the same relief and sea");
    check(!load_scene_or_die("mars.toml").bodies[0].relief.texture.empty(), "the Mars scene has its relief");
    const Scene kerbol = load_scene_or_die("kerbol.toml");
    int with_relief = 0;
    int with_sea = 0;
    for (const Body& b : kerbol.bodies) {
        if (b.relief.texture.empty()) {
            continue;
        }
        ++with_relief;
        const TestImage map = load_test_image(b.relief.texture.c_str());
        const TestImage albedo = load_test_image(b.texture.c_str());
        check(map.width > 0 && map.width == albedo.width && map.height == albedo.height,
              (b.name + ": the relief map matches the albedo map").c_str());
        if (b.sea.wind_m_s >= 0.0) {
            ++with_sea;
        }
    }
    check(with_relief == 15 && with_sea == 3, "every Kerbol world but Jool has relief; Kerbin, Eve and Laythe have seas",
          with_relief);
    check(kerbol.bodies[static_cast<size_t>(kerbol.find("Jool"))].relief.texture.empty() &&
              kerbol.bodies[static_cast<size_t>(kerbol.find("Kerbin"))].sea.wind_m_s >= 0.0 &&
              kerbol.bodies[static_cast<size_t>(kerbol.find("Mun"))].sea.wind_m_s < 0.0,
          "Jool has none, Kerbin a sea, the Mun none");

    // The baked maps against known places. Slopes in the map's sense: red up to the east,
    // green up to the north.
    auto slope = [](const TestImage& map, double lat, double lon, int c) { return relief_decode(map.at(lat, lon, c)); };
    if (mercury.width > 0) {
        // Debussy (34.0 S, 347.5 E): the east-west profile crosses both inner walls.
        check(slope(mercury, -33.9, 347.27, 0) < -0.02 && slope(mercury, -33.9, 347.93, 0) > 0.02,
              "Debussy's inner walls face the crater floor");
    }
    if (venus.width > 0) {
        check(slope(venus, 65.0, 0.5, 0) > 0.06 && slope(venus, 65.0, 5.5, 0) < -0.02,
              "Maxwell Montes rises from Lakshmi Planum and drops to the east");
        check(std::abs(slope(venus, 89.0, 0.0, 0)) < 0.001 && std::abs(slope(venus, 89.0, 0.0, 1)) < 0.001,
              "the missing Venus polar data leaves the surface normal smooth");
    }
    if (enceladus.width > 0) {
        // TestImage takes longitude relative to a 180 W left edge. Enceladus starts at 0 E.
        check(slope(enceladus, 56.8, 345.5 - 180.0, 0) < -0.1 && slope(enceladus, 56.8, 348.5 - 180.0, 0) > 0.1,
              "the trough on the east side of Ali Baba slopes toward its floor");
    }
    if (mimas.width > 0) {
        // Herschel (1.7 N, 111.8 W). The central peak lies between these inner walls.
        check(slope(mimas, 1.7, -130.8, 0) < -0.3 && slope(mimas, 1.7, -100.8, 0) > 0.25,
              "Herschel's inner walls face the crater floor");
    }
    if (tethys.width > 0) {
        check(slope(tethys, 32.8, -154.9, 0) < -0.1 && slope(tethys, 32.8, -106.9, 0) > 0.05,
              "the west and east inner slopes of Odysseus face the basin");
    }
    if (dione.width > 0) {
        // Evander lies south of the equatorial radius tile.
        check(slope(dione, -57.0, -145.0, 0) < -0.18 && slope(dione, -57.0, -141.0, 0) > 0.05,
              "opposing slopes near Evander are present in the south polar tile");
    }
    if (rhea.width > 0) {
        check(slope(rhea, 34.2, -154.7, 0) < -0.1 && slope(rhea, 34.2, -148.7, 0) > 0.06,
              "the floor of Tirawa lies between opposing slopes");
    }
    // Moon: Tycho (43.3 S, 11.2 W; 85 km across, 4.7 km deep), its inner walls face the center.
    const TestImage moon = load_test_image("textures/moon_relief.jpg");
    check(moon.width == 2048 && moon.height == 1024, "the Moon's relief map");
    if (moon.width > 0) {
        check(slope(moon, -43.3, -10.0, 0) > 0.1 && slope(moon, -43.3, -12.4, 0) < -0.1,
              "Tycho's east wall rises to the east, its west wall to the west", slope(moon, -43.3, -10.0, 0));
        check(slope(moon, -42.2, -11.2, 1) > 0.1 && slope(moon, -44.4, -11.2, 1) < -0.1,
              "its north wall rises to the north, its south wall to the south", slope(moon, -42.2, -11.2, 1));
        // Mare Imbrium is nearly level.
        check(std::abs(slope(moon, 35.0, -17.0, 0)) < 0.02 && std::abs(slope(moon, 35.0, -17.0, 1)) < 0.02,
              "Mare Imbrium is level", slope(moon, 35.0, -17.0, 0));
        check(moon.at(-43.3, -11.2, 2) < 0.02, "no water on the Moon");
    }
    // Mars: Olympus Mons (18.65 N, 226.2 E), 600 km across: its flanks rise toward the summit;
    // the north wall of Valles Marineris at Melas (9 S, 287 E) drops to the south; Hellas'
    // floor is level.
    const TestImage mars = load_test_image("textures/mars_relief.jpg");
    check(mars.width == 2048 && mars.height == 1024, "Mars' relief map");
    if (mars.width > 0) {
        check(slope(mars, 18.65, 223.5, 0) > 0.03 && slope(mars, 18.65, 229.0, 0) < -0.03,
              "Olympus Mons' flanks rise toward the summit, east and west", slope(mars, 18.65, 223.5, 0));
        check(slope(mars, 16.0, 226.2, 1) > 0.03 && slope(mars, 20.0, 226.2, 1) < -0.03,
              "and north and south", slope(mars, 16.0, 226.2, 1));
        check(std::abs(slope(mars, -42.4, 70.5, 0)) < 0.01, "Hellas' floor is level", slope(mars, -42.4, 70.5, 0));
        double steepest = 0.0;
        for (double lat = -8.0; lat >= -10.0; lat -= 0.1) {
            steepest = std::max(steepest, slope(mars, lat, 287.0, 1));
        }
        check(steepest > 0.08, "the north wall of Valles Marineris rises to the north", steepest);
    }
    // Earth: water where the sea and the great lakes are, not below-sea-level land, ice
    // shelves or sea ice; the Himalaya's front rises to the north.
    const TestImage globe = load_test_image("textures/earth_relief.jpg");
    check(globe.width == 2048 && globe.height == 1024, "the Earth's relief map");
    if (globe.width > 0) {
        check(globe.at(0.0, -150.0, 2) > 0.98 && globe.at(40.0, -40.0, 2) > 0.98 && globe.at(-40.0, 80.0, 2) > 0.98,
              "the open oceans are water");
        check(globe.at(23.0, 10.0, 2) < 0.02 && globe.at(50.0, 90.0, 2) < 0.02 && globe.at(-25.0, 135.0, 2) < 0.02,
              "the Sahara, Siberia and Australia are not");
        check(globe.at(47.7, -87.5, 2) > 0.9 && globe.at(42.0, 51.0, 2) > 0.9 && globe.at(-1.0, 33.0, 2) > 0.9,
              "Lake Superior, the Caspian Sea and Lake Victoria are water", globe.at(47.7, -87.5, 2));
        check(globe.at(29.9, 27.4, 2) < 0.05, "the Qattara depression, below sea level, is dry", globe.at(29.9, 27.4, 2));
        check(globe.at(-81.0, -175.0, 2) < 0.05 && globe.at(-78.0, -50.0, 2) < 0.05,
              "the Ross and Ronne ice shelves are not open water", globe.at(-81.0, -175.0, 2));
        check(globe.at(0.0, -150.0, 0) == globe.at(40.0, -40.0, 0) && std::abs(slope(globe, 0.0, -150.0, 0)) < 1e-4,
              "the sea is level");
        double steepest = 0.0;
        for (double lat = 26.5; lat <= 28.5; lat += 0.1) {
            steepest = std::max(steepest, slope(globe, lat, 86.0, 1));
        }
        check(steepest > 0.02, "the Himalaya's front rises to the north", steepest);
        check(std::abs(slope(globe, 33.0, 88.0, 1)) < 0.02, "Tibet is a plateau", slope(globe, 33.0, 88.0, 1));
    }
    // Kerbin: the invented sea covers about two thirds, and is level; the land is not.
    const TestImage kerbin = load_test_image("textures/kerbin_relief.jpg");
    if (kerbin.width > 0) {
        double water = 0.0, weight = 0.0, rough = 0.0, land = 0.0;
        for (int row = 0; row < 180; ++row) {
            const double lat = 89.5 - row;
            const double w = std::cos(glm::radians(lat));
            for (int col = 0; col < 360; ++col) {
                const double share = kerbin.at(lat, col - 179.5, 2);
                water += w * share;
                weight += w;
                if (share < 0.02) {
                    rough += w * std::hypot(slope(kerbin, lat, col - 179.5, 0), slope(kerbin, lat, col - 179.5, 1));
                    land += w;
                }
            }
        }
        check(water / weight > 0.55 && water / weight < 0.68, "Kerbin's open water", water / weight);
        check(rough / land > 0.02 && rough / land < 0.2, "Kerbin's land has slopes of a few degrees", rough / land);
    }

    // The Earth map shows a clear Sahara and a cloudy Southern Ocean.
    check(!earth.clouds.texture.empty() && earth.clouds.albedo > 0.5 && earth.clouds.albedo <= 0.9,
          "the Earth has clouds", earth.clouds.albedo);
    check(earth_moon.bodies[static_cast<size_t>(earth_moon.find("Earth"))].clouds.texture == earth.clouds.texture,
          "in the Earth-Moon scene too");
    check(solar.bodies[static_cast<size_t>(solar.find("Mars"))].clouds.texture.empty(), "Mars has no cloud map");
    for (const Body& b : kerbol.bodies) {
        const bool cloudy = b.name == "Kerbin" || b.name == "Eve" || b.name == "Laythe";
        check(cloudy == !b.clouds.texture.empty(), (b.name + ": procedural cloud assignment").c_str());
        if (!cloudy) {
            continue;
        }
        const TestImage clouds = load_test_image(b.clouds.texture.c_str());
        const TestImage albedo = load_test_image(b.texture.c_str());
        check(b.atmosphere.enabled && b.atmosphere.deck_altitude_km == 0.0 && b.clouds.albedo > 0.0 &&
                  b.clouds.tropical_drift_deg_per_day != 0.0,
              (b.name + ": clouds use the transparent atmosphere and drift").c_str());
        check(clouds.width > 0 && clouds.width == albedo.width && clouds.height == albedo.height,
              (b.name + ": the cloud map matches the albedo map").c_str());
        const auto [clear, overcast] = std::minmax_element(clouds.rgb.begin(), clouds.rgb.end());
        check(clear != clouds.rgb.end() && *clear < 32 && *overcast > 223,
              (b.name + ": the map contains clear sky and thick cloud").c_str());
    }
    const TestImage cloud_map = load_test_image(earth.clouds.texture.c_str());
    check(cloud_map.width == globe.width && cloud_map.height == globe.height && cloud_map.width > 0,
          "the cloud map matches the Earth's other maps");
    if (cloud_map.width > 0) {
        auto mean_cover = [&](double lat0, double lat1, double lon0, double lon1) {
            double sum = 0.0;
            int count = 0;
            for (double lat = lat0; lat <= lat1; lat += 0.5) {
                for (double lon = lon0; lon <= lon1; lon += 0.5) {
                    sum += cloud_map.at(lat, lon, 0);
                    ++count;
                }
            }
            return sum / count;
        };
        const double sahara = mean_cover(18.0, 28.0, 0.0, 25.0);
        const double southern_ocean = mean_cover(-60.0, -50.0, -179.5, 179.5);
        check(sahara < 0.03, "the Sahara is clear", sahara);
        check(southern_ocean > 0.3 && southern_ocean < 0.6, "the Southern Ocean is cloudy", southern_ocean);
    }
    check(cloud_glint_share(0.0) == 1.0 && cloud_glint_share(1.0) == 0.0 && cloud_lights_share(1.0) == 0.0 &&
              std::abs(cloud_glint_share(0.5) - 0.25) < 1e-15 && cloud_lights_share(0.5) == 0.5,
          "half a cloud cover leaves a quarter of the glint and half of the city lights");

    // The drift: west in the tropics, east in the middle latitudes, none at the poles; the two
    // copies of the picture always add up to one, and each is unseen when it starts over.
    check(earth.clouds.tropical_drift_deg_per_day < 0.0 && earth.clouds.midlatitude_drift_deg_per_day > 0.0,
          "the Earth's clouds drift");
    check(cloud_drift_deg_per_day(earth.clouds, 0.0) == earth.clouds.tropical_drift_deg_per_day &&
              std::abs(cloud_drift_deg_per_day(earth.clouds, glm::radians(45.0)) -
                       (earth.clouds.midlatitude_drift_deg_per_day + earth.clouds.tropical_drift_deg_per_day / 8.0)) < 1e-12 &&
              std::abs(cloud_drift_deg_per_day(earth.clouds, glm::radians(90.0))) < 1e-12 &&
              cloud_drift_deg_per_day(earth.clouds, glm::radians(-45.0)) == cloud_drift_deg_per_day(earth.clouds, glm::radians(45.0)),
          "the drift by latitude");
    double largest_shift = 0.0;
    bool weights_ok = true;
    for (int k = 0; k <= 200; ++k) {
        const double phase = k / 200.0;
        const CloudCopy a = cloud_copy(earth.clouds, glm::radians(45.0), phase, 0);
        const CloudCopy b = cloud_copy(earth.clouds, glm::radians(45.0), phase, 1);
        weights_ok = weights_ok && std::abs(a.weight + b.weight - 1.0) < 1e-12 && a.weight >= 0.0 && b.weight >= 0.0;
        largest_shift = std::max({largest_shift, std::abs(a.shift_deg), std::abs(b.shift_deg)});
    }
    check(weights_ok, "the two copies of the clouds add up to one");
    check(cloud_copy(earth.clouds, 0.7, 0.0, 0).weight == 0.0 && cloud_copy(earth.clouds, 0.7, 0.5, 1).weight == 0.0 &&
              cloud_copy(earth.clouds, 0.7, 0.5, 0).weight == 1.0 && cloud_copy(earth.clouds, 0.7, 0.5, 0).shift_deg == 0.0,
          "a copy starts over unseen, and is the still picture at its fullest");
    check(largest_shift < 20.0, "no copy is sheared far from the picture", largest_shift);
    // A copy moves at the drift's rate.
    const CloudCopy early = cloud_copy(earth.clouds, glm::radians(45.0), 0.4, 0);
    const CloudCopy late = cloud_copy(earth.clouds, glm::radians(45.0), 0.6, 0);
    check(std::abs((late.shift_deg - early.shift_deg) / (0.2 * kCloudCycleDays) -
                   cloud_drift_deg_per_day(earth.clouds, glm::radians(45.0))) < 1e-9,
          "a copy moves at the drift's rate");
    {
        // The cycle follows the clock, but slowly on screen whatever the warp; not while paused,
        // and backward when time runs backward.
        Simulation sim;
        std::string load_error;
        check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/earth_moon.toml", &load_error), "the Earth-Moon scene loads");
        sim.clock().paused = false;
        sim.clock().reverse = false;
        sim.clock().warp = 3600.0;
        const double start = sim.cloud_phase();
        for (int k = 0; k < 10; ++k) {
            sim.update(0.1);
        }
        const double slow = sim.cloud_phase() - start;
        check(std::abs(slow - 3600.0 / 86400.0 / kCloudCycleDays) < 1e-9, "an hour of drift per second at an hour per second", slow);
        sim.clock().warp = 1e7;
        const double before = sim.cloud_phase();
        for (int k = 0; k < 10; ++k) {
            sim.update(0.1);
        }
        const double fast = sim.cloud_phase() - before;
        check(std::abs(fast - kCloudMaxDaysPerSecond / kCloudCycleDays) < 1e-9, "capped at high warp", fast);
        sim.clock().reverse = true;
        const double turned = sim.cloud_phase();
        sim.update(0.1);
        check(sim.cloud_phase() < turned, "backward when time runs backward");
        sim.clock().paused = true;
        const double paused = sim.cloud_phase();
        sim.update(0.1);
        check(sim.cloud_phase() == paused, "still while paused");
    }

    std::string err;
    Scene bad;
    const char* cloudy = "name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [1.0]\ntexture = 'a.jpg'\n";
    const char* air = "[bodies.atmosphere]\nheight_km = 0.01\nrayleigh_depth = 0.1\nrayleigh_scale_height_km = 0.001\n";
    check(!load_scene_string("name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [1.0]\n"
                             "[bodies.clouds]\ntexture = 'c.jpg'\nalbedo = 0.8\n",
                             "bad_clouds", bad, &err),
          "clouds need the body's texture");
    check(!load_scene_string((std::string(cloudy) + air + "[bodies.clouds]\ntexture = 'c.jpg'\nalbedo = 1.5\n").c_str(),
                             "bad_cloud_albedo", bad, &err),
          "rejects a cloud albedo above 1");
    check(load_scene_string((std::string(cloudy) + air + "[bodies.clouds]\ntexture = 'c.jpg'\nalbedo = 0.7\n").c_str(),
                            "good_clouds", bad, &err) &&
              bad.bodies[0].clouds.albedo == 0.7,
          "reads a cloud layer");
    check(bad.bodies[0].clouds.tropical_drift_deg_per_day == 0.0 && bad.bodies[0].clouds.midlatitude_drift_deg_per_day == 0.0,
          "which stays still unless told to drift");
    check(!load_scene_string((std::string(cloudy) + "[bodies.clouds]\ntexture = 'c.jpg'\nalbedo = 0.7\n"
                                                    "midlatitude_drift_deg_per_day = 90.0\n").c_str(),
                             "bad_cloud_drift", bad, &err),
          "rejects a drift that would tear the picture apart");
    check(!load_scene_string("name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [1.0]\n"
                             "[bodies.relief]\ntexture = 'a.jpg'\n",
                             "bad_relief", bad, &err),
          "a relief map needs the body's texture");
    check(!load_scene_string("name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [1.0]\ntexture = 'a.jpg'\n"
                             "[bodies.sea]\nwind_m_s = 7.0\nrefractive_index = 1.333\n",
                             "bad_sea", bad, &err),
          "a sea needs the relief map");
    check(!load_scene_string("name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [1.0]\ntexture = 'a.jpg'\n"
                             "[bodies.relief]\ntexture = 'b.jpg'\n"
                             "[bodies.sea]\nwind_m_s = -1.0\nrefractive_index = 1.333\n",
                             "bad_wind", bad, &err),
          "rejects a negative wind");
    check(!load_scene_string("name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [1.0]\ntexture = 'a.jpg'\n"
                             "[bodies.relief]\ntexture = 'b.jpg'\nexaggeration = 0.0\n",
                             "bad_exaggeration", bad, &err),
          "rejects a zero exaggeration");
    check(load_scene_string("name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [1.0]\ntexture = 'a.jpg'\n"
                            "[bodies.relief]\ntexture = 'b.jpg'\nexaggeration = 2.0\n"
                            "[bodies.sea]\nwind_m_s = 3.0\nrefractive_index = 1.3\n",
                            "good_relief", bad, &err) &&
              bad.bodies[0].relief.exaggeration == 2.0 && bad.bodies[0].sea.wind_m_s == 3.0,
          "reads a relief map and a sea");
}

void test_refraction()
{
    const Scene solar = load_scene_or_die("solar_system.toml");
    const Body& earth = solar.bodies[static_cast<size_t>(solar.find("Earth"))];
    const Atmosphere& a = earth.atmosphere;
    const double radius = earth.equatorial_radius_km;
    check(a.refraction != nullptr, "Earth template builds the refraction table");
    if (!a.refraction) return;
    const double sigma2 = 1.0 / (0.55 * 0.55);
    const double ciddor = 1e-8 * (5792105.0 / (238.0185 - sigma2) + 167917.0 / (57.362 - sigma2));
    check(std::abs(a.refractivity - ciddor) < 1e-13, "Ciddor standard air at 550 nm (NIST A22/A32)");

    // Independent eikonal ODE in Cartesian coordinates, integrated in arc
    // length from the tangent point. The product integrates in radius instead.
    using State = std::array<double, 7>; // position xy, direction xy, optical depth RGB
    auto derivative = [&](const State& s) {
        const glm::dvec2 p(s[0], s[1]), v(s[2], s[3]);
        const double r = glm::length(p), h = r - radius;
        const double gas = std::exp(-h / a.rayleigh_scale_height_km);
        const double nu = a.refractivity * gas;
        const glm::dvec2 gradient = -p / r * (nu / ((1.0 + nu) * a.rayleigh_scale_height_km));
        const glm::dvec2 dv = gradient - v * glm::dot(v, gradient);
        const glm::dvec3 extinction = a.rayleigh_depth * (gas / a.rayleigh_scale_height_km) +
            a.haze_depth * (std::exp(-h / a.haze_scale_height_km) / a.haze_scale_height_km);
        return State{v.x, v.y, dv.x, dv.y, extinction.x, extinction.y, extinction.z};
    };
    auto advance = [](State s, const State& delta, double step) {
        for (size_t i = 0; i < s.size(); ++i) s[i] += step * delta[i];
        return s;
    };
    for (double height : {0.0, 3.0, 12.0, 40.0}) {
        State s{radius + height, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
        for (int step = 0; step < 3000; ++step) {
            const auto k1 = derivative(s);
            const auto k2 = derivative(advance(s, k1, 0.5));
            const auto k3 = derivative(advance(s, k2, 0.5));
            const auto k4 = derivative(advance(s, k3, 1.0));
            for (size_t i = 0; i < s.size(); ++i) s[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
        const auto ray = atmosphere_refracted_ray(a, radius, height);
        check(std::abs(ray.bending_rad + 2.0 * std::asin(s[2])) < 2e-8,
              "radial bending integral agrees with independent Cartesian ray propagation");
        check(std::abs(ray.impact_km - (s[0] * s[3] - s[1] * s[2])) < 1e-6,
              "curved ray preserves asymptotic impact parameter");
        for (int c = 0; c < 3; ++c) {
            check(std::abs(ray.transmission[c] - std::exp(-2.0 * s[static_cast<size_t>(4 + c)])) < 2e-6,
                  "curved-ray extinction agrees with independent path integral");
        }
    }
    const auto low = atmosphere_refracted_ray(a, radius, 3.0);
    check(low.transmission.r > 5.0 * low.transmission.g && low.transmission.g > 5.0 * low.transmission.b,
          "low tangent rays lose blue and green more strongly than red");
    check(low.bending_rad > 0.01 && low.bending_rad < 0.025, "low rays bend enough to reach the lunar umbra");
    Atmosphere vacuum = a;
    vacuum.refractivity = 0.0;
    vacuum.rayleigh_depth = vacuum.haze_depth = glm::dvec3(0.0);
    const auto straight = atmosphere_refracted_ray(vacuum, radius, 0.0);
    check(straight.bending_rad == 0.0 && straight.transmission == glm::dvec3(1.0),
          "vacuum leaves a ray straight and unattenuated");
    for (size_t i = 1; i < a.refraction->samples.size(); ++i) {
        const auto before = a.refraction->samples[i - 1], after = a.refraction->samples[i];
        check(after.w < before.w && glm::all(glm::greaterThanEqual(glm::vec3(after), glm::vec3(before))),
              "higher rays bend less and transmit more light");
    }

    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/solar_eclipse.toml", &error), "lunar eclipse scene loads");
    sim.stop_tour();
    sim.jump_to_event(2);
    sim.update(0.0);
    OutputView output;
    sim.compute_view(1280, 720, 1.0, output);
    const Scene& scene = sim.scene();
    const auto& event = scene.events[2];
    const glm::dvec3 moon = scene.bodies[static_cast<size_t>(scene.find("Moon"))].world_position;
    const glm::dvec3 observer = output.camera.position;
    check(event.observer == scene.find("Moon") && event.focus == scene.find("Earth"),
          "lunar eclipse observer looks from the Moon toward Earth");
    check(glm::length(observer - moon - event.observer_offset_km) < 1e-8 &&
          std::abs(glm::length(observer - moon) - 1937.4) < 1e-6,
          "lunar viewpoint is 200 km above the Moon at the reference time");
    const glm::dvec3 sun = scene.sun_position() - observer;
    const double alpha = std::asin(kSunRadiusKm / glm::length(sun));
    const double beta = std::asin(radius / glm::length(observer));
    const double separation = std::atan2(glm::length(glm::cross(glm::normalize(sun), glm::normalize(-observer))),
                                        glm::dot(glm::normalize(sun), glm::normalize(-observer)));
    check(beta > separation + alpha, "the Moon-side observer is in Earth's umbra");
    check(std::abs(glm::length(moon) - 390651.871445) < 0.01, "lunar distance agrees with Horizons bake");
    check(std::abs(output.camera.focal_y - 1.0 / std::tan(1.75 * kDegToRad)) < 1e-5,
          "lunar eclipse applies its 3.5 degree field of view");
    Scene bad;
    const std::string base = "name='x'\n[[bodies]]\nname='Earth'\ntemplate='bodies/earth.toml'\n[bodies.atmosphere]\n";
    for (const char* invalid : {"refractivity=-1", "refractivity=nan", "refractivity=0.01", "refraction_gain=inf"}) {
        check(!load_scene_string(base + invalid, "bad_refraction", bad, &error, ASTRAXIS_ASSET_DIR), "rejects invalid or ducting refraction optics");
        check(error.find("refraction requires") != std::string::npos || error.find("must be finite and non-negative") != std::string::npos,
              "invalid refraction fails in the optics validation");
    }
    check(load_scene_string(base + "refractivity=0", "no_refraction", bad, &error, ASTRAXIS_ASSET_DIR) && !bad.bodies[0].atmosphere.refraction,
          "zero refractivity disables the table through a template override");
}

void test_corona()
{
    auto utc = [](int y, int m, int d, int h) { return tdb_from_utc_jd(jd_from_calendar({y, m, d, h, 0, 0})); };
    check(std::abs(corona_electron_density(10.0) - 3790.0) < 1e-9,
          "Leblanc density at 10 solar radii, cm^-3");
    check(corona_electron_density(0.9) == 0.0, "no coronal density inside photosphere");
    const double column = 6.6524587051e-29 * kSunRadiusKm * 1e9;
    const glm::dvec2 surface = corona_scattering(1.0, 0.0);
    const double isotropic = 0.5 * column * corona_electron_density(1.0);
    check(std::abs(surface.x / isotropic - 1.0) < 1e-5 && std::abs(surface.y / isotropic) < 1e-5,
          "uniform hemispheric illumination gives isotropic Thomson brightness");
    const glm::dvec2 distant = corona_scattering(100.0);
    const double point = 3.0 / 16.0 * column * corona_electron_density(100.0) / 10000.0;
    check(std::abs(distant.x / point - 1.0) < 2e-4 && std::abs(distant.y / point - 1.0) < 2e-4,
          "finite disk approaches the point-source Thomson kernel");
    for (double r : {1.0, 1.01, 1.5, 2.0, 10.0, 32.0}) {
        const auto c = corona_scattering(r);
        check(std::isfinite(c.x) && c.x > 0.0 && c.x + c.y > 0.0,
              "finite positive coronal scattering at the limb and throughout the volume");
    }
    check(corona_visible_fraction(1.0, 1.1, 0.0) == 0.0, "total eclipse hides all photospheric light");
    check(corona_visible_fraction(1.0, 1.0, 2.1) == 1.0, "separate disks leave all photospheric light");
    check(corona_visible_fraction(1.0, 0.9, 0.0) > 0.0, "annular eclipse keeps a bright photospheric ring");
    double previous = 0.0;
    for (int i = 0; i <= 300; ++i) {
        const double visible = corona_visible_fraction(1.0, 1.1, i / 100.0);
        check(visible >= previous - 1e-10 && visible >= 0.0 && visible <= 1.0,
              "remaining photospheric flux increases continuously through egress");
        previous = visible;
    }
    check(corona_activity(utc(2019, 12, 1, 0)) == 0.0, "corona minimum anchor");
    check(corona_activity(utc(2024, 10, 15, 0)) == 1.0, "corona maximum-phase anchor");
    check(corona_activity(utc(2017, 8, 21, 0)) < corona_activity(utc(2024, 4, 8, 0)),
          "2017 descending activity is below 2024");

    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/solar_eclipse.toml", &error), "eclipse scene loads");
    sim.stop_tour();
    check(sim.scene().active_event() == 1, "eclipse scene opens its named 2024 event");
    for (size_t i = 0; i < 2; ++i) {
        sim.jump_to_event(i);
        sim.update(0.0);
        OutputView output;
        sim.compute_view(1280, 720, 1.0, output);
        const Scene& scene = sim.scene();
        const SceneEvent& event = scene.events[i];
        const glm::dvec3 p = output.camera.position;
        check(glm::length(p - event.observer_offset_km) < 1e-8,
              "observer is at its designed Earth-relative ICRF offset");
        check(std::abs(glm::length(p) - 7378.1366) < 1e-5, "eclipse viewpoint is 1000 km above Earth");
        const glm::dvec3 sun = scene.sun_position() - p;
        const glm::dvec3 moon = scene.bodies[static_cast<size_t>(scene.find("Moon"))].world_position - p;
        const double alpha = std::asin(kSunRadiusKm / glm::length(sun));
        const double beta = std::asin(1737.4 / glm::length(moon));
        const double sep = std::atan2(glm::length(glm::cross(glm::normalize(sun), glm::normalize(moon))),
                                      glm::dot(glm::normalize(sun), glm::normalize(moon)));
        check(beta > alpha + sep, "designed observer is in the lunar umbra");
        const double expected = i == 0 ? 1.0345665 : 1.0606412;
        check(std::abs(beta / alpha - expected) < 2e-6, "apparent disk sizes agree with independent bake");
        const double focal = 1.0 / std::tan(event.fov_deg * kDegToRad * 0.5);
        check(std::abs(output.camera.focal_y / focal - 1.0) < 1e-6, "eclipse event applies its narrow field of view");
        sim.camera().zoom(1.0);
        sim.update(0.0);
        sim.compute_view(1280, 720, 1.0, output);
        check(glm::length(output.camera.position - p) < 1e-8 && output.camera.focal_y > focal,
              "observer zoom changes optics without leaving the umbra");
    }
    sim.set_focus(sim.scene().find("Earth"));
    sim.update(3.0);
    OutputView normal;
    sim.compute_view(1280, 720, 1.0, normal);
    check(std::abs(normal.camera.focal_y - 1.0 / std::tan(22.5 * kDegToRad)) < 1e-6,
          "ordinary orbit shots restore the default field of view");

    Scene hidden;
    hidden.bodies.resize(2);
    hidden.bodies[0].equatorial_radius_km = 1.0;
    hidden.bodies[1].equatorial_radius_km = 0.001;
    hidden.bodies[1].world_position = glm::dvec3(10.0, 0.0, 0.0);
    OrbitCamera observer_camera;
    CameraShot shot;
    shot.target = 0;
    shot.observer = 1;
    shot.fov_deg = 10.0;
    observer_camera.fly_to(shot, 0.0, hidden);
    observer_camera.update(0.0, hidden);
    hidden.bodies[1].visible = false;
    observer_camera.update(0.0, hidden);
    const CameraView fallback = observer_camera.view(1280, 720);
    check(std::isfinite(glm::length(fallback.position)) &&
          std::abs(fallback.focal_y - 1.0 / std::tan(22.5 * kDegToRad)) < 1e-6,
          "unavailable observer returns to a finite ordinary orbit view");

    Scene bad;
    const std::string base = "name='x'\n[[bodies]]\nname='Earth'\nradii_km=[1.0]\n"
                             "[[bodies]]\nname='Moon'\nradii_km=[0.1]\n";
    check(!load_scene_string(base + "[[events]]\nname='bad'\ntime='2024-04-08'\nfocus='Moon'\n"
                            "observer='Earth'\nobserver_offset_km=[0,0,0]\n", "bad_observer", bad, &error),
          "rejects an observer inside a planet");
    check(!load_scene_string(base + "[[events]]\nname='bad'\ntime='2024-04-08'\nfov_deg=nan\n",
                            "bad_fov", bad, &error), "rejects non-finite field of view");
    check(!load_scene_string(base + "[corona]\ndisplay_gain=-1\n", "bad_corona", bad, &error),
          "rejects negative corona gain");
}

void test_jovian_aurora()
{
    Scene scene = load_scene_or_die("jupiter.toml");
    const Body& jupiter = scene.bodies[static_cast<size_t>(scene.find("Jupiter"))];
    check(bool(jupiter.aurora.jovian) && jupiter.aurora.model == AuroraModel::Jupiter && !jupiter.atmosphere.enabled,
          "Jupiter UV oval loads independently of Earth Kp and scattering atmosphere");
    if (!jupiter.aurora.jovian) return;
    const JovianAurora& oval = *jupiter.aurora.jovian;
    // Original Zenodo 12067692 CSV rows 0, 2500, 5000, 7500 (zero based,
    // excluding the header). System III WEST longitude, planetocentric latitude.
    const glm::dvec2 original[2][4] = {
        {{143.4605214673, 70.0201788883}, {182.7042086665, 57.6216824431},
         {220.9960067981, 68.7569755327}, {231.2970437547, 88.2575270876}},
        {{70.7399992225, -69.7600020669}, {146.1569465194, -82.6114458176},
         {301.3481340070, -77.2112209631}, {12.3009156519, -70.2662567914}}
    };
    for (int h = 0; h < 2; ++h) {
        const glm::dmat3 frame = JovianAurora::frame(h != 0);
        for (const auto& point : original[h]) {
            const double lon = -point.x * kDegToRad;
            const double lat = point.y * kDegToRad;
            const glm::dvec3 p{std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat)};
            const glm::dvec3 local = glm::transpose(frame) * p;
            const double angle = std::atan2(local.y, local.x);
            const double radius = oval.radius(h != 0, angle);
            const double error_km = std::abs(radius - std::acos(std::clamp(local.z, -1.0, 1.0))) * 71492.0;
            check(error_km < 100.0, "UV oval matches independent published lon/lat samples (km)", error_km);
            check(std::abs(radius - oval.radius(h != 0, angle + kTwoPi)) < 1e-12,
                  "UV oval wraps across the angular seam");
        }
        check(std::abs(oval.radius(h != 0, -1e-9) - oval.radius(h != 0, 1e-9)) < 1e-7,
              "UV oval is continuous through zero azimuth");
    }
    Scene solar = load_scene_or_die("solar_system.toml");
    const auto& other = solar.bodies[static_cast<size_t>(solar.find("Jupiter"))].aurora.jovian;
    check(other && other->radii == oval.radii, "Jupiter reference ovals agree across scene templates");

    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error), "UV event simulation loads");
    sim.stop_tour();
    int uv_events = 0;
    for (size_t i = 0; i < sim.scene().events.size(); ++i) {
        const SceneEvent& event = sim.scene().events[i];
        if (!event.aurora_ultraviolet || event.focus != sim.scene().find("Jupiter")) continue;
        ++uv_events;
        sim.jump_to_event(i);
        sim.update(0.0);
        sim.camera().update(3.0, sim.scene());
        OutputView view;
        sim.compute_view(1280, 720, 1.0, view);
        const Body& body = sim.scene().bodies[static_cast<size_t>(event.focus)];
        const double latitude = std::asin(glm::dot(glm::normalize(view.camera.position - body.world_position),
                                                   glm::normalize(body.orientation[2]))) / kDegToRad;
        check(std::abs(latitude - event.elevation_deg) < 0.01, "UV event looks at the intended pole", latitude);
        check(sim.last_event() == static_cast<int>(i), "UV event activates its overlay");
    }
    check(uv_events == 4, "both hemispheres have UV close-up events, and the moons' lines have two more", static_cast<double>(uv_events));
    sim.set_frame(0);
    check(sim.scene().active_event() == -1, "leaving the UV event clears its automatic display selection");
    sim.clock().paused = false;
    sim.clock().warp = 1e6;
    const double time = sim.aurora_time();
    sim.update(0.1);
    check(std::abs(sim.aurora_time() - time - 3.0) < 1e-9, "Jupiter illustrative animation caps at 30x");
    sim.clock().paused = true;
    const double paused = sim.aurora_time();
    sim.update(0.1);
    check(sim.aurora_time() == paused, "paused UV animation stays still");

    Scene bad;
    check(!load_scene_string("name='x'\n[[bodies]]\nname='Jupiter'\nradii_km=[71492,66854]\n"
                            "[bodies.aurora]\nkp='geomag/kp.bin'\njovian_oval='geomag/jupiter_oval.txt'\n",
                            "mixed_aurora", bad, &error, ASTRAXIS_ASSET_DIR),
          "rejects mixed terrestrial and Jovian aurora models");
}

void test_saturn_aurora()
{
    // An independent DFT recovers the published coefficients, including
    // the phase multiplier for each harmonic (Carbary 2012 equation 2).
    for (int h = 0; h < 2; ++h) {
        const double amplitudes[2][4] = {{15.1, 1.4, 0.2, 0.1}, {15.9, 2.4, 0.8, 0.3}};
        const double phases[2][3] = {{4.0, 20.7, 1.9}, {4.3, 21.9, 2.8}};
        for (int n = 0; n <= 3; ++n) {
            double cosine = 0.0, sine = 0.0;
            for (int k = 0; k < 240; ++k) {
                const double hour = k / 10.0;
                const double theta = SaturnAurora::profile(hour)[h] / kDegToRad;
                cosine += theta * std::cos(n * hour * kPi / 12.0) / 240.0;
                sine += theta * std::sin(n * hour * kPi / 12.0) / 240.0;
            }
            const double amplitude = n == 0 ? cosine : 2.0 * std::hypot(cosine, sine);
            check(std::abs(amplitude - amplitudes[h][n]) < 1e-10, "Saturn harmonic amplitude matches table 1");
            if (n > 0) {
                const double phase = n * phases[h][n - 1] * kPi / 12.0;
                check(std::abs(2.0 * cosine - amplitude * std::cos(phase)) < 1e-10 &&
                      std::abs(2.0 * sine - amplitude * std::sin(phase)) < 1e-10,
                      "Saturn local-time harmonic phase matches table 1");
            }
        }
    }
    for (double hour : {-36.0, -0.01, 0.0, 5.0, 12.0, 23.99}) {
        check(glm::length(SaturnAurora::profile(hour) - SaturnAurora::profile(hour + 24.0)) < 1e-12,
              "Saturn UV profile is periodic across midnight");
    }
    check(std::abs(SaturnAurora::profile(6).z * 2.354820045 / kDegToRad - 4.0) < 1e-10 &&
          std::abs(SaturnAurora::profile(18).z * 2.354820045 / kDegToRad - 8.0) < 1e-10,
          "Saturn dawn/dusk FWHM interpolates the observed widths");
    check(SaturnAurora::profile(5).w == 1.0 && SaturnAurora::profile(18).w < 0.003,
          "Saturn column brightness reaches the dawn peak and dusk baseline");
    check(glm::length(SaturnAurora::noon({0, 2, 3}) - glm::dvec3(0, 1, 0)) < 1e-12,
          "Saturn noon projects the Sun onto the spin equator");
    check(glm::length(SaturnAurora::noon({0, 0, 1})) == 1.0, "polar Sun has a finite noon fallback");

    Scene scene = load_scene_or_die("saturn.toml");
    const Body& saturn = scene.bodies[static_cast<size_t>(scene.find("Saturn"))];
    check(saturn.aurora.model == AuroraModel::Saturn && !saturn.atmosphere.enabled,
          "Saturn UV loads independently of Kp, Jupiter and atmospheric scattering");
    Scene solar = load_scene_or_die("solar_system.toml");
    check(solar.bodies[static_cast<size_t>(solar.find("Saturn"))].aurora.model == AuroraModel::Saturn,
          "solar-system Saturn inherits the same UV model");
    // Local-time geometry stays fixed relative to the Sun as System III turns.
    glm::dvec3 previous{};
    const double epoch = (jd_from_calendar({2008, 8, 26, 12, 0, 0}) - kJ2000Jd) * kSecondsPerDay;
    for (int k = 0; k < 2; ++k) {
        scene.update(epoch + k * 10000.0);
        const glm::dvec3 sun = scene.sun_position() - saturn.world_position;
        const glm::dvec3 noon = saturn.orientation * SaturnAurora::noon(glm::transpose(saturn.orientation) * sun);
        if (k) check(glm::dot(previous, noon) > 0.99999, "Saturn UV noon does not corotate with System III");
        previous = noon;
    }
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/saturn.toml", &error), "Saturn UV simulation loads");
    sim.stop_tour();
    int count = 0;
    for (size_t i = 0; i < sim.scene().events.size(); ++i) {
        const SceneEvent& event = sim.scene().events[i];
        if (!event.aurora_ultraviolet) continue;
        ++count;
        sim.jump_to_event(i);
        sim.update(0.0);
        sim.camera().update(3.0, sim.scene());
        OutputView view;
        sim.compute_view(1280, 720, 1.0, view);
        const Body& body = sim.scene().bodies[static_cast<size_t>(event.focus)];
        const double latitude = std::asin(glm::dot(glm::normalize(view.camera.position - body.world_position),
                                                   glm::normalize(body.orientation[2]))) / kDegToRad;
        check(std::abs(latitude - event.elevation_deg) < 0.01, "Saturn UV event views its intended pole", latitude);
        check(sim.scene().active_event() == static_cast<int>(i), "Saturn event activates UV overlay");
    }
    check(count == 2, "Saturn has both polar UV tours");

    // Dyudina 2015 Fig. 8: the boxes' clear-filter brightness, edge-on, and the color of the two layers.
    const AuroraProfile shell = aurora_profile(AuroraModel::Saturn, false, 3.0);
    const auto chapman = [](double h, double peak, double scale) {
        const double z = (h - peak) / scale;
        return std::exp(-z - std::exp(-z)) / scale;
    };
    const auto edge_on = [&](double tangent, double peak, double scale) {
        constexpr double kRadius = 54900.0; // 1 bar near 70 deg latitude
        double sum = 0.0;
        for (double s = 2.5; s < 12000.0; s += 5.0) {
            sum += 2.0 * chapman(std::hypot(kRadius + tangent, s) - kRadius, peak, scale) * 5.0;
        }
        return sum;
    };
    const auto clear_at = [&](double tangent) {
        return edge_on(tangent, shell.layers.x, shell.layers.y) +
               SaturnAurora::kVisibleHighColumn * edge_on(tangent, shell.layers.z, shell.layers.w);
    };
    double low_peak = 0.0;
    for (double h = 200.0; h < 700.0; h += 10.0) {
        low_peak = std::max(low_peak, clear_at(h));
    }
    check(std::abs(clear_at(shell.layers.z) / low_peak - 11.0 / 150.0) < 0.004,
          "Saturn's high layer is 11 kR against the low layer's 150 kR", clear_at(shell.layers.z) / low_peak);
    const auto ratio = [](const glm::dvec3& a, const glm::dvec3& b) { return glm::length(a / b - glm::dvec3(a.r / b.r)); };
    check(ratio(shell.low_rgb, SaturnAurora::kVisibleLowRgb) < 1e-12 && ratio(shell.high_rgb, SaturnAurora::kVisibleHighRgb) < 1e-12,
          "the shell keeps the measured colors");
    check(shell.high_rgb.r / shell.low_rgb.r < 0.1, "the high layer's column is an order of magnitude less");
    const glm::dvec3 low = SaturnAurora::kVisibleLowRgb;
    const glm::dvec3 high = SaturnAurora::kVisibleHighRgb;
    check(low.r > low.b && low.b > low.g && low.g > 0.3, "low layer: a pink, not a saturated magenta");
    check(high.b > high.r && high.r > high.g && high.g > 0.2 * high.b, "high layer: a purple with green in it");
    const auto visible = std::find_if(sim.scene().events.begin(), sim.scene().events.end(),
        [](const SceneEvent& e) { return e.name == "Saturn's visible aurora"; });
    check(visible != sim.scene().events.end(), "Saturn has a visible aurora tour");
    if (visible != sim.scene().events.end()) {
        check(visible->aurora_visible && !visible->aurora_ultraviolet, "visible aurora tour selects visible emission");
        sim.jump_to_event(static_cast<size_t>(visible - sim.scene().events.begin()));
        sim.update(0.0);
        sim.camera().update(3.0, sim.scene());
        OutputView view;
        sim.compute_view(1280, 720, 1.0, view);
        const Body& body = sim.scene().bodies[static_cast<size_t>(visible->focus)];
        const auto sun = sim.scene().sun_position();
        const auto direction = glm::normalize(view.camera.position - body.world_position);
        // Steep over the south pole so the oval clears the caption: mostly night, a sunlit crescent below.
        check(glm::dot(direction, glm::normalize(sun - body.world_position)) < -0.3,
              "visible aurora tour views the nightside", glm::dot(direction, glm::normalize(sun - body.world_position)));
        check(glm::dot(direction, glm::normalize(body.orientation[2])) < -0.3,
              "visible aurora tour views the southern curtain");
    }
    sim.set_frame(0);
    check(sim.scene().active_event() == -1, "leaving Saturn UV event clears automatic overlay");
    Scene bad;
    const std::string base = "name='x'\n[[bodies]]\nname='Saturn'\nradii_km=[60268,54364]\n[bodies.aurora]\n";
    for (const char* conflict : {"kp='geomag/kp.bin'", "jovian_oval='geomag/jupiter_oval.txt'",
                                 "footprints='geomag/jupiter_footprints.txt'"}) {
        check(!load_scene_string(base + "saturn_oval=true\n" + conflict, "mixed_saturn", bad, &error, ASTRAXIS_ASSET_DIR),
              "Saturn model rejects incompatible aurora fields");
    }
    check(!load_scene_string(base + "saturn_oval='true'", "bad_saturn", bad, &error, ASTRAXIS_ASSET_DIR),
          "Saturn aurora rejects wrong flag type");
    check(!load_scene_string("name='x'\n[[bodies]]\nname='Earth'\nradii_km=[6378]\n"
                            "[bodies.aurora]\nsaturn_oval=true", "bad_radii", bad, &error, ASTRAXIS_ASSET_DIR),
          "Saturn UV cannot silently use another planet's radii");
    check(load_scene_string(base + "saturn_oval=false", "disabled_saturn", bad, &error, ASTRAXIS_ASSET_DIR) &&
          bad.bodies[0].aurora.model == AuroraModel::None, "Saturn UV can be disabled explicitly");
}

void test_mars_aurora()
{
    const auto scene = load_scene_or_die("mars.toml");
    const auto solar = load_scene_or_die("solar_system.toml");
    check(scene.bodies[scene.find("Mars")].aurora.model == AuroraModel::Mars && solar.bodies[solar.find("Mars")].aurora.model == AuroraModel::Mars,
          "Mars aurora template loads in both scenes");
    int events = 0;
    for (const auto& e : scene.events) {
        if (e.name.find("Mars visible aurora") == std::string::npos) continue;
        ++events;
        check(e.aurora_visible && !e.aurora_ultraviolet && e.focus == scene.find("Mars"), "Mars visible tour settings");
    }
    check(events == 2, "Mars has the 2017 and 2024 visible aurora events");
    // Independent dense tangent-ray quadrature checks the orbit-limb scaling.
    double maximum = 0.0;
    for (int h = 20; h <= 140; ++h) {
        const double b = 3396.19 + h;
        const double end = std::sqrt(std::pow(3396.19 + 140.0, 2) - b*b);
        double integral = 0.0;
        for (int j = 0; j < 8000; ++j) {
            const double x = end * (j + 0.5) / 8000.0;
            const double altitude = std::sqrt(b*b + x*x) - 3396.19;
            const double sigma = altitude < 63.0 ? 12.0 : 17.0;
            integral += std::exp(-0.5 * std::pow((altitude - 63.0) / sigma, 2)) * 2.0 * end / 8000.0;
        }
        maximum = std::max(maximum, integral);
    }
    check(std::abs(maximum / MarsAurora::limb_integral() - 1.0) < 0.001, "Mars reference limb normalization", maximum);

    // Seen edge-on the layer has the shape of Schneider 2018 Fig. 2's CO2+ UVD (20 kR at its peak).
    const auto edge_on = [](double tangent) {
        constexpr double kRadius = 3396.19;
        const double b = kRadius + tangent;
        const double end = std::sqrt(std::pow(kRadius + MarsAurora::top_km, 2) - b * b);
        double sum = 0.0;
        for (int j = 0; j < 8000; ++j) {
            sum += MarsAurora::profile(std::hypot(b, end * (j + 0.5) / 8000.0) - kRadius) * 2.0 * end / 8000.0;
        }
        return sum;
    };
    const struct { double km, kr; } uvd[] = {{23, 11.2}, {35, 14.6}, {45, 17.3}, {57.5, 20.0}, {65, 16.2},
                                            {75, 10.3}, {85, 4.7}, {95, 1.8}};
    double worst = 0.0;
    for (const auto& point : uvd) {
        worst = std::max(worst, std::abs(20.0 * edge_on(point.km) / MarsAurora::limb_integral() - point.kr));
    }
    check(worst < 1.0, "Mars's layer edge-on follows the UVD limb profile (kR)", worst);
    // Lines of sight below the layer still cross it: edge-on, the peak is under the emission's.
    double limb_peak_km = 30.0;
    for (double h = 30.0; h < 90.0; h += 0.5) {
        if (edge_on(h) > edge_on(limb_peak_km)) limb_peak_km = h;
    }
    check(limb_peak_km > 50.0 && limb_peak_km < MarsAurora::peak_km, "the limb's peak is under the emission's", limb_peak_km);

    // Brightness by date. Days of September 2017.
    const auto september = [](double day) {
        return (jd_from_calendar({2017, 9, 1, 0, 0, 0}) - 1.0 + day - kJ2000Jd) * kSecondsPerDay;
    };
    const double peak = MarsAurora::limb_rayleigh_at(september(13.12));
    check(std::abs(peak - 21000.0) < 1.0, "the 2017 peak is ten times Knutsen's 2.1 kR", peak);
    double first_episode = 0.0;
    double later_episodes = 0.0;
    for (double day = 11.0; day < 22.0; day += 0.01) {
        const double r = MarsAurora::limb_rayleigh_at(september(day));
        if (day < 14.0) first_episode = std::max(first_episode, r);
        if (day >= 17.0) later_episodes = std::max(later_episodes, r);
    }
    check(first_episode > 10.0 * later_episodes, "the first episode is more than 10 times the others (Fig. 3)", first_episode / later_episodes);
    check(MarsAurora::limb_rayleigh_at(september(21.5)) > 3.0 * MarsAurora::limb_rayleigh_at(september(19.5)),
          "the third episode stands out of the lull");
    check(MarsAurora::limb_rayleigh_at(september(12.3)) < 0.05 * peak, "barely there the day before the peak");
    check(MarsAurora::limb_rayleigh_at(september(0.0)) == MarsAurora::limb_rayleigh &&
              MarsAurora::limb_rayleigh_at(september(40.0)) == MarsAurora::limb_rayleigh,
          "other dates use the reference observation");
    double largest_step = 1.0;
    for (double day = 5.0; day < 28.0; day += 0.001) {
        const double a = MarsAurora::limb_rayleigh_at(september(day));
        const double b = MarsAurora::limb_rayleigh_at(september(day + 0.001));
        largest_step = std::max(largest_step, std::max(a / b, b / a));
    }
    check(largest_step < 1.05, "no jumps, at either end of the record either", largest_step);

    // The simulation hands it to the shell.
    Simulation sim;
    std::string load_error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/mars.toml", &load_error), "Mars simulation loads");
    sim.stop_tour();
    const int mars = sim.scene().find("Mars");
    sim.clock().t_tdb = september(13.12);
    check(std::abs(sim.moon_aurora(mars).mars_limb_rayleigh - 21000.0) < 1.0, "the simulation reads the date");
    const AuroraProfile bright = aurora_profile(AuroraModel::Mars, false, 3.0, sim.moon_aurora(mars));
    const AuroraProfile reference = aurora_profile(AuroraModel::Mars, false, 3.0);
    check(std::abs(bright.low_rgb.g / reference.low_rgb.g - 10.0) < 1e-6, "the shell brightens with it");
    Scene bad;
    std::string error;
    const std::string base = "name='x'\n[[bodies]]\nname='Mars'\nradii_km=[3396.19,3376.20]\n[bodies.aurora]\n";
    for (const char* conflict : {"kp='x'", "europa_visible=true", "io_visible=true", "callisto_visible=true", "saturn_oval=true", "ganymede_map='x'", "jovian_oval='x'", "footprints='x'"})
        check(!load_scene_string(base + "mars_visible=true\n" + conflict, "bad", bad, &error, ASTRAXIS_ASSET_DIR), "Mars rejects mixed aurora models");
    check(!load_scene_string(base + "mars_visible='true'", "bad", bad, &error, ASTRAXIS_ASSET_DIR), "Mars flag requires boolean");
    check(load_scene_string(base + "mars_visible=false", "off", bad, &error, ASTRAXIS_ASSET_DIR) && bad.bodies[0].aurora.model == AuroraModel::None, "Mars can disable aurora");
    check(!load_scene_string("name='x'\n[[bodies]]\nname='x'\nradii_km=[6378]\n[bodies.aurora]\nmars_visible=true", "bad", bad, &error, ASTRAXIS_ASSET_DIR), "Mars validates radii");
}

void test_io_callisto_aurora()
{
    auto scene = load_scene_or_die("jupiter.toml");
    const auto solar = load_scene_or_die("solar_system.toml");
    for (const char* name : {"Io", "Callisto"}) {
        const int index = scene.find(name);
        const auto& moon = scene.bodies[index];
        const bool io = std::string(name) == "Io";
        const AuroraModel expected = io ? AuroraModel::Io : AuroraModel::Callisto;
        check(moon.aurora.model == expected, "visible moon aurora loads");
        const auto& other = solar.bodies[solar.find(name)];
        check(other.aurora.model == expected, "visible moon aurora template shared across scenes");
        const auto& model = io ? io_visible_aurora() : callisto_visible_aurora();
        int count = 0;
        for (const auto& event : scene.events) {
            if (event.focus != index || event.name.find("visible aurora") == std::string::npos) continue;
            ++count;
            scene.update(event.t_tdb);
            const auto& parent = scene.bodies[moon.parent];
            const auto sun = scene.sun_position() - moon.world_position;
            const auto direction = glm::normalize(sun);
            const auto oc = parent.world_position - moon.world_position;
            const double along = glm::dot(oc, direction);
            const double miss = glm::length(glm::cross(oc, direction));
            check(along > 0.0 && miss + moon.equatorial_radius_km + model.top_km +
                  along * 695700.0 / glm::length(sun) < parent.polar_radius_km,
                  "Io/Callisto event shell is within Jupiter umbra", miss);
            check(event.aurora_visible && !event.aurora_ultraviolet, "moon tour selects visible aurora");
        }
        check(count == 1, "moon has one visible aurora tour");
        // Independent view-ray integral checks both assumed layer heights.
        for (double H : {model.scale_km, model.sodium_scale_km}) {
            const double R = moon.equatorial_radius_km, top = R + model.top_km;
            double flux = 0.0;
            constexpr int radial = 4000, depth = 400;
            for (int i = 0; i < radial; ++i) {
                const double b2 = top * top * (i + 0.5) / radial;
                const double z1 = std::sqrt(top * top - b2);
                const double z0 = b2 < R * R ? std::sqrt(R * R - b2) : -z1;
                for (int j = 0; j < depth; ++j) {
                    const double z = z0 + (z1 - z0) * (j + 0.5) / depth;
                    const double h = std::sqrt(b2 + z * z) - R;
                    flux += std::exp(-h / H) / (H * (1.0 - std::exp(-model.top_km / H))) *
                            (z1 - z0) / depth * top * top / (R * R * radial);
                }
            }
            check(std::abs(flux - aurora_disk_factor(R, H, model.top_km)) < 0.004,
                  "moon column conversion matches independent aperture integral", flux);
        }
        Scene bad;
        std::string error;
        const std::string key = io ? "io_visible" : "callisto_visible";
        const std::string base = "name='x'\n[[bodies]]\nname='moon'\nradii_km=[" +
                                std::to_string(moon.equatorial_radius_km) + "]\n[bodies.aurora]\n";
        for (const char* conflict : {"europa_visible=true", "ganymede_map='x'", "kp='x'", "saturn_oval=true", "jovian_oval='x'", "footprints='x'"})
            check(!load_scene_string(base + key + "=true\n" + conflict, "mixed_moon", bad, &error, ASTRAXIS_ASSET_DIR),
                  "moon rejects incompatible aurora fields");
        check(!load_scene_string(base + key + "='true'", "bad_moon", bad, &error, ASTRAXIS_ASSET_DIR),
              "moon aurora requires boolean");
        check(!load_scene_string("name='x'\n[[bodies]]\nname='x'\nradii_km=[6378]\n[bodies.aurora]\n" +
                                key + "=true", "wrong_moon", bad, &error, ASTRAXIS_ASSET_DIR), "moon radius is validated");
    }
    check(kCallistoOxygenLines.size() == 2 && kCallistoOxygenLines[0].rayleigh == 19.0 &&
              std::abs(kCallistoOxygenLines[1].rayleigh / kCallistoOxygenLines[0].rayleigh - 1.0 / 3.09) < 1e-12,
          "Callisto uses detected red lines, not upper limits");
    check(kIoSodiumLines[0].rayleigh + kIoSodiumLines[1].rayleigh == 17170.0,
          "Io sodium brightness retains the measured eclipse average");
    check(io_visible_aurora().sodium_column_rgb.g > 0.0 && callisto_visible_aurora().column_rgb.g == 0.0,
          "Io sodium and Callisto oxygen have distinct spectra");
}

// The plasma sheet: Milby 2026's tables, by the longitudes he lists and in the scene.
void test_plasma_sheet()
{
    // Table 1 of Milby et al. 2025 and table 3.1 of the thesis: distance from Jupiter (RJ),
    // System III east longitude, height above the centrifugal equator.
    const struct { double r, east, height; } rows[] = {
        {5.913, 88.9, 0.232}, {5.913, 69.2, 0.003}, {5.877, 359.6, -0.633}, {5.917, 168.6, 0.671}, // Io
        {9.403, 69.1, 0.003}, {9.403, 33.3, -0.683}, {9.439, 312.4, -1.048}, {9.328, 165.4, 1.146}, // Europa
    };
    for (const auto& row : rows) {
        const double east = row.east * kDegToRad;
        const double height = centrifugal_equator_height_rj(row.r * 71492.0 * glm::dvec3(std::cos(east), std::sin(east), 0.0));
        check(std::abs(height - row.height) < 0.025, "the centrifugal equator's height matches Milby's tables", height - row.height);
    }

    Scene jupiter = load_scene_or_die("jupiter.toml");
    Scene solar = load_scene_or_die("solar_system.toml");
    // Times here are seen from Earth, as the observations were: light left the moon earlier.
    const auto at_moon = [&](const char* name, const char* utc) {
        double t = 0.0;
        check(parse_utc(utc, &t), "the observation's date parses");
        const int planet = solar.find("Jupiter"), earth = solar.find("Earth");
        const double light = glm::length(solar.icrf_state_at(planet, t).position - solar.icrf_state_at(earth, t).position) / kSpeedOfLightKmS;
        jupiter.update(t - light);
        return centrifugal_equator_height_rj(jupiter.bodies[static_cast<size_t>(jupiter.find(name))].in_parent_axes);
    };
    // Thesis tables 3.1 and 5.2. The tables give the first spectrum's start; the second
    // Europa row is the last one's, 5 min of its time are within the tolerance.
    const double europa_start = at_moon("Europa", "2021-05-20 13:40:00");
    std::printf("  info: Europa on 2021-05-20 at 13:40 UT: %.3f RJ from the sheet (Milby: 0.003 at his first spectrum)\n", europa_start);
    check(std::abs(europa_start) < 0.15, "Europa is at the sheet when the 2021-05-20 eclipse starts", europa_start);
    const struct { const char* utc; double height; } callisto[] = {
        {"2021-07-04 14:36:22", 2.52}, {"2021-07-04 15:07:33", 1.25}, {"2025-11-29 09:40:14", 4.23}, {"2025-11-29 11:23:55", 0.62}};
    for (const auto& row : callisto) {
        const double height = at_moon("Callisto", row.utc);
        std::printf("  info: Callisto at %s: %.2f RJ from the sheet (Milby: %.2f)\n", row.utc, height, row.height);
        // 2021 agrees to 0.03 RJ. 2025 is off by 0.3 RJ, eight minutes of Jupiter's turning;
        // the reason is not known (a slip of the clock, or of the scene's Callisto).
        check(std::abs(std::abs(height) - row.height) < 0.4, "Callisto's distance from the sheet matches Milby's table 5.2", height);
    }

    // Callisto's fit against his exposures (table 5.2, 2025-11-29): the mean of the fit over
    // the same times is close to the mean of the ones with a detection and non-detections.
    constexpr std::pair<double, double> kCallistoExposures[] = {
        {4.23, 9.0}, {4.07, 9.0}, {3.89, -3.0}, {3.66, -3.0}, {3.43, 9.0}, {3.18, -6.0}, {2.91, 7.0},
        {2.63, 12.0}, {2.31, 13.0}, {2.00, -6.0}, {1.64, 11.0}, {1.31, 25.0}, {0.97, 19.0}, {0.62, 26.0}};
    double observed = 0.0, model = 0.0;
    for (const auto& [height, brightness] : kCallistoExposures) {
        observed += brightness / std::size(kCallistoExposures);
        model += kCallistoOxygenLines[0].rayleigh * plasma_sheet_density(height, kCallistoPlasmaScaleRj) / std::size(kCallistoExposures);
    }
    std::printf("  info: Callisto's 630.0 nm over 2025-11-29: observed %.1f R (Milby's co-added 8.4), the fit %.1f R\n", observed, model);
    check(std::abs(observed - 8.4) < 1.0, "the exposures average to his co-added brightness", observed);
    check(model > observed && model < 2.0 * observed, "the fit is near them: it runs through the detections", model);

    // Europa: the reference is that of the eclipse; away from the sheet it is about half.
    const double reference = europa_reference_density();
    check(reference > 0.93 && reference < 0.97, "Europa's reference eclipse was near the sheet", reference);
    jupiter.update(0.0);
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error), "the sheet's simulation loads");
    sim.stop_tour();
    const int europa_index = sim.scene().find("Europa"), callisto_index = sim.scene().find("Callisto");
    sim.clock().warp = 60.0;
    double t0 = 0.0;
    check(parse_utc("2021-05-20 13:00:00", &t0), "the date parses");
    double brightest = 0.0, faintest = 1e9, callisto_brightest = 0.0, callisto_faintest = 1e9, mean = 0.0;
    constexpr int kSamples = 200;
    for (int k = 0; k < kSamples; ++k) {
        sim.clock().t_tdb = t0 + k * 600.0;
        sim.update(0.0);
        const double e = sim.moon_aurora(europa_index).plasma, c = sim.moon_aurora(callisto_index).plasma;
        brightest = std::max(brightest, e);
        faintest = std::min(faintest, e);
        callisto_brightest = std::max(callisto_brightest, c);
        callisto_faintest = std::min(callisto_faintest, c);
        mean += e / kSamples;
    }
    std::printf("  info: Europa's plasma factor %.2f .. %.2f, mean %.2f; Callisto's %.2f .. %.2f\n", faintest, brightest, mean, callisto_faintest, callisto_brightest);
    // At his 1.7 RJ the thesis' Figure 3.5 goes down about to 0.6 away from the equator.
    check(brightest > 1.0 && faintest / brightest > 0.5 && faintest / brightest < 0.7, "Europa dims by about 40% away from the sheet", faintest / brightest);
    check(callisto_faintest / callisto_brightest < 0.35, "Callisto, much farther out, swings more", callisto_faintest / callisto_brightest);
    sim.clock().warp = 86400.0;
    sim.update(0.0);
    check(std::abs(sim.moon_aurora(europa_index).plasma - mean) < 0.02, "at a day per second Europa takes its long-term mean", sim.moon_aurora(europa_index).plasma);
    MoonAuroraState half_plasma;
    half_plasma.plasma = 0.5;
    check(aurora_profile(AuroraModel::Europa, false, 3.0, half_plasma).low_rgb == aurora_profile(AuroraModel::Europa, false, 3.0).low_rgb * 0.5,
          "the plasma scales Europa's column");
}

// Europa's pattern: Roth et al. 2016 and 2021.
void test_europa_pattern()
{
    // At System III west longitude: the field line, the sign of the sheet, as the paper's text.
    const auto at_west = [](double west, double radius_rj = 9.4) {
        const double east = -west * kDegToRad;
        return radius_rj * 71492.0 * glm::dvec3(std::cos(east), std::sin(east), 0.0);
    };
    const auto sheet_at = [&](double west) { return -jovian_magnetic_latitude_deg(at_west(west)) / kJovianDipoleTiltDeg; };
    check(sheet_at(16.0) > 0.95, "farthest below the plasma sheet at 16 W (Roth section 3.2.1)", sheet_at(16.0));
    check(sheet_at(196.0) < -0.95, "farthest above it at 196 W", sheet_at(196.0));
    // Roth's Figure 10: above (181 W), the line's north end leans to the sub-Jovian side;
    // the projected tilt reaches about 20 degrees in Figure 7a.
    const glm::dvec3 above = jovian_field_line_at_moon(at_west(196.0));
    check(above.x > 0.0, "above the equator the line's north end leans toward Jupiter", above.x);
    const double tilt = std::atan2(above.x, above.z) / kDegToRad;
    check(std::abs(tilt - 19.0) < 4.0, "by about twice the magnetic latitude", tilt);
    const glm::dvec3 below = jovian_field_line_at_moon(at_west(16.0));
    check(below.x < 0.0 && std::abs(below.y) < 0.05 && std::abs(above.y) < 0.05, "and away from it below", below.x);

    // Means over the sphere by regions, in the body axes.
    constexpr int kSamples = 20000;
    const auto mean_over = [&](const glm::dvec3& field, double sheet, const glm::dvec3& sun, auto&& keep) {
        double sum = 0.0;
        int count = 0;
        for (int k = 0; k < kSamples; ++k) {
            const double z = 1.0 - 2.0 * (k + 0.5) / kSamples;
            const double rho = std::sqrt(1.0 - z * z);
            const double angle = k * 2.39996322972865332;
            const glm::dvec3 n(rho * std::cos(angle), rho * std::sin(angle), z);
            if (!keep(n)) continue;
            sum += europa_pattern(n, field, sheet, sun);
            ++count;
        }
        return sum / count;
    };
    const glm::dvec3 up(0.0, 0.0, 1.0), trailing(0.0, 1.0, 0.0), leading(0.0, -1.0, 0.0), to_jupiter(1.0, 0.0, 0.0);
    const auto all = [](const glm::dvec3&) { return true; };
    for (const auto& [field, sheet] : {std::pair{up, 0.0}, std::pair{above, -1.0}, std::pair{below, 1.0}}) {
        // Without the afternoon term's small correlation with the rest, a few percent.
        check(std::abs(mean_over(field, sheet, to_jupiter, all) / europa_pattern_mean(field, sheet) - 1.0) < 0.05,
              "the pattern's mean is what the shell divides by", mean_over(field, sheet, to_jupiter, all) / europa_pattern_mean(field, sheet));
    }
    // The afternoon side, seen as HST sees it (the sun behind the observer).
    for (const glm::dvec3& sun : {trailing, leading, to_jupiter}) {
        const glm::dvec3 afternoon = glm::cross(up, sun);
        const auto side = [&](double sign) {
            return mean_over(up, 0.0, sun, [&](const glm::dvec3& n) { return glm::dot(n, sun) > 0.0 && sign * glm::dot(n, afternoon) > 0.0; });
        };
        check(std::abs(side(1.0) / side(-1.0) - 1.59) < 0.1, "1.59 between the afternoon and morning halves (Figure 8a)", side(1.0) / side(-1.0));
    }
    // In eclipse the sun is behind Jupiter, which makes the afternoon side the trailing side,
    // as in de Kleer & Brown 2018.
    const auto hemisphere = [&](const glm::dvec3& field, double sheet, const glm::dvec3& view, double north) {
        return mean_over(field, sheet, view, [&](const glm::dvec3& n) { return glm::dot(n, view) > 0.0 && north * n.z > 0.5; });
    };
    const auto north_over_south = [&](const glm::dvec3& field, double sheet, const glm::dvec3& view) {
        return hemisphere(field, sheet, view, 1.0) / hemisphere(field, sheet, view, -1.0);
    };
    const double upstream_below = north_over_south(below, 1.0, trailing), upstream_above = north_over_south(above, -1.0, trailing);
    const double downstream_below = north_over_south(below, 1.0, leading), downstream_above = north_over_south(above, -1.0, leading);
    std::printf("  info: Europa's north over south: trailing %.2f, %.2f (Figure 7c: about 3.6, 0.26); leading %.2f, %.2f (Figure 7d: about 1.8, 0.7)\n",
                upstream_below, upstream_above, downstream_below, downstream_above);
    check(upstream_below > 2.0 && upstream_above < 0.5, "the side facing the plasma sheet is brighter, from upstream", upstream_below);
    check(downstream_below > 1.2 && downstream_above < 0.85, "and less so from downstream", downstream_below);
    check(upstream_below > downstream_below && upstream_above < downstream_above, "more on the trailing side");
    // Figure 7a: the quadrants with the magnetic poles are brighter than the others.
    const auto quadrants = [&](const glm::dvec3& field, double sheet, const glm::dvec3& view, double sign) {
        const glm::dvec3 right = glm::cross(up, view);
        const auto in_quadrant = [&](const glm::dvec3& n, double s) { return glm::dot(n, view) > 0.0 && s * glm::dot(n, right) * n.z > 0.0; };
        return mean_over(field, sheet, view, [&](const glm::dvec3& n) { return in_quadrant(n, sign); }) /
               mean_over(field, sheet, view, [&](const glm::dvec3& n) { return in_quadrant(n, -sign); });
    };
    // Trailing view: right is anti-Jovian, so above the equator the line runs from the upper left.
    check(quadrants(above, -1.0, trailing, -1.0) > 1.2, "the light gathers along the tilted field line", quadrants(above, -1.0, trailing, -1.0));
    check(quadrants(below, 1.0, leading, -1.0) > 1.2, "on the other side too", quadrants(below, 1.0, leading, -1.0));

    // Roth 2021 section 4.1, visit 17: out of eclipse at 146-161 W, into sunlight at 194-209 W,
    // 135.6 nm falls from 40.2 +- 3.4 to 36.3 +- 3.5 R, 10 +- 11%, moving away from the sheet.
    const double in_eclipse = plasma_sheet_density(centrifugal_equator_height_rj(at_west(153.5)), kEuropaPlasmaScaleRj);
    const double after_egress = plasma_sheet_density(centrifugal_equator_height_rj(at_west(201.5)), kEuropaPlasmaScaleRj);
    std::printf("  info: Europa's plasma from visit 17's eclipse to its egress: %.2f (Roth 2021: 0.90 +- 0.11)\n", after_egress / in_eclipse);
    check(std::abs(after_egress / in_eclipse - 0.90) < 1.5 * 0.11, "the plasma sheet's fall is within 1.5 sigma of the brightness across visit 17", after_egress / in_eclipse);
    // Roth 2016 Eq. 2: a fitted scale of 1.8 +- 0.1 RJ at 135.6 nm, 1.6 +- 0.1 at 130.4 nm.
    check(kEuropaPlasmaScaleRj >= 1.6 && kEuropaPlasmaScaleRj <= 1.8, "1.7 RJ lies between Roth's two lines");

    // The simulation: the state follows the moon, the pattern's total does not change.
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error), "Europa's simulation loads");
    sim.stop_tour();
    sim.clock().warp = 60.0;
    const int index = sim.scene().find("Europa");
    double t0 = 0.0;
    check(parse_utc("2021-05-20 13:00:00", &t0), "the date parses");
    double most_north = -2.0, most_south = 2.0;
    for (int k = 0; k < 40; ++k) {
        sim.clock().t_tdb = t0 + k * 900.0;
        sim.update(0.0);
        const MoonAuroraState state = sim.moon_aurora(index);
        most_north = std::max(most_north, state.sheet);
        most_south = std::min(most_south, state.sheet);
        const Body& moon = sim.scene().bodies[static_cast<size_t>(index)];
        check(std::abs(glm::dot(state.field, jovian_field_line_at_moon(moon.in_parent_axes)) - 1.0) < 1e-9, "the state's field line is the moon's");
    }
    check(most_north > 0.9 && most_south < -0.9, "Europa goes through the sheet", most_north);
    // Its body axes are the pattern's: x toward Jupiter, y trailing.
    const Body& moon = sim.scene().bodies[static_cast<size_t>(index)];
    const Body& parent = sim.scene().bodies[static_cast<size_t>(moon.parent)];
    const glm::dvec3 toward = glm::normalize(glm::transpose(moon.orientation) * (parent.world_position - moon.world_position));
    check(toward.x > 0.98, "Europa's x axis points at Jupiter", toward.x);
    sim.clock().warp = 86400.0;
    sim.update(0.0);
    check(sim.moon_aurora(index).sheet == 0.0 && sim.moon_aurora(index).field == glm::dvec3(0.0, 0.0, 1.0), "at a day per second it evens out");
    MoonAuroraState tilted;
    tilted.field = above;
    tilted.sheet = -1.0;
    const AuroraProfile shell = aurora_profile(AuroraModel::Europa, false, 3.0, tilted);
    check(glm::length(shell.low_rgb * europa_pattern_mean(above, -1.0) - aurora_profile(AuroraModel::Europa, false, 3.0).low_rgb) <
              1e-6 * glm::length(shell.low_rgb),
          "the pattern moves the light, the plasma sets how much");
    check(shell.map_sun.x == -kEuropaSheet && shell.map_sun.y == -kEuropaUpstream && shell.map_sun.z == kEuropaFieldTilt && shell.map_sun.w == kEuropaAfternoon,
          "the shader gets the pattern's numbers");
}

// Io through an eclipse: Schmidt et al. 2023, de Pater et al. 2020, Milby et al. 2025.
void test_io_eclipse_levels()
{
    const MoonEclipse eclipse{0.0, 8000.0};
    check(std::abs(io_sodium_level(eclipse, 0.0) - 1.0) < 1e-12, "sodium starts from its level before the eclipse");
    double first = 0.0;
    for (int k = 0; k < 600; ++k) first += io_sodium_level(eclipse, k + 0.5) / 600.0;
    check(std::abs(io_sodium_level(eclipse, 7000.0) / first - 0.27) < 0.005, "sodium's steady level is 27% of its first 10 minutes (Schmidt)", io_sodium_level(eclipse, 7000.0) / first);
    const double steady = io_sodium_level(eclipse, 7000.0);
    check(std::abs((io_sodium_level(eclipse, 600.0) - steady) / (1.0 - steady) - std::exp(-1.0)) < 1e-4, "it falls with an e-folding time of 10 minutes");
    check(io_sodium_level(eclipse, 1800.0) < steady * 1.25, "it is about steady within 30 minutes");
    check(std::abs(io_sodium_level(eclipse, 8000.0 + 7200.0) - (1.0 - (1.0 - steady) / std::exp(1.0))) < 1e-3, "it recovers with an e-folding time of 2 hours");
    check(std::abs(io_so2_level(eclipse, 7000.0) - 0.4) < 1e-3, "SO2 falls to 2.5 times less than in sunlight (de Pater)");
    check(std::abs(io_so2_level(eclipse, 8000.0 + 300.0) - 0.7) < 1e-3 && io_so2_level(eclipse, 8000.0 + 600.0) == 1.0, "and returns linearly within 10 minutes");

    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error), "Io's simulation loads");
    sim.stop_tour();
    sim.clock().warp = 60.0;
    const int io = sim.scene().find("Io");
    // The exposures of 2023-08-09 (Milby et al. 2025, Table 7: KOA file names give the start
    // in seconds of the day, five minutes each). Light left Io earlier.
    Scene solar = load_scene_or_die("solar_system.toml");
    double day = 0.0;
    check(parse_utc("2023-08-09 00:00:00", &day), "the date parses");
    const double light = glm::length(solar.icrf_state_at(solar.find("Jupiter"), day).position - solar.icrf_state_at(solar.find("Earth"), day).position) / kSpeedOfLightKmS;
    double mean_level = 0.0;
    MoonEclipse found;
    constexpr double kStarts[] = {49839.52, 50269.45, 50708.50, 51108.91, 51548.20, 51977.44, 52379.83};
    for (const double start : kStarts) {
        sim.clock().t_tdb = day + start + 150.0 - light;
        sim.update(0.0);
        check(sim.io_eclipse(found) && found.ingress <= sim.clock().t_tdb && sim.clock().t_tdb < found.egress,
              "each exposure of 2023-08-09 is in the eclipse");
        mean_level += io_sodium_level(found, sim.clock().t_tdb) / std::size(kStarts);
    }
    std::printf("  info: Io's ingress at %.0f s of the day, light time %.0f s, exposures %.1f .. %.1f min after it, mean sodium level %.4f\n",
                found.ingress - day, light, (day + kStarts[0] + 150.0 - light - found.ingress) / 60.0,
                (day + kStarts[6] + 150.0 - light - found.ingress) / 60.0, mean_level);
    check(std::abs(mean_level / kIoSodiumReferenceLevel - 1.0) < 0.01, "the sodium's reference level is the mean over Milby's exposures", mean_level);

    // A whole cycle, by the search and against the umbra as the scene says.
    const double eclipse_length = found.egress - found.ingress;
    check(eclipse_length > 1.9 * 3600.0 && eclipse_length < 2.4 * 3600.0, "Io's eclipse lasts about 2.2 hours", eclipse_length / 3600.0);
    sim.clock().t_tdb = found.ingress + 60.0;
    sim.update(0.0);
    const MoonAuroraState early = sim.moon_aurora(io);
    sim.clock().t_tdb = found.ingress + 3600.0;
    sim.update(0.0);
    const MoonAuroraState late = sim.moon_aurora(io);
    sim.clock().t_tdb = found.egress + 8.0 * 3600.0;
    sim.update(0.0);
    const MoonAuroraState after = sim.moon_aurora(io);
    check(early.sodium > 2.5 * late.sodium && after.sodium > early.sodium, "Io's sodium fades in the shadow and recovers after it");
    check(std::abs(late.sodium - 1.0) < 0.25, "the eclipse's reference is about its steady level", late.sodium);
    check(early.so2 > 1.5 * late.so2 && after.so2 > early.so2, "Io's SO2 too");
    check(late.plasma == 1.0, "Io's oxygen does not respond: Schmidt, Milby");
    // Backward in time, and at a warp too fast to follow: the level between eclipses.
    sim.clock().t_tdb = found.ingress - 3600.0;
    sim.update(0.0);
    check(std::abs(sim.moon_aurora(io).sodium * kIoSodiumReferenceLevel - 1.0) < 0.01, "before the ingress, the level between eclipses");
    sim.clock().t_tdb = found.ingress + 3600.0;
    sim.clock().warp = 86400.0;
    sim.update(0.0);
    check(std::abs(sim.moon_aurora(io).sodium * kIoSodiumReferenceLevel - 1.0) < 1e-9, "too fast to follow, also");
    MoonAuroraState io_state;
    io_state.sodium = 0.5;
    io_state.so2 = 0.25;
    const AuroraProfile io_shell = aurora_profile(AuroraModel::Io, false, 3.0, io_state);
    const auto& io_model = io_visible_aurora();
    check(io_shell.high_rgb == io_model.sodium_column_rgb * 0.5 * io_model.display_gain, "sodium scales the corona");
    check(glm::length(io_shell.low_rgb - (io_model.column_rgb + io_model.so2_column_rgb * 0.25) * io_model.display_gain) < 1e-9 * glm::length(io_shell.low_rgb),
          "SO2 scales its part of the low layer");
    check(io_model.so2_column_rgb.b > 10.0 * io_model.so2_column_rgb.g && io_model.so2_column_rgb.r == 0.0, "SO2's band is blue");

    // The glows are sub- and anti-Jovian: the body's x axis points at Jupiter.
    const Scene& scene = sim.scene();
    const Body& moon = scene.bodies[static_cast<size_t>(io)];
    const glm::dvec3 to_jupiter = glm::normalize(glm::transpose(moon.orientation) *
                                                 (scene.bodies[static_cast<size_t>(moon.parent)].world_position - moon.world_position));
    check(to_jupiter.x > 0.98, "Io's body x axis points at Jupiter", to_jupiter.x);
}

void test_aurora_profiles()
{
    using Model = AuroraModel;
    using Pattern = AuroraPattern;
    using Vertical = AuroraVertical;
    using Lighting = AuroraLighting;
    const struct {
        Model model;
        bool ultraviolet;
        Pattern pattern;
        Vertical vertical;
        Lighting lighting;
    } shells[] = {
        {Model::Earth, false, Pattern::EarthOval, Vertical::EarthChapman, Lighting::NightSide},
        {Model::Jupiter, true, Pattern::JovianOval, Vertical::ChapmanPair, Lighting::Always},
        {Model::Saturn, true, Pattern::SaturnOval, Vertical::ChapmanPair, Lighting::Always},
        {Model::Saturn, false, Pattern::SaturnOval, Vertical::ChapmanPair, Lighting::NightSide},
        {Model::Ganymede, true, Pattern::BodyMapShares, Vertical::UniformColumn, Lighting::Always},
        {Model::Uranus, false, Pattern::BodyMapSeasons, Vertical::Table, Lighting::Always},
        {Model::Neptune, false, Pattern::BodyMapDusk, Vertical::ChapmanPair, Lighting::Always},
        {Model::Europa, false, Pattern::Europa, Vertical::Exponential, Lighting::NightOrEclipsed},
        {Model::Io, false, Pattern::IoFlanks, Vertical::Exponential, Lighting::NightOrEclipsed},
        {Model::Callisto, false, Pattern::Uniform, Vertical::Exponential, Lighting::NightOrEclipsed},
        {Model::Mars, false, Pattern::Uniform, Vertical::Gaussian, Lighting::NightSide},
    };
    for (const auto& s : shells) {
        const AuroraProfile p = aurora_profile(s.model, s.ultraviolet, 3.0);
        check(p.pattern == s.pattern && p.vertical == s.vertical && p.lighting == s.lighting,
              "each aurora model has its shader choices");
        check(p.top_km > p.bottom_km && p.bottom_km >= 0.0, "aurora shell has thickness");
        check(glm::dot(p.low_rgb, glm::dvec3(1.0)) > 0.0, "aurora shell emits");
        check(p.layers.x > 0.0 && (s.vertical == Vertical::UniformColumn || (p.layers.y > 0.0 && p.layers.z >= 0.0)),
              "aurora layers have positive sizes");
        if (s.vertical == Vertical::Exponential) {
            check(p.layers.z > 0.0 && p.layers.y == p.top_km, "exponential layers have both scales and the shell's top");
        }
    }
    const AuroraProfile spots = jovian_footprint_profile();
    check(spots.pattern == Pattern::JovianFootprints && spots.sigma_km > 0.0 && spots.layers.x == 900.0 &&
              spots.top_km > aurora_profile(Model::Jupiter, true, 3.0).top_km,
          "Jupiter's footprints have a higher shell of their own");

    // Which models each view shows.
    for (const Model m : {Model::Earth, Model::Saturn, Model::Europa, Model::Io, Model::Callisto, Model::Mars}) {
        check(aurora_shown(m, true, false, false) && aurora_shown(m, true, true, false) && aurora_shown(m, true, false, true) && !aurora_shown(m, false, false, false),
              "visible aurora shows in both views");
    }
    for (const Model m : {Model::Jupiter, Model::Ganymede}) {
        check(!aurora_shown(m, true, false, true) && aurora_shown(m, true, true, false), "UV data shows only in the UV view");
    }
    for (const Model m : {Model::Uranus, Model::Neptune}) {
        check(!aurora_shown(m, true, true, false) && aurora_shown(m, true, false, true), "IR data shows only in the IR view");
    }
    check(!aurora_shown(Model::None, true, true, true), "no model, no aurora");
    check(aurora_surface_gain(Model::Uranus, true) < 0.1 && aurora_surface_gain(Model::Neptune, true) < 0.1 &&
              aurora_surface_gain(Model::Uranus, false) == 1.0 && aurora_surface_gain(Model::Earth, true) == 1.0,
          "only the infrared view darkens only the ice giants");

    // Saturn's two views are different shells; the oval is the same.
    const AuroraProfile saturn_uv = aurora_profile(Model::Saturn, true, 3.0);
    const AuroraProfile saturn_visible = aurora_profile(Model::Saturn, false, 3.0);
    check(saturn_uv.high_rgb == glm::dvec3(0.0) && saturn_visible.high_rgb.b > saturn_visible.high_rgb.g &&
              saturn_visible.low_rgb.r > saturn_visible.low_rgb.b,
          "Saturn is false color in the UV view, pink and purple in the visible");
    check(saturn_uv.oval[10] == saturn_visible.oval[10] && saturn_uv.oval[10].x > 0.0f, "Saturn's oval in both views");

    // Earth's activity reaches the profile.
    const AuroraProfile quiet = aurora_profile(Model::Earth, false, 1.0);
    const AuroraProfile storm = aurora_profile(Model::Earth, false, 9.0);
    check(storm.blue_kr > quiet.blue_kr && storm.oval[0].y > quiet.oval[0].y, "Kp sets Earth's column and oval");

    // The moons' spectra: the column colors, with their exposure.
    const glm::dvec3 europa = (kRedPerKr * (0.573 + 0.190 * 630.0 / 636.4) + kGreenPerKr * 0.037) * 3.0e7 /
                              aurora_disk_factor(1560.8, 150.0, 600.0);
    check(glm::length(aurora_profile(Model::Europa, false, 3.0).low_rgb - europa) < 1e-9 * glm::length(europa),
          "Europa's column is its disk-average oxygen lines");
    check(aurora_profile(Model::Europa, false, 3.0).high_rgb == glm::dvec3(0.0) &&
              aurora_profile(Model::Callisto, false, 3.0).high_rgb == glm::dvec3(0.0) &&
              aurora_profile(Model::Io, false, 3.0).high_rgb.g > 0.0,
          "only Io has a sodium layer");
    const glm::dvec3 mars = kGreenPerKr * (MarsAurora::limb_rayleigh / 1000.0 * MarsAurora::display_gain / MarsAurora::limb_integral());
    check(aurora_profile(Model::Mars, false, 3.0).low_rgb == mars && mars.g > mars.r && mars.b == 0.0,
          "Mars glows in the green line");
}

// The ice giants' H3+ in the infrared: Melin et al. 2025 (Uranus's map and its Table 1),
// Tiranti et al. 2026 (its vertical profile), Melin et al. 2025 (Neptune).
struct SubPoint {
    double east_deg, latitude_deg, light_s;
};

// Where a line of sight from the Earth at t_tdb meets the body, with the light time.
SubPoint sub_earth_point(Scene& scene, int index, double t_tdb)
{
    const int earth_index = scene.find("Earth");
    scene.update(t_tdb);
    const glm::dvec3 earth = scene.bodies[static_cast<size_t>(earth_index)].icrf_position;
    double light = 0.0;
    for (int i = 0; i < 3; ++i) {
        scene.update(t_tdb - light);
        light = glm::length(scene.bodies[static_cast<size_t>(index)].icrf_position - earth) / 299792.458;
    }
    scene.update(t_tdb - light);
    const Body& body = scene.bodies[static_cast<size_t>(index)];
    const glm::dvec3 n = glm::normalize(glm::transpose(scene.body_to_icrf(body, t_tdb - light)) * (earth - body.icrf_position));
    return {std::atan2(n.y, n.x) / kDegToRad, std::asin(n.z) / kDegToRad, light};
}

double angle_difference_deg(double a, double b)
{
    return std::remainder(a - b, 360.0);
}

// Neptune's dark spots against Wong et al. 2022 (Tables 2-5), and their events.
void test_dark_spots()
{
    Scene scene = load_scene_or_die("neptune.toml");
    const int index = scene.find("Neptune");
    Body& neptune = scene.bodies[static_cast<size_t>(index)];
    const auto find = [&](const char* name) -> const DarkSpot* {
        for (const DarkSpot& s : neptune.dark_spots) {
            if (s.name == name) return &s;
        }
        check(false, name);
        return nullptr;
    };
    const DarkSpot* gds = find("GDS");
    const DarkSpot* nds = find("NDS-2018");
    const DarkSpot* second = find("Secondary dark feature");
    if (!gds || !nds || !second) return;
    const double radius = neptune.equatorial_radius_km;
    const auto year = [](double y) { return tdb_from_julian_year(y); };
    // The measured contrast, without what is added for the eye.
    const auto state = [&](const DarkSpot& s, double y) {
        DarkSpotState out = dark_spot_state(s, radius, year(y));
        out.contrast /= s.display_gain;
        return out;
    };
    const auto planetographic = [&](double lat_rad) {
        const double ratio = neptune.equatorial_radius_km / neptune.polar_radius_km;
        return std::atan(ratio * ratio * std::tan(lat_rad)) / kDegToRad;
    };

    // Table 3, 4 and 5, at their epochs.
    const DarkSpotState peak = state(*nds, 2020.48);
    check(std::abs(peak.contrast.z + 0.187) < 1e-6, "NDS-2018 at its darkest, -18.7 percent in F467M", peak.contrast.z);
    check(std::abs(planetographic(peak.lat_rad) - 17.56) < 1e-6, "its latitude in 2020.48 (planetographic)", planetographic(peak.lat_rad));
    check(std::abs(peak.width_ew_rad * radius - 6656.0) < 1e-6 && std::abs(peak.width_ns_rad * radius - 4019.0) < 1e-6, "its size in 2020.48");
    check(std::abs(state(*nds, 2020.02).contrast.y + 0.126) < 1e-6, "F547M in 2020.02", state(*nds, 2020.02).contrast.y);
    const double southmost = planetographic(state(*nds, 2020.02).lat_rad);
    check(planetographic(state(*nds, 2020.48).lat_rad) > southmost, "it turned back in early 2020", southmost);
    const double migration = (planetographic(state(*nds, 2018.85).lat_rad) - southmost) / (2020.02 - 2018.85);
    check(std::abs(migration - 3.34) < 0.3, "3.34 deg/yr toward the equator until then", migration);
    check(std::abs(state(*second, 2020.02).width_ew_rad * radius - 7248.0) < 1e-6, "the secondary feature's size in 2020.02");

    // They come and go.
    check(state(*gds, 1987.0).contrast == glm::dvec3(0.0), "no Great Dark Spot before Voyager's time");
    check(std::abs(state(*gds, 1989.65).contrast.z + 0.12) < 1e-6, "12 percent in blue at the flyby", state(*gds, 1989.65).contrast.z);
    check(std::abs(planetographic(state(*gds, 1989.65).lat_rad) + 20.0) < 1e-6, "20 S at the flyby", planetographic(state(*gds, 1989.65).lat_rad));
    const double monthly = (planetographic(state(*gds, 1989.65).lat_rad) - planetographic(state(*gds, 1989.0).lat_rad)) / (0.65 * 12.0);
    check(std::abs(monthly - 1.3) < 0.05, "1.3 deg a month toward the equator", monthly);
    check(state(*gds, 1994.8).contrast == glm::dvec3(0.0), "gone when Hubble could see it in 1994");
    check(state(*nds, 2017.5).contrast == glm::dvec3(0.0) && state(*nds, 2023.0).contrast == glm::dvec3(0.0), "NDS-2018 only in its lifetime");
    check(state(*second, 2019.0).contrast == glm::dvec3(0.0) && state(*second, 2021.0).contrast == glm::dvec3(0.0), "the secondary feature in 2019-2020");
    check(std::abs(state(*nds, 2020.95).contrast.z / state(*nds, 2020.95).contrast.y - 0.166 / 0.119) < 1e-6, "bluer than green");

    // The shape: half depth half a width away, along both axes.
    const DarkSpotState s = state(*nds, 2020.02);
    const auto at = [](double lon, double lat) { return glm::dvec3(std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat)); };
    check(std::abs(dark_spot_weight(s, at(s.east_lon_rad, s.lat_rad)) - 1.0) < 1e-12, "the spot's center");
    check(std::abs(dark_spot_weight(s, at(s.east_lon_rad, s.lat_rad + 0.5 * s.width_ns_rad)) - 0.5) < 1e-6, "half depth north");
    check(std::abs(dark_spot_weight(s, at(s.east_lon_rad + 0.5 * s.width_ew_rad / std::cos(s.lat_rad), s.lat_rad)) - 0.5) < 1e-6, "half depth east");
    check(std::abs(dark_spot_weight(s, at(s.east_lon_rad + kTwoPi - 0.5 * s.width_ew_rad / std::cos(s.lat_rad), s.lat_rad)) - 0.5) < 1e-6, "across the seam");

    // It goes round westward at the measured speed.
    const double day = 86400.0;
    const DarkSpotState later = dark_spot_state(*nds, radius, nds->anchor_tdb + day);
    const double moved = std::remainder(later.east_lon_rad - nds->anchor_east_lon_rad, kTwoPi) * radius * std::cos(dark_spot_state(*nds, radius, nds->anchor_tdb).lat_rad);
    check(std::abs(moved * 1000.0 / day + 333.0) < 0.01, "333 m/s westward", moved * 1000.0 / day);

    // Each event opens with its spot on the sunlit side.
    const struct { const char* event; const DarkSpot* spot; double west_of_center_deg; } shots[] = {
        {"The Great Dark Spot", gds, 25.0},
        {"A dark spot that turned back", nds, 55.0},
    };
    for (const auto& shot : shots) {
        const auto e = std::find_if(scene.events.begin(), scene.events.end(), [&](const SceneEvent& v) { return v.name.find(shot.event) != std::string::npos; });
        check(e != scene.events.end(), shot.event);
        if (e == scene.events.end()) continue;
        scene.update(e->t_tdb);
        const glm::dvec3 to_sun = glm::transpose(neptune.orientation) * (scene.sun_position() - neptune.world_position);
        const double sub_solar = std::atan2(to_sun.y, to_sun.x) / kDegToRad;
        const double east = dark_spot_state(*shot.spot, radius, e->t_tdb).east_lon_rad / kDegToRad;
        const double off = std::remainder(sub_solar - shot.west_of_center_deg - east, 360.0);
        std::printf("  info: %s: sub-solar longitude %.2f; the spot should be at %.2f, is off by %.2f deg\n", e->name.c_str(), sub_solar, sub_solar - shot.west_of_center_deg, off);
        check(std::abs(off) < 2.0, "the spot opens on the sunlit side", off);
    }
}

void test_ice_giant_aurora()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    const int uranus = scene.find("Uranus");
    const int neptune = scene.find("Neptune");
    const Body& u = scene.bodies[static_cast<size_t>(uranus)];
    const Body& n = scene.bodies[static_cast<size_t>(neptune)];
    check(u.aurora.model == AuroraModel::Uranus && u.aurora.map && n.aurora.model == AuroraModel::Neptune && n.aurora.map,
          "the ice giants have H3+ maps");
    if (!u.aurora.map || !n.aurora.map) return;

    // Melin 2025 Table 1: the sub-JWST longitude (IAU new, east), from the Earth here
    // (JWST is 0.03 deg away). This fixes the rotation's rate, phase and direction.
    constexpr double kJan8 = 2459952.5;
    const struct { double day, hours, east; } table1[] = {
        {0, 15 + 33 / 60.0 + 21 / 3600.0, 194.3}, {0, 15 + 44 / 60.0 + 16 / 3600.0, 198.0},
        {0, 15 + 55 / 60.0 + 11 / 3600.0, 201.8}, {0, 16 + 6 / 60.0 + 6 / 3600.0, 205.6},
        {0, 21 + 2 / 60.0 + 10 / 3600.0, 308.6}, {0, 21 + 13 / 60.0 + 5 / 3600.0, 312.4},
        {0, 21 + 24 / 60.0, 316.2}, {0, 21 + 35 / 60.0, 320.0},
        {1, 2 + 47 / 60.0 + 59 / 3600.0, 68.9}, {1, 2 + 58 / 60.0 + 54 / 3600.0, 72.7},
        {1, 3 + 9 / 60.0 + 49 / 3600.0, 76.5}, {1, 3 + 20 / 60.0 + 44 / 3600.0, 80.3},
    };
    double worst = 0.0;
    SubPoint first{};
    for (const auto& row : table1) {
        const SubPoint p = sub_earth_point(scene, uranus, tdb_from_utc_jd(kJan8 + row.day + row.hours / 24.0));
        worst = std::max(worst, std::abs(angle_difference_deg(p.east_deg, row.east)));
        if (row.east == 194.3) first = p;
    }
    check(worst < 1.0, "Uranus's sub-observer longitude follows Melin 2025 Table 1", worst);
    check(std::abs(first.latitude_deg - 55.1) < 1.5, "sub-observer latitude: the IAU north pole faces us", first.latitude_deg);
    // Tiranti 2026 section 4: the sub-solar latitude on 2025-01-19 is in the same sense.
    const double tiranti_t = tdb_from_utc_jd(2460694.5);
    scene.update(tiranti_t);
    const Body& uranus_body = scene.bodies[static_cast<size_t>(uranus)];
    const glm::dvec3 to_sun = glm::transpose(scene.body_to_icrf(uranus_body, tiranti_t)) * (scene.sun_position() - uranus_body.world_position);
    const double subsolar_deg = std::asin(glm::normalize(to_sun).z) / kDegToRad;
    check(std::abs(subsolar_deg - 67.7) < 1.5, "the sub-solar latitude is Tiranti's", subsolar_deg);
    // The polar night in the map (tools/geomag/ice_giant_h3p_bake.py) is for Melin's dates.
    const double melin_t = tdb_from_utc_jd(kJan8 + 0.875);
    scene.update(melin_t);
    const glm::dvec3 to_sun_2023 = glm::transpose(scene.body_to_icrf(uranus_body, melin_t)) * (scene.sun_position() - uranus_body.world_position);
    const double subsolar_2023_deg = std::asin(glm::normalize(to_sun_2023).z) / kDegToRad;
    check(std::abs(subsolar_2023_deg - 58.83) < 0.05, "the map's polar night has the sub-solar latitude of Melin's dates", subsolar_2023_deg);

    // The map, in ULS west longitude (IAU new east + 50.6) and IAU latitude.
    const AuroraMap& map = *u.aurora.map;
    const AuroraProfile uranus_profile = aurora_profile(AuroraModel::Uranus, false, 0.0);
    const auto uranus_with_sun = [&](double east_deg, double latitude_deg, const glm::dvec3& sun) {
        return aurora_map_pattern(uranus_profile, map, east_deg * kDegToRad, latitude_deg * kDegToRad, sun);
    };
    const glm::dvec3 sun_2023 = glm::normalize(to_sun_2023);
    const auto uranus_2023 = [&](double east_deg, double latitude_deg) { return uranus_with_sun(east_deg, latitude_deg, sun_2023); };
    const auto uranus_at = [&](double uls_west_deg, double latitude_deg) { return uranus_2023(uls_west_deg - 50.6, latitude_deg); };
    // Melin's southern aurora (ULS south is IAU north), and the dark spot just poleward of it.
    check(uranus_at(232.0, 37.5) > 1.3 * uranus_at(217.0, 62.5), "the southern aurora outshines the dark spot beside it",
          uranus_at(232.0, 37.5) / uranus_at(217.0, 62.5));
    // Tiranti's limb: bright near 85 W, an emission dip near 195 W.
    check(uranus_at(85.0, -15.0) > 1.2 * uranus_at(195.0, -15.0), "the limb's emission dip", uranus_at(85.0, -15.0) / uranus_at(195.0, -15.0));
    double peak = 0.0;
    for (int lon = 0; lon < 360; lon += 5) {
        for (int lat = 10; lat <= 85; lat += 5) peak = std::max(peak, uranus_at(lon, lat));
    }
    check(peak > 0.4 && peak < 1.0, "Uranus's total emission is of order 0.3 uW m-2 sr-1", peak);

    // South of 25 S nothing was observed: a model from the observed profiles in latitude
    // and in the magnetic shells of Melin's Q3 tables. Check where its structure lies with
    // Q3's own coefficients, which the tables do not go through (Connerney et al. 1987,
    // the file gsfcq3.dat of libinternalfield, see assets/geomag/SOURCES.md; ULS, like AH5).
    double brightest_south = 0.0, brightest_night = 0.0, brightest_pole = 0.0;
    double brightest_row_share = 0.0;
    glm::dvec3 brightest_direction(0.0);
    MagneticField q3;
    q3.set_coefficient(false, 1, 0, 11893.0);
    q3.set_coefficient(false, 1, 1, 11579.0);
    q3.set_coefficient(true, 1, 1, -15684.0);
    q3.set_coefficient(false, 2, 0, -6030.0);
    q3.set_coefficient(false, 2, 1, -12587.0);
    q3.set_coefficient(true, 2, 1, 6116.0);
    q3.set_coefficient(false, 2, 2, 196.0);
    q3.set_coefficient(true, 2, 2, 4759.0);
    const glm::dmat3 uls_to_body = magnetic_model_to_body(u, 0.0);
    const auto body_direction = [](double east_deg, double latitude_deg) {
        const double east = east_deg * kDegToRad, lat = latitude_deg * kDegToRad;
        return glm::dvec3(std::cos(lat) * std::cos(east), std::cos(lat) * std::sin(east), std::sin(lat));
    };
    for (double lat = -87.5; lat <= -37.5; lat += 2.5) {
        double row_mean = 0.0;
        for (int east = 0; east < 360; east += 5) row_mean += uranus_2023(east, lat) / 72.0;
        for (int east = 0; east < 360; east += 5) {
            const double value = uranus_2023(east, lat);
            brightest_south = std::max(brightest_south, value);
            if (lat <= -62.5) brightest_night = std::max(brightest_night, value);
            if (lat <= -85.0) brightest_pole = std::max(brightest_pole, value);
            // North of this the limb's longitude curve still shows through the smoothing.
            if (lat > -65.0 && lat <= -42.5 && value / row_mean > brightest_row_share) {
                brightest_row_share = value / row_mean;
                brightest_direction = body_direction(east, lat);
            }
        }
    }
    check(brightest_south < peak, "the unobserved south is no brighter than what was observed", brightest_south / peak);
    // The polar night. Non-auroral H3+ comes from the sun and lives about 1 h (Melin 2025):
    // the far south keeps only what particles ionize, a shell-mean far below the bright spots.
    const auto zonal_mean = [&](double latitude_deg) {
        double mean = 0.0;
        for (int east = 0; east < 360; east += 5) mean += uranus_2023(east, latitude_deg) / 72.0;
        return mean;
    };
    check(zonal_mean(-15.0) > zonal_mean(-30.0) && zonal_mean(-30.0) > zonal_mean(-45.0) && zonal_mean(-45.0) > zonal_mean(-60.0),
          "the south darkens away from the sun");
    // At the terminator at noon, against the same distance from the equator in the north.
    constexpr double kNoonTerminatorLatitude = 90.0 - 58.83;
    const double terminator_share = zonal_mean(-kNoonTerminatorLatitude) / zonal_mean(kNoonTerminatorLatitude);
    check(terminator_share > 0.2 && terminator_share < 0.45, "about Tiranti's column near the terminator over Melin's sunlit one (0.29)", terminator_share);
    check(brightest_pole < 0.005, "nothing at the pole of the polar night", brightest_pole);
    check(brightest_night > 0.01 && brightest_night < 0.15, "the aurora alone, weak", brightest_night);
    // Melin 2025 calls L < 3 non-auroral: the map's excess should lie on Q3's shells around
    // L = 3 and up, and not on them if the field is turned the other way.
    FieldLineSettings shell_settings;
    shell_settings.max_radius = 30.0;
    const auto shell_of = [&](const glm::dvec3& position) {
        double shell = 0.0;
        for (const glm::dvec4& point : trace_field_line(q3, position * 1.001, shell_settings)) {
            shell = std::max(shell, glm::length(glm::dvec3(point)));
        }
        return shell;
    };
    const auto shell_contrast = [&](double turn_deg) {
        const double c = std::cos(turn_deg * kDegToRad), s = std::sin(turn_deg * kDegToRad);
        double outer = 0.0, inner = 0.0;
        int outer_cells = 0, inner_cells = 0;
        for (double lat = -62.5; lat <= -42.5; lat += 2.5) {
            double row_mean = 0.0;
            for (int east = 0; east < 360; east += 5) row_mean += uranus_2023(east, lat) / 72.0;
            for (int east = 0; east < 360; east += 10) {
                const glm::dvec3 p = glm::transpose(uls_to_body) * body_direction(east, lat);
                const glm::dvec3 position(c * p.x - s * p.y, s * p.x + c * p.y, p.z);
                if (glm::dot(q3.field_nt(position), position) <= 0.0) continue; // around ULS north only
                const double shell = shell_of(position);
                const double value = uranus_2023(east, lat) / row_mean;
                if (shell >= 2.0 && shell < 6.0) {
                    outer += value;
                    ++outer_cells;
                } else if (shell >= 1.1 && shell < 1.6) {
                    inner += value;
                    ++inner_cells;
                }
            }
        }
        return outer_cells > 10 && inner_cells > 10 ? (outer / outer_cells) / (inner / inner_cells) : 0.0;
    };
    const double contrast = shell_contrast(0.0);
    const double turned_contrast = shell_contrast(180.0);
    const double brightest_shell = shell_of(glm::transpose(uls_to_body) * brightest_direction);
    std::printf("  info: Uranus's modeled south: up to %.3f, in the night %.3f, noon terminator share %.2f; L 2-6 over L 1.1-1.6 %.3f (field turned 180 deg: %.3f); brightest for its row %.2f, on L = %.1f\n",
                brightest_south, brightest_night, terminator_share, contrast, turned_contrast, brightest_row_share, brightest_shell);
    check(contrast > 1.05, "the modeled south is brighter on Q3's auroral shells", contrast);
    check(turned_contrast < 1.0, "and not with the field turned away", turned_contrast);
    check(brightest_shell > 2.0 && brightest_shell < 6.0, "the brightest of the modeled south is on those shells", brightest_shell);
    // Nothing of a seam where the model takes over: no steeper than the observed map itself.
    const auto steepest_step = [&](double from_lat, double to_lat) {
        double steepest = 0.0;
        for (int east = 0; east < 360; east += 5) {
            for (double lat = from_lat; lat < to_lat; lat += 2.5) {
                steepest = std::max(steepest, std::abs(uranus_2023(east, lat) -
                                                       uranus_2023(east, lat + 2.5)));
            }
        }
        return steepest;
    };
    check(steepest_step(-45.0, -15.0) < steepest_step(10.0, 80.0), "no seam where the model takes over", steepest_step(-45.0, -15.0));

    // The seasons: the same two layers by the sun of the moment, from the scene's own sun.
    const auto sun_at = [&](double t) {
        scene.update(t);
        return glm::dvec3(glm::normalize(glm::transpose(scene.body_to_icrf(uranus_body, t)) * (scene.sun_position() - uranus_body.world_position)));
    };
    const auto zonal_with_sun = [&](double latitude_deg, const glm::dvec3& sun) {
        double mean = 0.0;
        for (int east = 0; east < 360; east += 5) mean += uranus_with_sun(east, latitude_deg, sun) / 72.0;
        return mean;
    };
    check(std::abs(uranus_profile.map_sun.z / kDegToRad - 120.2) < 0.1, "no sun at 4000 km beyond 120.2 deg of zenith angle", uranus_profile.map_sun.z / kDegToRad);
    // Voyager 2, 1986-01-24: the sun over ULS north, the IAU's south. The other way round.
    const glm::dvec3 sun_1986 = sun_at(tdb_from_utc_jd(2446454.5));
    const double subsolar_1986_deg = std::asin(sun_1986.z) / kDegToRad;
    check(subsolar_1986_deg < -70.0, "the sun was over the IAU's south at Voyager 2", subsolar_1986_deg);
    check(zonal_with_sun(-60.0, sun_1986) > 0.25, "the far south is lit in 1986", zonal_with_sun(-60.0, sun_1986));
    check(zonal_with_sun(60.0, sun_1986) < 0.05, "the far north is in its polar night in 1986", zonal_with_sun(60.0, sun_1986));
    // What was seen of the aurora stays in the night, and only that.
    double aurora_night = 0.0;
    for (double west = 215.0; west <= 250.0; west += 2.5) {
        for (double lat = 25.0; lat <= 50.0; lat += 2.5) aurora_night = std::max(aurora_night, uranus_with_sun(west - 50.6, lat, sun_1986));
    }
    const double beside_aurora = uranus_with_sun(310.0 - 50.6, 60.0, sun_1986);
    check(aurora_night > 0.08 && aurora_night < peak, "Melin's southern aurora stays in the polar night", aurora_night);
    check(beside_aurora < 0.005, "with nothing around it", beside_aurora);
    // The equinox of December 2007: both hemispheres lit alike, by the same rule.
    const glm::dvec3 sun_2007 = sun_at(tdb_from_utc_jd(2454441.5));
    check(std::abs(std::asin(sun_2007.z) / kDegToRad) < 2.0, "the sun on Uranus's equator in December 2007", std::asin(sun_2007.z) / kDegToRad);
    for (const double lat : {30.0, 45.0}) {
        const double ratio = zonal_with_sun(lat, sun_2007) / zonal_with_sun(-lat, sun_2007);
        check(ratio > 0.75 && ratio < 1.35, "the hemispheres alike at the equinox", ratio);
    }
    check(zonal_with_sun(85.0, sun_2007) > 0.05 && zonal_with_sun(-85.0, sun_2007) > 0.05, "the poles grazed by the sun at the equinox");
    std::printf("  info: Uranus's seasons: sub-solar %.1f deg at Voyager 2, aurora in the night %.3f; at the equinox 60 N %.3f, 60 S %.3f\n",
                subsolar_1986_deg, aurora_night, zonal_with_sun(60.0, sun_2007), zonal_with_sun(-60.0, sun_2007));

    // The vertical profile holds a whole column; its emission peaks well below the temperature's.
    const AuroraProfile& profile = uranus_profile;
    double column = 0.0;
    size_t top = 0;
    for (size_t k = 0; k + 1 < profile.table.size(); ++k) {
        column += 0.5 * (profile.table[k] + profile.table[k + 1]) * profile.layers.y;
        if (profile.table[k] > profile.table[top]) top = k;
    }
    check(std::abs(column - 1.0) < 0.05, "the table is a share of the column per km", column);
    const double top_km = profile.layers.x + profile.layers.y * static_cast<double>(top);
    check(top_km > 1000.0 && top_km < 2000.0, "emission peaks near the density peak", top_km);
    check(profile.top_km > 5000.0 && profile.bottom_km < 475.0, "the shell holds the observed altitudes");

    // The event is at the observation: the longitude Melin's first dither saw.
    const Scene uranus_scene = load_scene_or_die("uranus.toml");
    const Scene neptune_scene = load_scene_or_die("neptune.toml");
    const auto find_event = [](const Scene& where, const char* text) {
        return std::find_if(where.events.begin(), where.events.end(),
                            [&](const SceneEvent& e) { return e.name.find(text) != std::string::npos; });
    };
    const auto uranus_event = find_event(uranus_scene, "Uranus's infrared");
    const auto neptune_event = find_event(neptune_scene, "Neptune's first infrared");
    check(uranus_event != uranus_scene.events.end() && neptune_event != neptune_scene.events.end(), "the ice giants' events");
    if (uranus_event == uranus_scene.events.end() || neptune_event == neptune_scene.events.end()) return;
    for (const auto& e : {*uranus_event, *neptune_event}) {
        check(e.aurora_infrared && !e.aurora_ultraviolet && !e.aurora_visible, "the ice giants' events are infrared");
    }
    // These scenes have no Earth: the cameras are written as offsets toward it.
    const auto check_earth_view = [&](const SceneEvent& e, int planet) {
        scene.update(e.t_tdb);
        const glm::dvec3 to_earth = scene.bodies[static_cast<size_t>(scene.find("Earth"))].icrf_position - scene.bodies[static_cast<size_t>(planet)].icrf_position;
        const glm::dvec3 wanted = glm::normalize(to_earth) * 150000.0;
        const double off_deg = std::acos(std::min(1.0, glm::dot(glm::normalize(wanted), glm::normalize(e.observer_offset_km)))) / kDegToRad;
        std::printf("  info: %s: toward the Earth [%.1f, %.1f, %.1f], off by %.4f deg\n", e.name.c_str(), wanted.x, wanted.y, wanted.z, off_deg);
        check(off_deg < 0.01 && std::abs(glm::length(e.observer_offset_km) - 150000.0) < 1.0, "the event's camera is toward the Earth", off_deg);
    };
    check_earth_view(*uranus_event, uranus);
    check_earth_view(*neptune_event, neptune);
    const double uranus_observed = tdb_from_utc_jd(kJan8 + (15 + 33 / 60.0 + 21 / 3600.0) / 24.0);
    check(std::abs(uranus_event->t_tdb - (uranus_observed - first.light_s)) < 300.0, "Uranus's event less the light time",
          uranus_event->t_tdb - (uranus_observed - first.light_s));

    // Neptune: Melin 2025 gives the IAU longitude at the middle of the second view,
    // whose phase is otherwise unknown; the map is aligned to the scene's rotation there.
    constexpr double kNeptuneView2Jd = 2460117.5 + (14 + 1 / 60.0 + 13 / 3600.0) / 24.0;
    constexpr double kNeptuneView2West = 303.2;
    const SubPoint view2 = sub_earth_point(scene, neptune, tdb_from_utc_jd(kNeptuneView2Jd));
    const double scene_minus_paper = std::remainder(-view2.east_deg - kNeptuneView2West, 360.0);
    std::printf("info: Neptune's scene west longitude minus Melin's %.2f deg\n", scene_minus_paper);
    check(std::abs(neptune_event->t_tdb - (tdb_from_utc_jd(kNeptuneView2Jd) - view2.light_s)) < 300.0, "Neptune's event less the light time",
          neptune_event->t_tdb - (tdb_from_utc_jd(kNeptuneView2Jd) - view2.light_s));
    const AuroraMap& spots = *n.aurora.map;
    double brightest = 0.0, brightest_west = 0.0, brightest_lat = 0.0;
    for (int j = 0; j < spots.height; ++j) {
        for (int i = 0; i < spots.width; ++i) {
            const double v = spots.values[static_cast<size_t>(j * spots.width + i)];
            if (v > brightest) {
                brightest = v;
                brightest_west = (i + 0.5) * 360.0 / spots.width;
                brightest_lat = -90.0 + (j + 0.5) * 180.0 / spots.height;
            }
        }
    }
    const double paper_west = std::fmod(brightest_west - scene_minus_paper + 720.0, 360.0);
    check(paper_west > 200.0 && paper_west < 280.0 && brightest_lat > -60.0 && brightest_lat < -30.0,
          "Neptune's patch is where Melin 2025 finds it", paper_west);
    check(brightest > 1.4 && brightest < 2.2, "about twice as bright as its surroundings", brightest);
    check(std::abs(spots.brightness(0.0, 40.0 * kDegToRad) - 1.0) < 0.01, "no patch in the north, where none was seen");

    // The background rises from dawn to dusk (Melin 2025), fixed to the sun; the patch is not.
    const AuroraProfile neptune_profile = aurora_profile(AuroraModel::Neptune, false, 0.0);
    const double dusk_gradient = neptune_profile.map_sun.x;
    check(dusk_gradient > 0.2 && dusk_gradient < 0.27, "Neptune's gradient as fitted to each view alone (0.21, 0.27)", dusk_gradient);
    const auto neptune_with_sun = [&](double east_deg, double latitude_deg, double sun_east_deg) {
        const glm::dvec3 sun(std::cos(sun_east_deg * kDegToRad), std::sin(sun_east_deg * kDegToRad), 0.0);
        return aurora_map_pattern(neptune_profile, spots, east_deg * kDegToRad, latitude_deg * kDegToRad, sun);
    };
    check(std::abs(neptune_with_sun(30.0, 30.0, 30.0) - 1.0) < 0.01, "no gradient at noon");
    check(std::abs(neptune_with_sun(120.0, 0.0, 30.0) - (1.0 + dusk_gradient)) < 0.01, "brightest at the dusk limb, east of noon");
    check(std::abs(neptune_with_sun(-60.0, 0.0, 30.0) - (1.0 - dusk_gradient)) < 0.01, "dimmest at the dawn limb, west of noon");
    check(std::abs(neptune_with_sun(120.0, 0.0, 120.0) - 1.0) < 0.01, "the gradient goes with the sun, not the planet");
    // The patch adds the same amount wherever the sun is.
    const double patch_east = -brightest_west, patch_lat = brightest_lat;
    const double patch_excess = spots.brightness(patch_east * kDegToRad, patch_lat * kDegToRad) - 1.0;
    check(std::abs(neptune_with_sun(patch_east, patch_lat, patch_east) - (1.0 + patch_excess)) < 0.01, "the patch at noon");
    check(neptune_with_sun(patch_east, patch_lat, patch_east - 90.0) > 1.0 + patch_excess, "the patch on the dusk side");
    // East of noon is afternoon only if the sun moves west over the body: from the scene.
    const auto neptune_sun_east = [&](double t) {
        scene.update(t);
        const Body& neptune_body = scene.bodies[static_cast<size_t>(neptune)];
        const glm::dvec3 sun = glm::transpose(scene.body_to_icrf(neptune_body, t)) * (scene.sun_position() - neptune_body.world_position);
        return std::atan2(sun.y, sun.x) / kDegToRad;
    };
    const double sun_motion = angle_difference_deg(neptune_sun_east(neptune_event->t_tdb + 3600.0), neptune_sun_east(neptune_event->t_tdb));
    check(sun_motion < -10.0, "the sun moves west over Neptune", sun_motion);
}

void test_europa_aurora()
{
    auto scene = load_scene_or_die("jupiter.toml");
    const int index = scene.find("Europa");
    check(index >= 0 && scene.bodies[index].aurora.model == AuroraModel::Europa, "Europa visible aurora loads");
    const auto solar = load_scene_or_die("solar_system.toml");
    check(solar.bodies[solar.find("Europa")].aurora.model == AuroraModel::Europa, "Europa template shared across scenes");
    // Independent line-of-sight integration, including the off-limb shell.
    constexpr double R = 1560.8, top = R + 600.0;
    double flux = 0.0;
    constexpr int radial = 4000, depth = 400;
    for (int i = 0; i < radial; ++i) {
        const double b2 = top * top * (i + 0.5) / radial;
        const double z1 = std::sqrt(top * top - b2);
        const double z0 = b2 < R * R ? std::sqrt(R * R - b2) : -z1;
        double column = 0.0;
        for (int j = 0; j < depth; ++j) {
            const double z = z0 + (z1 - z0) * (j + 0.5) / depth;
            const double h = std::sqrt(b2 + z * z) - R;
            column += std::exp(-h / 150.0) / (150.0 * (1.0 - std::exp(-4.0))) * (z1 - z0) / depth;
        }
        flux += column * top * top / (R * R * radial);
    }
    check(std::abs(flux - aurora_disk_factor(R, 150.0, 600.0)) < 0.003,
          "Europa aperture brightness agrees with independent ray integral", flux);
    int count = 0;
    for (const auto& event : scene.events) {
        if (event.focus != index || event.name.find("oxygen aurora") == std::string::npos) continue;
        ++count;
        scene.update(event.t_tdb);
        const auto& moon = scene.bodies[index];
        const auto& parent = scene.bodies[moon.parent];
        const auto light = glm::normalize(scene.sun_position() - moon.world_position);
        const auto oc = parent.world_position - moon.world_position;
        const double along = glm::dot(oc, light);
        const double miss = glm::length(oc - along * light);
        check(along > 0.0 && miss + R + 600.0 + along * 695700.0 / glm::length(scene.sun_position() - moon.world_position) < parent.polar_radius_km,
              "Europa aurora tour lies inside Jupiter shadow", miss);
        check(event.aurora_visible && !event.aurora_ultraviolet, "Europa tour uses visible spectrum");
        scene.update(event.t_tdb + 3.0 * 3600.0);
        const auto day_light = glm::normalize(scene.sun_position() - moon.world_position);
        const auto day_oc = parent.world_position - moon.world_position;
        check(glm::length(glm::cross(day_oc, day_light)) > parent.equatorial_radius_km + R + 600.0,
              "Europa leaves Jupiter shadow three hours after tour");
    }
    check(count == 1, "Europa has an eclipse tour");
    Scene bad;
    std::string error;
    const std::string base = "name='x'\n[[bodies]]\nname='Europa'\nradii_km=[1560.8]\n[bodies.aurora]\n";
    for (const char* conflict : {"kp='geomag/kp.bin'", "saturn_oval=true", "ganymede_map='x'", "jovian_oval='x'"})
        check(!load_scene_string(base + "europa_visible=true\n" + conflict, "mixed_europa", bad, &error, ASTRAXIS_ASSET_DIR),
              "Europa rejects conflicting aurora models");
    check(!load_scene_string(base + "europa_visible='true'", "bad_europa", bad, &error, ASTRAXIS_ASSET_DIR),
          "Europa validates boolean type");
}

void test_ganymede_aurora()
{
    Scene scene = load_scene_or_die("jupiter.toml");
    const auto& moon = scene.bodies[static_cast<size_t>(scene.find("Ganymede"))];
    check(bool(moon.aurora.map) && moon.aurora.model == AuroraModel::Ganymede,
          "Ganymede loads its own UV map independently of Jupiter's footprint");
    if (!moon.aurora.map) return;
    const auto& map = *moon.aurora.map;
    Scene solar = load_scene_or_die("solar_system.toml");
    const auto& other = solar.bodies[static_cast<size_t>(solar.find("Ganymede"))].aurora.map;
    check(other && other->values == map.values, "Ganymede UV map is shared across scenes");
    // Independent published morphology: Marzok 2022 section 3.3, Fig. 8.
    // Leading (90 W) crescents lie closer to the equator than trailing (270 W).
    for (int sign : {-1, 1}) {
        double peaks[2]{};
        for (int h = 0; h < 2; ++h) {
            double maximum = -1.0;
            for (double lat = 5.0; lat <= 65.0; lat += 0.25) {
                const double b = map.brightness(-(90.0 + h * 180.0) * kDegToRad, sign * lat * kDegToRad);
                if (b > maximum) { maximum = b; peaks[h] = lat; }
            }
        }
        check(peaks[0] < 30.0 && peaks[1] > 30.0 && peaks[1] > peaks[0] + 10.0,
              "HST leading crescents reach lower latitude than trailing crescents", peaks[1] - peaks[0]);
    }
    for (double lat : {-90.0, -42.75, 0.0, 20.25, 90.0}) {
        check(std::abs(map.brightness(-1e-9, lat * kDegToRad) - map.brightness(1e-9, lat * kDegToRad)) < 1e-5,
              "Ganymede UV map is continuous across the prime meridian");
        check(std::abs(map.brightness(0.3, lat * kDegToRad) - map.brightness(0.3 + kTwoPi, lat * kDegToRad)) < 1e-10,
              "Ganymede east-to-west lookup wraps longitude");
    }
    const auto pixels = map.texture_rgba();
    for (int y = 0; y < map.height; ++y) {
        for (int x = 0; x < map.width; ++x) {
            const size_t k = static_cast<size_t>(y * map.width + x);
            const double b = map.brightness(-(2.25 + 4.5 * x) * kDegToRad, (-87.75 + 4.5 * y) * kDegToRad);
            check(std::abs(b - map.values[k]) < 1e-10, "UV cell centers preserve west longitude and latitude");
            check(std::abs(pixels[4 * k] * map.range / 255.0 - b) < 0.52,
                  "GPU UV texture quantization is about half a Rayleigh");
            check(std::abs(pixels[4 * k + 3] * map.range / 255.0 - map.more_layers[2][k]) < 0.52,
                  "the fourth layer is in the alpha channel");
        }
    }
    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error), "Ganymede UV simulation loads");
    sim.stop_tour();
    int count = 0;
    for (size_t i = 0; i < sim.scene().events.size(); ++i) {
        const auto& event = sim.scene().events[i];
        // The event with the field lines looks at them, not at one hemisphere's aurora.
        if (!event.aurora_ultraviolet || event.focus != sim.scene().find("Ganymede") || event.magnetic_lines) continue;
        ++count;
        sim.jump_to_event(i);
        sim.update(0.0);
        sim.camera().update(3.0, sim.scene());
        OutputView view;
        sim.compute_view(1280, 720, 1.0, view);
        const auto& body = sim.scene().bodies[static_cast<size_t>(event.focus)];
        const glm::dvec3 direction = glm::transpose(body.orientation) * (view.camera.position - body.world_position);
        const double east = std::atan2(direction.y, direction.x) / kDegToRad;
        check(std::abs(std::remainder(east - event.phase_deg, 360.0)) < 3.0,
              "Ganymede UV tours face their intended leading or trailing hemisphere", east);
        check(sim.scene().active_event() == static_cast<int>(i), "Ganymede event activates UV overlay");
    }
    check(count == 2, "Ganymede has leading and trailing UV tours");

    // Marzok et al. 2022 Table 1: the current sheet's dipole, System III west longitude and
    // magnetic latitude of a few exposures, across its whole range.
    constexpr std::pair<double, double> kMarzokRows[] = {{229.7, 8.32}, {291.8, -0.16}, {331.0, -6.13}, {26.49, -9.45},
                                                         {197.2, 9.48}, {145.2, 5.36}, {263.8, 4.31}, {345.6, -7.77}};
    for (const auto& [west, latitude] : kMarzokRows) {
        const double east = -west * kDegToRad;
        check(std::abs(jovian_magnetic_latitude_deg(glm::dvec3(std::cos(east), std::sin(east), 0.0)) - latitude) < 0.08,
              "the dipole gives Marzok's magnetic latitudes", latitude);
    }
    // The layers: the file holds Figure 5 and the three groups of Figure 6.
    check(map.more_layers.size() == 3, "Ganymede's map has the current-sheet groups");
    const auto layer_value = [&](size_t layer, double west, double latitude) {
        return map.layer_brightness(layer, -west * kDegToRad, latitude * kDegToRad);
    };
    // The time Ganymede spends in each group: the groups average back to the long-term map,
    // and a cell that was not observed stays dark in all of them.
    std::array<double, 3> time_shares{};
    constexpr int kLongitudes = 3600;
    for (int k = 0; k < kLongitudes; ++k) {
        const double east = kTwoPi * (k + 0.5) / kLongitudes;
        const auto shares = ganymede_sheet_shares(jovian_magnetic_latitude_deg(glm::dvec3(std::cos(east), std::sin(east), 0.0)));
        for (size_t g = 0; g < 3; ++g) time_shares[g] += shares[g] / kLongitudes;
    }
    std::printf("  info: Ganymede's time above, inside and below the sheet: %.3f %.3f %.3f (HST's exposures 0.304 0.413 0.283)\n",
                time_shares[0], time_shares[1], time_shares[2]);
    check(std::abs(time_shares[0] - 14.0 / 46.0) < 0.03 && std::abs(time_shares[1] - 19.0 / 46.0) < 0.03,
          "Ganymede's time in each group is about HST's share of exposures");
    double worst_mean = 0.0;
    bool dark_stays_dark = true;
    for (size_t k = 0; k < map.values.size(); ++k) {
        double mean = 0.0;
        for (size_t g = 0; g < 3; ++g) {
            mean += time_shares[g] * map.more_layers[g][k];
            dark_stays_dark = dark_stays_dark && (map.values[k] > 0.0 || map.more_layers[g][k] == 0.0);
        }
        worst_mean = std::max(worst_mean, std::abs(mean - map.values[k]));
    }
    check(worst_mean < 2.0, "the groups' time average is the long-term map to 2 R", worst_mean);
    check(dark_stays_dark, "cells without data stay dark in every group");
    // Marzok section 3.3: inside the sheet the downstream (leading) aurora brightens by
    // 1.3 +- 0.31 and the upstream (trailing) one dims by 0.78 +- 0.19, in the ovals' bands.
    const std::tuple<double, double, double> kOvals[] = {{32.3, -16.9, 1.5}, {-29.4, 11.6, 8.8}};
    const auto band_mean = [&](size_t layer, double west_low, double west_high) {
        double sum = 0.0;
        int count = 0;
        for (double west = west_low; west <= west_high; west += 4.5) {
            for (const auto& [mean, amplitude, phase] : kOvals) {
                const double center = mean + amplitude * std::sin((west + phase) * kDegToRad);
                for (int step = -3; step <= 3; ++step) {
                    const double b = layer_value(layer, west, center + step * 4.5);
                    if (b > 0.0) {
                        sum += b;
                        ++count;
                    }
                }
            }
        }
        return sum / count;
    };
    const auto inside_over_outside = [&](double west_low, double west_high) {
        const double above = band_mean(1, west_low, west_high), below = band_mean(3, west_low, west_high);
        return band_mean(2, west_low, west_high) / ((14.0 * above + 13.0 * below) / 27.0);
    };
    const double downstream = inside_over_outside(40.0, 150.0), upstream = inside_over_outside(240.0, 340.0);
    std::printf("  info: Ganymede inside/outside the sheet: downstream %.2f, upstream %.2f\n", downstream, upstream);
    check(std::abs(downstream - 1.3) < 0.31 && downstream > 1.1, "the leading aurora brightens inside the sheet", downstream);
    check(std::abs(upstream - 0.78) < 0.19 && upstream < 0.95, "the trailing aurora dims inside the sheet", upstream);
    // Between 180 and 210 W nearly only the exposures above the sheet looked, so the others
    // are the long-term map (which is that group) there, without a seam around the gap.
    double gap_ratio = 0.0;
    for (double latitude = -40.0; latitude <= 40.0; latitude += 4.5) {
        const double main_value = layer_value(0, 195.0, latitude);
        if (main_value > 5.0) gap_ratio = std::max(gap_ratio, std::abs(layer_value(2, 195.0, latitude) / main_value - 1.0));
    }
    check(gap_ratio < 0.45, "a group's gap falls back on the long-term map", gap_ratio);
    double step = 0.0;
    for (size_t g = 1; g <= 3; ++g) {
        for (double latitude = -60.0; latitude <= 60.0; latitude += 4.5) {
            for (double west = 2.25; west < 360.0; west += 4.5) {
                const double a = layer_value(0, west, latitude), b = layer_value(0, west + 4.5, latitude);
                if (a < 20.0 || b < 20.0) continue;
                const double ratio = layer_value(g, west, latitude) / a, next = layer_value(g, west + 4.5, latitude) / b;
                step = std::max(step, std::abs(next / ratio - 1.0));
            }
        }
    }
    check(step < 0.35, "a group's ratio to the long-term map changes smoothly from cell to cell", step);

    // The view: which layer a moment shows, and the mean at a high warp.
    sim.stop_tour();
    const int moon_index = sim.scene().find("Ganymede");
    sim.clock().warp = 30.0;
    double t = sim.clock().t_tdb;
    bool saw[3] = {false, false, false};
    for (int k = 0; k < 40; ++k) {
        sim.clock().t_tdb = t + k * 1800.0;
        sim.update(0.0);
        const MoonAuroraState state = sim.moon_aurora(moon_index);
        const double latitude = jovian_magnetic_latitude_deg(sim.scene().bodies[static_cast<size_t>(moon_index)].in_parent_axes);
        check(std::abs(state.map_shares.x + state.map_shares.y + state.map_shares.z + state.map_shares.w - 1.0) < 1e-12,
              "the shares sum to one");
        if (latitude > 7.0) { saw[0] = true; check(state.map_shares.y == 1.0, "above the sheet: that map"); }
        if (std::abs(latitude) < 5.0) { saw[1] = true; check(state.map_shares.z == 1.0, "inside the sheet: that map"); }
        if (latitude < -7.0) { saw[2] = true; check(state.map_shares.w == 1.0, "below the sheet: that map"); }
    }
    check(saw[0] && saw[1] && saw[2], "Ganymede goes through all three groups");
    sim.clock().warp = 86400.0;
    sim.update(0.0);
    check(sim.moon_aurora(moon_index).map_shares == glm::dvec4(1.0, 0.0, 0.0, 0.0),
          "at a day per second the long-term map is drawn, not flicker");
    MoonAuroraState inside_sheet;
    inside_sheet.map_shares = glm::dvec4(0.0, 0.0, 1.0, 0.0);
    const AuroraProfile shares_profile = aurora_profile(AuroraModel::Ganymede, true, 0.0, inside_sheet);
    check(std::abs(aurora_map_pattern(shares_profile, map, -90.0 * kDegToRad, 20.0 * kDegToRad, glm::dvec3(1.0, 0.0, 0.0)) -
                   layer_value(2, 90.0, 20.0) * shares_profile.map_gain) < 1e-12,
          "the pattern reads the layers by their shares");
    Scene bad;
    const std::string base = "name='x'\n[[bodies]]\nname='Ganymede'\nradii_km=[2631.2]\n[bodies.aurora]\n";
    const std::string field = "ganymede_map='geomag/ganymede_aurora.txt'\n";
    for (const char* conflict : {"kp='geomag/kp.bin'", "jovian_oval='geomag/jupiter_oval.txt'",
                                "saturn_oval=false", "footprints='geomag/jupiter_footprints.txt'"})
        check(!load_scene_string(base + field + conflict, "mixed_ganymede", bad, &error, ASTRAXIS_ASSET_DIR),
              "Ganymede map rejects conflicting aurora fields");
    check(!load_scene_string(base + "ganymede_map=true", "bad_map_type", bad, &error, ASTRAXIS_ASSET_DIR),
          "Ganymede map requires a file path");
    check(!load_scene_string(base + "ganymede_map='geomag/jupiter_oval.txt'", "bad_map", bad, &error, ASTRAXIS_ASSET_DIR),
          "Ganymede map rejects wrong resource format");
    check(!load_scene_string("name='x'\n[[bodies]]\nname='Earth'\nradii_km=[6378]\n[bodies.aurora]\n" + field,
                            "wrong_moon", bad, &error, ASTRAXIS_ASSET_DIR), "Ganymede map checks the body's radius");
}

void test_magnetic_field()
{
    // The expansion itself: a pure dipole, then the harmonic field of a high degree.
    MagneticField dipole;
    dipole.set_coefficient(false, 1, 0, -30000.0);
    const glm::dvec3 equator = dipole.field_nt({1.0, 0.0, 0.0});
    check(std::abs(equator.z - 30000.0) < 1e-6 && std::abs(equator.x) < 1e-6, "a g10 < 0 dipole points north at the equator");
    check(std::abs(dipole.field_nt({0.0, 0.0, 1.0}).z + 60000.0) < 1e-3, "and is twice as strong at the pole");
    check(std::abs(dipole.field_nt({2.0, 0.0, 0.0}).z - 30000.0 / 8.0) < 1e-6, "and falls as r^-3");

    // A dipole's line: r = L cos^2(latitude), so it lands at cos^2(latitude) = 1/L.
    const FieldLineSettings plain;
    const std::vector<glm::dvec4> shell = trace_field_line(dipole, {4.0, 0.0, 0.0}, plain);
    check(shell.size() > 10 && shell.size() < 400, "a line is a few dozen points", static_cast<double>(shell.size()));
    double worst_shape = 0.0;
    for (const glm::dvec4& p : shell) {
        const double r = glm::length(glm::dvec3(p));
        const double cos_latitude2 = 1.0 - p.z * p.z / (r * r);
        worst_shape = std::max(worst_shape, std::abs(r - 4.0 * cos_latitude2));
    }
    check(worst_shape < 1e-3, "the line follows r = L cos^2(latitude)", worst_shape);
    const glm::dvec3 end_a(shell.front());
    const glm::dvec3 end_b(shell.back());
    check(std::abs(glm::length(end_a) - 1.0) < 1e-2 && std::abs(glm::length(end_b) - 1.0) < 1e-2 && end_a.z * end_b.z < 0.0,
          "it ends on the surface in both hemispheres");
    check(std::abs(std::abs(end_a.z) - std::sqrt(0.75)) < 2e-3, "at latitude 60 for L = 4", std::abs(end_a.z));

    Scene solar = load_scene_or_die("solar_system.toml");
    Scene saturn_scene = load_scene_or_die("saturn.toml");
    Scene earth_scene = load_scene_or_die("earth_moon.toml");
    Scene jupiter_scene = load_scene_or_die("jupiter.toml");
    const Body& uranus = solar.bodies[static_cast<size_t>(solar.find("Uranus"))];
    const Body& neptune = solar.bodies[static_cast<size_t>(solar.find("Neptune"))];
    const Body& saturn = saturn_scene.bodies[static_cast<size_t>(saturn_scene.find("Saturn"))];
    const Body& earth = earth_scene.bodies[static_cast<size_t>(earth_scene.find("Earth"))];
    const Body& jupiter = jupiter_scene.bodies[static_cast<size_t>(jupiter_scene.find("Jupiter"))];
    for (const Body* body : {&uranus, &neptune, &saturn, &earth, &jupiter}) {
        check(body->magnetic_field.model != nullptr, "the planet has a field model");
    }
    if (!uranus.magnetic_field.model || !neptune.magnetic_field.model || !saturn.magnetic_field.model ||
        !earth.magnetic_field.model || !jupiter.magnetic_field.model) {
        return;
    }
    const MagneticField& uranus_model = *uranus.magnetic_field.model;
    const MagneticField& neptune_model = *neptune.magnetic_field.model;
    const MagneticField& saturn_model = *saturn.magnetic_field.model;
    const MagneticField& earth_model = *earth.magnetic_field.model;
    const MagneticField& jupiter_model = *jupiter.magnetic_field.model;

    // B = -grad V with V harmonic: no divergence and no curl, in every degree and order.
    for (const MagneticField* model : {&jupiter_model, &earth_model, &uranus_model}) {
        const glm::dvec3 x(1.3, 0.7, -0.6);
        constexpr double kStep = 1e-3;
        glm::dvec3 d[3];
        for (int axis = 0; axis < 3; ++axis) {
            glm::dvec3 step(0.0);
            step[axis] = kStep;
            d[axis] = (model->field_nt(x + step) - model->field_nt(x - step)) / (2.0 * kStep);
        }
        const double size = glm::length(model->field_nt(x));
        const double divergence = d[0].x + d[1].y + d[2].z;
        const double curl = glm::length(glm::dvec3(d[1].z - d[2].y, d[2].x - d[0].z, d[0].y - d[1].x));
        check(std::abs(divergence) < 1e-3 * size, "the field has no divergence", divergence / size);
        check(curl < 1e-3 * size, "and no curl", curl / size);
    }

    // The models against their papers.
    check(std::abs(neptune_model.dipole_tilt_deg() - 46.9) < 0.1, "O8: dipole tilted 46.9 deg (A&A 2024, OTD-O8)", neptune_model.dipole_tilt_deg());
    check(std::abs(neptune_model.dipole_nt() - 14200.0) < 100.0, "O8: dipole 0.142 G R^3 (Connerney 1991)", neptune_model.dipole_nt());
    const glm::dvec3 neptune_offset = neptune_model.dipole_offset();
    check(glm::length(neptune_offset - glm::dvec3(0.05, 0.48, 0.00)) < 0.01, "O8: its center is 0.48 R off (A&A 2024, OTD-O8)", glm::length(neptune_offset));
    check(std::abs(uranus_model.dipole_tilt_deg() - 59.85) < 0.05, "AH5: dipole tilted about 59 deg", uranus_model.dipole_tilt_deg());
    check(std::abs(jupiter_model.dipole_tilt_deg() - 10.25) < 0.05, "JRM33: dipole tilted about 10 deg", jupiter_model.dipole_tilt_deg());
    check(saturn_model.dipole_tilt_deg() < 0.007, "Cassini 11: no tilt (Cao 2020: below 0.007 deg)", saturn_model.dipole_tilt_deg());
    check(saturn_model.dipole_offset().z * saturn.magnetic_field.reference_radius_km > 2000.0,
          "Cassini 11: the magnetic equator is north of the planet's", saturn_model.dipole_offset().z * saturn.magnetic_field.reference_radius_km);

    // Longitudes. Neptune: NLS' prime meridian is defined by Voyager 2 being at 167.7 deg W at
    // 1989-08-25 03:56:00, its closest approach (the PDS catalog, VG2-N-POS-5-SUMM-NLSCOORDS).
    // The scene's rough trajectory gets there 7 minutes late, so its direction at its own closest
    // approach goes with the planet's rotation at the real one. NLS' longitudes are west.
    {
        const double closest_real = tdb_from_utc_jd(2447763.5 + (3.0 + 56.0 / 60.0) / 24.0);
        const int voyager = solar.find("Voyager 2");
        const int planet = solar.find("Neptune");
        glm::dvec3 relative(1e30);
        double closest_scene = closest_real;
        for (double t = closest_real - 1800.0; t < closest_real + 1800.0; t += 5.0) {
            const glm::dvec3 r = solar.icrf_state_at(voyager, t).position - solar.icrf_state_at(planet, t).position;
            if (glm::length(r) < glm::length(relative)) {
                relative = r;
                closest_scene = t;
            }
        }
        const glm::dmat3 model_to_icrf = solar.body_to_icrf(neptune, closest_real) * magnetic_model_to_body(neptune, closest_real);
        const glm::dvec3 local = glm::transpose(model_to_icrf) * relative;
        const double west = -std::atan2(local.y, local.x) / kDegToRad;
        std::printf("  info: Voyager 2 at Neptune: %.2f deg W in NLS (PDS 167.7), %.0f km away, %.1f min after the real closest approach\n",
                    wrap_two_pi(west * kDegToRad) / kDegToRad, glm::length(relative), (closest_scene - closest_real) / 60.0);
        check(std::abs(angle_difference_deg(west, 167.7)) < 3.0, "Voyager 2's NLS longitude at Neptune", angle_difference_deg(west, 167.7));
    }
    // O8's pole leans toward 72 deg W (Connerney 1991), in the model's frame east is positive.
    check(std::abs(angle_difference_deg(std::atan2(neptune_model.dipole_axis().y, neptune_model.dipole_axis().x) / kDegToRad, -72.0)) < 0.5,
          "O8's pole leans toward 72 deg W");

    // Uranus: the field enters the planet where Melin 2025 sees the southern aurora, near
    // 232 W and IAU latitude 37.5 (ULS south is IAU north), so ULS is placed right.
    {
        const glm::dmat3 to_body = magnetic_model_to_body(uranus, 0.0);
        check(glm::length(to_body * glm::dvec3(0.0, 0.0, 1.0) - glm::dvec3(0.0, 0.0, -1.0)) < 1e-12, "ULS north is the IAU's south");
        const glm::dvec3 pole = to_body * -uranus_model.dipole_axis();
        const double latitude = std::asin(pole.z) / kDegToRad;
        const double uls_west = std::atan2(pole.y, pole.x) / kDegToRad + 50.6;
        const double along = std::abs(angle_difference_deg(uls_west, 232.0)) * std::cos(37.5 * kDegToRad);
        check(std::hypot(along, latitude - 37.5) < 10.0, "AH5's southern magnetic pole is by Melin's southern aurora", std::hypot(along, latitude - 37.5));
        std::printf("  info: AH5's southern pole at %.1f deg W (ULS), IAU latitude %.1f; Melin's aurora 232, 37.5\n",
                    wrap_two_pi(uls_west * kDegToRad) / kDegToRad, latitude);
    }

    // Earth: the same pole as the aurora's, from the same IGRF-14 terms.
    {
        const double t = tdb_from_julian_year(2025.0);
        const glm::dvec3 pole = geomagnetic_north_pole(t);
        const double angle = std::acos(std::clamp(glm::dot(pole, -earth_model.dipole_axis()), -1.0, 1.0)) / kDegToRad;
        check(angle < 0.01, "IGRF-14's dipole is the aurora's geomagnetic pole", angle);
    }

    // Jupiter: the internal field alone lands near the table with the magnetodisc: Io's line
    // at its equatorial longitude led eastward, as test_jovian_footprints' oracle has it.
    {
        const auto& oval = *jupiter.aurora.jovian;
        const double radius = 0.5 * (oval.footpaths[0].min_radius + oval.footpaths[0].max_radius);
        const double latitude = 0.37 * kDegToRad;
        for (int h = 0; h < 2; ++h) {
            glm::dvec3 table;
            check(oval.footprint(0, h != 0, {radius * 71492.0 * std::cos(latitude), 0.0, radius * 71492.0 * std::sin(latitude)}, table),
                  "the table has Io");
            const double lead = JovianAurora::lead_angle(0, h != 0, 0.0);
            FieldLineSettings settings;
            settings.surface_equatorial = (71492.0 + 900.0) / 71492.0;
            settings.surface_polar = (66854.0 + 900.0) / 71492.0;
            const std::vector<glm::dvec4> line = trace_field_line(
                jupiter_model,
                radius * glm::dvec3(std::cos(latitude) * std::cos(lead), std::cos(latitude) * std::sin(lead), std::sin(latitude)),
                settings);
            const glm::dvec3 foot = glm::normalize(glm::dvec3(line.front().z * (h != 0 ? -1.0 : 1.0) > 0.0 ? line.front() : line.back()));
            const double angle = std::acos(std::clamp(glm::dot(foot, table), -1.0, 1.0)) / kDegToRad;
            check(angle < 3.0, "Io's line lands near the table's footprint", angle);
        }
    }

    // What the lines show: from a fixed point 5 radii out, how far the top of the line
    // through it swings as the planet turns. Almost nothing on Saturn, a lot on the ice giants.
    const auto apex_swing_deg = [](const MagneticField& model) {
        double low = 90.0;
        double high = -90.0;
        for (int degree = 0; degree < 360; degree += 15) {
            const double a = degree * kDegToRad;
            const std::vector<glm::dvec4> line = trace_field_line(model, {5.0 * std::cos(a), 5.0 * std::sin(a), 0.0}, FieldLineSettings{});
            const auto apex = *std::max_element(line.begin(), line.end(), [](const glm::dvec4& p, const glm::dvec4& q) {
                return glm::length(glm::dvec3(p)) < glm::length(glm::dvec3(q));
            });
            const double latitude = std::asin(apex.z / glm::length(glm::dvec3(apex))) / kDegToRad;
            low = std::min(low, latitude);
            high = std::max(high, latitude);
        }
        return high - low;
    };
    check(apex_swing_deg(uranus_model) > 50.0, "Uranus's lines swing widely", apex_swing_deg(uranus_model));
    check(apex_swing_deg(neptune_model) > 60.0, "and Neptune's", apex_swing_deg(neptune_model));
    check(apex_swing_deg(jupiter_model) < 30.0, "Jupiter's less", apex_swing_deg(jupiter_model));
    check(apex_swing_deg(saturn_model) < 0.5, "and Saturn's not at all", apex_swing_deg(saturn_model));

    // Mercury: Anderson's offset dipole, whose field is much stronger in the north, inside
    // Winslow's magnetopause at 1.45 radii.
    Scene mercury_scene = load_scene_or_die("mercury.toml");
    const Body& mercury = mercury_scene.bodies[static_cast<size_t>(mercury_scene.find("Mercury"))];
    check(mercury.magnetic_field.model != nullptr, "Mercury has a field model");
    if (!mercury.magnetic_field.model) return;
    const MagneticField& mercury_model = *mercury.magnetic_field.model;
    const double mercury_offset_km = mercury_model.dipole_offset().z * mercury.magnetic_field.reference_radius_km;
    check(std::abs(mercury_offset_km - 484.0) < 15.0, "Mercury's dipole is 484 km north of its center (Anderson 2011)", mercury_offset_km);
    check(mercury_model.dipole_tilt_deg() < 2.5, "and tilted less than 2.5 deg", mercury_model.dipole_tilt_deg());
    const double north_pole = mercury_model.field_nt({0.0, 0.0, 1.0}).z;
    const double south_pole = mercury_model.field_nt({0.0, 0.0, -1.0}).z;
    check(std::abs(north_pole / south_pole) > 3.0, "the field at the north pole is over three times that at the south", north_pole / south_pole);
    for (const auto& line : mercury.magnetic_field.lines) {
        double top = 0.0;
        for (const glm::dvec4& p : line) top = std::max(top, glm::length(glm::dvec3(p)));
        check(top < 1.45, "Mercury's lines stay inside the magnetopause", top);
        const glm::dvec3 a(line.front());
        const glm::dvec3 b(line.back());
        const double north_foot = std::asin((a.z > b.z ? a : b).z / glm::length(a.z > b.z ? a : b)) / kDegToRad;
        const double south_foot = std::asin((a.z > b.z ? b : a).z / glm::length(a.z > b.z ? b : a)) / kDegToRad;
        check(north_foot + south_foot > 10.0, "their feet are displaced to the north", north_foot + south_foot);
    }

    // Ganymede: Kivelson's pole, in the southern hemisphere, 24 deg from the Jupiter-facing
    // meridian toward the trailing side. So the model's frame is the body's: x faces Jupiter
    // and y is the trailing side, as here (east of the sub-Jupiter point).
    const Body& ganymede = jupiter_scene.bodies[static_cast<size_t>(jupiter_scene.find("Ganymede"))];
    check(ganymede.magnetic_field.model != nullptr, "Ganymede has a field model");
    if (!ganymede.magnetic_field.model) return;
    const MagneticField& ganymede_model = *ganymede.magnetic_field.model;
    check(std::abs(ganymede_model.dipole_nt() - 719.0) < 1.0, "Ganymede's dipole is 719 nT (Kivelson 2002)", ganymede_model.dipole_nt());
    const glm::dvec3 ganymede_pole = ganymede_model.dipole_axis();
    check(ganymede_pole.z < 0.0 && std::abs(std::atan2(ganymede_pole.y, ganymede_pole.x) / kDegToRad - 24.0) < 3.0,
          "its pole leans toward Jupiter and the trailing side", std::atan2(ganymede_pole.y, ganymede_pole.x) / kDegToRad);
    {
        const double t = tdb_from_julian_year(2021.5);
        jupiter_scene.update(t);
        const int moon_index = jupiter_scene.find("Ganymede");
        const int planet_index = jupiter_scene.find("Jupiter");
        const glm::dvec3 away = jupiter_scene.icrf_state_at(planet_index, t).position - jupiter_scene.icrf_state_at(moon_index, t).position;
        const glm::dvec3 velocity = jupiter_scene.icrf_state_at(moon_index, t).velocity - jupiter_scene.icrf_state_at(planet_index, t).velocity;
        const glm::dmat3 axes = jupiter_scene.body_to_icrf(ganymede, t);
        check(glm::dot(axes[0], glm::normalize(away)) > 0.99, "the body's x axis faces Jupiter", glm::dot(axes[0], glm::normalize(away)));
        check(glm::dot(axes[1], glm::normalize(velocity)) < -0.99, "and its y axis the trailing side", glm::dot(axes[1], glm::normalize(velocity)));
    }
    // Ganymede stands in Jupiter's field: the scene adds it, uniform, to the dipole. The lines
    // near the equator close; the rest leave along Jupiter's field, where the line from the moon
    // to its footprint goes on.
    // The scene follows Jupiter's field only while someone draws the lines.
    const auto follow_field = [&](double t) {
        ganymede.magnetic_field.lines_in_view = true;
        jupiter_scene.update(t);
    };
    {
        const BodyMagneticField& field = ganymede.magnetic_field;
        const Body& planet = jupiter_scene.bodies[static_cast<size_t>(jupiter_scene.find("Jupiter"))];
        const double t = tdb_from_julian_year(2021.5);
        check(field.ambient_nt == glm::dvec3(0.0) && field.lines.size() == 72, "with its lines off, Ganymede's field is left alone");
        follow_field(t);
        check(!field.lines_in_view, "the request lasts one update");
        check(field.ambient_from_parent && glm::length(field.ambient_nt) > 60.0 && field.model->ambient_nt() == field.ambient_nt,
              "Ganymede's model has Jupiter's field around it", glm::length(field.ambient_nt));
        // The same vector as Jupiter's table gives, turned into the moon's axes.
        glm::dvec3 table_field;
        check(planet.aurora.jovian->field_at_moon(2, planet.jovian_moons[2], table_field), "Jupiter's table has the field at Ganymede");
        const glm::dvec3 expected = glm::transpose(jupiter_scene.body_to_icrf(ganymede, t)) * (jupiter_scene.body_to_icrf(planet, t) * table_field);
        check(glm::length(expected - field.ambient_nt) < 1e-9 * glm::length(expected), "it is the field of Jupiter's table there");
        check(field.ambient_nt.z < -60.0, "mostly southward, against the dipole's field at the equator", field.ambient_nt.z);
        // Closed lines end on the moon at both ends; open ones at max_radius, along the ambient field.
        const auto on_surface = [](const glm::dvec4& p) { return glm::length(glm::dvec3(p)) < 1.01; };
        int closed = 0, open = 0;
        double farthest_closed = 0.0, worst_parallel = 1.0;
        const glm::dvec3 ambient = glm::normalize(field.ambient_nt);
        for (const auto& line : field.lines) {
            if (on_surface(line.front()) && on_surface(line.back())) {
                ++closed;
                for (const glm::dvec4& p : line) farthest_closed = std::max(farthest_closed, glm::length(glm::dvec3(p)));
            } else {
                ++open;
                const bool starts_out = !on_surface(line.front());
                check(on_surface(line.front()) != on_surface(line.back()), "an open line has one foot on the moon");
                const glm::dvec3 end = glm::normalize(starts_out ? glm::dvec3(line[0]) - glm::dvec3(line[1])
                                                                  : glm::dvec3(line.back()) - glm::dvec3(line[line.size() - 2]));
                worst_parallel = std::min(worst_parallel, std::abs(glm::dot(end, ambient)));
            }
        }
        std::printf("  info: Ganymede in %.0f nT: %d closed lines out to %.2f radii, %d open ones, within %.1f deg of Jupiter's field "
                    "at their ends\n", glm::length(field.ambient_nt), closed, farthest_closed, open,
                    std::acos(std::min(worst_parallel, 1.0)) / kDegToRad);
        check(open == 48 && closed >= 16 && closed <= 24, "the polar lines leave, the equatorial ones close", closed);
        check(farthest_closed > 1.5 && farthest_closed < 2.3, "the closed lines stay within about 2 radii", farthest_closed);
        check(worst_parallel > std::cos(12.0 * kDegToRad), "the open lines leave along Jupiter's field", worst_parallel);
        // Jupiter's own line through the moon (the flux tube's) goes the same way.
        std::vector<glm::dvec3> tube;
        check(planet.aurora.jovian->flux_tube(2, false, planet.jovian_moons[2], false, tube), "Jupiter's line through Ganymede");
        if (tube.size() > 1) {
            const glm::dvec3 along = glm::transpose(jupiter_scene.body_to_icrf(ganymede, t)) *
                                     (jupiter_scene.body_to_icrf(planet, t) * glm::normalize(tube[1] - tube[0]));
            check(std::abs(glm::dot(along, ambient)) > std::cos(4.0 * kDegToRad), "runs on from the moon's open lines", glm::dot(along, ambient));
        }
        // The field swings as Jupiter turns: traced again after an hour, not after ten seconds.
        const auto before = field.model;
        follow_field(t + 10.0);
        check(field.model == before, "ten seconds on, the same lines");
        follow_field(t + 3600.0);
        check(field.model != before && glm::length(field.ambient_nt - before->ambient_nt()) > 5.0, "an hour on, Jupiter's field has swung",
              glm::length(field.ambient_nt - before->ambient_nt()));
        // Over a turn of Jupiter the closed shells stay: none blinks out.
        int fewest_closed = 100;
        for (int step = 0; step < 42; ++step) {
            follow_field(t + step * 900.0);
            int count = 0;
            for (const auto& line : field.lines) count += on_surface(line.front()) && on_surface(line.back()) ? 1 : 0;
            fewest_closed = std::min(fewest_closed, count);
        }
        check(fewest_closed == 24, "the three closed shells are closed all the way round Jupiter's turn", fewest_closed);
        // Tracing them again is cheap: it happens every few minutes of the moon's time.
        const auto started = std::chrono::steady_clock::now();
        size_t points = 0;
        for (int i = 0; i < 20; ++i) points += trace_field_lines(*field.model, field.line_settings).size();
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started).count() / 20.0;
        std::printf("  info: tracing Ganymede's %zu lines takes %.2f ms\n", points / 20, ms);
        follow_field(t);
    }
    // The aurora lies at the edge of the closed lines (Duling 2022), which Jupiter's plasma pushes
    // toward the poles upstream (the trailing side) and toward the equator downstream (Marzok
    // 2022). Without plasma the edge is the same on both sides: between the two crescents.
    if (ganymede.aurora.map) {
        const auto crescent = [&](double east_deg, int sign) {
            double best = -1.0, at = 0.0;
            for (double latitude = 5.0; latitude <= 70.0; latitude += 0.5) {
                const double b = ganymede.aurora.map->brightness(east_deg * kDegToRad, sign * latitude * kDegToRad);
                if (b > best) {
                    best = b;
                    at = latitude;
                }
            }
            return at;
        };
        // The latitude where lines stop closing, on a meridian of the moon (body axes are the model's).
        FieldLineSettings settings = ganymede.magnetic_field.line_settings;
        const auto edge = [&](double east_deg, int sign) {
            double low = 5.0, high = 85.0;
            for (int i = 0; i < 30; ++i) {
                const double mid = 0.5 * (low + high);
                const double latitude = sign * mid * kDegToRad;
                const glm::dvec3 seed = 1.002 * glm::dvec3(std::cos(latitude) * std::cos(east_deg * kDegToRad),
                                                           std::cos(latitude) * std::sin(east_deg * kDegToRad), std::sin(latitude));
                const auto line = trace_field_line(*ganymede.magnetic_field.model, seed, settings);
                const bool closes = glm::length(glm::dvec3(line.front())) < 1.01 && glm::length(glm::dvec3(line.back())) < 1.01;
                (closes ? low : high) = mid;
            }
            return low;
        };
        for (int sign : {-1, 1}) {
            const double trailing = crescent(90.0, sign);
            const double leading = crescent(-90.0, sign);
            // Averaged over Jupiter's turn, as the map is over years.
            double trailing_edge = 0.0, leading_edge = 0.0;
            constexpr int kPhases = 21;
            for (int step = 0; step < kPhases; ++step) {
                follow_field(tdb_from_julian_year(2021.5) + step * 1800.0);
                trailing_edge += edge(90.0, sign) / kPhases;
                leading_edge += edge(-90.0, sign) / kPhases;
            }
            std::printf("  info: Ganymede %s: lines stop closing at %.1f deg (trailing), %.1f (leading); the aurora's crescents at %.1f, %.1f\n",
                        sign > 0 ? "north" : "south", trailing_edge, leading_edge, trailing, leading);
            for (double edge_latitude : {trailing_edge, leading_edge}) {
                check(edge_latitude > leading - 2.0 && edge_latitude < trailing + 2.0,
                      "the edge of the closed lines is between the leading and trailing crescents", edge_latitude);
            }
        }
        follow_field(tdb_from_julian_year(2021.5));
    }

    // The scene's lines: within their limits, cheap to draw and to make.
    const struct { const Body* body; double max_radius; } limits[] = {
        {&uranus, 12.0}, {&neptune, 12.0}, {&saturn, 12.0}, {&earth, 8.0}, {&jupiter, 10.0}, {&mercury, 1.7}, {&ganymede, 4.0},
    };
    for (const auto& limit : limits) {
        const BodyMagneticField& field = limit.body->magnetic_field;
        size_t points = 0;
        double farthest = 0.0;
        double fade_at_end = 1.0;
        for (const auto& line : field.lines) {
            points += line.size();
            for (const glm::dvec4& p : line) {
                farthest = std::max(farthest, glm::length(glm::dvec3(p)));
                if (glm::length(glm::dvec3(p)) > limit.max_radius - 1e-6) fade_at_end = std::min(fade_at_end, p.w);
            }
        }
        // Ganymede's also has the lines that leave: three latitudes in each hemisphere.
        check(field.lines.size() == (limit.body == &ganymede ? 72u : 24u), "three shells of eight lines", static_cast<double>(field.lines.size()));
        check(points > 200 && points < 3000, "a line costs little", static_cast<double>(points));
        check(farthest <= limit.max_radius + 1e-6, "lines stop at max_radius", farthest);
        check(fade_at_end > 0.999, "where they have faded out", fade_at_end);
    }
    // Lines end on the oblate surface, not the reference sphere: Jupiter's polar radius is 6% less.
    double shortest = 1e9;
    for (const auto& line : jupiter.magnetic_field.lines) {
        shortest = std::min({shortest, glm::length(glm::dvec3(line.front())), glm::length(glm::dvec3(line.back()))});
    }
    check(shortest < 0.99, "Jupiter's lines reach the flattened surface", shortest);

    // Neptune's model turns at its own rate, 4.8 deg/day slower than the scene's IAU rate.
    const glm::dvec3 prime_a = magnetic_model_to_body(neptune, 0.0)[0];
    const glm::dvec3 prime_b = magnetic_model_to_body(neptune, kSecondsPerDay)[0];
    check(std::abs(std::acos(glm::dot(prime_a, prime_b)) / kDegToRad - 4.8269265) < 1e-5, "Neptune's field keeps NLS' rate");
    check(magnetic_rotation_deg_per_day(neptune) == 536.3128492, "which is how fast it turns");
    check(magnetic_rotation_deg_per_day(uranus) == uranus.pm_rate_deg_per_day, "Uranus's turns with Lamy's rate");

    // Display limits: readable while a turn takes more than 4 s; a planet of a few pixels shows none.
    check(field_line_visibility(500.9, 1800.0) == 1.0, "1,800 x warp shows them");
    check(field_line_visibility(500.9, 3.0e5) == 0.0, "a turn a second does not");
    check(field_line_size_fade(1.0) == 0.0 && field_line_size_fade(10.0) == 1.0, "nor a dot of a planet");

    // Off by default; only these events turn them on, and with something to show.
    Scene uranus_scene = load_scene_or_die("uranus.toml");
    Scene neptune_scene = load_scene_or_die("neptune.toml");
    for (const Scene* scene : {&solar,&uranus_scene, &neptune_scene, &saturn_scene, &earth_scene, &jupiter_scene, &mercury_scene}) {
        for (const SceneEvent& event : scene->events) {
            if (event.magnetic_lines) {
                check(event.focus >= 0 && scene->bodies[static_cast<size_t>(event.focus)].magnetic_field.model != nullptr,
                      "an event's field lines are on its focus");
                check(event.focus_orbits_only, "and no distant orbits cross them");
            }
        }
    }
    const auto shown = [](const Scene& scene) {
        return std::count_if(scene.events.begin(), scene.events.end(), [](const SceneEvent& e) { return e.magnetic_lines; });
    };
    check(shown(uranus_scene) == 1 && shown(neptune_scene) == 1 && shown(saturn_scene) == 1 && shown(mercury_scene) == 1 && shown(jupiter_scene) == 1,
          "the ice giants, Saturn, Mercury and Ganymede each have their event");
    check(shown(solar) == 0, "the Solar System's tour leaves them to the planets' own scenes");

    // A wrong description is refused.
    const auto loads = [](const std::string& field) {
        const std::string text = "name = 'x'\n[[bodies]]\nname = 'Planet'\nradii_km = [25600.0]\n[bodies.magnetic_field]\n" + field;
        Scene out;
        std::string error;
        return load_scene_string(text, "test", out, &error, ASTRAXIS_ASSET_DIR);
    };
    const std::string good = "coefficients = 'geomag/field_uranus_ah5.txt'\nreference_radius_km = 25600.0\n";
    check(loads(good), "a field model loads", 0.0);
    check(!loads("reference_radius_km = 25600.0\n"), "a field needs coefficients");
    check(!loads("coefficients = 'geomag/field_uranus_ah5.txt'\n"), "and a reference radius");
    check(!loads(good + "longitude_offset_deg = 1.0\nprime_meridian_deg = [0.0, 1.0]\n"), "and only one way to place it");
    check(!loads(good + "l_shells = [0.5]\n"), "and shells outside the planet");
}

// Magnetopauses and bow shocks (scene/magnetosphere.hpp; assets/geomag/SOURCES.md).
void test_magnetospheres()
{
    const auto utc = [](const char* text) {
        double t = 0.0;
        check(parse_utc(text, &t), "a time parses");
        return t;
    };
    // A body in another's magnetosphere axes, in the model's radii.
    const auto in_axes = [](Scene& scene, const char* planet, const char* craft, double t, const SolarWind& wind) {
        scene.update(t);
        const int p = scene.find(planet);
        const glm::dvec3 offset = scene.icrf_state_at(scene.find(craft), t).position - scene.icrf_state_at(p, t).position;
        const BodyMagnetosphere& sphere = scene.bodies[static_cast<size_t>(p)].magnetosphere;
        return glm::transpose(scene.magnetosphere_axes(p, t, wind)) * offset / sphere.radius_km - sphere.center_offset_radii;
    };

    // --- The Earth: Shue et al. 1998 -------------------------------------------------------
    Scene earth_scene = load_scene_or_die("earth_moon.toml");
    const BodyMagnetosphere& earth = earth_scene.bodies[static_cast<size_t>(earth_scene.find("Earth"))].magnetosphere;
    check(earth.model == MagnetosphereModel::Shue1998, "the Earth uses Shue 1998");
    const BoundarySurface quiet = magnetopause(earth, earth.typical);
    std::printf("Earth: quiet standoff %.2f Re, flaring %.3f, shock nose %.2f Re\n", quiet.standoff, quiet.flaring,
                bow_shock(earth, earth.typical).nose());
    check(std::abs(quiet.standoff - 10.25) < 0.05, "a 2 nPa wind stands off at 10.25 radii", quiet.standoff);
    check(std::abs(quiet.radius({0.0, 1.0, 0.0}) - quiet.standoff * std::pow(2.0, quiet.flaring)) < 1e-9,
          "the terminator is 2^alpha standoffs out");
    // Slavin et al. 1985: the tail is 60 +/- 5 radii across at 130-225 radii, no longer flaring.
    for (double x : {130.0, 180.0, 225.0}) {
        const double half = quiet.tail_half_width(x);
        check(std::abs(half - 30.0) < 2.5, "the far tail is 60 +/- 5 radii across (Slavin 1985)", half);
    }
    check(std::abs(quiet.tail_half_width(20.0) / quiet.tail_half_width(120.0) - 1.0 / 1.3) < 0.2,
          "it widens by about 30 percent from 20 to 120 radii (Slavin 1983)",
          quiet.tail_half_width(120.0) / quiet.tail_half_width(20.0));
    // OMNI, 2024-05-11 02-03 UT: 48.76 nPa, Bz -20.5 nT. Inside geosynchronous orbit.
    const SceneEvent* storm = nullptr;
    for (const SceneEvent& e : earth_scene.events) {
        if (!e.magnetosphere.empty()) storm = &e;
    }
    check(storm != nullptr, "the Earth-Moon scene has a magnetosphere event");
    if (storm) {
        // The event follows the hourly record of the wind.
        const MagnetosphereState& live = storm->magnetosphere.front();
        check(live.recorded && earth.record != nullptr, "the storm's magnetopause follows the record of the wind");
        if (!live.recorded || !earth.record) return;
        // The middle of an hour is that hour's mean (to the record's steps: 0.01 nPa, 0.5 nT, 8 km/s).
        const SolarWind peak = earth_scene.magnetosphere_wind(live, utc("2024-05-11 02:30"), 1.0);
        check(std::abs(peak.pressure_npa - 48.76) < 0.006 && std::abs(peak.bz_nt + 20.5) < 0.26 && std::abs(peak.speed_km_s - 738.0) < 4.1,
              "OMNI's hour of 2024-05-11 02 UT: 48.76 nPa, Bz -20.5 nT, 738 km/s", peak.pressure_npa);
        // Between two hours' middles, on the line between them (03 UT: 23.68 nPa, -22.8 nT).
        const SolarWind between = earth_scene.magnetosphere_wind(live, utc("2024-05-11 03:00"), 1.0);
        check(std::abs(between.pressure_npa - 0.5 * (48.76 + 23.68)) < 0.006 && std::abs(between.bz_nt + 0.5 * (20.5 + 22.8)) < 0.26,
              "and half way to the next hour's", between.pressure_npa);
        const BoundarySurface squeezed = magnetopause(earth, peak);
        std::printf("Earth: storm standoff %.2f Re, flaring %.3f\n", squeezed.standoff, squeezed.flaring);
        check(squeezed.standoff < 6.6 && squeezed.standoff > 4.0, "the storm's magnetopause is inside geosynchronous orbit",
              squeezed.standoff);
        check(squeezed.flaring > quiet.flaring, "and a southward field flares the tail", squeezed.flaring);
        check(storm->magnetosphere.size() == 2 && storm->magnetosphere.back().faint && !storm->magnetosphere.back().recorded,
              "with the quiet one for comparison");
        // The event opens before the shock arrives (17 UT on the 10th) and sees the edge come in.
        const double opening = magnetopause(earth, earth_scene.magnetosphere_wind(live, storm->t_tdb, storm->warp)).standoff;
        double least = 1e9, least_hours = 0.0, inside_hours = 0.0;
        for (double hours = 0.0; hours <= 36.0; hours += 0.25) {
            const double standoff = magnetopause(earth, earth_scene.magnetosphere_wind(live, storm->t_tdb + 3600.0 * hours, storm->warp)).standoff;
            inside_hours += standoff < 6.6 ? 0.25 : 0.0;
            if (standoff < least) {
                least = standoff;
                least_hours = hours;
            }
        }
        std::printf("Earth: the storm event opens at %.2f Re; least %.2f Re %.2f h in; inside 6.6 Re for %.2f h of 36\n", opening,
                    least, least_hours, inside_hours);
        check(opening > 9.0 && storm->t_tdb < utc("2024-05-10 17:00"), "it opens on a quiet magnetosphere", opening);
        check(least < 5.5 && inside_hours > 8.0, "which is then squeezed inside geosynchronous orbit for hours", least);
        // Hours without a measurement, and times outside the record, are the typical wind's.
        const SolarWindRecord& record = *earth.record;
        const double gap = utc("1983-10-27 04:30");
        check(record.coverage(gap - 1800.0, gap + 1800.0) == 0.0, "OMNI has nothing for 1983-10-27");
        const SolarWind in_gap = earth_scene.magnetosphere_wind(live, gap, 1.0);
        check(in_gap.pressure_npa == earth.typical.pressure_npa && in_gap.bz_nt == earth.typical.bz_nt &&
                  in_gap.speed_km_s == earth.typical.speed_km_s,
              "a gap is filled with the typical wind", in_gap.pressure_npa);
        check(utc("1963-01-01 00:00") == record.start_tdb() && record.end_tdb() > utc("2026-10-01 00:00"), "the record runs from 1963 to 2026");
        for (double t : {record.start_tdb() - 86400.0, record.end_tdb() + 86400.0, utc("2060-01-01 00:00")}) {
            check(earth_scene.magnetosphere_wind(live, t, 1.0).pressure_npa == earth.typical.pressure_npa, "outside the record: the typical wind");
        }
        // Coverage: about a third of the hours of the 1960s, nearly all since 2000.
        const double early = record.coverage(utc("1964-01-01 00:00"), utc("1970-01-01 00:00"));
        const double late = record.coverage(utc("2010-01-01 00:00"), utc("2020-01-01 00:00"));
        std::printf("Earth: the wind is measured in %.1f%% of 1964-69 and %.1f%% of 2010-19\n", 100.0 * early, 100.0 * late);
        check(early > 0.3 && early < 0.6 && late > 0.98, "the record's coverage", late);
        // At a high warp the window is long: the storm is averaged away, a year is near the median.
        const SolarWind blurred = earth_scene.magnetosphere_wind(live, utc("2024-05-11 02:30"), 86400.0 * 20.0);
        check(blurred.pressure_npa > 1.2 && blurred.pressure_npa < 4.0, "a month's window averages the storm away", blurred.pressure_npa);
        // A state of its own wind is not touched, nor is the quiet one.
        const MagnetosphereState& faint = storm->magnetosphere.back();
        check(earth_scene.magnetosphere_wind(faint, utc("2024-05-11 02:30"), 1.0).pressure_npa == earth.typical.pressure_npa,
              "the comparison stays at the typical wind");

        Scene bad;
        std::string error;
        const std::string head = "name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [6378.0]\n[bodies.magnetosphere]\nmodel = 'shue1998'\n";
        const std::string event = "[[events]]\nname = 'e'\ntime_tdb = '2024-05-11 02:00'\nfocus = 'P'\ndistance_km = 1e5\n";
        check(!load_scene_string((head + event + "magnetosphere = [{ body = 'P', recorded = true }]\n").c_str(), "no_record", bad, &error,
                                 ASTRAXIS_ASSET_DIR),
              "recorded needs a record");
        check(!load_scene_string((head + "solar_wind = 'geomag/solar_wind.bin'\n" + event +
                                  "magnetosphere = [{ body = 'P', recorded = true, pressure_npa = 3.0 }]\n").c_str(),
                                 "both_winds", bad, &error, ASTRAXIS_ASSET_DIR),
              "recorded excludes a wind of the event's own");
        check(load_scene_string((head + "solar_wind = 'geomag/solar_wind.bin'\n" + event + "magnetosphere = [{ body = 'P', recorded = true }]\n").c_str(),
                                "recorded", bad, &error, ASTRAXIS_ASSET_DIR),
              "a recorded wind loads");
        check(!load_scene_string((head + "solar_wind = 'geomag/kp.bin'\n").c_str(), "not_a_record", bad, &error, ASTRAXIS_ASSET_DIR),
              "another file is not a solar wind record");
    }
    // The axis leans into the Earth's motion: about 4 degrees at 430 km/s (Fairfield 1971).
    {
        const double t = utc("1983-02-08 00:00");
        earth_scene.update(t);
        const int e = earth_scene.find("Earth");
        const glm::dmat3 axes = earth_scene.magnetosphere_axes(e, t, earth.typical);
        const glm::dvec3 to_sun = glm::normalize(earth_scene.sun_icrf_state_at(t).position - earth_scene.icrf_state_at(e, t).position);
        const glm::dvec3 velocity = earth_scene.icrf_state_at(e, t).velocity - earth_scene.sun_icrf_state_at(t).velocity;
        const double lean = std::acos(glm::dot(axes[0], to_sun)) / kDegToRad;
        check(std::abs(lean - 4.0) < 0.4, "the Earth's tail is aberrated by about 4 deg", lean);
        check(glm::dot(axes[0], velocity) > 0.0, "upstream leans toward the Earth's motion");
        check(glm::dot(axes[1], velocity) < 0.0, "y points to dusk, against the motion");
    }

    // ISEE-3 in the deep tail: where it is in the average tail on its two long passes.
    {
        Scene isee = load_scene_or_die("isee3.toml");
        const BodyMagnetosphere& sphere = isee.bodies[static_cast<size_t>(isee.find("Earth"))].magnetosphere;
        const BoundarySurface pause = magnetopause(sphere, sphere.typical);
        int inside = 0;
        int count = 0;
        double farthest = 0.0;
        double farthest_off_axis = 0.0;
        for (double t = utc("1982-10-15 00:00"); t < utc("1983-12-20 00:00"); t += 6.0 * 3600.0) {
            const glm::dvec3 p = in_axes(isee, "Earth", "ISEE-3", t, sphere.typical);
            if (p.x > -60.0) continue;
            ++count;
            inside += pause.contains(p) ? 1 : 0;
            if (-p.x > farthest) {
                farthest = -p.x;
                farthest_off_axis = std::hypot(p.y, p.z);
            }
        }
        std::printf("ISEE-3 beyond 60 Re downtail: %d of %d six-hour samples inside the average tail; farthest %.1f Re, %.1f Re off the axis\n",
                    inside, count, farthest, farthest_off_axis);
        check(farthest > 215.0 && farthest < 240.0, "ISEE-3 went about 220-237 radii down the tail", farthest);
        check(count > 0 && inside > 0.4 * count && inside < 0.75 * count, "and was inside the average tail about half of that time", count > 0 ? double(inside) / count : 0.0);
        bool shown = false;
        for (const SceneEvent& e : isee.events) {
            for (const MagnetosphereState& s : e.magnetosphere) shown = shown || s.tail_radii >= 236.0;
        }
        check(shown, "an ISEE-3 event draws the tail as far as it went");
    }

    // --- Jupiter: Joy et al. 2002 ----------------------------------------------------------
    Scene jupiter_scene = load_jupiter();
    const BodyMagnetosphere& jupiter = jupiter_scene.bodies[static_cast<size_t>(jupiter_scene.find("Jupiter"))].magnetosphere;
    check(jupiter.model == MagnetosphereModel::Joy2002, "Jupiter uses Joy 2002");
    // Rutala et al. 2025, Table: Joy's standoffs at the 10th, 50th and 90th percentile pressures.
    const double joy_table[3][3] = {{0.02, 101.0, 132.0}, {0.071, 83.0, 108.0}, {0.22, 68.0, 82.0}};
    for (const auto& row : joy_table) {
        SolarWind wind;
        wind.pressure_npa = row[0];
        const double pause = magnetopause(jupiter, wind).nose();
        const double shock = bow_shock(jupiter, wind).nose();
        std::printf("Joy at %.3f nPa: magnetopause %.1f, bow shock %.1f Rj\n", row[0], pause, shock);
        check(std::abs(pause - row[1]) < 0.8, "Joy's magnetopause standoff", pause);
        check(std::abs(shock - row[2]) < 0.6, "Joy's bow shock standoff", shock);
    }
    // Joy's two most probable magnetopauses, 63 and 92 radii out, and the pressures the paper
    // gives for them; the mean bow shock, 84 radii, is at its median pressure (Table 2).
    {
        SolarWind wind;
        wind.pressure_npa = 0.306;
        const double compressed = magnetopause(jupiter, wind).nose();
        wind.pressure_npa = 0.039;
        const double expanded = magnetopause(jupiter, wind).nose();
        wind.pressure_npa = 0.258;
        const double shock = bow_shock(jupiter, wind).nose();
        std::printf("Joy: compressed %.1f, expanded %.1f Rj; median bow shock %.1f Rj\n", compressed, expanded, shock);
        check(std::abs(compressed - 63.0) < 1.0, "Joy's compressed magnetopause is 63 Rj out", compressed);
        check(std::abs(expanded - 92.0) < 1.0, "the expanded one 92 Rj", expanded);
        check(std::abs(shock - 84.0) < 8.0, "the bow shock about 84 Rj at its median pressure", shock);
    }
    {
        const BoundarySurface pause = magnetopause(jupiter, jupiter.typical);
        const double dusk = pause.radius({0.0, 1.0, 0.0});
        const double dawn = pause.radius({0.0, -1.0, 0.0});
        const double pole = pause.radius({0.0, 0.0, 1.0});
        std::printf("Joy typical: nose %.1f, dusk %.1f, dawn %.1f, pole %.1f Rj\n", pause.nose(), dusk, dawn, pole);
        check(pole < 0.95 * 0.5 * (dusk + dawn), "Jupiter's magnetopause is flattened at the poles", pole / dusk);
        // The surface satisfies Joy's equation where radius() puts it.
        const glm::dvec3 d = glm::normalize(glm::dvec3(-0.3, 0.6, 0.5));
        const glm::dvec3 p = pause.radius(d) * d / pause.scale;
        const double* c = pause.coefficients;
        check(std::abs(p.z * p.z - (c[0] + c[1] * p.x + c[2] * p.x * p.x + c[3] * p.y + c[4] * p.y * p.y + c[5] * p.x * p.y)) < 1e-12,
              "radius() solves Joy's quadric");
    }
    // Kurth et al. 2002: on 2001-01-10 Cassini left the magnetosphere at 20:35 and Galileo at
    // 20:52 (at the spacecraft), over 100 Rj apart in x. Galileo's crossing is on Joy's "25
    // percent" surface, Cassini's beyond the "10 percent" one: no one surface holds both.
    {
        const auto pressure_through = [&](const glm::dvec3& p) {
            double low = 1e-4;  // large
            double high = 2.0;  // small
            for (int i = 0; i < 60; ++i) {
                SolarWind wind;
                wind.pressure_npa = std::sqrt(low * high);
                (magnetopause(jupiter, wind).contains(p) ? low : high) = wind.pressure_npa;
            }
            return std::sqrt(low * high);
        };
        const glm::dvec3 cassini = in_axes(jupiter_scene, "Jupiter", "Cassini", utc("2001-01-10 20:35"), jupiter.typical);
        const glm::dvec3 galileo = in_axes(jupiter_scene, "Jupiter", "Galileo", utc("2001-01-10 20:52"), jupiter.typical);
        const double p_cassini = pressure_through(cassini);
        const double p_galileo = pressure_through(galileo);
        SolarWind wind;
        wind.pressure_npa = p_cassini;
        const double nose_cassini = magnetopause(jupiter, wind).nose();
        wind.pressure_npa = p_galileo;
        const double nose_galileo = magnetopause(jupiter, wind).nose();
        std::printf("Cassini at (%.1f, %.1f, %.1f) Rj: Joy pressure %.4f nPa, nose %.1f Rj\n", cassini.x, cassini.y, cassini.z,
                    p_cassini, nose_cassini);
        std::printf("Galileo at (%.1f, %.1f, %.1f) Rj: Joy pressure %.4f nPa, nose %.1f Rj\n", galileo.x, galileo.y, galileo.z,
                    p_galileo, nose_galileo);
        check(galileo.x - cassini.x > 100.0, "Galileo was over 100 Rj sunward of Cassini (Kurth 2002)", galileo.x - cassini.x);
        check(cassini.y > 0.0 && galileo.y > 0.0, "both on the dusk side");
        check(p_cassini < p_galileo, "Cassini's crossing needs a larger magnetosphere than Galileo's", p_cassini / p_galileo);
        // Joy's Table 2: 0.048 nPa for the 25 percent surface, 0.030 nPa for the 10 percent one.
        check(std::abs(p_galileo / 0.048 - 1.0) < 0.15, "Galileo's crossing is on Joy's 25 percent surface (Kurth 2002)", p_galileo);
        check(p_cassini < 0.030, "Cassini's is beyond the 10 percent one", p_cassini);
        // z of the axes is toward Jupiter's spin axis.
        {
            const double t = utc("2001-01-10 20:52");
            const int j = jupiter_scene.find("Jupiter");
            const glm::dmat3 axes = jupiter_scene.magnetosphere_axes(j, t, jupiter.typical);
            const glm::dvec3 spin = jupiter_scene.body_to_icrf(jupiter_scene.bodies[static_cast<size_t>(j)], t)[2];
            check(std::abs(glm::dot(axes[1], spin)) < 1e-12 && glm::dot(axes[2], spin) > 0.99, "Joy's z is toward the spin axis",
                  glm::dot(axes[2], spin));
        }
        const SceneEvent* event = nullptr;
        for (const SceneEvent& e : jupiter_scene.events) {
            if (e.magnetosphere.size() == 2) event = &e;
        }
        check(event != nullptr, "the Jupiter scene shows the two magnetopauses");
        if (event) {
            const double p_now = event->magnetosphere[0].wind.pressure_npa;
            const double p_before = event->magnetosphere[1].wind.pressure_npa;
            check(std::abs(p_now / p_galileo - 1.0) < 0.03, "one passes through Galileo's crossing", p_now);
            check(std::abs(p_before / p_cassini - 1.0) < 0.03, "the other, faint, through Cassini's", p_before);
            check(event->magnetosphere[1].faint && !event->magnetosphere[0].faint, "the earlier one is the faint one");
        }
    }
    // Cassini first met the bow shock on 2000-12-28 at 04:19, close to Joy's surface for
    // 0.02 nPa; the wind's pressure from the plasma frequency then was 0.018 nPa (Kurth 2002).
    {
        const glm::dvec3 cassini = in_axes(jupiter_scene, "Jupiter", "Cassini", utc("2000-12-28 04:19"), jupiter.typical);
        SolarWind wind;
        wind.pressure_npa = 0.02;
        const BoundarySurface shock = bow_shock(jupiter, wind);
        const double miss = glm::length(cassini) / shock.radius(glm::normalize(cassini));
        std::printf("Cassini's first bow shock: %.1f Rj out at (%.1f, %.1f, %.1f), %.3f of Joy's 0.02 nPa shock\n",
                    glm::length(cassini), cassini.x, cassini.y, cassini.z, miss);
        check(std::abs(glm::length(cassini) - 140.2) < 1.5, "Cassini met the shock 140.2 Rj out (Kurth 2002, Fig. 2)", glm::length(cassini));
        check(std::abs(miss - 1.0) < 0.1, "close to Joy's bow shock for 0.02 nPa", miss);
    }

    // --- Uranus: Jasinski et al. 2024 ------------------------------------------------------
    Scene uranus_scene = load_scene_or_die("uranus.toml");
    const BodyMagnetosphere& uranus = uranus_scene.bodies[static_cast<size_t>(uranus_scene.find("Uranus"))].magnetosphere;
    check(uranus.model == MagnetosphereModel::PressureBalance, "Uranus uses a pressure balance");
    check(std::abs(magnetopause(uranus, uranus.typical).standoff - 22.2) < 1e-6, "its mean standoff is 22.2 radii");
    {
        // Their pressures and standoffs: 0.0184 nPa at the flyby, 0.031 and 0.00043 nPa the extremes.
        const double table[3][2] = {{0.0184, 17.3}, {0.031, 15.8}, {0.00043, 33.5}};
        for (const auto& row : table) {
            SolarWind wind;
            wind.pressure_npa = row[0];
            const double standoff = magnetopause(uranus, wind).standoff;
            check(std::abs(standoff - row[1]) < 0.1, "Jasinski's standoff for a pressure", standoff);
            check(std::abs(uranus.pressure_for_standoff(row[1]) / row[0] - 1.0) < 0.04, "and back");
        }
        // Voyager 2 crossed the magnetopause at 10:07 at (15.784, 8.684, -0.632) radii, 29 deg
        // from the nose, which their shape puts 17.3 radii out.
        SolarWind flyby;
        flyby.pressure_npa = uranus.pressure_for_standoff(17.3);
        const glm::dvec3 voyager = in_axes(uranus_scene, "Uranus", "Voyager 2", utc("1986-01-24 10:07"), flyby);
        const BoundarySurface pause = magnetopause(uranus, flyby);
        const double ratio = glm::length(voyager) / pause.radius(glm::normalize(voyager));
        std::printf("Voyager 2 at the magnetopause: (%.3f, %.3f, %.3f) Ru, %.4f of the model's distance\n", voyager.x, voyager.y,
                    voyager.z, ratio);
        check(glm::length(voyager - glm::dvec3(15.784, 8.684, -0.632)) < 0.15, "Voyager 2 where Jasinski has it at 10:07",
              glm::length(voyager - glm::dvec3(15.784, 8.684, -0.632)));
        check(std::abs(ratio - 1.0) < 0.01, "on the 17.3-radii magnetopause", ratio);
        check(magnetopause(uranus, uranus.typical).contains(voyager), "and well inside the mean one");
        // The bow shock at 07:28, 23.5 radii out (Bagenal et al. 1987 as quoted by Treumann &
        // Jaroschek 2008): the assumed parabola 1.275 standoffs out is within a radius of it.
        const glm::dvec3 at_shock = in_axes(uranus_scene, "Uranus", "Voyager 2", utc("1986-01-24 07:28"), flyby);
        const double shock = bow_shock(uranus, flyby).radius(glm::normalize(at_shock));
        std::printf("Voyager 2 at the bow shock: %.2f Ru out, the model's shock %.2f Ru\n", glm::length(at_shock), shock);
        check(std::abs(glm::length(at_shock) - 23.5) < 0.3, "Voyager 2 met the shock 23.5 radii out", glm::length(at_shock));
        check(std::abs(shock - glm::length(at_shock)) < 1.0, "where the model's shock is", shock);
        bool shown = false;
        for (const SceneEvent& e : uranus_scene.events) {
            if (e.magnetosphere.size() == 2) {
                shown = std::abs(magnetopause(uranus, e.magnetosphere[0].wind).standoff - 17.3) < 1e-6 &&
                        std::abs(magnetopause(uranus, e.magnetosphere[1].wind).standoff - 22.2) < 1e-6 && e.magnetosphere[1].faint;
            }
        }
        check(shown, "the Uranus scene compares the squeezed magnetopause with the mean one");
    }

    // --- Saturn: Kanani et al. 2010 and Went et al. 2011 -----------------------------------
    Scene saturn_scene = load_scene_or_die("saturn.toml");
    const BodyMagnetosphere& saturn = saturn_scene.bodies[static_cast<size_t>(saturn_scene.find("Saturn"))].magnetosphere;
    {
        check(saturn.model == MagnetosphereModel::PowerLaw, "Saturn uses a power law");
        SolarWind wind;
        wind.pressure_npa = 0.04;
        const BoundarySurface shock = bow_shock(saturn, wind);
        const BoundarySurface pause = magnetopause(saturn, wind);
        std::printf("Saturn at 0.04 nPa: magnetopause %.2f Rs (flaring %.3f), bow shock %.2f Rs; typical %.2f and %.2f Rs\n",
                    pause.nose(), pause.flaring, shock.nose(), magnetopause(saturn, saturn.typical).nose(),
                    bow_shock(saturn, saturn.typical).nose());
        // A subsolar shock 27 radii out for about 0.04 nPa (Sulaiman et al. 2016 on Went's model).
        check(std::abs(shock.nose() - 27.0) < 0.5, "Went's bow shock is 27 Rs out at 0.04 nPa", shock.nose());
        check(std::abs(shock.radius({0.0, 1.0, 0.0}) / shock.nose() - 1.84) < 1e-9, "with eccentricity 0.84");
        check(std::abs(pause.nose() - 10.3 * std::pow(0.04, -0.2)) < 1e-9, "Kanani's standoff is 10.3 P^-0.2");
        check(std::abs(pause.flaring - (0.73 + 0.4 * 0.04)) < 1e-12, "and the flaring 0.73 + 0.4 P");
        // The event: Titan, 20.3 radii out, inside the typical magnetopause and outside the
        // one for 2.5 times the pressure when it is on the day side.
        const SceneEvent* event = nullptr;
        for (const SceneEvent& e : saturn_scene.events) {
            if (e.magnetosphere.size() == 2) event = &e;
        }
        check(event != nullptr, "the Saturn scene has a magnetosphere event");
        if (event) {
            const glm::dvec3 titan = in_axes(saturn_scene, "Saturn", "Titan", event->t_tdb, saturn.typical);
            const double from_noon = std::acos(titan.x / glm::length(titan)) / kDegToRad;
            std::printf("Titan at the event: %.2f Rs out, %.1f deg from the Sun's direction\n", glm::length(titan), from_noon);
            check(std::abs(glm::length(titan) - 20.3) < 0.7, "Titan is about 20.3 Rs out", glm::length(titan));
            check(from_noon < 40.0, "and on the day side at the event", from_noon);
            check(magnetopause(saturn, event->magnetosphere[0].wind).contains(titan), "inside the typical magnetopause");
            check(!magnetopause(saturn, event->magnetosphere[1].wind).contains(titan) && event->magnetosphere[1].faint,
                  "outside the one for 2.5 times the pressure");
            check(bow_shock(saturn, event->magnetosphere[1].wind).contains(titan), "but still behind the shock");
        }
    }

    // --- Mercury: Winslow et al. 2013 ------------------------------------------------------
    Scene mercury_scene = load_scene_or_die("mercury.toml");
    const BodyMagnetosphere& mercury = mercury_scene.bodies[static_cast<size_t>(mercury_scene.find("Mercury"))].magnetosphere;
    {
        check(mercury.model == MagnetosphereModel::PowerLaw, "Mercury uses a power law");
        // The abstract's standoffs: 1.45 radii at 14.3 nPa, 1.55 to 1.35 for 8.8 to 21.6 nPa.
        const double table[3][2] = {{14.3, 1.45}, {8.8, 1.55}, {21.6, 1.35}};
        for (const auto& row : table) {
            SolarWind wind;
            wind.pressure_npa = row[0];
            check(std::abs(magnetopause(mercury, wind).nose() - row[1]) < 0.02, "Winslow's magnetopause standoff", magnetopause(mercury, wind).nose());
        }
        const BoundarySurface pause = magnetopause(mercury, mercury.typical);
        const BoundarySurface shock = bow_shock(mercury, mercury.typical);
        check(std::abs(shock.nose() - 1.96) < 0.005, "the bow shock is 1.96 radii out", shock.nose());
        // Winslow's equation 3 with the fit corrected for the Mach number (Table 2): focus 0.5
        // radii up the axis, p = 2.90 radii, e = 1.02, so |q| (1 + e cos) = p e about the focus.
        {
            const glm::dvec3 d = glm::normalize(glm::dvec3(-0.4, 0.7, 0.3));
            const glm::dvec3 q = shock.radius(d) * d - glm::dvec3(0.5, 0.0, 0.0);
            check(std::abs(glm::length(q) + 1.02 * q.x - 2.90 * 1.02) < 0.005, "the shock is Winslow's conic about its focus",
                  glm::length(q) + 1.02 * q.x);
            check(shock.radius({0.0, 1.0, 0.0}) > pause.radius({0.0, 1.0, 0.0}), "outside the magnetopause abeam");
        }
        check(std::abs(mercury.center_offset_radii.z - 0.196) < 1e-9, "centered 0.196 radii north (Winslow's Zd)");
        // About 7 deg of aberration on average, 3.5 to 10.2 over the orbit and the wind's speeds.
        {
            const int m = mercury_scene.find("Mercury");
            double low = 90.0;
            double high = 0.0;
            for (int day = 0; day < 88; day += 4) {
                const double t = utc("2013-01-01 00:00") + 86400.0 * day;
                mercury_scene.update(t);
                const glm::dvec3 to_sun = glm::normalize(mercury_scene.sun_icrf_state_at(t).position - mercury_scene.icrf_state_at(m, t).position);
                const double lean = std::acos(std::min(1.0, glm::dot(mercury_scene.magnetosphere_axes(m, t, mercury.typical)[0], to_sun))) / kDegToRad;
                low = std::min(low, lean);
                high = std::max(high, lean);
            }
            std::printf("Mercury's aberration over an orbit: %.1f to %.1f deg\n", low, high);
            check(low > 5.0 && high < 9.0 && 0.5 * (low + high) > 6.0 && 0.5 * (low + high) < 8.0, "Mercury's aberration is about 7 deg", 0.5 * (low + high));
        }
        check(pause.nose() - 1.0 + mercury.center_offset_radii.z < 0.7, "the magnetopause is under 0.7 radii above the ground at noon");
        // MESSENGER over the event's day: inside the magnetopause at its lowest, beyond the shock at its highest.
        const SceneEvent* event = nullptr;
        for (const SceneEvent& e : mercury_scene.events) {
            if (!e.magnetosphere.empty()) event = &e;
        }
        check(event != nullptr, "the Mercury scene has a magnetosphere event");
        if (event) {
            int pause_crossings = 0;
            int shock_crossings = 0;
            bool was_in_pause = false;
            bool was_in_shock = false;
            for (int i = 0; i <= 24 * 60; ++i) {
                const glm::dvec3 p = in_axes(mercury_scene, "Mercury", "MESSENGER", event->t_tdb + 60.0 * i, mercury.typical);
                const bool in_pause = pause.contains(p);
                const bool in_shock = shock.contains(p);
                pause_crossings += i > 0 && in_pause != was_in_pause ? 1 : 0;
                shock_crossings += i > 0 && in_shock != was_in_shock ? 1 : 0;
                was_in_pause = in_pause;
                was_in_shock = in_shock;
            }
            std::printf("MESSENGER over a day: %d magnetopause and %d bow shock crossings of the average surfaces\n",
                        pause_crossings, shock_crossings);
            check(pause_crossings >= 4 && shock_crossings >= 4, "MESSENGER crosses both surfaces on every orbit", shock_crossings);
        }
    }

    // --- Neptune: Voyager 2's crossings ----------------------------------------------------
    Scene neptune_scene = load_scene_or_die("neptune.toml");
    const BodyMagnetosphere& neptune = neptune_scene.bodies[static_cast<size_t>(neptune_scene.find("Neptune"))].magnetosphere;
    {
        check(neptune.model == MagnetosphereModel::PowerLaw, "Neptune's surfaces are of a given size");
        const BoundarySurface pause = magnetopause(neptune, neptune.typical);
        const BoundarySurface shock = bow_shock(neptune, neptune.typical);
        check(std::abs(pause.nose() - 26.5) < 1e-9 && std::abs(shock.nose() - 34.9) < 0.01, "26.5 and 34.9 radii out", shock.nose());
        SolarWind gust;
        gust.pressure_npa = 10.0 * neptune.typical.pressure_npa;
        check(magnetopause(neptune, gust).nose() == pause.nose(), "whatever the pressure");
        // Bow shock 14:38 UT, magnetopause 18:00 UT (Lepping et al. 1992), 34.9 and 26.5 radii
        // out, near the line to the Sun (Ness et al. 1989).
        const glm::dvec3 at_shock = in_axes(neptune_scene, "Neptune", "Voyager 2", utc("1989-08-24 14:38"), neptune.typical);
        const glm::dvec3 at_pause = in_axes(neptune_scene, "Neptune", "Voyager 2", utc("1989-08-24 18:00"), neptune.typical);
        std::printf("Voyager 2 at Neptune: shock %.2f Rn out, %.1f deg from the Sun's line (model %.2f); magnetopause %.2f Rn, %.1f deg (model %.2f)\n",
                    glm::length(at_shock), std::acos(at_shock.x / glm::length(at_shock)) / kDegToRad, shock.radius(glm::normalize(at_shock)),
                    glm::length(at_pause), std::acos(at_pause.x / glm::length(at_pause)) / kDegToRad, pause.radius(glm::normalize(at_pause)));
        check(std::abs(glm::length(at_shock) - 34.9) < 0.4, "Voyager 2 met Neptune's shock 34.9 radii out", glm::length(at_shock));
        check(std::abs(glm::length(at_pause) - 26.5) < 0.4, "and the magnetopause 26.5 radii out", glm::length(at_pause));
        check(std::acos(at_pause.x / glm::length(at_pause)) / kDegToRad < 30.0, "near the line to the Sun");
        check(std::abs(shock.radius(glm::normalize(at_shock)) / glm::length(at_shock) - 1.0) < 0.05, "the model's shock is there", shock.radius(glm::normalize(at_shock)));
        check(std::abs(pause.radius(glm::normalize(at_pause)) / glm::length(at_pause) - 1.0) < 0.05, "and its magnetopause", pause.radius(glm::normalize(at_pause)));
        bool shown = false;
        for (const SceneEvent& e : neptune_scene.events) shown = shown || !e.magnetosphere.empty();
        check(shown, "the Neptune scene has a magnetosphere event");

        // The shock follows the planet's turn (Mejnertsen et al. 2016, Figure 4, read by
        // tools/geomag/neptune_boundaries_fig4.py: on the axis 32.75 radii pole-on and 31.25
        // Earth-like in both planes): its nose is 1.5 radii farther with the dipole toward the Sun
        // than half a turn later, and 34.9 radii out at the phase Voyager 2 met it. The
        // magnetopause stays.
        const int n = neptune_scene.find("Neptune");
        const Body& planet = neptune_scene.bodies[static_cast<size_t>(n)];
        MagnetosphereState state;
        state.body = n;
        state.wind = neptune.typical;
        const auto nose_at = [&](double t, double warp) {
            neptune_scene.update(t);
            const SolarWind wind = neptune_scene.magnetosphere_wind(state, t, warp);
            check(magnetopause(neptune, wind).nose() == pause.nose(), "the magnetopause does not move");
            return bow_shock(neptune, wind).nose();
        };
        const auto alignment = [&](double nose) { return (nose - shock.nose()) / neptune.shock_dipole_radii + neptune.shock_dipole_reference; };
        const double t_shock = utc("1989-08-24 14:38");
        const double t_pause = utc("1989-08-24 18:00");
        const double t_closest = utc("1989-08-25 03:56");
        const double turn = 360.0 / magnetic_rotation_deg_per_day(planet) * 86400.0;
        double low = 1e9, high = 0.0, mean = 0.0;
        int high_phase = 0;
        constexpr int kPhases = 480;
        for (int k = 0; k < kPhases; ++k) {
            const double nose = nose_at(t_shock + turn * k / kPhases, 1.0);
            low = std::min(low, nose);
            mean += nose / kPhases;
            if (nose > high) {
                high = nose;
                high_phase = k;
            }
        }
        // Earth-like: half a turn from pole-on, the dipole then nearest the orbit's pole.
        const double earth_like = nose_at(t_shock + turn * (high_phase + kPhases / 2) / kPhases, 1.0);
        std::printf("Neptune's shock over a turn (%.2f h): %.2f .. %.2f Rn, mean %.2f, %.2f half a turn from pole-on; cos^2 of the "
                    "dipole to the Sun %.3f .. %.3f, %.3f half a turn from pole-on; %.3f at the shock crossing, %.3f at the "
                    "magnetopause's, %.3f at closest approach\n",
                    turn / 3600.0, low, high, mean, earth_like, alignment(low), alignment(high), alignment(earth_like),
                    alignment(nose_at(t_shock, 1.0)), alignment(nose_at(t_pause, 1.0)), alignment(nose_at(t_closest, 1.0)));
        check(std::abs(nose_at(t_shock, 1.0) - 34.9) < 0.02, "the shock is 34.9 radii out when Voyager 2 met it", nose_at(t_shock, 1.0));
        check(std::abs(high - earth_like - (32.75 - 31.25)) < 0.03, "pole-on, its nose is 1.5 radii farther than Earth-like", high - earth_like);
        // Voyager 2 met the shock with the dipole across the Sun's line, where it is nearest. (Before
        // that Cairns et al. 1991 saw its nose move between 38.5 and 34.2 radii: more than this.)
        check(nose_at(t_shock, 1.0) < low + 0.05 && high > 36.5 && high < 37.0, "from 34.9 radii out to 36.7", high);
        // Voyager 2 entered through the cusp, near pole-on, and passed closest in the Earth-like
        // configuration (Mejnertsen et al. 2016, section 3.5).
        check(alignment(nose_at(t_pause, 1.0)) > 0.6, "Voyager 2 crossed the magnetopause near pole-on", alignment(nose_at(t_pause, 1.0)));
        check(alignment(nose_at(t_closest, 1.0)) < 0.25, "and passed closest in the Earth-like configuration",
              alignment(nose_at(t_closest, 1.0)));
        // Turning too fast to follow: the mean over a turn, from the tilt and the season.
        check(std::abs(nose_at(t_shock, 3.0e5) - mean) < 0.02, "at a high warp the shock stays at its mean", nose_at(t_shock, 3.0e5));
        check(std::abs(mean_dipole_alignment(0.0, 0.3) - std::pow(std::cos(0.3), 2.0)) < 1e-12 &&
                  std::abs(mean_dipole_alignment(0.5 * kPi, 0.5 * kPi) - 0.5) < 1e-12,
              "the mean alignment of an aligned dipole and of one turning in the plane of the wind");
        // A shock given by an event is left alone.
        MagnetosphereState given = state;
        given.wind.shock_radii = 38.5;
        check(bow_shock(neptune, neptune_scene.magnetosphere_wind(given, t_pause, 1.0)).nose() == 38.5, "a given shock stays where it is given");
    }

    // --- The Sun: the heliopause and the termination shock ---------------------------------
    {
        Scene solar = load_scene_or_die("solar_system.toml");
        const int sun = solar.find("Sun");
        const BodyMagnetosphere& helio = solar.bodies[static_cast<size_t>(sun)].magnetosphere;
        check(helio.model == MagnetosphereModel::Heliosphere, "the Sun has a heliosphere");
        const double au = helio.radius_km;
        check(std::abs(au - kAuKm) < 1.0, "in au");
        // The nose: ecliptic longitude 255.7, latitude 5.1 (McComas et al. 2015), i.e.
        // opposite to where the interstellar helium flows, 75.7 and -5.1.
        const glm::dmat3 axes = solar.magnetosphere_axes(sun, 0.0, helio.typical);
        const glm::dvec3 pole = ecliptic_pole_icrf();
        const double nose_lat = std::asin(glm::dot(axes[0], pole)) / kDegToRad;
        const glm::dvec3 in_plane = axes[0] - glm::dot(axes[0], pole) * pole;
        const double nose_lon = wrap_two_pi(std::atan2(glm::dot(in_plane, glm::cross(pole, glm::dvec3(1.0, 0.0, 0.0))), in_plane.x)) / kDegToRad;
        check(std::abs(nose_lat - 5.1) < 1e-6 && std::abs(nose_lon - 255.7) < 1e-6, "the nose is toward 255.7, +5.1 deg", nose_lon);
        check(glm::dot(axes[2], pole) > 0.99, "z is toward the ecliptic's north");
        // Parker's shape: r = R / cos(theta / 2); twice the nose's distance wide far down the tail.
        const BoundarySurface typical = magnetopause(helio, helio.typical);
        check(std::abs(typical.radius({0.0, 1.0, 0.0}) - 120.0 * std::sqrt(2.0)) < 1e-9, "Parker's heliopause is sqrt 2 noses out abeam");
        check(std::abs(typical.tail_half_width(1e5) / 240.0 - 1.0) < 0.01, "and two noses wide far downstream", typical.tail_half_width(1e5));
        check(bow_shock(helio, helio.typical).nose() < typical.nose(), "the termination shock is inside the heliopause");

        // Where the Voyagers were at their crossings, against the published distances:
        // termination shock 94.0 au (V1, 2004-12-16; Stone et al. 2005) and 83.7 au (V2,
        // 2007-08-30; Stone et al. 2008); heliopause 121.6 au (V1, 2012-08-25; Stone et al.
        // 2013) and 119.0 au (V2, 2018-11-05; Stone et al. 2019).
        struct Crossing {
            const char* craft;
            const char* shock_date;
            double shock_au;
            const char* pause_date;
            double pause_au;
        };
        const Crossing crossings[2] = {{"Voyager 1", "2004-12-16 00:00", 94.0, "2012-08-25 00:00", 121.6},
                                       {"Voyager 2", "2007-08-30 00:00", 83.7, "2018-11-05 00:00", 119.0}};
        double noses[2] = {0.0, 0.0};
        for (int i = 0; i < 2; ++i) {
            const Crossing& c = crossings[i];
            const glm::dvec3 at_shock = in_axes(solar, "Sun", c.craft, utc(c.shock_date), helio.typical);
            const glm::dvec3 at_pause = in_axes(solar, "Sun", c.craft, utc(c.pause_date), helio.typical);
            const double theta = std::acos(at_pause.x / glm::length(at_pause));
            noses[i] = glm::length(at_pause) * std::cos(0.5 * theta);
            std::printf("%s: termination shock %.2f au, heliopause %.2f au, %.1f deg from the nose (%.1f deg north); Parker nose %.2f au\n",
                        c.craft, glm::length(at_shock), glm::length(at_pause), theta / kDegToRad,
                        std::asin(at_pause.z / glm::length(at_pause)) / kDegToRad, noses[i]);
            check(std::abs(glm::length(at_shock) - c.shock_au) < 0.3, "a Voyager at its termination shock crossing", glm::length(at_shock));
            check(std::abs(glm::length(at_pause) - c.pause_au) < 0.3, "a Voyager at its heliopause crossing", glm::length(at_pause));
            // The event draws the heliopause through the craft and the shock where it met it.
            const SceneEvent* event = nullptr;
            for (const SceneEvent& e : solar.events) {
                if (!e.magnetosphere.empty() && e.name.find(c.craft) != std::string::npos) event = &e;
            }
            check(event != nullptr, "a Voyager's heliopause event draws the heliosphere");
            if (!event) continue;
            const glm::dvec3 then = in_axes(solar, "Sun", c.craft, event->t_tdb, helio.typical);
            const SolarWind& wind = event->magnetosphere.front().wind;
            const double miss = glm::length(then) / magnetopause(helio, wind).radius(glm::normalize(then));
            check(std::abs(miss - 1.0) < 0.003, "its heliopause passes through the craft", miss);
            check(std::abs(bow_shock(helio, wind).nose() - c.shock_au) < 1e-9, "its termination shock is where the craft met it");
            check(std::abs(event->t_tdb - utc(c.pause_date)) < 86400.0, "it opens on the day of the crossing");
        }
        // One Parker surface cannot hold both: Voyager 2, farther from the nose, was nearer.
        std::printf("Voyager 2's crossing is at %.3f of the surface through Voyager 1's\n", noses[1] / noses[0]);
        check(noses[1] < noses[0], "Voyager 2's heliopause is the smaller of the two", noses[1] / noses[0]);
    }

    // --- The wire surfaces -----------------------------------------------------------------
    for (const BodyMagnetosphere* sphere : {&earth, &jupiter, &uranus, &saturn, &mercury, &neptune}) {
        for (int which = 0; which < 2; ++which) {
            const BoundarySurface surface = which == 0 ? magnetopause(*sphere, sphere->typical) : bow_shock(*sphere, sphere->typical);
            std::vector<std::vector<glm::dvec4>> lines;
            boundary_wireframe(surface, sphere->tail_radii, 12, lines);
            check(lines.size() > 14, "a wire surface has meridians and rings", static_cast<double>(lines.size()));
            double worst = 0.0;
            double longest = 0.0;
            double deepest = 0.0;
            size_t points = 0;
            for (const auto& line : lines) {
                points += line.size();
                for (size_t i = 0; i < line.size(); ++i) {
                    const glm::dvec3 p(line[i]);
                    worst = std::max(worst, std::abs(glm::length(p) / surface.radius(glm::normalize(p)) - 1.0));
                    deepest = std::max(deepest, -p.x);
                    check(line[i].w >= 0.0 && line[i].w <= 1.0, "fades are within 0..1");
                    if (i > 0) longest = std::max(longest, glm::length(p - glm::dvec3(line[i - 1])));
                }
            }
            check(worst < 1e-6, "its points are on the surface", worst);
            check(std::abs(deepest / sphere->tail_radii - 1.0) < 1e-3, "it reaches the cut and no farther", deepest);
            check(longest < 0.2 * (2.0 * surface.nose() + sphere->tail_radii), "without long chords", longest);
            check(points < 4000, "and not too many points", static_cast<double>(points));
        }
    }
}

void test_flux_tubes()
{
    Scene scene = load_scene_or_die("jupiter.toml");
    const int index = scene.find("Jupiter");
    const Body& jupiter = scene.bodies[static_cast<size_t>(index)];
    check(jupiter.aurora.jovian != nullptr, "Jupiter has the moons' lines");
    if (!jupiter.aurora.jovian) return;
    const JovianAurora& model = *jupiter.aurora.jovian;

    // The table: from the nominal place of the moon, the time spent rising from 0 to 1, and
    // ending at the footprint table's foot of the same field.
    const double nominal[] = {421800.0, 671100.0, 1070400.0};
    for (int moon = 0; moon < 3; ++moon) {
        const auto& tube = model.flux_tubes[static_cast<size_t>(moon)];
        check(tube.size() == 2u * JovianAurora::kTubeLongitudes * JovianAurora::kTubePoints, "a line for each longitude and hemisphere", static_cast<double>(tube.size()));
        double worst_foot = 0.0, worst_start = 0.0;
        bool rising = true;
        for (int longitude = 0; longitude < JovianAurora::kTubeLongitudes; ++longitude) {
            const double east = longitude * kTwoPi / JovianAurora::kTubeLongitudes;
            const glm::dvec3 seed(nominal[moon] * std::cos(east), nominal[moon] * std::sin(east), 0.0);
            for (int south = 0; south < 2; ++south) {
                const glm::vec4* points = &tube[(static_cast<size_t>(longitude) * 2 + static_cast<size_t>(south)) * JovianAurora::kTubePoints];
                worst_start = std::max(worst_start, glm::length(glm::dvec3(points[0]) - glm::dvec3(nominal[moon] / 71492.0, 0.0, 0.0)) * 71492.0);
                rising = rising && points[0].w == 0.0f && points[JovianAurora::kTubePoints - 1].w == 1.0f;
                for (int k = 1; k < JovianAurora::kTubePoints; ++k) {
                    rising = rising && points[k].w >= points[k - 1].w;
                }
                const glm::dvec3 foot_turned(points[JovianAurora::kTubePoints - 1]);
                const glm::dvec3 foot(std::cos(east) * foot_turned.x - std::sin(east) * foot_turned.y,
                                      std::sin(east) * foot_turned.x + std::cos(east) * foot_turned.y, foot_turned.z);
                glm::dvec3 direction;
                check(model.instant_footprint(moon, south != 0, seed, direction), "the nominal moon is in the footprints' table");
                worst_foot = std::max(worst_foot, glm::length(glm::normalize(foot) - direction) * 72000.0);
            }
        }
        std::printf("  info: %s: lines' feet within %.1f km of the footprint table\n", JovianAurora::kMoons[static_cast<size_t>(moon)], worst_foot);
        check(worst_start < 1.0, "the lines start at the moon's nominal place", worst_start);
        check(rising, "the wave's time rises from nothing to all of it");
        check(worst_foot < 5.0, "and they end where the footprint table says", worst_foot);
    }

    // The lines of the real moons, over a few turns of Jupiter.
    const SceneEvent* wide = nullptr;
    for (const SceneEvent& event : scene.events) {
        if (event.flux_tubes && !wide) wide = &event;
    }
    check(wide != nullptr, "an event shows them");
    if (!wide) return;
    double smallest_gap[3] = {1e30, 1e30, 1e30}, largest_gap[3] = {0.0, 0.0, 0.0}, largest_jump[3] = {0.0, 0.0, 0.0};
    bool ends_right = true, leads_east = true;
    std::vector<glm::dvec3> line, wing;
    for (int step = 0; step < 300; ++step) {
        scene.update(wide->t_tdb + step * 1800.0);
        for (int moon = 0; moon < 3; ++moon) {
            const glm::dvec3& p = jupiter.jovian_moons[static_cast<size_t>(moon)];
            for (int south = 0; south < 2; ++south) {
                glm::dvec3 instant, spot;
                const bool have = model.flux_tube(moon, south != 0, p, false, line) &&
                                  model.flux_tube(moon, south != 0, p, true, wing) &&
                                  model.instant_footprint(moon, south != 0, p, instant) &&
                                  model.footprint(moon, south != 0, p, spot);
                check(have, "the moons stay inside the tables");
                if (!have) return;
                ends_right = ends_right && line.front() == p && wing.front() == p &&
                             glm::length(glm::normalize(line.back()) - instant) < 1e-9 &&
                             glm::length(glm::normalize(wing.back()) - spot) < 1e-9 &&
                             glm::length(spot - jupiter.jovian_footprints[static_cast<size_t>(2 * moon + south)]) < 1e-9;
                const double gap = glm::length(line.back() - wing.back());
                smallest_gap[moon] = std::min(smallest_gap[moon], gap);
                largest_gap[moon] = std::max(largest_gap[moon], gap);
                const double lead = JovianAurora::lead_angle(moon, south != 0, -std::atan2(p.y, p.x));
                // Halfway down the wave is east of the line, by no more than the whole lead and a bit.
                const glm::dvec3& a = line[line.size() / 2];
                const glm::dvec3& b = wing[wing.size() / 2];
                const double ahead = std::remainder(std::atan2(b.y, b.x) - std::atan2(a.y, a.x), kTwoPi);
                leads_east = leads_east && ahead > 0.0 && ahead < lead * 1.3 + 2.0 * kDegToRad;
                for (const std::vector<glm::dvec3>* path : {&line, &wing}) {
                    for (size_t k = 1; k < path->size(); ++k) {
                        largest_jump[moon] = std::max(largest_jump[moon], glm::length((*path)[k] - (*path)[k - 1]));
                    }
                }
            }
        }
    }
    check(ends_right, "both lines start at the moon, the field line ends at the instant footprint and the wave at the spot");
    check(leads_east, "the wave runs east of the field line, as far as the lead angle");
    for (int moon = 0; moon < 3; ++moon) {
        std::printf("  info: %s: the spot is %.0f .. %.0f km from the end of the field line, longest step %.0f km\n",
                    JovianAurora::kMoons[static_cast<size_t>(moon)], smallest_gap[moon], largest_gap[moon], largest_jump[moon]);
        check(smallest_gap[moon] > 300.0, "the spot is well clear of the line's end, larger than the spot", smallest_gap[moon]);
        check(largest_jump[moon] < 3.0e5, "and the lines are smooth", largest_jump[moon]);
    }
    check(largest_gap[2] > largest_gap[0], "Ganymede's gap is the largest", largest_gap[2] / largest_gap[0]);

    // The field at the moons, from the same model (JRM33 + CON2020): along the table's lines
    // where they start. Its strength at Ganymede is 72 to 117 nT, of which 70 to 79 southward
    // (compared with Galileo's passes below); at Io about 2,000 nT.
    for (int moon = 0; moon < 3; ++moon) {
        double weakest = 1e30, strongest = 0.0, least_bz = 1e30, most_bz = -1e30, worst_alignment = 1.0;
        for (double east = 0.0; east < 360.0; east += 2.5) {
            const glm::dvec3 p = nominal[moon] * glm::dvec3(std::cos(east * kDegToRad), std::sin(east * kDegToRad), 0.0);
            glm::dvec3 field;
            check(model.field_at_moon(moon, p, field), "the table has the field at the moon");
            weakest = std::min(weakest, glm::length(field));
            strongest = std::max(strongest, glm::length(field));
            least_bz = std::min(least_bz, field.z);
            most_bz = std::max(most_bz, field.z);
            // Jupiter's field points south at the equator: the northern line leaves against it.
            if (model.flux_tube(moon, false, p, false, line)) {
                worst_alignment = std::min(worst_alignment, glm::dot(glm::normalize(line[1] - line[0]), -glm::normalize(field)));
            }
        }
        std::printf("  info: %s stands in %.0f .. %.0f nT, %.0f .. %.0f along the spin axis; lines leave within %.2f deg of it\n",
                    JovianAurora::kMoons[static_cast<size_t>(moon)], weakest, strongest, least_bz, most_bz,
                    std::acos(std::min(worst_alignment, 1.0)) / kDegToRad);
        // By the end of the first step of the curve, a third of the way to the table's next point.
        check(worst_alignment > std::cos(8.0 * kDegToRad), "the lines leave the moon along the field", worst_alignment);
        if (moon == 0) check(weakest > 1700.0 && strongest < 2200.0, "Io stands in about 2,000 nT", strongest);
        if (moon == 2) {
            check(weakest > 60.0 && weakest < 80.0 && strongest > 105.0 && strongest < 125.0, "Ganymede in 70 to 115 nT", strongest);
            check(least_bz > -90.0 && most_bz < -65.0, "of which 70 to 80 nT southward", least_bz);
        }
    }
    // Kivelson, Khurana & Volwerk 2002 (Icarus 157, 507), Table I: Jupiter's field around
    // Ganymede at Galileo's six passes, from the data before and after each, in GphiO (X along
    // the corotating flow, Y from Ganymede to Jupiter, Z along the spin axis; nT).
    {
        const struct { const char* pass; const char* utc; glm::dvec3 field; } passes[] = {
            {"G1", "1996-06-27 06:29:07", {6.0, -79.0, -79.0}},  {"G2", "1996-09-06 18:59:34", {17.0, -73.0, -85.0}},
            {"G7", "1997-04-05 07:09:58", {-3.0, 84.0, -76.0}},  {"G8", "1997-05-07 15:56:09", {-11.0, 11.0, -77.0}},
            {"G28", "2000-05-20 10:10:10", {-7.0, 78.0, -76.0}}, {"G29", "2000-12-28 08:58:00", {-9.0, -83.0, -79.0}},
        };
        double worst = 0.0, sum = 0.0;
        for (const auto& pass : passes) {
            double t = 0.0;
            check(parse_utc(pass.utc, &t), "a time parses");
            scene.update(t);
            const glm::dvec3& p = jupiter.jovian_moons[2];
            glm::dvec3 field;
            check(model.field_at_moon(2, p, field), "the table has the field at Ganymede");
            const glm::dvec3 toward_jupiter = -glm::normalize(glm::dvec3(p.x, p.y, 0.0));
            const glm::dvec3 spin(0.0, 0.0, 1.0);
            const glm::dvec3 in_gphio(glm::dot(field, glm::cross(toward_jupiter, spin)), glm::dot(field, toward_jupiter), field.z);
            const double miss = glm::length(in_gphio - pass.field);
            std::printf("  info: %s: Ganymede at %.1f deg W, the model's field (%.0f, %.0f, %.0f) nT, Galileo's (%.0f, %.0f, %.0f): %.0f nT apart\n",
                        pass.pass, wrap_two_pi(-std::atan2(p.y, p.x)) / kDegToRad, in_gphio.x, in_gphio.y, in_gphio.z, pass.field.x,
                        pass.field.y, pass.field.z, miss);
            worst = std::max(worst, miss);
            sum += miss * miss;
            check(in_gphio.y * pass.field.y > 0.0 || std::abs(pass.field.y) < 20.0, "the field leans toward or away from Jupiter as measured");
        }
        std::printf("  info: the model's field at Ganymede is within %.0f nT of Galileo's at the six passes, rms %.0f\n", worst,
                    std::sqrt(sum / 6.0));
        check(worst < 25.0, "Jupiter's field at Ganymede is as Galileo found it at its six passes", worst);
        scene.update(wide->t_tdb);
    }
    glm::dvec3 none;
    check(!model.field_at_moon(3,glm::dvec3(1.0e6, 0.0, 0.0), none) && none == glm::dvec3(0.0), "no field for another moon");

    // Outside the table, nothing: a trace is not extrapolated.
    const glm::dvec3 too_far(1.5 * nominal[0], 0.0, 0.0);
    check(!model.flux_tube(0, false, too_far, false, line), "no line outside the table's domain");
    check(!model.flux_tube(3, false, too_far, false, line), "nor for another moon");
}

void test_magnetic_maps()
{
    Scene mars_scene = load_scene_or_die("mars.toml");
    Scene jupiter_scene = load_scene_or_die("jupiter.toml");
    const Body& mars = mars_scene.bodies[static_cast<size_t>(mars_scene.find("Mars"))];
    const Body& jupiter = jupiter_scene.bodies[static_cast<size_t>(jupiter_scene.find("Jupiter"))];
    check(mars.magnetic_map.surface && !mars.magnetic_map.deep, "Mars has a map of its crustal field");
    check(jupiter.magnetic_map.surface && jupiter.magnetic_map.deep, "Jupiter has two");
    if (!mars.magnetic_map.surface || !jupiter.magnetic_map.surface || !jupiter.magnetic_map.deep) return;

    const auto cell_latitude = [](const FieldMap& m, int row) { return (0.5 - (row + 0.5) / m.height) * kPi; };
    const auto cell_longitude = [](const FieldMap& m, int column) { return ((column + 0.5) / m.width - 0.5) * kTwoPi; };

    // Mars: the field is where Acuna 1999 found it, and not in the great basins.
    const FieldMap& crust = *mars.magnetic_map.surface;
    double peak = 0.0, peak_latitude = 0.0, peak_east = 0.0, north = 0.0, south = 0.0;
    for (int row = 0; row < crust.height; ++row) {
        const double latitude = cell_latitude(crust, row);
        for (int column = 0; column < crust.width; ++column) {
            const double east = cell_longitude(crust, column);
            const double br = crust.value_nt(east, latitude);
            if (std::abs(br) > peak) {
                peak = std::abs(br);
                peak_latitude = latitude;
                peak_east = east;
            }
            (latitude > 0.0 ? north : south) += br * br * std::cos(latitude);
        }
    }
    check(peak > 1000.0 && peak < 1500.0, "Mars's strongest radial field 150 km up is near 1,200 nT (Langlais 2019)", peak);
    check(std::abs(peak_latitude / kDegToRad + 45.0) < 15.0 && std::abs(angle_difference_deg(peak_east / kDegToRad, 180.0)) < 30.0,
          "in Terra Cimmeria and Terra Sirenum", peak_east / kDegToRad);
    check(south > 5.0 * north, "the southern highlands hold most of it", south / north);
    const auto basin = [&](double latitude_deg, double east_deg, double radius_deg) {
        double strongest = 0.0;
        for (int a = 0; a < 16; ++a) {
            for (double r = 0.0; r <= radius_deg; r += radius_deg / 4.0) {
                const double b = a * kTwoPi / 16.0;
                const double latitude = latitude_deg + r * std::sin(b);
                const double east = east_deg + r * std::cos(b) / std::cos(latitude * kDegToRad);
                strongest = std::max(strongest, std::abs(crust.value_nt(east * kDegToRad, latitude * kDegToRad)));
            }
        }
        return strongest;
    };
    check(basin(-42.4, 70.5, 12.0) < 0.1 * peak, "Hellas shows almost none", basin(-42.4, 70.5, 12.0) / peak);
    check(basin(-49.7, 316.0, 8.0) < 0.1 * peak, "nor Argyre", basin(-49.7, 316.0, 8.0) / peak);

    // Jupiter: both maps are JRM33 where the body's mesh puts each cell.
    const MagneticField& jrm33 = *jupiter.magnetic_field.model;
    const auto check_map = [&](const FieldMap& map, bool deep, const char* what) {
        double worst = 0.0;
        for (int row = 3; row < map.height; row += 17) {
            const double beta = cell_latitude(map, row);
            const double cos_latitude = std::cos(beta);
            const double sin_latitude = std::sin(beta);
            const double r = deep ? 0.81 : 1.0 / std::hypot(cos_latitude, sin_latitude * 71492.0 / 66854.0);
            for (int column = 2; column < map.width; column += 23) {
                const double east = cell_longitude(map, column);
                const glm::dvec3 direction(cos_latitude * std::cos(east), cos_latitude * std::sin(east), sin_latitude);
                const double expected = glm::dot(jrm33.field_nt(r * direction), direction);
                const double value = map.value_nt(east, beta);
                // A byte's step of the square root, in nT, near this value.
                const double step = 2.0 * std::sqrt(std::abs(value) * map.scale_nt) / 127.0 + map.scale_nt / (127.0 * 127.0);
                worst = std::max(worst, std::abs(value - expected) / step);
            }
        }
        check(worst <= 1.0, what, worst);
    };
    check_map(*jupiter.magnetic_map.surface, false, "the cloud-top map is JRM33 on the flattened surface");
    check_map(*jupiter.magnetic_map.deep, true, "the deep map is JRM33 at 0.81 radii");

    // The Great Blue Spot: near the equator, near 90 W (Moore 2018), and it stands out deep down.
    const auto spot = [&](const FieldMap& map, double& strength, double& latitude, double& west) {
        strength = 0.0;
        double sum = 0.0;
        for (int row = 0; row < map.height; ++row) {
            const double beta = cell_latitude(map, row);
            for (int column = 0; column < map.width; ++column) {
                const double east = cell_longitude(map, column);
                const double br = map.value_nt(east, beta);
                sum += br * br * std::cos(beta);
                if (std::abs(beta) < 15.0 * kDegToRad && br < strength) {
                    strength = br;
                    latitude = beta / kDegToRad;
                    west = -east / kDegToRad;
                }
            }
        }
        return std::sqrt(sum / (map.width * map.height * 2.0 / kPi));
    };
    double strength = 0.0, latitude = 0.0, west = 0.0;
    const double top_rms = spot(*jupiter.magnetic_map.surface, strength, latitude, west);
    check(std::abs(strength + 1.28e6) < 0.05e6, "the Great Blue Spot is about 13 gauss at the cloud tops", strength);
    check(std::abs(latitude) < 12.0 && std::abs(angle_difference_deg(west, 90.0)) < 10.0, "near the equator, near 90 W", west);
    const double top_contrast = -strength / top_rms;
    const double deep_rms = spot(*jupiter.magnetic_map.deep, strength, latitude, west);
    check(std::abs(strength + 7.9e6) < 0.2e6, "and about 80 gauss at 0.81 radii", strength);
    check(std::abs(angle_difference_deg(west, 84.0)) < 5.0, "in the same place", west);
    check(-strength / deep_rms > 2.0 * top_contrast, "where it stands out far more", -strength / deep_rms / top_contrast);

    // The Moon: the strongest field by Gerasimovich (22.9 S, 122.6 W), on the far side, a patch
    // under Reiner Gamma (7.5 N, 59 W), almost nothing in Mare Crisium (17 N, 59 E), and all of
    // it far weaker than any planet's.
    Scene earth_scene = load_scene_or_die("earth_moon.toml");
    const Body& moon = earth_scene.bodies[static_cast<size_t>(earth_scene.find("Moon"))];
    const Body& earth = earth_scene.bodies[static_cast<size_t>(earth_scene.find("Earth"))];
    check(moon.magnetic_map.surface && earth.magnetic_map.surface, "the Moon and the Earth have maps");
    if (!moon.magnetic_map.surface || !earth.magnetic_map.surface) return;
    const FieldMap& lunar = *moon.magnetic_map.surface;
    check(!lunar.magnitude, "the Moon's is of the radial field");
    double lunar_peak = 0.0, lunar_latitude = 0.0, lunar_east = 0.0, near_side = 0.0, far_side = 0.0;
    double near_area = 0.0, far_area = 0.0;
    for (int row = 0; row < lunar.height; ++row) {
        const double beta = cell_latitude(lunar, row);
        for (int column = 0; column < lunar.width; ++column) {
            const double east = cell_longitude(lunar, column);
            const double br = lunar.value_nt(east, beta);
            if (std::abs(br) > lunar_peak) {
                lunar_peak = std::abs(br);
                lunar_latitude = beta / kDegToRad;
                lunar_east = east / kDegToRad;
            }
            const bool near = std::abs(east) < kPi / 2.0;
            (near ? near_side : far_side) += br * br * std::cos(beta);
            (near ? near_area : far_area) += std::cos(beta);
        }
    }
    check(lunar_peak > 15.0 && lunar_peak < 40.0, "the Moon's strongest radial field 30 km up is about 25 nT", lunar_peak);
    check(std::abs(lunar_latitude + 22.9) < 3.0 && std::abs(angle_difference_deg(lunar_east, -122.6)) < 3.0,
          "by the crater Gerasimovich", lunar_east);
    check(far_side / far_area > 2.0 * near_side / near_area, "the far side's field is the stronger", (far_side / far_area) / (near_side / near_area));
    const auto lunar_strongest = [&](double latitude_deg, double east_deg, double radius_deg) {
        double most = 0.0;
        for (double a = -radius_deg; a <= radius_deg; a += 0.5) {
            for (double b = -radius_deg; b <= radius_deg; b += 0.5) {
                most = std::max(most, std::abs(lunar.value_nt((east_deg + a) * kDegToRad, (latitude_deg + b) * kDegToRad)));
            }
        }
        return most;
    };
    const double reiner = lunar_strongest(7.5, -59.0, 2.0);
    const double crisium = lunar_strongest(17.0, 59.0, 2.0);
    check(reiner > 4.0, "Reiner Gamma's patch is several nT", reiner);
    check(reiner > 5.0 * crisium, "and Mare Crisium has little", reiner / crisium);

    // The Earth: the strength of IGRF-14 cell by cell, weakest in the South Atlantic Anomaly.
    const FieldMap& strength_map = *earth.magnetic_map.surface;
    const MagneticField& igrf = *earth.magnetic_field.model;
    check(strength_map.magnitude, "the Earth's is of the strength");
    double worst_strength = 0.0, weakest = 1e30, weakest_latitude = 0.0, weakest_east = 0.0;
    int weakest_cells = 0;
    for (int row = 0; row < strength_map.height; ++row) {
        const double beta = cell_latitude(strength_map, row);
        const double r = 1.0 / std::hypot(std::cos(beta) / 6378.1366, std::sin(beta) / 6356.7519) / 6371.2;
        for (int column = 0; column < strength_map.width; ++column) {
            const double east = cell_longitude(strength_map, column);
            const double value = strength_map.value_nt(east, beta);
            // Several cells share the lowest byte: their middle.
            if (value < weakest) {
                weakest = value;
                weakest_latitude = 0.0;
                weakest_east = 0.0;
                weakest_cells = 0;
            }
            if (value == weakest) {
                weakest_latitude += beta / kDegToRad;
                weakest_east += east / kDegToRad;
                ++weakest_cells;
            }
            if (row % 11 == 0 && column % 13 == 0) {
                const glm::dvec3 direction(std::cos(beta) * std::cos(east), std::cos(beta) * std::sin(east), std::sin(beta));
                const double expected_nt = glm::length(igrf.field_nt(r * direction));
                worst_strength = std::max(worst_strength, std::abs(value - expected_nt));
            }
        }
    }
    weakest_latitude /= weakest_cells;
    weakest_east /= weakest_cells;
    const double byte_nt = (strength_map.scale_nt - strength_map.min_nt) / 254.0;
    check(worst_strength <= byte_nt, "the Earth's map is IGRF-14 on the flattened surface", worst_strength);
    check(std::abs(weakest - 22100.0) < 300.0 && std::abs(weakest_latitude + 27.0) < 3.0 && std::abs(angle_difference_deg(weakest_east, -59.0)) < 4.0,
          "weakest over the South Atlantic, about 22,000 nT", weakest);
    check(strength_map.scale_nt > 2.9 * weakest, "a third of the strongest field", strength_map.scale_nt / weakest);

    // The field of other years (IGRF-14, 1900 to 2030): the model, its lines and this map follow
    // the clock by whole years.
    {
        Body& body = earth_scene.bodies[static_cast<size_t>(earth_scene.find("Earth"))];
        const BodyMagneticField& field = body.magnetic_field;
        check(field.series != nullptr && field.model_year == 2025.0, "the Earth's field is a series, at 2025 until the clock says otherwise");
        if (!field.series) return;
        check(field.series->first_year() == 1900.0 && field.series->last_epoch() == 2025.0 && field.series->last_year() == 2030.0,
              "IGRF-14 spans 1900 to 2025, with rates to 2030");
        const std::vector<signed char> baked = strength_map.cells;
        const double baked_low = strength_map.min_nt, baked_high = strength_map.scale_nt;
        const std::vector<glm::dvec4> line_2025 = field.lines.front();
        // Where the field is weakest, from the model itself (the map's bytes are 193 nT apart).
        const auto weakest_of = [&](double& latitude, double& east) {
            double least = 1e30;
            for (int row = 0; row < 180; ++row) {
                const double beta = (89.5 - row) * kDegToRad;
                const double r = 1.0 / std::hypot(std::cos(beta) / 6378.1366, std::sin(beta) / 6356.7519) / 6371.2;
                for (int column = 0; column < 360; ++column) {
                    const double lambda = (column - 179.5) * kDegToRad;
                    const double b = glm::length(field.model->field_nt(
                        r * glm::dvec3(std::cos(beta) * std::cos(lambda), std::cos(beta) * std::sin(lambda), std::sin(beta))));
                    if (b < least) {
                        least = b;
                        latitude = beta / kDegToRad;
                        east = lambda / kDegToRad;
                    }
                }
            }
            return least;
        };
        check(!follow_magnetic_field_year(body, tdb_from_julian_year(2025.3)), "nothing changes within the year");
        double latitude_2025 = 0.0, east_2025 = 0.0;
        const double weakest_2025 = weakest_of(latitude_2025, east_2025);

        // The same year from another: the map comes back to the baked bytes.
        check(follow_magnetic_field_year(body, tdb_from_julian_year(1900.0)) && field.model_year == 1900.0, "1900 is another field");
        check(field.model->g(1, 0) == -31543.0 && field.model->h(1, 1) == 5922.0 && field.model->degree() == 13,
              "with the coefficients of its epoch", field.model->g(1, 0));
        double latitude_1900 = 0.0, east_1900 = 0.0;
        const double weakest_1900 = weakest_of(latitude_1900, east_1900);
        const FieldMap& map_1900 = *body.magnetic_map.surface;
        check(map_1900.magnitude && map_1900.cells != baked && map_1900.min_nt == baked_low && map_1900.scale_nt == baked_high,
              "the map is redrawn in the same range");
        check(std::abs(map_1900.value_nt(east_1900 * kDegToRad, latitude_1900 * kDegToRad) - weakest_1900) < byte_nt,
              "and holds the model's weakest field", weakest_1900);
        std::printf("  info: weakest field %.0f nT at %.1f, %.1f E in 1900; %.0f nT at %.1f, %.1f E in 2025\n", weakest_1900,
                    latitude_1900, east_1900, weakest_2025, latitude_2025, east_2025);
        // The anomaly has deepened and drifted west (about 0.2 deg a year).
        check(weakest_1900 > 25000.0 && weakest_1900 < 26000.0 && weakest_1900 > 1.13 * weakest_2025,
              "the South Atlantic Anomaly was a seventh stronger in 1900", weakest_1900 / weakest_2025);
        const double drift = angle_difference_deg(east_1900, east_2025);
        check(drift > 18.0 && drift < 28.0 && latitude_1900 > latitude_2025, "and some 23 deg farther east and nearer the equator", drift);
        check(glm::length(glm::dvec3(field.lines.front().front()) - glm::dvec3(line_2025.front())) > 0.01,
              "the lines are traced again", glm::length(glm::dvec3(field.lines.front().front()) - glm::dvec3(line_2025.front())));
        // The dipole of each year is the aurora's (its own table of the same terms).
        double worst_pole = 0.0;
        for (double year : {1900.0, 1933.0, 1962.0, 1987.0, 2004.0, 2023.0, 2028.0, 2030.0}) {
            follow_magnetic_field_year(body, tdb_from_julian_year(year));
            const glm::dvec3 pole = geomagnetic_north_pole(tdb_from_julian_year(year));
            // By the cross product: the arc cosine of a dot product within a rounding error of 1
            // is anywhere up to 1e-6 deg, whatever the vectors.
            worst_pole = std::max(worst_pole, std::asin(std::min(glm::length(glm::cross(pole, -field.model->dipole_axis())), 1.0)) / kDegToRad);
        }
        check(worst_pole < 1e-6, "between the epochs and past the last, the dipole is the aurora's geomagnetic pole", worst_pole);
        check(std::abs(field.model->g(1, 0) - (-29350.0 + 5.0 * 12.6)) < 1e-9, "2030 is 2025 with five years of its rates", field.model->g(1, 0));
        // Held beyond both ends.
        follow_magnetic_field_year(body, tdb_from_julian_year(2300.0));
        check(field.model_year == 2030.0, "later years keep the field of 2030", field.model_year);
        follow_magnetic_field_year(body, tdb_from_julian_year(1700.0));
        check(field.model_year == 1900.0, "earlier ones that of 1900", field.model_year);
        // Scene::update does it, and 2025 gives the baked map back byte for byte.
        earth_scene.update(tdb_from_julian_year(2025.2));
        check(field.model_year == 2025.0 && body.magnetic_map.surface->cells == baked,
              "the map of 2025 is the baked one, cell for cell");
        check(field.lines.front().size() == line_2025.size() &&
                  glm::length(glm::dvec3(field.lines.front().front()) - glm::dvec3(line_2025.front())) < 1e-12,
              "and its lines those of the loaded scene");

        Scene bad;
        std::string error;
        const std::string head = "name = 'x'\n[[bodies]]\nname = 'P'\nradii_km = [6378.0]\n[bodies.magnetic_field]\nreference_radius_km = 6371.2\n";
        check(load_scene_string((head + "coefficients_by_year = 'geomag/field_earth_igrf14.txt'\n").c_str(), "by_year", bad, &error, ASTRAXIS_ASSET_DIR),
              "a field by year loads");
        check(!load_scene_string((head + "coefficients_by_year = 'geomag/field_earth_igrf14.txt'\ncoefficients = 'geomag/field_neptune_o8.txt'\n").c_str(),
                                 "both", bad, &error, ASTRAXIS_ASSET_DIR),
              "coefficients and coefficients_by_year exclude each other");
        check(!load_scene_string((head + "coefficients_by_year = 'geomag/field_neptune_o8.txt'\n").c_str(), "not_series", bad, &error, ASTRAXIS_ASSET_DIR),
              "a single model is not a series");
    }

    // The Earth's crust: its own map, at most a tenth of the main field and usually far less.
    check(earth.magnetic_map.crust != nullptr, "the Earth has a map of its crust");
    check(!moon.magnetic_map.crust && !jupiter.magnetic_map.crust, "and nothing else does");
    if (!earth.magnetic_map.crust) return;
    const FieldMap& crust_map = *earth.magnetic_map.crust;
    check(!crust_map.magnitude && crust_map.width == 2048 && crust_map.height == 1024, "the crust's is of the radial field, about 20 km a cell");
    check(crust_map.scale_nt == 600.0, "saturating at 600 nT", crust_map.scale_nt);
    const auto crust_at = [&](double latitude_deg, double east_deg) {
        return crust_map.value_nt(east_deg * kDegToRad, latitude_deg * kDegToRad);
    };
    const auto crust_rms_near = [&](double latitude_deg, double east_deg, double radius_deg, double* most_in = nullptr) {
        double sum = 0.0;
        int count = 0;
        for (double a = -radius_deg; a <= radius_deg; a += 0.25) {
            for (double b = -radius_deg; b <= radius_deg; b += 0.25) {
                const double value = crust_at(latitude_deg + b, east_deg + a);
                if (most_in) *most_in = std::max(*most_in, -value);
                sum += value * value;
                ++count;
            }
        }
        return std::sqrt(sum / count);
    };
    // Kursk at 51 N, 37 E; Bangui around 6 N, 18 E. The South Pacific, far from the continents
    // and the ridges, has little.
    double kursk_in = 0.0;
    check(crust_rms_near(51.0, 37.0, 2.0, &kursk_in) > 300.0, "the Kursk anomaly stands out", crust_rms_near(51.0, 37.0, 2.0));
    check(kursk_in >= 599.0, "and goes in, past the top of the scale", kursk_in);
    const double bangui = crust_rms_near(6.0, 18.0, 6.0);
    const double pacific = crust_rms_near(-30.0, -160.0, 6.0);
    check(bangui > 2.0 * pacific, "the Bangui anomaly is far stronger than the open Pacific", bangui / pacific);
    double crust_square = 0.0, crust_area = 0.0;
    for (int row = 0; row < crust_map.height; row += 7) {
        const double beta = cell_latitude(crust_map, row);
        for (int column = 0; column < crust_map.width; column += 7) {
            const double value = crust_map.value_nt(cell_longitude(crust_map, column), beta);
            crust_square += value * value * std::cos(beta);
            crust_area += std::cos(beta);
        }
    }
    const double crust_rms = std::sqrt(crust_square / crust_area);
    check(crust_rms > 50.0 && crust_rms < 200.0, "about a hundred nT over the whole sphere, below a third of a percent of the main field", crust_rms);

    // Uranus and Neptune: the models where the mesh puts each cell, also as Neptune's model turns
    // against the body, and the strongest field on the side the off-center dipole is nearest.
    Scene giants = load_scene_or_die("solar_system.toml");
    const Body& uranus = giants.bodies[static_cast<size_t>(giants.find("Uranus"))];
    const Body& neptune = giants.bodies[static_cast<size_t>(giants.find("Neptune"))];
    check(uranus.magnetic_map.surface && neptune.magnetic_map.surface, "Uranus and Neptune have maps");
    if (!uranus.magnetic_map.surface || !neptune.magnetic_map.surface) return;
    const auto check_giant = [&](const Body& body, double t, double reference_km, const char* what) {
        const FieldMap& map = *body.magnetic_map.surface;
        const MagneticField& model = *body.magnetic_field.model;
        const glm::dmat3 model_to_body = magnetic_model_to_body(body, t);
        const double turn = magnetic_map_turn_deg(body, t) * kDegToRad;
        double worst = 0.0;
        for (int row = 3; row < map.height; row += 7) {
            const double latitude = cell_latitude(map, row);
            const double r = 1.0 / std::hypot(std::cos(latitude) / body.equatorial_radius_km, std::sin(latitude) / body.polar_radius_km);
            for (int column = 2; column < map.width; column += 11) {
                const double east = cell_longitude(map, column) + turn;
                const glm::dvec3 in_body(std::cos(latitude) * std::cos(east), std::cos(latitude) * std::sin(east), std::sin(latitude));
                const glm::dvec3 in_model = glm::transpose(model_to_body) * in_body;
                const double expected = glm::dot(model.field_nt(r / reference_km * in_model), in_model);
                const double value = map.value_nt(east - turn, latitude);
                const double step = 2.0 * std::sqrt(std::abs(value) * map.scale_nt) / 127.0 + map.scale_nt / (127.0 * 127.0);
                worst = std::max(worst, std::abs(value - expected) / step);
            }
        }
        check(worst <= 1.0, what, worst);
    };
    check_giant(uranus, 0.0, 25600.0, "the map of Uranus is AH5 in the body's axes, wherever its pole points");
    check_giant(neptune, 0.0, 24765.0, "the map of Neptune is O8");
    check_giant(neptune, -3.3e8, 24765.0, "also in 1989, when the model had turned against the body");
    const auto check_lopsided = [&](const Body& body, double t, double most_nt, const char* what) {
        const FieldMap& map = *body.magnetic_map.surface;
        const MagneticField& model = *body.magnetic_field.model;
        double peak_value = 0.0, peak_latitude = 0.0, peak_east = 0.0, other_side = 0.0;
        for (int row = 0; row < map.height; ++row) {
            for (int column = 0; column < map.width; ++column) {
                const double value = map.value_nt(cell_longitude(map, column), cell_latitude(map, row));
                if (std::abs(value) > std::abs(peak_value)) {
                    peak_value = value;
                    peak_latitude = cell_latitude(map, row);
                    peak_east = cell_longitude(map, column);
                }
            }
        }
        for (int row = 0; row < map.height; ++row) {
            for (int column = 0; column < map.width; ++column) {
                const double value = map.value_nt(cell_longitude(map, column), cell_latitude(map, row));
                if (value * peak_value < 0.0) other_side = std::max(other_side, std::abs(value));
            }
        }
        const double east = peak_east + magnetic_map_turn_deg(body, t) * kDegToRad;
        const glm::dvec3 peak_direction(std::cos(peak_latitude) * std::cos(east), std::cos(peak_latitude) * std::sin(east), std::sin(peak_latitude));
        const glm::dvec3 offset = glm::normalize(magnetic_model_to_body(body, t) * model.dipole_offset());
        const double from_offset = std::acos(glm::dot(peak_direction, offset)) / kDegToRad;
        std::printf("  info: %s: strongest %.0f nT, opposite sign at most %.0f, %.0f deg from the dipole's offset\n", what, peak_value, other_side, from_offset);
        check(std::abs(std::abs(peak_value) - most_nt) < 3000.0, what, peak_value);
        check(from_offset < 50.0, "on the side the dipole's center is displaced to", from_offset);
        check(std::abs(peak_value) > 3.0 * model.dipole_nt(), "more than the 2 g of a centered dipole's pole", std::abs(peak_value) / model.dipole_nt());
        check(std::abs(peak_value) > 1.4 * other_side, "and well above the other sign's", std::abs(peak_value) / other_side);
    };
    check_lopsided(uranus, 0.0, 149000.0, "Uranus's strongest radial field at the surface");
    check_lopsided(neptune, -3.3e8, 89000.0, "Neptune's");

    // The events open with the feature coming round into the sunlight, or for the Moon, which
    // turns slowly, with its two patches on either side of the sub-solar point.
    Scene uranus_scene = load_scene_or_die("uranus.toml");
    Scene neptune_scene = load_scene_or_die("neptune.toml");
    const auto sub_solar_east = [](Scene& scene, const Body& body, double t) {
        scene.update(t);
        const glm::dvec3 to_sun = glm::transpose(body.orientation) * (scene.sun_position() - body.world_position);
        return std::atan2(to_sun.y, to_sun.x) / kDegToRad;
    };
    // Uranus turns the other way in longitude: the Sun comes to a feature from its west.
    // The Earth's strength map has two events: 2025 and 1900, the anomaly then 23 deg farther east.
    const struct { Scene* scene; const Body* body; double opening_east; FieldMapView view; int events = 1; } expected[] = {
        {&mars_scene, &mars, 215.0, FieldMapView::Surface},
        {&jupiter_scene, &jupiter, 315.5, FieldMapView::Surface},
        {&jupiter_scene, &jupiter, 316.5, FieldMapView::Deep},
        {&earth_scene, &moon, 270.0, FieldMapView::Surface},
        {&earth_scene, &earth, 325.5, FieldMapView::Surface, 2},
        {&earth_scene, &earth, 50.0, FieldMapView::Crust},
        {&uranus_scene, &uranus_scene.bodies[static_cast<size_t>(uranus_scene.find("Uranus"))], 155.5, FieldMapView::Surface},
        {&neptune_scene, &neptune_scene.bodies[static_cast<size_t>(neptune_scene.find("Neptune"))], 117.5, FieldMapView::Surface},
    };
    for (const auto& e : expected) {
        int found = 0;
        for (const SceneEvent& event : e.scene->events) {
            if (event.magnetic_map != e.view || event.focus != e.scene->find(e.body->name.c_str())) continue;
            ++found;
            check(event.focus_orbits_only || event.hide_orbits, "the event looks at the body, without distant orbits across it");
            // Neptune's map is given in its model's longitudes.
            const double opening_east = e.opening_east + magnetic_map_turn_deg(*e.body, event.t_tdb);
            const double offset = angle_difference_deg(sub_solar_east(*e.scene, *e.body, event.t_tdb), opening_east);
            std::printf("  info: %s: sub-solar longitude %+.1f deg from where it should open\n", event.name.c_str(), offset);
            check(std::abs(offset) < 5.0, "and opens with the feature coming round", offset);
        }
        check(found == e.events, "the body has its map event", static_cast<double>(found));
    }

    // The texture holds a value's sign about the middle, for interpolation.
    FieldMap small;
    small.width = 3;
    small.height = 1;
    small.scale_nt = 100.0;
    small.cells = {-127, 0, 127};
    const std::vector<unsigned char> pixels = small.texture_rgba();
    check(pixels[0] == 0 && pixels[4] == 128 && pixels[8] == 255, "the texture holds -1 .. 1 as 0 .. 255");
    check(small.value_nt(-kPi + 0.1, 0.0) == -100.0 && small.value_nt(kPi - 0.1, 0.0) == 100.0, "and the map's sign and scale");

    // A wrong description is refused.
    const auto loads = [](const std::string& extra) {
        const std::string text = "name = 'x'\n[[bodies]]\nname = 'Planet'\nradii_km = [3396.0]\n[bodies.magnetic_map]\n" + extra;
        Scene out;
        std::string error;
        return load_scene_string(text, "test", out, &error, ASTRAXIS_ASSET_DIR);
    };
    check(loads("surface = 'geomag/mars_br.fieldmap'\n"), "a map loads");
    check(!loads("deep = 'geomag/mars_br.fieldmap'\n"), "a body needs a surface map");
    check(!loads("surface = 'geomag/jupiter_oval.txt'\n"), "and a file of the right kind");
    const std::string event_head = "name = 'x'\n[[bodies]]\nname = 'Planet'\nradii_km = [1.0]\n[[events]]\nname = 'e'\ntime_tdb = '2000-01-01 12:00'\nfocus = 'Planet'\n";
    Scene bad;
    std::string error;
    check(!load_scene_string(event_head + "magnetic_map = 'nowhere'\n", "bad", bad, &error), "an event's map is surface or deep");
}

void test_jovian_footprints()
{
    Scene scene = load_scene_or_die("jupiter.toml");
    const int index = scene.find("Jupiter");
    const auto& model = *scene.bodies[static_cast<size_t>(index)].aurora.jovian;
    // Independent scipy DOP853 integration with a terminal spheroid event;
    // tools/geomag/jupiter_footprint_check.py. Input longitude 0 E, latitude
    // 0.37 N, nominal radius; each hemisphere has its own eastward lead.
    const glm::dvec3 oracle[3][2] = {
        {{.1818101235137252, .1240225903079459, .9754811510640435},
         {.4723908679903595, -.0355644655884435, -.8806713556297395}},
        {{.0885374619738735, .1132414547584144, .9896148193875346},
         {.4034604622620909, -.0292434238759410, -.9145296482625765}},
        {{.0302769013630110, .1401158610503802, .9896720945474642},
         {.3607216865484141, -.0081370681424426, -.9326380074689820}}
    };
    for (int m = 0; m < 3; ++m) {
        const auto& path = model.footpaths[static_cast<size_t>(m)];
        check(path.points.size() == 9 * JovianAurora::kLongitudes, "satellite magnetic table loads");
        for (int h = 0; h < 2; ++h) {
            const double left = JovianAurora::lead_angle(m, h != 0, -1e-8);
            const double right = JovianAurora::lead_angle(m, h != 0, 1e-8);
            check(std::abs(left-right) < 1e-8, "empirical lead is continuous across System III seam");
            for (int degree = -360; degree <= 720; degree += 7) {
                const double x = degree * kDegToRad;
                const double lead = JovianAurora::lead_angle(m, h != 0, x);
                check(lead > 0 && lead < 22*kDegToRad &&
                      std::abs(lead-JovianAurora::lead_angle(m, h != 0, x+kTwoPi)) < 1e-12,
                      "lead remains positive and periodic through multiple rotations");
            }
            glm::dvec3 a, b;
            const double radius = .5*(path.min_radius+path.max_radius)*71492;
            const double latitude = .37*kDegToRad;
            check(model.footprint(m, h != 0, {radius*std::cos(latitude), 0, radius*std::sin(latitude)}, a) &&
                  glm::length(a-oracle[m][h])*71492 < 30,
                  "footprints agree with independent integration after hemisphere-specific eastward lead");
            check(model.footprint(m, h != 0, {radius, -1e-4, 0}, a) &&
                  model.footprint(m, h != 0, {radius, 1e-4, 0}, b) && glm::length(a-b) < 1e-8,
                  "magnetic lookup and lead correction share a continuous seam");
            check(!model.footprint(m, h != 0, {1e7, 0, 0}, a) && a == glm::dvec3(0),
                  "out-of-domain radius suppresses footprint");
            check(!model.footprint(m, h != 0, {0, 0, radius}, a), "out-of-domain latitude suppresses footprint");
        }
    }
    check(std::abs(JovianAurora::lead_angle(0, false, kPi/2)/kDegToRad - 4.56) < 1e-12 &&
          std::abs(JovianAurora::lead_angle(0, true, kPi/2)/kDegToRad - 3.22) < 1e-12,
          "Io lead agrees with Hue 2023 table 1 at 90 degrees west");

    // Cover historical flybys, the UV events and ephemeris fallback dates.
    // The 3-D domain must contain the actual scene positions, not just a
    // nominal equatorial circular orbit.
    const double start = scene.events.back().t_tdb;
    for (int sample = 0; sample < 120; ++sample) {
        const double t = start + (sample - 60) * 337.17 * 86400.0;
        scene.update(t);
        const auto& p = scene.bodies[static_cast<size_t>(index)].jovian_footprints;
        for (size_t k = 0; k < p.size(); ++k) {
            check(std::abs(glm::length(p[k])-1.0) < 1e-10 && (k%2 == 0 ? p[k].z > .5 : p[k].z < -.5),
                  "all six footprints cover scene ephemerides and fallback orbits");
        }
    }
    scene.update(start);
    const auto original = scene.bodies[static_cast<size_t>(index)].jovian_footprints;
    scene.update(start + 3600);
    const auto advanced = scene.bodies[static_cast<size_t>(index)].jovian_footprints;
    check(glm::length(original[0]-advanced[0]) > .01, "footprints follow the satellite synodic motion");
    scene.update(start);
    check(scene.bodies[static_cast<size_t>(index)].jovian_footprints == original,
          "backward jump recovers the same footprints without animation history");
    scene.set_active_frame(0);
    scene.update(start);
    check(scene.bodies[static_cast<size_t>(index)].jovian_footprints == original,
          "display reference frame does not change magnetic coordinates");

    Simulation sim;
    std::string error;
    check(sim.load_scene(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", &error), "footprint simulation loads");
    sim.stop_tour();
    sim.clock().t_tdb = start;
    sim.clock().paused = true;
    sim.update(0);
    const auto current = sim.scene().bodies[static_cast<size_t>(index)].jovian_footprints;
    OutputView small, large;
    sim.compute_view(640, 480, 1, small);
    sim.compute_view(2560, 1440, 2, large);
    sim.update(.2);
    check(sim.scene().bodies[static_cast<size_t>(index)].jovian_footprints == current,
          "different output sizes and pause preserve the shared footprint positions");
    sim.clock().paused = false;
    sim.clock().warp = 3600;
    sim.clock().reverse = true;
    sim.update(1);
    scene.update(start-3600);
    check(sim.scene().bodies[static_cast<size_t>(index)].jovian_footprints ==
          scene.bodies[static_cast<size_t>(index)].jovian_footprints,
          "reverse playback uses current ephemerides without the aurora animation speed cap");
    Scene bad;
    check(!load_scene_string("name='x'\n[[bodies]]\nname='Earth'\nradii_km=[6378]\n"
                            "[bodies.aurora]\nkp='geomag/kp.bin'\nfootprints='geomag/jupiter_footprints.txt'\n",
                            "bad_footprints", bad, &error, ASTRAXIS_ASSET_DIR),
          "terrestrial Kp cannot acquire Jovian footprints");
    // --- The secondary spots (Bonfond et al. 2008, 2009, 2013) -------------------------------
    using Spot = JovianAurora::Spot;
    const auto lead_deg = [](int moon, bool south, double west_deg, Spot spot) {
        return JovianAurora::spot_lead_angle(moon, south, west_deg * kDegToRad, spot) / kDegToRad;
    };
    // The beam spot is on the other hemisphere's main spot's field line, so the two hemispheres'
    // gaps between main and beam spots are the same on the equator; the reflected wing crosses
    // the other half of the sheet twice more.
    bool mirrored = true;
    double widest_gap = 0.0, farthest_reflected[2] = {0.0, 0.0};
    for (double west = 0.0; west < 360.0; west += 1.0) {
        for (int moon = 0; moon < 3; ++moon) {
            const double north_gap = lead_deg(moon, false, west, Spot::Beam) - lead_deg(moon, false, west, Spot::Main);
            const double south_gap = lead_deg(moon, true, west, Spot::Beam) - lead_deg(moon, true, west, Spot::Main);
            mirrored = mirrored && std::abs(north_gap + south_gap) < 1e-12 &&
                       lead_deg(moon, false, west, Spot::Beam) == lead_deg(moon, true, west, Spot::Main) &&
                       std::abs(lead_deg(moon, false, west, Spot::Reflected) - lead_deg(moon, false, west, Spot::Main) -
                                2.0 * lead_deg(moon, true, west, Spot::Main)) < 1e-12;
            if (moon == 0) widest_gap = std::max(widest_gap, std::abs(north_gap));
        }
        for (int south = 0; south < 2; ++south) {
            farthest_reflected[south] = std::max(farthest_reflected[south], lead_deg(0, south != 0, west, Spot::Reflected) -
                                                                                lead_deg(0, south != 0, west, Spot::Main));
        }
    }
    check(mirrored, "the beam spot's gap is equal and opposite in the two hemispheres");
    std::printf("  info: Io's main and beam spots are up to %.2f deg apart on the equator; the reflected spot up to %.1f (north), "
                "%.1f (south) behind the main one\n", widest_gap, farthest_reflected[0], farthest_reflected[1]);
    // Bonfond 2008: the first two spots are at most about 4 deg apart, the third at most 12 deg
    // behind the main one.
    check(widest_gap > 3.5 && widest_gap < 5.5, "Io's first two spots are at most about 4 deg apart", widest_gap);
    check(farthest_reflected[0] > 10.0 && farthest_reflected[0] < 15.0 && farthest_reflected[1] > 10.0 && farthest_reflected[1] < 15.0,
          "its third spot at most about 12 deg behind the main one", farthest_reflected[0]);
    // Bonfond 2008: a leading spot in the north with Io at System III longitudes 0 to 100 deg (Io
    // south of the torus), in the south at 130 to 300 deg; elsewhere the beam spot trails.
    bool leads_north = true, leads_south = true;
    for (double west = 0.0; west <= 100.0; west += 5.0) {
        leads_north = leads_north && lead_deg(0, false, west, Spot::Beam) < lead_deg(0, false, west, Spot::Main);
    }
    for (double west = 130.0; west <= 280.0; west += 5.0) {
        leads_south = leads_south && lead_deg(0, true, west, Spot::Beam) < lead_deg(0, true, west, Spot::Main);
    }
    check(leads_north, "Io's northern beam spot leads the main one at longitudes 0 to 100 deg");
    check(leads_south, "the southern one at 130 to 280 deg");
    // Bonfond 2008, paragraph 10: on 2007-03-05 the southern leading spot was 3.1 deg ahead of
    // the main one at 09:02 UT, the northern secondary spot 3.4 deg behind at 09:10.
    // Io then, from JPL Horizons (501 from 599, ICRF, km; the scene's Io has no ephemeris for
    // 2007): the light seen at 09:06 UT left 44.26 minutes earlier.
    {
        double t = 0.0;
        check(parse_utc("2007-03-05 08:21:44", &t), "a time parses");
        scene.update(t);
        const Body& planet = scene.bodies[static_cast<size_t>(index)];
        const glm::dvec3 io = glm::transpose(scene.body_to_icrf(planet, t)) * glm::dvec3(360617.0, -198069.0, -88436.0);
        const double west = -std::atan2(io.y, io.x) / kDegToRad;
        const double gap = lead_deg(0, false, west, Spot::Beam) - lead_deg(0, false, west, Spot::Main);
        std::printf("  info: 2007-03-05 09:06 UT: Io at %.1f deg W; the northern beam spot %.2f deg behind the main one "
                    "(Bonfond 2008: gaps of 3.1 and 3.4 deg)\n", wrap_two_pi(west * kDegToRad) / kDegToRad, gap);
        // By the paper's own figure 2 the leading spot at this longitude is the northern one, as here;
        // its example names the southern. Only the size of the gap is compared.
        check(std::abs(gap) > 2.0 && std::abs(gap) < 5.5, "the spots of 2007-03-05 are about as far apart as observed", gap);
        check(parse_utc("2021-06-07 18:00", &t), "a time parses");
        scene.update(t);
        // In the scene: unit vectors in their hemispheres, the beam spots apart from the main ones.
        for (size_t k = 0; k < planet.jovian_secondary_spots.size(); ++k) {
            const glm::dvec3& p = planet.jovian_secondary_spots[k];
            check(std::abs(glm::length(p) - 1.0) < 1e-10 && (k % 2 == 0 ? p.z > 0.5 : p.z < -0.5), "the secondary spots are placed");
        }
        const double apart_km = std::acos(glm::dot(planet.jovian_secondary_spots[0], planet.jovian_footprints[0])) * 71492.0;
        const double behind_km = std::acos(glm::dot(planet.jovian_secondary_spots[6], planet.jovian_footprints[0])) * 71492.0;
        std::printf("  info: 2021-06-07 18:00: Io's northern beam spot %.0f km from the main one, the reflected one %.0f km\n", apart_km, behind_km);
        check(apart_km < 5000.0 && behind_km > 1500.0 && behind_km < 12000.0, "Io's secondary spots are within a few thousand km of the main one", behind_km);
    }
    // Bonfond 2013: Ganymede's two spots are at most 4,300 km apart on the planet.
    {
        double widest_km = 0.0;
        for (double east = 0.0; east < 360.0; east += 2.0) {
            const glm::dvec3 moon = 1070400.0 * glm::dvec3(std::cos(east * kDegToRad), std::sin(east * kDegToRad), 0.0);
            for (int south = 0; south < 2; ++south) {
                glm::dvec3 main_spot, beam_spot;
                check(model.footprint(2, south != 0, moon, main_spot) && model.footprint(2, south != 0, moon, beam_spot, Spot::Beam),
                      "Ganymede's spots are in the table");
                widest_km = std::max(widest_km, std::acos(std::min(glm::dot(main_spot, beam_spot), 1.0)) * (71492.0 + 900.0));
            }
        }
        // Hue's two lead angles differ by up to 14 deg, more than the 10 deg Bonfond expected and
        // twice what Hubble saw in 2007, over part of the longitudes.
        std::printf("  info: Ganymede's main and beam spots are up to %.0f km apart (Bonfond 2013: 4,300 km seen)\n", widest_km);
        check(widest_km > 4300.0 && widest_km < 9000.0, "Ganymede's two spots are thousands of km apart at most", widest_km);
    }

    Scene alone;
    check(load_scene_string("name='x'\n[[bodies]]\nname='Jupiter'\ntemplate='bodies/jupiter.toml'\n",
                            "no_moons", alone, &error, ASTRAXIS_ASSET_DIR), "Jupiter template loads without satellites");
    if (!alone.bodies.empty()) {
        alone.update(start);
        for (const auto& p : alone.bodies[0].jovian_footprints)
            check(p == glm::dvec3(0), "missing satellites leave no stale footprints");
        for (const auto& p : alone.bodies[0].jovian_secondary_spots)
            check(p == glm::dvec3(0), "nor secondary spots");
    }
}

} // namespace

// Comet Shoemaker-Levy 9's impacts: Harrington & Deming 2001 [HD01], Hammel et al. 1995 and 1996 [HAM95],
// Beebe 1996 [BEE96], Sanchez-Lavega et al. 1998 [SL98].
void test_impact_sites()
{
    Scene scene = load_scene_or_die("jupiter.toml");
    const int jupiter = scene.find("Jupiter");
    const Body& planet = scene.bodies[static_cast<size_t>(jupiter)];
    check(planet.impacts.size() == 21, "Jupiter has 21 impact sites", static_cast<double>(planet.impacts.size()));
    if (planet.impacts.empty()) {
        return;
    }

    // Each site is where its fragment's table ends.
    double worst_time = 0.0;
    double worst_place = 0.0;
    for (const ImpactSite& site : planet.impacts) {
        const int f = scene.find(("SL9 " + site.name).c_str());
        check(f > 0, "a site's fragment exists");
        const MotionSource& motion = *scene.bodies[static_cast<size_t>(f)].motion;
        double lo = site.t_tdb - 30.0;
        double hi = site.t_tdb + 30.0;
        check(motion.valid_at(lo) && !motion.valid_at(hi), "the fragment ends at its site's time");
        for (int i = 0; i < 40; ++i) {
            const double mid = 0.5 * (lo + hi);
            (motion.valid_at(mid) ? lo : hi) = mid;
        }
        worst_time = std::max(worst_time, std::abs(lo - site.t_tdb));
        const glm::dvec3 r = scene.icrf_state_at(f, lo).position - scene.icrf_state_at(jupiter, lo).position;
        const glm::dvec3 fixed = glm::transpose(scene.body_to_icrf(planet, lo)) * r;
        const double lat = glm::degrees(std::asin(fixed.z / glm::length(fixed)));
        const double west = glm::degrees(-std::atan2(fixed.y, fixed.x));
        worst_place = std::max({worst_place, std::abs(lat - site.lat_deg), std::abs(std::remainder(west - site.west_lon_deg, 360.0))});
    }
    check(worst_time < 1.0, "the sites' times match the fragments' tables", worst_time);
    check(worst_place < 0.02, "the sites' places match the fragments' tables", worst_place);

    auto site = [&](const char* name) -> const ImpactSite& {
        for (const ImpactSite& s : planet.impacts) {
            if (s.name == name) {
                return s;
            }
        }
        return planet.impacts[0];
    };
    // [HAM95] Table 2, and its site near G's ring center: 43.65 S, 25.7 W.
    check(site("G").klass == ImpactClass::C1 && site("K").klass == ImpactClass::C1 && site("L").klass == ImpactClass::C1,
          "G, K and L are class 1");
    check(site("B").klass == ImpactClass::C3 && site("F").klass == ImpactClass::C4 && site("Q1").klass == ImpactClass::C2b,
          "other site classes");
    check(std::abs(site("G").lat_deg + 43.65) < 0.3 && std::abs(site("G").west_lon_deg - 25.7) < 0.6,
          "G against Hammel's ring center", site("G").west_lon_deg);
    // The solar system scene has the same Jupiter.
    Scene solar = load_scene_or_die("solar_system.toml");
    check(solar.bodies[static_cast<size_t>(solar.find("Jupiter"))].impacts.size() == 21, "the solar system's Jupiter has them too");

    // The plume [HD01] against the observations [HAM95]: 3000 km within 6-8 minutes or so, settled in 15 to 20.
    check(std::abs(impact_plume_vanguard_share() - 0.22) < 0.01, "the vanguard's mass share, Harrington & Deming", impact_plume_vanguard_share());
    const double peak_s = kPlumeVMaxKmS / kPlumeGravityKmS2;
    check(std::abs(impact_plume_top_km(peak_s) - 3000.0) < 20.0, "plume top 3000 km, Hammel", impact_plume_top_km(peak_s));
    check(peak_s > 6.0 * 60.0 && peak_s < 9.0 * 60.0, "plume at its top after 6-9 minutes", peak_s / 60.0);
    check(impact_plume_end_s() > 15.0 * 60.0 && impact_plume_end_s() < 20.0 * 60.0, "plume down within 15-20 minutes", impact_plume_end_s() / 60.0);
    // Its density integrates to the whole mass, at any time, while nothing has come down.
    const double age = 200.0;
    double mass = 0.0;
    const double step = 160.0;
    for (double x = -3200.0; x < 3200.0; x += step) {
        for (double y = -3200.0; y < 3200.0; y += step) {
            for (double z = 0.5 * step; z < 3200.0; z += step) {
                mass += impact_plume_density(glm::dvec3(x, y, z), age) * step * step * step;
            }
        }
    }
    check(std::abs(mass - 1.0) < 0.1, "the plume's density holds its mass", mass);
    // It leans southeast (Hammel: the fragments came from the south, 16 deg toward the east).
    check(impact_plume_density(glm::dvec3(500.0, -700.0, 1500.0), age) > 5.0 * impact_plume_density(glm::dvec3(-500.0, 700.0, 1500.0), age),
          "the plume leans southeast");
    check(impact_plume_density(glm::dvec3(0.0, 0.0, 100.0), -1.0) == 0.0, "no plume before the impact");

    // The fireball, with the swarm's top, clears the limb about when the second precursor starts (Chodas &
    // Yeomans Table 1: A +49 s, G +59-70 s, W +8 s), for material 750, 444 and 137 km up.
    auto clears = [](double km) {
        double t = 0.0;
        while (impact_plume_top_km(t) < km) {
            t += 0.5;
        }
        return t;
    };
    check(std::abs(clears(750.0) - 49.0) < 25.0, "A's fireball clears the limb", clears(750.0));
    check(std::abs(clears(444.0) - 65.0) < 30.0, "G's fireball clears the limb", clears(444.0));
    check(std::abs(clears(137.0) - 8.0) < 10.0, "W's fireball clears the limb", clears(137.0));
    // The flashes. [MAR95] and [CHA95]: up in 2 s; a plateau of about 10 s, over by 35-50 s, for the large ones;
    // K's width at half height 27 s; N's flash 6 s and a faint 9 s, at 5-10 s only 5-10% of K (peaking at 40% of K).
    const double large_share = impact_class_params(ImpactClass::C1).fireball_share;
    const double small_share = impact_class_params(ImpactClass::C3).fireball_share;
    check(impact_flash(0.0, large_share) == 0.0 && impact_flash(2.0, large_share) > 0.95, "the flash's 2 s rise", impact_flash(2.0, large_share));
    check(impact_flash(12.0, large_share) > 0.8, "the plateau of about 10 s", impact_flash(12.0, large_share));
    check(impact_flash(52.0, large_share) < 0.01, "over within K's 52 s", impact_flash(52.0, large_share));
    double half = 2.0;
    while (impact_flash(half, large_share) > 0.5) {
        half += 0.1;
    }
    check(std::abs(half - 1.0 - 27.0) < 3.0, "K's width at half height, Chapman", half - 1.0);
    const double n_at_7 = 0.4 * impact_flash(7.5, small_share) / impact_flash(7.5, large_share);
    check(n_at_7 > 0.05 && n_at_7 < 0.12, "N's faint glow against K, Chapman", n_at_7);
    check(impact_flash(16.0, small_share) < 0.1, "N is over in 15 s", impact_flash(16.0, small_share));
    // The peaks: the seven measured, relative to L 2.4 : 2.1 : 1.6 : 1.0 (L, G, H, Q1); W was 1.5% in green; the
    // seven differ by a factor under 6; K, of the largest class, is between Q1 and H.
    check(std::abs(site("L").flash_share - 0.06) < 1e-3, "L at 6% of Jupiter, Martin", site("L").flash_share);
    check(std::abs(site("L").flash_share / site("Q1").flash_share - 2.4) < 0.1, "L against Q1, Martin", site("L").flash_share / site("Q1").flash_share);
    check(std::abs(site("W").flash_share - 0.015) < 2e-3, "W in green, Chapman", site("W").flash_share);
    check(site("L").flash_share / site("N").flash_share < 6.0, "the bolides differ by under 6 times, Chapman");
    check(site("K").flash_share > site("Q1").flash_share && site("K").flash_share < site("H").flash_share, "K between Q1 and H, Chapman");
    check(impact_flash_share(site("C")) > 0.0 && impact_flash_share(site("T")) < 0.01, "the flashes nobody measured");
    // The color temperatures [MAR95, CHA95, CAR95].
    check(impact_flash_temperature_k(1.0) > 10000.0, "over 10,000 K at the peak, Martin", impact_flash_temperature_k(1.0));
    check(std::abs(impact_flash_temperature_k(5.0) - 7600.0) < 1.0, "7600 K in G's first seconds, Chapman");
    check(std::abs(impact_flash_temperature_k(21.0) - 2200.0) < 1.0, "2200 K, Carlson");
    check(std::abs(impact_flash_temperature_k(65.0) - 450.0) < 1.0, "450 K, Martin quoting Carlson");
    // The fireball [CAR95]: 40 km across about 21 s in, 1 km/s or less at first, 2-3 km/s after 20-30 s, and low.
    check(std::abs(2.0 * impact_fireball_sigma_km(21.0) - 40.0) < 5.0, "G's fireball is 40 km across", 2.0 * impact_fireball_sigma_km(21.0));
    const double rate = impact_fireball_height_km(25.0) - impact_fireball_height_km(24.0);
    check(rate > 2.0 && rate < 3.0, "it rises at 2-3 km/s, Carlson", rate);
    check(impact_fireball_height_km(5.0) - impact_fireball_height_km(4.0) <= 1.0, "1 km/s or less at first");
    check(impact_fireball_height_km(60.0) < 444.0, "it stays hidden from the Earth, 444 km for G, Hammel", impact_fireball_height_km(60.0));

    // The plume's brightness against the measurements. G from the side, 640 s in: [HAM95] I/F 0.090 +- 0.01. All of it,
    // [CHA95]: under a quarter of Io for K at 9 minutes, here by the optically thin limit's order.
    const glm::dvec3 view = glm::normalize(glm::dvec3(0.8, 0.6, -0.15));
    const glm::dvec3 across = glm::normalize(glm::cross(view, glm::dvec3(0.0, 0.0, 1.0)));
    const glm::dvec3 upward = glm::cross(across, view);
    const double g2 = kPlumeAsymmetry * kPlumeAsymmetry;
    const double back = (1.0 - g2) / std::pow(1.0 + g2 + 2.0 * kPlumeAsymmetry, 1.5);
    const auto brightness = [&](double age_s, double* total) {
        const double cell = 120.0;
        double best = 0.0;
        double sum = 0.0;
        for (double u = -7200.0; u < 7200.0; u += cell) {
            for (double w = -7200.0; w < 7200.0; w += cell) {
                double column = 0.0;
                for (double s = -7200.0; s < 7200.0; s += cell) {
                    column += impact_plume_density(u * across + w * upward + s * view, age_s) * cell;
                }
                const double i_f = 0.25 * back * kPlumeAlbedo * (1.0 - std::exp(-kPlumeExtinctionKm2 * column));
                best = std::max(best, i_f);
                sum += i_f * cell * cell;
            }
        }
        *total = sum;
        return best;
    };
    double total = 0.0;
    const double g_if = brightness(640.0, &total);
    check(std::abs(g_if - 0.09) < 0.02, "G's plume I/F, Hammel", g_if);
    const double quarter_io = 0.25 * 0.6 * kPi * 1821.0 * 1821.0;
    check(total < 1.2 * quarter_io, "the plume under a quarter of Io, Chapman", total / quarter_io);

    // Rings: 454 m/s from 586 km; at 1.5 h, 3,037 km.
    check(std::abs(impact_ring_radius_km(5400.0, false) - 3037.6) < 0.1, "the main ring", impact_ring_radius_km(5400.0, false));
    const ImpactClassParams& large = impact_class_params(ImpactClass::C1);
    auto along_axis = [&](double r, double age_s) {
        const double a = kPlumeAzimuthDeg * kDegToRad;
        return impact_scar_shape(large, r * std::sin(a), r * std::cos(a), age_s, 0.0);
    };
    // The ring is a dark line, opposite the crescent; and gone in 10 hours.
    const double a = (kPlumeAzimuthDeg + 180.0) * kDegToRad;
    const double r_ring = impact_ring_radius_km(5400.0, false);
    const double on_ring = impact_scar_shape(large, r_ring * std::sin(a), r_ring * std::cos(a), 5400.0, 0.0);
    const double off_ring = impact_scar_shape(large, (r_ring + 500.0) * std::sin(a), (r_ring + 500.0) * std::cos(a), 5400.0, 0.0);
    check(on_ring > off_ring + 0.03, "a ring shows at its radius", on_ring - off_ring);
    const double r_late = impact_ring_radius_km(11.0 * 3600.0, false);
    check(std::abs(impact_scar_shape(large, r_late * std::sin(a), r_late * std::cos(a), 11.0 * 3600.0, 0.0) -
                   impact_scar_shape(large, (r_late + 500.0) * std::sin(a), (r_late + 500.0) * std::cos(a), 11.0 * 3600.0, 0.0)) < 0.01,
          "and it is gone after 10 hours");

    // West 1996 Fig. 3: the crescent's optical depth at 890 nm is about 0.8 against the core's 3.4, at about 2 h.
    // G's crescent: outer edge 13,000 km to the southeast, and much less on the other side. Nothing
    // comes down in the first minutes.
    const double t_look = 5400.0;
    double crest = 0.0;
    double crest_r = 0.0;
    for (double r = 5000.0; r < 18000.0; r += 100.0) {
        if (along_axis(r, t_look) > crest) {
            crest = along_axis(r, t_look);
            crest_r = r;
        }
    }
    double edge_r = crest_r;
    while (along_axis(edge_r, t_look) > 0.5 * crest && edge_r < 30000.0) {
        edge_r += 50.0;
    }
    check(std::abs(edge_r - 13000.0) < 1500.0, "G's crescent reaches 13,000 km, Hammel", edge_r);
    check(along_axis(crest_r, t_look) > 3.0 * along_axis(-crest_r, t_look), "the crescent is to the southeast");
    const double crescent_depth = impact_core_depth(t_look) * crest;
    check(crescent_depth > 0.5 && crescent_depth < 1.1, "the crescent's optical depth, West Fig. 3", crescent_depth);
    const double core_depth = impact_core_depth(t_look) * along_axis(-300.0, t_look);
    check(std::abs(core_depth - 3.4) < 0.8, "the core's optical depth, Sanchez-Lavega Table IV", core_depth);
    // What it does to the clouds in visible light, by tools/sl9/debris_optics.py from West Fig. 3: G's core at 0.50 of
    // its surroundings (he measures 0.64 at 955 nm and 0.56 at 275 nm), the crescent at 0.74 (0.84 and 0.65). Not black,
    // slightly red; thin debris shows more toward the limb, and thick debris never goes below its own albedo.
    const double mu_west = 0.77;
    const glm::dvec3 core_factor = impact_debris_factor(3.4, mu_west, mu_west);
    check(std::abs(core_factor.g - 0.50) < 0.03, "a core's dimming, West Fig. 3", core_factor.g);
    check(core_factor.r > core_factor.b, "the debris is slightly red");
    const double thin = impact_debris_factor(0.8, mu_west, mu_west).g;
    check(std::abs(thin - 0.74) < 0.03, "a crescent's dimming, West Fig. 3", thin);
    check(impact_debris_factor(0.8, 0.2, mu_west).g < thin - 0.05, "thin debris shows more toward the limb");
    check(impact_debris_factor(100.0, mu_west, mu_west).g >= kDebrisRatioGreen - 1e-9, "thick debris shows its own albedo");
    check(along_axis(crest_r, 200.0) < 0.01 * crest, "no ejecta in the first minutes");
    // The small and the undetected.
    const ImpactClassParams& small = impact_class_params(ImpactClass::C3);
    check(impact_scar_shape(small, 5000.0, 0.0, t_look, 0.0) < 0.01 && impact_scar_shape(small, 0.0, 0.0, t_look, 0.0) > 0.9,
          "a class 3 site is only a dark patch under 3,000 km");
    check(impact_class_params(ImpactClass::C4).weight == 0.0, "a class 4 site leaves nothing");
    check(impact_site_fade(small, 90.0 * 3600.0) < 0.5, "B's site was gone before N arrived, Beebe", impact_site_fade(small, 90.0 * 3600.0));

    // The optical depths at 890 nm, [SL98] Table IV: 3.4 at first, 2.2 on Aug 28, 0.6 in December.
    check(std::abs(impact_core_depth(2.0 * kSecondsPerDay) - 3.4) < 0.01, "a core's depth at first", impact_core_depth(2.0 * kSecondsPerDay));
    check(std::abs(impact_core_depth(41.0 * kSecondsPerDay) - 2.2) < 0.01, "on Aug 28", impact_core_depth(41.0 * kSecondsPerDay));
    const double dec17 = (2449703.8 - 2449550.0); // days after the first impact (July 16.8)
    check(std::abs(impact_band_depth(dec17) - 0.6) < 0.05, "the band in December 1994", impact_band_depth(dec17));

    // The places: dark at G after 1.5 h, not before the impact; in the band, later, at the latitude of the
    // sites; none near the equator.
    const ImpactSite& g = site("G");
    check(impact_scar_depth(planet.impacts, g.lat_deg, g.west_lon_deg, g.t_tdb + 5400.0) > 2.0, "G is dark",
          impact_scar_depth(planet.impacts, g.lat_deg, g.west_lon_deg, g.t_tdb + 5400.0));
    check(impact_scar_depth(planet.impacts, g.lat_deg, g.west_lon_deg, planet.impacts[0].t_tdb - 3600.0) == 0.0, "nothing before");
    const double t_band = g.t_tdb + 150.0 * kSecondsPerDay;
    check(impact_scar_depth(planet.impacts, -43.0, 120.0, t_band) > 0.4, "in December a band at 43 S",
          impact_scar_depth(planet.impacts, -43.0, 120.0, t_band));
    check(impact_scar_depth(planet.impacts, 0.0, 120.0, t_band) < 0.01, "not at the equator");
    check(impact_scar_depth(planet.impacts, -43.0, 120.0, g.t_tdb + 4000.0 * kSecondsPerDay) < 1e-6, "gone long after");
    // East of a site's latitude, the wind [BEE96]: the debris moves east to the north, west to the south.
    check(impact_wind_km_s(2000.0) > 0.0 && impact_wind_km_s(-2000.0) < 0.0, "the shear's direction");
    check(impact_wind_km_s(1e5) == 0.020 && impact_wind_km_s(-1e5) == -0.010, "the shear's limits, Table III");
}

// The fragments' tails against Sekanina 1996, Sec. 3.4 and 4.4 and Fig. 5.
void test_bound_dust_tail()
{
    Scene scene = load_scene_or_die("jupiter.toml");
    const int k = scene.find("SL9 K");
    const Body& fragment = scene.bodies[static_cast<size_t>(k)];
    check(fragment.comet.bound_tail && fragment.comet.dust_only, "the fragments have tails released once");
    BoundDustTail tail;
    const auto started = std::chrono::steady_clock::now();
    tail.build(scene, k);
    const double build_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started).count();
    std::printf("info: a fragment's tail takes %.1f ms and %zu steps\n", build_ms, tail.row_count());
    check(tail.row_count() > 100 && tail.row_count() < 3000, "the tail's table is of a sensible size", static_cast<double>(tail.row_count()));

    // On the sky, from the Earth (Horizons, from Jupiter), with north the ICRF pole; east toward increasing RA.
    const auto on_sky = [](const glm::dvec3& d, const glm::dvec3& earth) {
        const glm::dvec3 to_earth = glm::normalize(earth);
        return d - to_earth * glm::dot(d, to_earth);
    };
    const auto position_angle = [](const glm::dvec3& d, const glm::dvec3& earth) {
        const glm::dvec3 to_earth = glm::normalize(earth);
        const glm::dvec3 north = glm::normalize(glm::dvec3(0.0, 0.0, 1.0) - to_earth * to_earth.z);
        const glm::dvec3 east = glm::cross(north, -to_earth);
        return std::fmod(glm::degrees(std::atan2(glm::dot(d, east), glm::dot(d, north))) + 360.0, 360.0);
    };
    const double jan94 = tdb_from_jd_tdb(2449377.5); // 1994-01-25
    const double jul94 = tdb_from_jd_tdb(2449547.5); // 1994-07-14
    const glm::dvec3 earth_jan(601460748.0, 520171339.9, 206808453.3);
    const glm::dvec3 earth_jul(622908828.5, 406956663.6, 159976989.5);
    const auto earth_from_k = [&](const glm::dvec3& earth, double t) {
        return earth - (scene.icrf_state_at(k, t).position - scene.icrf_state_at(0, t).position);
    };
    const glm::dvec3 none(0.0);
    const int first = 0;
    const int last = BoundDustTail::kReleases - 1;

    // [SEK] Sec. 4.4: in July 1994 grains 15,000 km from the fragment have beta 5e-5 if released in early July 1992,
    // 2e-4 if in late 1992.
    const glm::dvec3 from_earth_jul = earth_from_k(earth_jul, jul94);
    const double beta_early = 15000.0 / glm::length(on_sky(tail.offset(first, 1.0, none, jul94), from_earth_jul));
    const double beta_late = 15000.0 / glm::length(on_sky(tail.offset(last, 1.0, none, jul94), from_earth_jul));
    check(beta_early > 2.5e-5 && beta_early < 7.5e-5, "beta of the early grains, Sekanina", beta_early);
    check(beta_late > 1.3e-4 && beta_late < 3.0e-4, "beta of the late grains, Sekanina", beta_late);

    // [SEK] Fig. 5: the tails at position angles of about 270-295 deg in January 1994, 245-255 in mid-July. Sec. 3.4: 30 deg
    // from the train (itself at about 243 and 240).
    const glm::dvec3 from_earth_jan = earth_from_k(earth_jan, jan94);
    const double pa_jan = position_angle(tail.offset(last / 2, 1.0, none, jan94), from_earth_jan);
    const double pa_jul = position_angle(tail.offset(last / 2, 1.0, none, jul94), from_earth_jul);
    check(pa_jan > 265.0 && pa_jan < 295.0, "the tails' position angle in January 1994, Sekanina Fig. 5", pa_jan);
    check(pa_jul > 240.0 && pa_jul < 262.0, "the tails' position angle in July 1994, Sekanina Fig. 5", pa_jul);
    const int a_body = scene.find("SL9 A");
    const int w_body = scene.find("SL9 W");
    const glm::dvec3 train = scene.icrf_state_at(w_body, jan94).position - scene.icrf_state_at(a_body, jan94).position;
    const double pa_train = position_angle(train, from_earth_jan);
    check(std::abs(pa_train - 243.0) < 8.0, "the train's position angle in January 1994, Sekanina Fig. 5", pa_train);
    check(std::abs(std::remainder(pa_jan - pa_train, 360.0) - 30.0) < 12.0, "30 deg from the train, Sekanina", pa_jan - pa_train);

    // [SEK] Sec. 3.4: over 100,000 km long in January 1994 for the major ones, and 8,000 to 25,000 km wide from speeds
    // of 0.1 to 0.4 m/s (normal to the orbit).
    const double longest = glm::length(on_sky(tail.offset(last / 2, BoundDustTail::kBetaMax, none, jan94), from_earth_jan));
    check(longest > 100000.0, "the tails are over 100,000 km long, Sekanina", longest);
    const glm::dvec3 normal = glm::normalize(glm::cross(scene.icrf_state_at(k, jan94).position - scene.icrf_state_at(0, jan94).position,
                                                        scene.icrf_state_at(k, jan94).velocity - scene.icrf_state_at(0, jan94).velocity));
    const double across = 2.0 * glm::length(tail.offset(last / 2, 0.0, normal * BoundDustTail::kSpeedMaxKmS, jan94));
    check(across > 8000.0 && across < 30000.0, "their width, Sekanina", across);

    // The grains: nothing before the release or after the impact; the same for the same time; none far off.
    tail.update(scene, k, jan94);
    const size_t count = tail.grains().size();
    check(count > 150, "most grains are there in January 1994", static_cast<double>(count));
    double weights = 0.0;
    for (const DustGrain& g : tail.grains()) {
        weights += g.weight;
        check(g.age_s > 0.0, "a grain's age");
    }
    check(weights > 0.5 && weights <= 1.0 + 1e-9, "the weights add up", weights);
    tail.update(scene, k, tdb_from_jd_tdb(2448800.0));
    check(tail.grains().empty(), "no grains before the breakup");
    tail.update(scene, k, tdb_from_jd_tdb(2449600.0));
    check(tail.grains().empty(), "no grains after the impact");
}

void run_appearance_tests()
{
    test_impact_sites();
    test_relief();
    test_bound_dust_tail();
    test_comet_tails();
    test_dust_tail();
    test_atmospheres();
    test_plumes();
    test_night_and_aurora();
    test_jovian_aurora();
    test_saturn_aurora();
    test_aurora_profiles();
    test_ice_giant_aurora();
    test_dark_spots();
    test_magnetic_field();
    test_flux_tubes();
    test_magnetospheres();
    test_magnetic_maps();
    test_ganymede_aurora();
    test_europa_aurora();
    test_io_callisto_aurora();
    test_plasma_sheet();
    test_europa_pattern();
    test_io_eclipse_levels();
    test_mars_aurora();
    test_jovian_footprints();
    test_refraction();
    test_corona();
}
