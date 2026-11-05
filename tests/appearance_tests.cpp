// Tests for body appearance: comet tails, atmospheres, plumes, night lights and aurora.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "scene/aurora.hpp"
#include "scene/aurora_model.hpp"
#include "scene/satellite_aurora.hpp"
#include "scene/saturn_aurora.hpp"
#include "scene/mars_aurora.hpp"
#include "scene/magnetic_field.hpp"
#include "scene/comet.hpp"
#include "scene/corona.hpp"
#include "scene/scene_loader.hpp"
#include "scene/simulation.hpp"

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
        check(glm::dot(direction, glm::normalize(sun - body.world_position)) < -0.7,
              "visible aurora tour views the nightside");
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
    check(events == 1, "Mars has one visible aurora event");
    // Independent dense tangent-ray quadrature checks the orbit-limb scaling.
    double maximum = 0.0;
    for (int h = 20; h <= 140; ++h) {
        const double b = 3396.19 + h;
        const double end = std::sqrt(std::pow(3396.19 + 140.0, 2) - b*b);
        double integral = 0.0;
        for (int j = 0; j < 8000; ++j) {
            const double x = end * (j + 0.5) / 8000.0;
            const double altitude = std::sqrt(b*b + x*x) - 3396.19;
            integral += std::exp(-0.5 * std::pow((altitude - 60.0) / 15.0, 2)) * 2.0 * end / 8000.0;
        }
        maximum = std::max(maximum, integral);
    }
    check(std::abs(maximum / MarsAurora::limb_integral() - 1.0) < 0.001, "Mars reference limb normalization", maximum);
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
    check(kCallistoOxygenLines.size() == 2 && kCallistoOxygenLines[0].rayleigh == 6.1 &&
          kCallistoOxygenLines[1].rayleigh == 3.3, "Callisto uses detected red lines, not upper limits");
    check(kIoSodiumLines[0].rayleigh + kIoSodiumLines[1].rayleigh == 17170.0,
          "Io sodium brightness retains the measured eclipse average");
    check(io_visible_aurora().sodium_column_rgb.g > 0.0 && callisto_visible_aurora().column_rgb.g == 0.0,
          "Io sodium and Callisto oxygen have distinct spectra");
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
        {Model::Ganymede, true, Pattern::BodyMap, Vertical::UniformColumn, Lighting::Always},
        {Model::Uranus, false, Pattern::BodyMapSeasons, Vertical::Table, Lighting::Always},
        {Model::Neptune, false, Pattern::BodyMapDusk, Vertical::ChapmanPair, Lighting::Always},
        {Model::Europa, false, Pattern::EuropaTrailing, Vertical::Exponential, Lighting::NightOrEclipsed},
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
                              aurora_disk_factor(1560.8, 50.0, 300.0);
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
    const auto find_event = [&](const char* text) {
        return std::find_if(scene.events.begin(), scene.events.end(),
                            [&](const SceneEvent& e) { return e.name.find(text) != std::string::npos; });
    };
    const auto uranus_event = find_event("Uranus's infrared");
    const auto neptune_event = find_event("Neptune's first infrared");
    check(uranus_event != scene.events.end() && neptune_event != scene.events.end(), "the ice giants' events");
    if (uranus_event == scene.events.end() || neptune_event == scene.events.end()) return;
    for (const auto& e : {*uranus_event, *neptune_event}) {
        check(e.aurora_infrared && !e.aurora_ultraviolet && !e.aurora_visible, "the ice giants' events are infrared");
    }
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
    constexpr double R = 1560.8, top = R + 300.0;
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
            column += std::exp(-h / 50.0) / (50.0 * (1.0 - std::exp(-6.0))) * (z1 - z0) / depth;
        }
        flux += column * top * top / (R * R * radial);
    }
    check(std::abs(flux - aurora_disk_factor(R, 50.0, 300.0)) < 0.003,
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
        check(along > 0.0 && miss + R + 300.0 + along * 695700.0 / glm::length(scene.sun_position() - moon.world_position) < parent.polar_radius_km,
              "Europa aurora tour lies inside Jupiter shadow", miss);
        check(event.aurora_visible && !event.aurora_ultraviolet, "Europa tour uses visible spectrum");
        scene.update(event.t_tdb + 3.0 * 3600.0);
        const auto day_light = glm::normalize(scene.sun_position() - moon.world_position);
        const auto day_oc = parent.world_position - moon.world_position;
        check(glm::length(glm::cross(day_oc, day_light)) > parent.equatorial_radius_km + R + 300.0,
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
            check(std::abs(pixels[4 * k] * map.range / 255.0 - b) < 0.4,
                  "GPU UV texture quantization is below 0.4 Rayleigh");
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
    // The aurora lies at the edge of the closed lines (Duling 2022), but Jupiter's plasma moves
    // that edge: the leading side's crescents are far nearer the equator (Marzok 2022). The
    // dipole's outer shell can only match the trailing side's, and lies poleward of the leading's.
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
        for (int sign : {-1, 1}) {
            double foot_latitude = 0.0;
            for (const auto& line : ganymede.magnetic_field.lines) {
                for (const glm::dvec4& p : {line.front(), line.back()}) {
                    foot_latitude = std::max(foot_latitude, std::asin(p.z / glm::length(glm::dvec3(p))) / kDegToRad * sign);
                }
            }
            const double trailing = crescent(90.0, sign);
            const double leading = crescent(-90.0, sign);
            std::printf("  info: Ganymede %s: the lines' outermost feet at %.1f deg, the aurora's crescents at %.1f (trailing), %.1f (leading)\n",
                        sign > 0 ? "north" : "south", foot_latitude, trailing, leading);
            check(std::abs(foot_latitude - trailing) < 15.0, "the outer lines land near the trailing crescent", foot_latitude - trailing);
            check(foot_latitude > leading + 10.0, "and poleward of the leading one", foot_latitude - leading);
        }
    }

    // The scene's lines: within their limits, cheap to draw and to make.
    const struct { const Body* body; double max_radius; } limits[] = {
        {&uranus, 12.0}, {&neptune, 12.0}, {&saturn, 12.0}, {&earth, 8.0}, {&jupiter, 10.0}, {&mercury, 1.7}, {&ganymede, 2.5},
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
        check(field.lines.size() == 24, "three shells of eight lines", static_cast<double>(field.lines.size()));
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
    for (const Scene* scene : {&solar, &saturn_scene, &earth_scene, &jupiter_scene, &mercury_scene}) {
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
    check(shown(solar) == 2 && shown(saturn_scene) == 1 && shown(mercury_scene) == 1 && shown(jupiter_scene) == 1,
          "the ice giants, Saturn, Mercury and Ganymede each have their event");

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
    const auto sub_solar_east = [](Scene& scene, const Body& body, double t) {
        scene.update(t);
        const glm::dvec3 to_sun = glm::transpose(body.orientation) * (scene.sun_position() - body.world_position);
        return std::atan2(to_sun.y, to_sun.x) / kDegToRad;
    };
    // Uranus turns the other way in longitude: the Sun comes to a feature from its west.
    const struct { Scene* scene; const Body* body; double opening_east; FieldMapView view; } expected[] = {
        {&mars_scene, &mars, 215.0, FieldMapView::Surface},
        {&jupiter_scene, &jupiter, 315.5, FieldMapView::Surface},
        {&jupiter_scene, &jupiter, 316.5, FieldMapView::Deep},
        {&earth_scene, &moon, 270.0, FieldMapView::Surface},
        {&earth_scene, &earth, 325.5, FieldMapView::Surface},
        {&earth_scene, &earth, 50.0, FieldMapView::Crust},
        {&giants, &uranus, 155.5, FieldMapView::Surface},
        {&giants, &neptune, 117.5, FieldMapView::Surface},
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
        check(found == 1, "the body has its map event", static_cast<double>(found));
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
    Scene alone;
    check(load_scene_string("name='x'\n[[bodies]]\nname='Jupiter'\ntemplate='bodies/jupiter.toml'\n",
                            "no_moons", alone, &error, ASTRAXIS_ASSET_DIR), "Jupiter template loads without satellites");
    if (!alone.bodies.empty()) {
        alone.update(start);
        for (const auto& p : alone.bodies[0].jovian_footprints)
            check(p == glm::dvec3(0), "missing satellites leave no stale footprints");
    }
}

} // namespace

void run_appearance_tests()
{
    test_comet_tails();
    test_dust_tail();
    test_atmospheres();
    test_plumes();
    test_night_and_aurora();
    test_jovian_aurora();
    test_saturn_aurora();
    test_aurora_profiles();
    test_ice_giant_aurora();
    test_magnetic_field();
    test_flux_tubes();
    test_magnetic_maps();
    test_ganymede_aurora();
    test_europa_aurora();
    test_io_callisto_aurora();
    test_mars_aurora();
    test_jovian_footprints();
    test_refraction();
    test_corona();
}
