// Tests for body appearance: comet tails, atmospheres, plumes, night lights and aurora.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "scene/aurora.hpp"
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
    check(earth.aurora != nullptr && earth.night_lights.gain > 0.0 && !earth.night_lights.texture.empty(),
          "Earth has night lights and an aurora");
    for (const char* file : {"didymos.toml", "earth_moon.toml", "isee3.toml", "mercury.toml", "parker.toml",
                             "rosetta.toml", "space_telescopes.toml", "venus.toml"}) {
        const Scene other = load_scene_or_die(file);
        const Body& e = other.bodies[static_cast<size_t>(other.find("Earth"))];
        check(e.aurora != nullptr && e.night_lights.texture == earth.night_lights.texture &&
                  e.night_lights.color == earth.night_lights.color && e.night_lights.gain == earth.night_lights.gain,
              (std::string(file) + ": Earth's night lights and aurora as in solar_system.toml").c_str());
    }
    for (const char* name : {"Mars", "Venus", "Moon"}) {
        const Body& b = scene.bodies[static_cast<size_t>(scene.find(name))];
        check(!b.aurora && b.night_lights.texture.empty(), name);
    }

    // The Kp record: 1932 on, 3-hour UT intervals.
    const KpRecord& kp = *earth.aurora;
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

} // namespace

void run_appearance_tests()
{
    test_comet_tails();
    test_dust_tail();
    test_atmospheres();
    test_plumes();
    test_night_and_aurora();
    test_refraction();
    test_corona();
}
