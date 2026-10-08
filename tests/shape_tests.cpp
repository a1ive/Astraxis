// Tests for the shape models of small bodies.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

using namespace astraxis;

namespace {

// Signed volume of a closed mesh (positive: counter-clockwise seen from outside).
double signed_volume(const ShapeModel& shape)
{
    double volume = 0.0;
    for (size_t k = 0; k < shape.indices.size(); k += 3) {
        const glm::dvec3 a(shape.positions[shape.indices[k]]);
        const glm::dvec3 b(shape.positions[shape.indices[k + 1]]);
        const glm::dvec3 c(shape.positions[shape.indices[k + 2]]);
        volume += glm::dot(a, glm::cross(b, c)) / 6.0;
    }
    return volume;
}

// Largest difference (turns) between a mesh's map u and the east longitude of its
// vertex directions (seam vertices may be 0 or 1; poles skipped).
double map_u_error(const ShapeModel& shape)
{
    double worst = 0.0;
    for (size_t k = 0; k < shape.positions.size() && k < shape.map_u.size(); ++k) {
        const glm::vec3& p = shape.positions[k];
        if (std::hypot(p.x, p.y) < 1e-3 * glm::length(p)) {
            continue;
        }
        const double d = std::atan2(p.y, p.x) / (2.0 * kPi) - shape.map_u[k];
        worst = std::max(worst, std::abs(d - std::round(d)));
    }
    return worst;
}

// Inertia tensor per unit mass of a closed mesh of uniform density, about its centroid
// (summed over the tetrahedra from the origin to each triangle).
glm::dmat3 uniform_inertia(const ShapeModel& shape)
{
    double volume = 0.0;
    glm::dvec3 first(0.0);
    glm::dmat3 second(0.0); // integral of r r^T
    for (size_t k = 0; k < shape.indices.size(); k += 3) {
        const glm::dvec3 a(shape.positions[shape.indices[k]]);
        const glm::dvec3 b(shape.positions[shape.indices[k + 1]]);
        const glm::dvec3 c(shape.positions[shape.indices[k + 2]]);
        const double v = glm::dot(a, glm::cross(b, c)) / 6.0;
        const glm::dvec3 sum = a + b + c;
        volume += v;
        first += v * sum / 4.0;
        second += v / 20.0 *
                  (glm::outerProduct(a, a) + glm::outerProduct(b, b) + glm::outerProduct(c, c) +
                   glm::outerProduct(sum, sum));
    }
    const glm::dvec3 centroid = first / volume;
    const glm::dmat3 s = second / volume - glm::outerProduct(centroid, centroid);
    return (s[0][0] + s[1][1] + s[2][2]) * glm::dmat3(1.0) - s;
}

// Eigenvectors (columns, by ascending eigenvalue) of a symmetric matrix, cyclic Jacobi.
glm::dmat3 symmetric_eigen(glm::dmat3 m, glm::dvec3& values)
{
    glm::dmat3 vectors(1.0);
    for (int sweep = 0; sweep < 50; ++sweep) {
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(m[q][p]) <= 1e-15 * (std::abs(m[p][p]) + std::abs(m[q][q]))) {
                    continue;
                }
                // Rotation in the p-q plane that zeroes m[p][q].
                const double angle = 0.5 * std::atan2(2.0 * m[q][p], m[q][q] - m[p][p]);
                glm::dmat3 r(1.0);
                r[p][p] = r[q][q] = std::cos(angle);
                r[p][q] = -std::sin(angle);
                r[q][p] = std::sin(angle);
                m = glm::transpose(r) * m * r;
                vectors = vectors * r;
            }
        }
    }
    int order[3] = {0, 1, 2};
    std::sort(order, order + 3, [&](int a, int b) { return m[a][a] < m[b][b]; });
    glm::dmat3 sorted(1.0);
    for (int k = 0; k < 3; ++k) {
        sorted[k] = vectors[order[k]];
        values[k] = m[order[k]][order[k]];
    }
    return sorted;
}

// Arrokoth's shape model and pole against independent numbers in Porter et al. 2024,
// Table 2: the volume (equal-volume diameter 19.896 km) and the pole's inclination to
// the heliocentric orbit (100.39 deg) and to New Horizons' approach direction (41.096 deg).
void test_arrokoth_shape()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    const int arrokoth = scene.find("Arrokoth");
    const Body& body = scene.bodies[static_cast<size_t>(arrokoth)];
    check(body.shape != nullptr, "Arrokoth has a shape model");
    if (!body.shape) {
        return;
    }
    const ShapeModel& shape = *body.shape;
    glm::dvec3 lo(1e300), hi(-1e300);
    for (const glm::vec3& p : shape.positions) {
        lo = glm::min(lo, glm::dvec3(p));
        hi = glm::max(hi, glm::dvec3(p));
    }
    const glm::dvec3 size = hi - lo;
    check(std::abs(size.x - 34.546) < 0.01 && std::abs(size.y - 19.838) < 0.01 && std::abs(size.z - 13.822) < 0.01,
          "Arrokoth overall dimensions (km)", size.x);
    check(std::abs(size.x - 2.0 * body.equatorial_radius_km) < 0.01 &&
              std::abs(size.z - 2.0 * body.polar_radius_km) < 0.01,
          "Arrokoth radii are half its dimensions", body.equatorial_radius_km);

    const double d_equal = std::cbrt(6.0 * signed_volume(shape) / kPi);
    check(std::abs(d_equal - 19.896) < 0.01, "Arrokoth equal-volume diameter (km)", d_equal);
    const auto [min_albedo, max_albedo] = std::minmax_element(shape.albedo.begin(), shape.albedo.end());
    check(*min_albedo > 0.5f && *max_albedo < 2.5f, "Arrokoth relative albedo range", *max_albedo);
    check(shape.map_u.empty(), "Arrokoth's albedo is per vertex (no map coordinates)");

    auto angle_deg = [](const glm::dvec3& a, const glm::dvec3& b) {
        return std::acos(std::clamp(glm::dot(glm::normalize(a), glm::normalize(b)), -1.0, 1.0)) / kDegToRad;
    };
    const double ca = tdb_from_jd_tdb(2458484.73230); // 2019-01-01 05:34:31 TDB
    const State orbit = scene.icrf_state_at(arrokoth, ca);
    const double to_orbit = angle_deg(body.pole, glm::cross(orbit.position, orbit.velocity));
    check(std::abs(to_orbit - 100.39) < 0.05, "Arrokoth pole to its orbit normal (deg)", to_orbit);
    // Toward the approaching spacecraft, a day out (the approach is nearly straight).
    const glm::dvec3 approach =
        scene.icrf_state_at(scene.find("New Horizons"), ca - kSecondsPerDay).position -
        scene.icrf_state_at(arrokoth, ca - kSecondsPerDay).position;
    const double to_approach = angle_deg(body.pole, approach);
    check(std::abs(to_approach - (180.0 - 41.096)) < 0.1, "Arrokoth pole to the approach direction (deg)", to_approach);
    std::printf("info: Arrokoth pole to orbit normal %.2f deg, to the approaching spacecraft %.3f deg\n", to_orbit,
                to_approach);
}

// Phobos and Deimos: Thomas' shape models (west longitudes in the tables, east in the
// meshes) with map coordinates for their textures. Volumes against the mean radii in
// JPL SSD's satellite physical parameters (Archinal et al. 2018), and Stickney (49 W)
// as the check that the longitudes were turned around.
void test_mars_moon_shapes()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    struct Moon {
        const char* name;
        double mean_radius_km; // [SPHY]
        double tolerance_km;
    };
    for (const Moon& moon : {Moon{"Phobos", 11.08, 0.05}, Moon{"Deimos", 6.2, 0.1}}) {
        const Body& body = scene.bodies[static_cast<size_t>(scene.find(moon.name))];
        const std::string name = moon.name;
        check(body.shape != nullptr && !body.texture.empty(), (name + " has a shape model and a texture").c_str());
        if (!body.shape) {
            continue;
        }
        const ShapeModel& shape = *body.shape;
        check(shape.map_u.size() == shape.positions.size(), (name + " mesh has map coordinates").c_str());
        const double r_equal = std::cbrt(3.0 * signed_volume(shape) / (4.0 * kPi));
        check(std::abs(r_equal - moon.mean_radius_km) < moon.tolerance_km, (name + " equal-volume radius (km)").c_str(),
              r_equal);

        // Long axis toward Mars (the prime meridian), shortest along the pole, like the PCK ellipsoid.
        glm::dvec3 extent(0.0);
        for (const glm::vec3& p : shape.positions) {
            extent = glm::max(extent, glm::abs(glm::dvec3(p)));
        }
        check(extent.x > extent.y && extent.y > extent.z, (name + " axes ordered x > y > z").c_str(), extent.x);
        check(std::abs(extent.x - body.equatorial_radius_km) < 1.0, (name + " long axis near the PCK radius (km)").c_str(),
              extent.x);

        const double u_error = map_u_error(shape);
        check(u_error < 1e-6, (name + " map u matches the vertex longitude").c_str(), u_error);
    }

    // Stickney: the mean radius within 12 deg of its center is well below that of the
    // surrounding ring at 2 S 50 W, and not at the mirrored 2 S 50 E.
    const ShapeModel* phobos = scene.bodies[static_cast<size_t>(scene.find("Phobos"))].shape.get();
    if (!phobos) {
        return;
    }
    auto crater_depth = [&](double lon_east_deg) {
        const double lat = -2.0 * kDegToRad, lon = lon_east_deg * kDegToRad;
        const glm::dvec3 center(std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat));
        double inner = 0.0, ring = 0.0;
        int n_inner = 0, n_ring = 0;
        for (const glm::vec3& p : phobos->positions) {
            const double r = glm::length(glm::dvec3(p));
            const double angle = std::acos(std::clamp(glm::dot(glm::dvec3(p) / r, center), -1.0, 1.0)) / kDegToRad;
            if (angle < 12.0) {
                inner += r;
                ++n_inner;
            } else if (angle > 16.0 && angle < 24.0) {
                ring += r;
                ++n_ring;
            }
        }
        return ring / n_ring - inner / n_inner;
    };
    const double depth_west = crater_depth(-50.0);
    const double depth_east = crater_depth(50.0);
    check(depth_west > 0.5, "Stickney is a depression at 50 W (km)", depth_west);
    check(depth_west > depth_east + 0.5, "no Stickney at the mirrored 50 E (km)", depth_east);
    std::printf("info: Phobos mean radius within 12 deg of 2 S 50 W is %.2f km below the ring (50 E: %.2f km)\n",
                depth_west, depth_east);
}

// Vesta: the DLR Dawn HAMO DTM as a 1.5 deg grid, with the mosaic of the same release.
// Its volume against SBDB's equivalent diameter (522.77 +- 0.1 km, Park et al. 2025);
// the 2013 DTM and the grid's smoothing allow a few tenths of a km.
void test_vesta_shape()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    const Body& body = scene.bodies[static_cast<size_t>(scene.find("Vesta"))];
    check(body.shape != nullptr && !body.texture.empty(), "Vesta has a shape model and a texture");
    if (!body.shape) {
        return;
    }
    const ShapeModel& shape = *body.shape;
    const double d_equal = 2.0 * std::cbrt(3.0 * signed_volume(shape) / (4.0 * kPi));
    check(std::abs(d_equal - 522.77) < 0.6, "Vesta equal-volume diameter (km)", d_equal);
    glm::dvec3 lo(1e300), hi(-1e300);
    for (const glm::vec3& p : shape.positions) {
        lo = glm::min(lo, glm::dvec3(p));
        hi = glm::max(hi, glm::dvec3(p));
    }
    const glm::dvec3 size = hi - lo;
    check(size.x > size.y && size.y > size.z, "Vesta axes ordered x > y > z", size.x);
    const double u_error = map_u_error(shape);
    check(shape.map_u.size() == shape.positions.size() && u_error < 1e-6, "Vesta map u matches the vertex longitude",
          u_error);
}

// Pallas: the VLT/SPHERE MPCD shape model ([PALM] in solar_system.toml). Volume against
// Vernazza et al. 2021 (Table 1: D = 511 +- 4 km); the scene's IAU rotation against the
// spin solution of the ADAM model whose body frame the mesh shares (DAMIT model 4395,
// Marsset et al. 2020): lambda 42, beta -15 deg, P 7.81322 h, phi0 = 0 at JD 2433827.77154,
// body -> J2000 ecliptic = Rz(lambda) Ry(90 deg - beta) Rz(phi).
void test_pallas_shape()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    scene.set_active_frame(0); // inertial (ICRF axes)
    const int pallas = scene.find("Pallas");
    const Body& body = scene.bodies[static_cast<size_t>(pallas)];
    check(body.shape != nullptr && body.texture.empty(), "Pallas has a shape model (no texture)");
    if (!body.shape) {
        return;
    }
    const double d_equal = std::cbrt(6.0 * signed_volume(*body.shape) / kPi);
    check(std::abs(d_equal - 511.0) < 4.0, "Pallas equal-volume diameter (km)", d_equal);

    const double period_h = 360.0 * 24.0 / body.pm_rate_deg_per_day;
    check(std::abs(period_h - 7.81322) < 1e-5, "Pallas sidereal period (h)", period_h);
    // DAMIT's IAU conversion (from the pole rounded to RA 44, Dec 1): W0 = 41.7 deg.
    check(std::abs(body.pm_w0_deg - 41.7) < 0.15, "Pallas W0 near DAMIT's IAU value (deg)", body.pm_w0_deg);

    auto rotation_y = [](double a) {
        glm::dmat3 m(1.0); // columns
        m[0] = glm::dvec3(std::cos(a), 0.0, -std::sin(a));
        m[2] = glm::dvec3(std::sin(a), 0.0, std::cos(a));
        return m;
    };
    const double obliquity = 84381.448 / 3600.0 * kDegToRad; // IAU 1976
    const glm::dmat3 expected =
        rotation_x(obliquity) * rotation_z(42.0 * kDegToRad) * rotation_y((90.0 + 15.0) * kDegToRad);
    scene.update(tdb_from_jd_tdb(2433827.77154)); // phi = 0
    double worst = 0.0;
    for (int k : {0, 2}) {
        const double cos_angle = glm::dot(scene.bodies[static_cast<size_t>(pallas)].orientation[k], expected[k]);
        worst = std::max(worst, std::acos(std::clamp(cos_angle, -1.0, 1.0)) / kDegToRad);
    }
    check(worst < 0.05, "Pallas body axes match the DAMIT spin solution at t0 (deg)", worst);
}

// Amalthea and Thebe (Stooke, Jupiter scene) and Hyperion (Thomas, Cassini; Saturn
// scene): volumes against the mean radii in JPL SSD's satellite physical parameters
// (within their uncertainties). The Jovian moons' long axes point at Jupiter (+x).
void test_jupiter_saturn_small_moon_shapes()
{
    Scene jupiter = load_scene_or_die("jupiter.toml");
    Scene saturn = load_scene_or_die("saturn.toml");
    struct Moon {
        const Scene* scene;
        const char* name;
        double mean_radius_km; // [SPHY]
        double sigma_km;
        bool long_axis_x;
        bool textured;
    };
    for (const Moon& m : {Moon{&jupiter, "Amalthea", 83.5, 3.0, true, true}, Moon{&jupiter, "Thebe", 49.3, 4.0, true, false},
                          Moon{&saturn, "Hyperion", 135.0, 4.0, false, false}}) {
        const std::string name = m.name;
        const Body& body = m.scene->bodies[static_cast<size_t>(m.scene->find(m.name))];
        check(body.shape != nullptr, (name + " has a shape model").c_str());
        if (!body.shape) {
            continue;
        }
        const double r_equal = std::cbrt(3.0 * signed_volume(*body.shape) / (4.0 * kPi));
        check(std::abs(r_equal - m.mean_radius_km) < m.sigma_km, (name + " equal-volume radius (km)").c_str(), r_equal);
        check(body.texture.empty() != m.textured, (name + " texture as expected").c_str());
        if (m.textured) {
            check(map_u_error(*body.shape) < 1e-6, (name + " map u matches the vertex longitude").c_str());
        }
        if (m.long_axis_x) {
            glm::dvec3 extent(0.0);
            for (const glm::vec3& p : body.shape->positions) {
                extent = glm::max(extent, glm::abs(glm::dvec3(p)));
            }
            check(extent.x > extent.y && extent.x > extent.z, (name + " long axis toward Jupiter").c_str(), extent.x);
        }
    }
}

// Hyperion (Saturn scene) against Harbison, Thomas & Nicholson 2011 ([HSPN]). The mesh's
// uniform-density principal axes: Sect. 2 gives A/C = 0.58, B/C = 0.87 (+- 0.03), and
// Table 1 the 2005-09-25 spin, i.e. the mesh's z axis, as (0.902, 0.133, 0.411) along them.
// Table 1's Euler angles at that time (principal axes -> xyz = Rz(theta) Rx(phi) Rz(psi),
// z = Saturn's pole, x = Saturn -> pericenter) reproduce Table 2's omega_xyz, and give the
// orientation the scene should show when the osculating M is 303 deg.
void test_hyperion_spin_state()
{
    Scene scene = load_scene_or_die("saturn.toml");
    scene.set_active_frame(0); // inertial (ICRF axes)
    const int hyperion = scene.find("Hyperion");
    const Body& body = scene.bodies[static_cast<size_t>(hyperion)];
    if (!body.shape) {
        check(false, "Hyperion has a shape model");
        return;
    }

    glm::dvec3 moments(0.0);
    glm::dmat3 axes = symmetric_eigen(uniform_inertia(*body.shape), moments); // columns A, B, C
    check(std::abs(moments.x / moments.z - 0.58) < 0.03, "Hyperion A/C (uniform density)", moments.x / moments.z);
    check(std::abs(moments.y / moments.z - 0.87) < 0.03, "Hyperion B/C (uniform density)", moments.y / moments.z);
    for (int k = 0; k < 3; ++k) {
        if (axes[k].z < 0.0) { // the signs that make Table 1's spin components positive
            axes[k] = -axes[k];
        }
    }
    check(glm::dot(glm::cross(axes[0], axes[1]), axes[2]) > 0.999, "Hyperion principal axes right-handed");
    const glm::dvec3 spin_abc = glm::normalize(glm::dvec3(0.902, 0.133, 0.411));
    const glm::dvec3 z_abc(axes[0].z, axes[1].z, axes[2].z); // the mesh's z axis along A, B, C
    const double z_to_spin = std::acos(std::clamp(glm::dot(z_abc, spin_abc), -1.0, 1.0)) / kDegToRad;
    check(z_to_spin < 3.0, "Hyperion mesh z axis vs Table 1 spin direction (deg)", z_to_spin);

    // The Euler angle convention: Table 1's spin along the principal axes, rotated to xyz.
    const glm::dmat3 abc_to_xyz = rotation_z(2.989) * rotation_x(1.685) * rotation_z(1.641);
    const glm::dvec3 omega_xyz(1.151, 2.018, 3.565); // Table 2, units of n
    const double omega_error = glm::length(abc_to_xyz * (4.255 * spin_abc) - omega_xyz);
    check(omega_error < 0.02, "Hyperion Table 1 Euler angles reproduce Table 2 omega (n)", omega_error);

    // The epoch: osculating M = 303 deg (Horizons' GM for Saturn).
    scene.update(tdb_from_utc_jd(jd_from_calendar({2005, 9, 25, 17, 33, 22})));
    const Body& saturn = scene.bodies[static_cast<size_t>(scene.find("Saturn"))];
    const glm::dvec3 r = body.icrf_position - saturn.icrf_position;
    const glm::dvec3 v = body.icrf_velocity - saturn.icrf_velocity;
    const glm::dvec3 h = glm::cross(r, v);
    const glm::dvec3 ecc = glm::cross(v, h) / 3.7931206604853049e7 - glm::normalize(r);
    const double e = glm::length(ecc);
    const double f = std::atan2(glm::dot(glm::cross(ecc, r), glm::normalize(h)), glm::dot(ecc, r));
    const double ea = 2.0 * std::atan(std::sqrt((1.0 - e) / (1.0 + e)) * std::tan(0.5 * f));
    const double mean_anomaly = wrap_two_pi(ea - e * std::sin(ea)) / kDegToRad;
    check(std::abs(mean_anomaly - 303.0) < 0.05, "Hyperion osculating M at the Table 1 epoch (deg)", mean_anomaly);
    check(std::abs(e - 0.113) < 0.001, "Hyperion osculating e at the Table 1 epoch", e);

    const glm::dvec3 z = saturn.orientation[2]; // Saturn's pole
    const glm::dvec3 x = glm::normalize(ecc - glm::dot(ecc, z) * z);
    const glm::dmat3 xyz_to_icrf(x, glm::cross(z, x), z);
    const glm::dmat3 expected = xyz_to_icrf * abc_to_xyz * glm::transpose(axes); // mesh -> ICRF
    const glm::dmat3 relative = glm::transpose(body.orientation) * expected;
    const double trace = relative[0][0] + relative[1][1] + relative[2][2];
    const double angle = std::acos(std::clamp(0.5 * (trace - 1.0), -1.0, 1.0)) / kDegToRad;
    check(angle < 0.5, "Hyperion orientation at 2005-09-25 matches Table 1 (deg)", angle);
    const glm::dvec3 omega_icrf = glm::normalize(xyz_to_icrf * omega_xyz);
    const double pole_to_omega = std::acos(std::clamp(glm::dot(body.orientation[2], omega_icrf), -1.0, 1.0)) / kDegToRad;
    check(pole_to_omega < 3.0, "Hyperion pole vs Table 2 spin direction (deg)", pole_to_omega);
    const double spin_over_n = body.pm_rate_deg_per_day / 16.94; // n from Sect. 1
    check(std::abs(spin_over_n - 4.255) < 1e-4, "Hyperion spin rate / mean motion", spin_over_n);
    std::printf("info: Hyperion A/C %.3f, B/C %.3f, mesh z %.1f deg from the Table 1 spin; orientation within %.2f deg "
                "of Table 1, pole %.1f deg from Table 2\n",
                moments.x / moments.z, moments.y / moments.z, z_to_spin, angle, pole_to_omega);
}

// Halley (ISEE-3 scene): Stooke's shape model ([STKS]) against Szego 1991's independent
// Vega/Giotto model ([SZG]: a 15.3 x 7.2 x 7.22 km box, 365 km^3), and the long-axis mode
// of Belton et al. 1991 ([BLT]): from the long axis at 66.0 deg to M, 3.69 d around it and
// 7.1 d about itself, the paper derives a total spin of 2.84 d at 21.4 deg from M. Then the
// phases: the long axis at the encounters ([BLT], [SA91]) and Stooke's longitudes ([STKS]).
void test_halley_shape()
{
    Scene scene = load_scene_or_die("isee3.toml");
    scene.set_active_frame(0); // inertial (ICRF axes)
    const int halley = scene.find("Halley");
    const Body& body = scene.bodies[static_cast<size_t>(halley)];
    check(body.shape != nullptr && body.texture.empty(), "Halley has a shape model (no texture)");
    check(body.free_precession.enabled && body.pm_rate_deg_per_day == 0.0, "Halley precesses freely");
    if (!body.shape) {
        return;
    }
    glm::dvec3 lo(1e300), hi(-1e300);
    for (const glm::vec3& p : body.shape->positions) {
        lo = glm::min(lo, glm::dvec3(p));
        hi = glm::max(hi, glm::dvec3(p));
    }
    const glm::dvec3 size = hi - lo;
    check(std::abs(size.z - 15.3) < 0.5, "Halley long axis along z (km)", size.z);
    check(std::abs(size.x - 7.2) < 0.5 && std::abs(size.y - 7.2) < 0.5, "Halley short axes (km)",
          std::max(size.x, size.y));
    const double volume = signed_volume(*body.shape);
    check(std::abs(volume / 365.0 - 1.0) < 0.15, "Halley volume near Szego's model (km^3)", volume);

    // Angular velocity from the orientation a minute apart: R' R^T = [omega]x.
    const double t0 = tdb_from_utc_jd(jd_from_calendar({1986, 3, 14, 0, 0, 0}));
    const double dt = 60.0;
    double worst_nutation = 0.0;
    for (double days : {0.0, 1.3, 2.9}) {
        const double t = t0 + days * kSecondsPerDay;
        scene.update(t);
        const glm::dmat3 r0 = scene.bodies[static_cast<size_t>(halley)].orientation;
        scene.update(t + dt);
        const glm::dmat3 r1 = scene.bodies[static_cast<size_t>(halley)].orientation;
        const glm::dmat3 w = (r1 - r0) * glm::transpose(r0) / dt;
        const glm::dvec3 omega(w[1][2], w[2][0], w[0][1]);
        const double period_d = kTwoPi / glm::length(omega) / kSecondsPerDay;
        const double spin_to_m = std::acos(glm::dot(glm::normalize(omega), body.pole)) / kDegToRad;
        const double nutation = std::acos(glm::dot(r0[2], body.pole)) / kDegToRad;
        worst_nutation = std::max(worst_nutation, std::abs(nutation - 66.0));
        check(std::abs(period_d - 2.84) < 0.01, "Halley total spin period (d)", period_d);
        check(std::abs(spin_to_m - 21.4) < 0.1, "Halley spin vector to M (deg)", spin_to_m);
        if (days == 0.0) {
            std::printf("info: Halley mesh %.2f x %.2f x %.2f km, %.1f km^3; total spin %.3f d at %.2f deg from M\n",
                        size.x, size.y, size.z, volume, period_d, spin_to_m);
        }
    }
    check(worst_nutation < 1e-6, "Halley long axis stays 66 deg from M (deg)", worst_nutation);

    // The phases. B1950 directions to J2000 with SPICE's FK4 -> J2000 rotation (pxform), and
    // B1950 ecliptic to equator with the IAU 1980 obliquity at B1950.0.
    const glm::dmat3 fk4_to_j2000(0.9999256794956877, 0.0111814832391717, 0.0048590037723143,
                                  -0.0111814832204662, 0.9999374848933135, -0.0000271702937440,
                                  -0.0048590038153592, -0.0000271625947142, 0.9999881946023742);
    const double t_b1950 = (2433282.4235 - kJ2000Jd) / 36525.0;
    const double eps_b1950 = (84381.448 - 46.8150 * t_b1950) / 3600.0 * kDegToRad;
    const auto radec = [](double ra_deg, double dec_deg) {
        const double ra = ra_deg * kDegToRad;
        const double dec = dec_deg * kDegToRad;
        return glm::dvec3(std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec));
    };
    const auto big_end_at = [&](double t) {
        scene.update(t);
        return scene.bodies[static_cast<size_t>(halley)].orientation[2];
    };
    const auto angle_deg = [](const glm::dvec3& a, const glm::dvec3& b) {
        return std::acos(std::clamp(glm::dot(glm::normalize(a), glm::normalize(b)), -1.0, 1.0)) / kDegToRad;
    };
    // [BLT] sec. 2.2: at the Vega 2 encounter (JD 2446498.80556 UT, the closest approach of
    // [VEGA]) the long axis pierces the "big" end (Stooke's body +z) toward RA, Dec (B1950) =
    // 313.20, -7.52.
    const double vega2_ca = tdb_from_utc_jd(jd_from_calendar({1986, 3, 9, 7, 20, 0}));
    const double p_error = angle_deg(big_end_at(vega2_ca), fk4_to_j2000 * radec(313.20, -7.52));
    check(p_error < 0.1, "Halley big end at the Vega 2 encounter (Belton's P, deg)", p_error);
    // Independent: [SA91] Table II, combination (2,1,2), the big end in ecliptic (B1950)
    // longitude, latitude at the other two encounters (closest approaches from [VEGA], [GIO]).
    const struct {
        const char* name;
        double t;
        double lon, lat;
    } encounters[] = {
        {"Vega 1", tdb_from_utc_jd(jd_from_calendar({1986, 3, 6, 7, 20, 6})), 259.0, -15.0},
        {"Vega 2", vega2_ca, 310.0, 9.0},
        {"Giotto", tdb_from_utc_jd(jd_from_calendar({1986, 3, 14, 0, 3, 2})), 235.0, -32.0},
    };
    for (const auto& e : encounters) {
        const glm::dvec3 ecl = radec(e.lon, e.lat);
        const glm::dvec3 eq(ecl.x, ecl.y * std::cos(eps_b1950) - ecl.z * std::sin(eps_b1950),
                            ecl.y * std::sin(eps_b1950) + ecl.z * std::cos(eps_b1950));
        const double error = angle_deg(big_end_at(e.t), fk4_to_j2000 * eq);
        check(error < 5.0, (std::string("Halley big end at ") + e.name + " vs Samarasinha & A'Hearn (deg)").c_str(),
              error);
    }
    // [STKS]: 90 deg east (V1.0's 270 W) below Vega 2 in its image 1.5 s before the closest
    // approach; that image's header [VTV] has 8,032 km and a phase angle of 28.7 deg.
    {
        scene.update(vega2_ca - 1.5);
        const auto position = [&](const char* name) {
            return scene.bodies[static_cast<size_t>(scene.find(name))].icrf_position;
        };
        const glm::dvec3 to_vega2 = position("Vega 2") - position("Halley");
        const glm::dvec3 sub = glm::transpose(scene.bodies[static_cast<size_t>(halley)].orientation) *
                               glm::normalize(to_vega2);
        const double lon = std::atan2(sub.y, sub.x) / kDegToRad;
        const double phase = angle_deg(to_vega2, position("Sun") - position("Halley"));
        check(std::abs(glm::length(to_vega2) - 8032.0) < 20.0, "Vega 2 range in its closest image (km)",
              glm::length(to_vega2));
        check(std::abs(phase - 28.7) < 1.0, "Vega 2 phase angle in its closest image (deg)", phase);
        check(std::abs(lon - 90.0) < 0.1, "Halley longitude below Vega 2 (deg east)", lon);
    }
}

} // namespace

void run_shape_tests()
{
    test_arrokoth_shape();
    test_mars_moon_shapes();
    test_vesta_shape();
    test_pallas_shape();
    test_jupiter_saturn_small_moon_shapes();
    test_hyperion_spin_state();
    test_halley_shape();
}
