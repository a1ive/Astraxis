// Tests for the shape models of small bodies.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>

using namespace astraxis;

namespace {

// Plane intersection and oriented edge tests provide a geometric reference
// independent of the shader's barycentric triangle intersection.
double triangle_hit(const glm::dvec3& o, const glm::dvec3& d, const glm::dvec3& a,
                    const glm::dvec3& b, const glm::dvec3& c, double limit)
{
    const glm::dvec3 n = glm::cross(b - a, c - a);
    const double denominator = glm::dot(n, d);
    if (std::abs(denominator) < 1e-12) {
        return limit;
    }
    const double t = glm::dot(n, a - o) / denominator;
    const glm::dvec3 p = o + t * d;
    if (t > 0.0 && t < limit && glm::dot(n, glm::cross(b - a, p - a)) >= 0.0 &&
        glm::dot(n, glm::cross(c - b, p - b)) >= 0.0 && glm::dot(n, glm::cross(a - c, p - c)) >= 0.0) {
        return t;
    }
    return limit;
}

void test_shape_occlusion()
{
    check(build_shape_occlusion(ShapeModel{}).empty(), "empty mesh has no occlusion tree");
    for (const char* file : {"halley.mesh", "67p.mesh"}) {
        ShapeModel shape;
        std::string error;
        const bool loaded = load_shape_model(std::filesystem::path(ASTRAXIS_ASSET_DIR) / "shapes" / file, shape, &error);
        check(loaded, "load comet occlusion mesh");
        if (!loaded) {
            continue;
        }
        const auto nodes = build_shape_occlusion(shape);
        check(nodes.size() == 2 * shape.indices.size() / 3 - 1 && nodes.front().escape == nodes.size(),
              "occlusion tree covers all triangles");
        bool links_valid = true;
        for (size_t i = 0; i < nodes.size(); ++i) {
            const auto& n = nodes[i];
            links_valid = links_valid && n.escape > i && n.escape <= nodes.size();
            if (n.leaf) {
                links_valid = links_valid && n.escape == i + 1;
            } else {
                const size_t right = nodes[i + 1].escape;
                links_valid = links_valid && right < n.escape && nodes[right].escape == n.escape;
            }
        }
        check(links_valid, "occlusion subtree escape links");
        if (!links_valid) {
            continue;
        }
        double radius = 0.0;
        for (const auto& p : shape.positions) {
            radius = std::max(radius, glm::length(glm::dvec3(p)));
        }
        int hits = 0, misses = 0;
        double worst = 0.0;
        // Six perspective views include axis-parallel rays, the limb, and both
        // lobes of 67P. Limit each ray to a foreground plane and beyond the nucleus.
        for (int axis = 0; axis < 3; ++axis) {
            for (double side : {-1.0, 1.0}) {
                for (int x = -6; x <= 6; ++x) {
                    for (int y = -6; y <= 6; ++y) {
                        glm::dvec3 o(0.0), target(0.0);
                        o[axis] = side * 4.0 * radius;
                        target[(axis + 1) % 3] = x * radius / 5.0;
                        target[(axis + 2) % 3] = y * radius / 5.0;
                        const glm::dvec3 d = glm::normalize(target - o);
                        for (double limit : {3.5 * radius, 10.0 * radius}) {
                            double reference = limit;
                            for (size_t i = 0; i < shape.indices.size(); i += 3) {
                                reference = triangle_hit(o, d, glm::dvec3(shape.positions[shape.indices[i]]),
                                    glm::dvec3(shape.positions[shape.indices[i + 1]]),
                                    glm::dvec3(shape.positions[shape.indices[i + 2]]), reference);
                            }
                            double hit = limit;
                            size_t index = 0;
                            while (index < nodes.size()) {
                                const auto& n = nodes[index];
                                const glm::dvec3 a(n.a), b(n.b), c(n.c);
                                if (n.leaf) {
                                    hit = triangle_hit(o, d, a, b, c, hit);
                                    index = n.escape;
                                } else {
                                    double near_t = 0.0, far_t = hit;
                                    for (int k = 0; k < 3; ++k) {
                                        if (std::abs(d[k]) < 1e-8) {
                                            if (o[k] < a[k] || o[k] > b[k]) {
                                                far_t = -1.0;
                                            }
                                        } else {
                                            const double t0 = (a[k] - o[k]) / d[k];
                                            const double t1 = (b[k] - o[k]) / d[k];
                                            near_t = std::max(near_t, std::min(t0, t1));
                                            far_t = std::min(far_t, std::max(t0, t1));
                                        }
                                    }
                                    index = near_t <= far_t ? index + 1 : n.escape;
                                }
                            }
                            worst = std::max(worst, std::abs(hit - reference));
                            if (reference < limit) {
                                ++hits;
                            } else {
                                ++misses;
                            }
                        }
                    }
                }
            }
        }
        check(hits > 50 && misses > 50, "occlusion rays cover surface and sky");
        check(worst < 1e-4, "comet BVH hits match unaccelerated mesh (km)", worst);
    }
}

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

void test_proteus_shape()
{
    Scene scene = load_scene_or_die("solar_system.toml");
    const int id = scene.find("Proteus");
    check(id >= 0, "Proteus exists in the Solar System");
    if (id < 0) {
        return;
    }
    const Body& body = scene.bodies[static_cast<size_t>(id)];
    check(body.parent == scene.find("Neptune") && body.shape != nullptr, "Proteus mesh orbits Neptune");
    if (!body.shape) {
        return;
    }
    const ShapeModel& shape = *body.shape;
    const double radius = std::cbrt(3.0 * signed_volume(shape) / (4.0 * kPi));
    check(std::abs(radius - 201.0) < 0.1, "Proteus mesh has km units and outward winding", radius);
    // Independent PDS n8proteus.tab: west longitude 5 deg, latitude -50 deg,
    // radius 220 km. A reflected longitude convention misses this point.
    const double lon = -5.0 * kDegToRad, lat = -50.0 * kDegToRad;
    const glm::dvec3 reference = 220.0 * glm::dvec3(std::cos(lat) * std::cos(lon),
                                                  std::cos(lat) * std::sin(lon), std::sin(lat));
    double nearest = 1e300;
    for (const glm::vec3& p : shape.positions) {
        nearest = std::min(nearest, glm::length(glm::dvec3(p) - reference));
    }
    check(nearest < 0.00002, "Proteus orientation matches PDS west longitudes (km)", nearest);
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

// Amalthea and Thebe (Stooke, Jupiter scene), Hyperion (Thomas, Cassini; Saturn scene)
// and Phoebe (Gaskell, Cassini; Saturn scene): volumes against the mean radii in JPL
// SSD's satellite physical parameters (within their uncertainties). The Jovian moons'
// long axes point at Jupiter (+x); Phoebe's shortest axis is its spin axis (z), as in
// the PCK's triaxial radii.
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
        bool short_axis_z;
        bool textured;
    };
    for (const Moon& m : {Moon{&jupiter, "Amalthea", 83.5, 3.0, true, false, true},
                          Moon{&jupiter, "Thebe", 49.3, 4.0, true, false, false},
                          Moon{&saturn, "Hyperion", 135.0, 4.0, false, false, false},
                          Moon{&saturn, "Janus", 89.2, 0.8, true, true, false},
                          Moon{&saturn, "Epimetheus", 58.2, 1.2, true, true, false},
                          Moon{&saturn, "Phoebe", 106.5, 0.7, false, true, false}}) {
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
        glm::dvec3 extent(0.0);
        for (const glm::vec3& p : body.shape->positions) {
            extent = glm::max(extent, glm::abs(glm::dvec3(p)));
        }
        if (m.long_axis_x) {
            check(extent.x > extent.y && extent.x > extent.z, (name + " long axis toward the planet").c_str(), extent.x);
        }
        if (m.short_axis_z) {
            check(extent.z < extent.x && extent.z < extent.y, (name + " shortest axis along the spin axis").c_str(), extent.z);
        }
    }

    // The PCK triaxial ellipsoids give an independent scale check for the seven
    // Thomas meshes. Craters and equatorial ridges change their volume slightly.
    for (const char* name : {"Prometheus", "Pandora", "Pan", "Daphnis", "Telesto", "Calypso", "Helene"}) {
        const int id = saturn.find(name);
        check(id > 0, "Saturn has the ring and Trojan moon");
        if (id <= 0) {
            continue;
        }
        const Body& body = saturn.bodies[static_cast<size_t>(id)];
        check(body.shape != nullptr, "ring/Trojan moon has a Thomas shape");
        if (!body.shape) {
            continue;
        }
        const double r_mesh = std::cbrt(3.0 * signed_volume(*body.shape) / (4.0 * kPi));
        const double r_pck = std::cbrt(body.equatorial_radius_km * body.equatorial_radius_b_km * body.polar_radius_km);
        check(std::abs(r_mesh / r_pck - 1.0) < 0.1, "Thomas mesh scale agrees with PCK radii", r_mesh / r_pck);
        check(body.synchronous_with == saturn.find("Saturn") && body.texture.empty(),
              "small moons use synchronous, untextured shapes");
        check(body.shape->positions.size() > 13000, "Thomas mesh retains the published vertices");
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

// Didymos and Dimorphos (didymos.toml): the DART shape models (Didymos v003, Dimorphos
// v004). Dimorphos' volume against Daly et al. 2023 (Table 1: 0.00181 km^3 +- 10%), that
// of Didymos against the oblate spheroid of Naidu et al. 2024 (Table 1: 394 x 290 m, to
// 8% given its uncertainties); the radii are half the extents.
void test_didymos_shapes()
{
    const Scene scene = load_scene_or_die("didymos.toml");
    struct Expected {
        const char* name;
        double volume_km3;
        double tolerance;
    };
    const double spheroid = 4.0 / 3.0 * kPi * 0.394 * 0.394 * 0.290;
    for (const Expected& e : {Expected{"Didymos", spheroid, 0.10}, Expected{"Dimorphos", 0.00181, 0.10}}) {
        const Body& body = scene.bodies[static_cast<size_t>(scene.find(e.name))];
        const std::string name = e.name;
        check(body.shape != nullptr, (name + " has a shape model").c_str());
        if (!body.shape) {
            continue;
        }
        const double volume = signed_volume(*body.shape);
        check(std::abs(volume / e.volume_km3 - 1.0) < e.tolerance, (name + " volume (km^3)").c_str(), volume);
        glm::dvec3 lo(1e300), hi(-1e300);
        for (const glm::vec3& p : body.shape->positions) {
            lo = glm::min(lo, glm::dvec3(p));
            hi = glm::max(hi, glm::dvec3(p));
        }
        const glm::dvec3 half = 0.5 * (hi - lo);
        check(std::abs(half.x - body.equatorial_radius_km) < 0.001 && std::abs(half.y - body.equatorial_radius_b_km) < 0.001 &&
                  std::abs(half.z - body.polar_radius_km) < 0.001,
              (name + " radii are half its extents (km)").c_str(), half.x);
        std::printf("info: %s volume %.5f km^3 (expected %.5f), extents %.0f x %.0f x %.0f m\n", e.name, volume,
                    e.volume_km3, 2e3 * half.x, 2e3 * half.y, 2e3 * half.z);
    }
}

// Endurance (gargantua.toml): the tribute model of tools/shapes/make_endurance.py, a
// fixed hull and a turning ring, built with hard edges.
void test_endurance_shape()
{
    Scene scene = load_scene_or_die("gargantua.toml");
    const int craft = scene.find("Endurance");
    check(craft > 0, "Endurance in the Gargantua scene");
    if (craft <= 0) {
        return;
    }
    const Body& body = scene.bodies[static_cast<size_t>(craft)];
    check(body.kind == BodyKind::Spacecraft && body.shape && body.rotor.shape, "Endurance has a hull and a rotor");
    if (!body.shape || !body.rotor.shape) {
        return;
    }
    const ShapeModel& hull = *body.shape;
    const ShapeModel& ring = *body.rotor.shape;
    check(signed_volume(hull) > 0.0 && signed_volume(ring) > 0.0, "Endurance meshes are closed and face outward");

    // Hard edges: the flat tops of the hull's end rings (z = 10.2 m) face straight up,
    // and their round walls straight out; merged normals would be 45 deg between.
    int up = 0, out = 0, bent = 0;
    for (size_t k = 0; k < hull.positions.size(); ++k) {
        if (std::abs(hull.positions[k].z - 0.0102f) > 1e-6f) {
            continue;
        }
        const float nz = hull.normals[k].z;
        (std::abs(nz - 1.0f) < 1e-5f ? up : std::abs(nz) < 1e-5f ? out : bent) += 1;
    }
    check(up > 0 && out > 0 && bent == 0, "Endurance hull edges stay hard", static_cast<double>(bent));

    // The ring: 70 m across about the z axis it turns on (the pallets' corners a little
    // more); the modules reach 33 m out.
    double r_max = 0.0, r_module = 0.0;
    for (const glm::vec3& p : ring.positions) {
        const double r = std::hypot(p.x, p.y) * 1000.0;
        r_max = std::max(r_max, r);
        if (std::abs(p.z) * 1000.0 > 3.0) { // above the tube and pallets: the modules
            r_module = std::max(r_module, r);
        }
    }
    check(std::abs(r_max - 35.0) < 0.2, "Endurance ring 70 m across", 2.0 * r_max);
    check(r_module > 32.5 && r_module < 34.0, "Endurance modules reach 33 m out", r_module);
    check(body.equatorial_radius_km * 1000.0 >= r_max, "Endurance radius encloses the ring",
          body.equatorial_radius_km);

    // One turn every 11.44 s gives 1 g on the floor, 32.5 m from the axis.
    const double omega = kTwoPi / body.rotor.period_s;
    const double g = omega * omega * 32.5 / 9.80665;
    check(std::abs(g - 1.0) < 1e-3, "Endurance floor gravity (g)", g);
}

} // namespace

// The share of the mesh's vertices across the Sun line in a region that its bake left in the light,
// and how many there are. Their normals can point either way: a sheet is drawn from both sides.
template <typename Region>
std::pair<double, int> lit_share(const ShapeModel& shape, Region region)
{
    int lit = 0, count = 0;
    for (size_t k = 0; k < shape.positions.size(); ++k) {
        const glm::dvec3 p_m = glm::dvec3(shape.positions[k]) * 1000.0;
        if (std::abs(shape.normals[k].z) > 0.3f && region(p_m)) {
            ++count;
            lit += shape.albedo[k] > 0.5f ? 1 : 0;
        }
    }
    return {count > 0 ? static_cast<double>(lit) / count : -1.0, count};
}

const ShapeModel* telescopes_scene_model(const Scene& scene, const char* name)
{
    return scene.bodies[static_cast<size_t>(scene.find(name))].shape.get();
}

glm::dvec3 scene_direction(const Scene& scene, int from, int to, double t_tdb)
{
    return glm::normalize(scene.icrf_state_at(to, t_tdb).position - scene.icrf_state_at(from, t_tdb).position);
}

void test_spacecraft_models()
{
    struct Case {
        const char* file;
        const char* craft;
        const char* target;
        double t_tdb;
    };
    const double t2026 = tdb_from_julian_year(2026.84);
    const double perihelion = tdb_from_julian_year(2024.98);
    const Case cases[] = {{"solar_system.toml", "Voyager 1", "Earth", t2026},
                          {"solar_system.toml", "Voyager 2", "Earth", t2026},
                          {"solar_system.toml", "New Horizons", "Earth", t2026},
                          {"solar_system.toml", "Pioneer 10", "Earth", tdb_from_julian_year(1973.92)},
                          {"space_telescopes.toml", "Hubble", "Sun", t2026},
                          {"space_telescopes.toml", "Chandra", "Sun", tdb_from_julian_year(2010.0)},
                          {"space_telescopes.toml", "Roman", "Sun", tdb_from_julian_year(2027.5)},
                          {"space_telescopes.toml", "Spitzer", "Sun", tdb_from_julian_year(2010.0)},
                          {"space_telescopes.toml", "Kepler", "Sun", tdb_from_julian_year(2012.0)},
                          {"space_telescopes.toml", "TESS", "Sun", tdb_from_julian_year(2020.0)},
                          {"space_telescopes.toml", "WMAP", "Sun", tdb_from_julian_year(2005.0)},
                          {"mars.toml", "MAVEN", "Sun", tdb_from_julian_year(2016.0)},
                          {"mars.toml", "MRO", "Sun", tdb_from_julian_year(2025.3)},
                          {"mars.toml", "Odyssey", "Sun", tdb_from_julian_year(2025.3)},
                          {"mercury.toml", "MESSENGER", "Sun", tdb_from_julian_year(2012.0)},
                          {"venus.toml", "MESSENGER", "Sun", tdb_from_julian_year(2007.4)},
                          {"earth_moon.toml", "LADEE", "Sun", tdb_from_julian_year(2013.9)},
                          {"earth_moon.toml", "ARTEMIS-P1", "Sun", tdb_from_julian_year(2014.0)},
                          {"earth_moon.toml", "Artemis II", "Sun", tdb_from_julian_year(2026.26)},
                          {"space_telescopes.toml", "Gaia", "Sun", tdb_from_julian_year(2020.0)},
                          {"space_telescopes.toml", "Euclid", "Sun", tdb_from_julian_year(2026.4)},
                          {"mars.toml", "TGO", "Sun", tdb_from_julian_year(2025.3)},
                          {"mars.toml", "Europa Clipper", "Sun", tdb_from_julian_year(2025.165)},
                          {"didymos.toml", "DART", "Sun", tdb_from_julian_year(2022.6)},
                          {"rosetta.toml", "Rosetta", "Sun", tdb_from_julian_year(2015.6)},
                          {"solar_system.toml", "Pioneer 11", "Earth", tdb_from_julian_year(1979.6)},
                          {"saturn.toml", "Cassini", "Sun", tdb_from_julian_year(2007.3)},
                          {"jupiter.toml", "Cassini", "Sun", tdb_from_julian_year(2000.99)},
                          {"jupiter.toml", "Galileo", "Sun", tdb_from_julian_year(1997.0)},
                          {"jupiter.toml", "Juno", "Sun", tdb_from_julian_year(2024.8)},
                          {"space_telescopes.toml", "JWST", "Sun", t2026},
                          {"parker.toml", "Parker Solar Probe", "Sun", perihelion}};
    for (const Case& c : cases) {
        const std::string name = c.craft;
        Scene scene = load_scene_or_die(c.file);
        const int craft = scene.find(c.craft);
        check(craft > 0, (name + " is in its scene").c_str());
        if (craft <= 0) {
            continue;
        }
        const Body& body = scene.bodies[static_cast<size_t>(craft)];
        check(body.kind == BodyKind::Spacecraft && body.shape && !body.texture.empty(), (name + " has a model").c_str());
        if (!body.shape) {
            continue;
        }
        const ShapeModel& shape = *body.shape;
        check(shape.map_v.size() == shape.positions.size() && shape.map_u.size() == shape.positions.size(),
              (name + " has texture coordinates at every vertex").c_str());
        bool in_atlas = true, bake_is_binary = true;
        int shadowed = 0;
        for (size_t k = 0; k < shape.positions.size(); ++k) {
            in_atlas = in_atlas && shape.map_u[k] >= 0.0f && shape.map_u[k] <= 1.0f && shape.map_v[k] >= 0.0f &&
                       shape.map_v[k] <= 1.0f;
            bake_is_binary = bake_is_binary && (shape.albedo[k] == 0.0f || shape.albedo[k] == 1.0f);
            shadowed += shape.albedo[k] == 0.0f ? 1 : 0;
        }
        check(in_atlas, (name + " stays inside its atlas").c_str());
        check(bake_is_binary && shadowed > 0, (name + " has baked sunlight").c_str());

        // The attitude: +z at the target, +y toward the ecliptic north, a right-handed frame.
        const glm::dmat3 axes = scene.body_to_icrf(body, c.t_tdb);
        const glm::dvec3 self = scene.icrf_state_at(craft, c.t_tdb).position;
        const int target = scene.find(c.target);
        const glm::dvec3 to = target >= 0 ? scene.icrf_state_at(target, c.t_tdb).position
                                          : scene.sun_icrf_state_at(c.t_tdb).position;
        check(glm::dot(axes[2], glm::normalize(to - self)) > 1.0 - 1e-12, (name + " points +z at its target").c_str());
        const glm::dvec3 north = ecliptic_pole_icrf();
        check(std::abs(glm::dot(axes[0], north)) < 1e-12 && glm::dot(axes[1], north) >= 0.0,
              (name + " keeps +y toward the ecliptic north").c_str());
        check(std::abs(glm::determinant(axes) - 1.0) < 1e-12, (name + " has a right-handed frame").c_str());
    }

    // Scale and bake, in metres. Voyager: the dish, 3.7 m across about +z, its rim 0.6 to 0.9 m
    // above the origin; the only parts there further out are past 2 m.
    Scene solar = load_scene_or_die("solar_system.toml");
    const ShapeModel& voyager = *solar.bodies[static_cast<size_t>(solar.find("Voyager 1"))].shape;
    double dish = 0.0, reach = 0.0;
    for (const glm::vec3& p : voyager.positions) {
        const glm::dvec3 m = glm::dvec3(p) * 1000.0;
        const double r = std::hypot(m.x, m.y);
        if (m.z > 0.6 && m.z < 0.9 && r < 2.0) {
            dish = std::max(dish, r);
        }
        reach = std::max(reach, glm::length(m));
    }
    check(std::abs(2.0 * dish - 3.7) < 0.08, "the Voyager dish is 3.7 m across", 2.0 * dish);
    check(reach > 13.0 && reach < 16.0, "the Voyager boom reaches 13 m past the bus", reach);
    const auto bus = lit_share(voyager, [](const glm::dvec3& m) { return m.z < -0.2 && std::hypot(m.x, m.y) < 0.8; });
    check(bus.second > 0 && bus.first < 0.05, "the dish shades the Voyager bus", bus.first);

    // New Horizons: the dish, 2.1 m across, its rim above everything else; the generator reaches
    // out of its shade, about 2.1 m from the axis (the sunward end: its fins shade their roots),
    // and the bus does not.
    const ShapeModel& horizons = *solar.bodies[static_cast<size_t>(solar.find("New Horizons"))].shape;
    double horizons_dish = 0.0;
    for (const glm::vec3& p : horizons.positions) {
        const glm::dvec3 m = glm::dvec3(p) * 1000.0;
        if (m.z > 0.75) {
            horizons_dish = std::max(horizons_dish, std::hypot(m.x, m.y));
        }
    }
    check(std::abs(2.0 * horizons_dish - 2.1) < 0.08, "the New Horizons dish is 2.1 m across", 2.0 * horizons_dish);
    const auto horizons_bus =
        lit_share(horizons, [](const glm::dvec3& m) { return m.z < 0.2 && std::hypot(m.x, m.y) < 0.8; });
    check(horizons_bus.second > 0 && horizons_bus.first < 0.05, "the dish shades the New Horizons bus", horizons_bus.first);
    const auto generator = lit_share(horizons, [](const glm::dvec3& m) { return m.y > 1.5 && m.z > -0.12; });
    check(generator.second > 0 && generator.first > 0.5, "the New Horizons generator reaches out of the shade", generator.first);

    // Pioneer 10: the dish, 2.74 m across [PTM], shades the bus behind its own apex (0.53 m back);
    // the generators stand clear.
    const ShapeModel& pioneer = *solar.bodies[static_cast<size_t>(solar.find("Pioneer 10"))].shape;
    double pioneer_dish = 0.0;
    for (const glm::vec3& p : pioneer.positions) {
        const glm::dvec3 m = glm::dvec3(p) * 1000.0;
        if (m.z > -0.6 && m.z < 1.1 && std::abs(m.x) < 1.7) {
            pioneer_dish = std::max(pioneer_dish, std::abs(m.x));
        }
    }
    check(std::abs(2.0 * pioneer_dish - 2.74) < 0.08, "the Pioneer 10 dish is 2.74 m across", 2.0 * pioneer_dish);
    const auto pioneer_bus =
        lit_share(pioneer, [](const glm::dvec3& m) { return m.z < -0.6 && std::hypot(m.x, m.y) < 0.9; });
    check(pioneer_bus.second > 0 && pioneer_bus.first < 0.05, "the dish shades the Pioneer 10 bus", pioneer_bus.first);
    const auto pioneer_rtg = lit_share(pioneer, [](const glm::dvec3& m) { return std::abs(m.x) > 1.6; });
    check(pioneer_rtg.second > 0 && pioneer_rtg.first > 0.3, "the Pioneer 10 generators reach out of the shade", pioneer_rtg.first);

    // Cassini: 4 m wide across the dish and about 6.7 m high [CLK], the dish shading the bus.
    Scene saturn_scene = load_scene_or_die("saturn.toml");
    const ShapeModel& cassini = *saturn_scene.bodies[static_cast<size_t>(saturn_scene.find("Cassini"))].shape;
    double cassini_dish = 0.0, cassini_top = -1e9;
    for (const glm::vec3& p : cassini.positions) {
        const glm::dvec3 m = glm::dvec3(p) * 1000.0;
        cassini_top = std::max(cassini_top, m.z);
        if (const double r = std::hypot(m.x, m.y); m.z > 1.6 && r < 2.5) {   // not the long antennas
            cassini_dish = std::max(cassini_dish, r);
        }
    }
    check(std::abs(2.0 * cassini_dish - 4.0) < 0.1, "the Cassini dish is 4 m across", 2.0 * cassini_dish);
    check(std::abs(cassini_top - 3.39) < 0.1, "the Cassini model is to scale", cassini_top);
    const auto cassini_bus =
        lit_share(cassini, [](const glm::dvec3& m) { return m.z < 0.0 && std::hypot(m.x, m.y) < 1.2; });
    check(cassini_bus.second > 0 && cassini_bus.first < 0.05, "the dish shades the Cassini bus", cassini_bus.first);

    // Jupiter's scene: Juno's panels span about 20 m and face the Sun, with no shade but their
    // own; Galileo is 5.3 m from the top of its antenna to the bottom of its probe [GLF].
    Scene jupiter_scene = load_scene_or_die("jupiter.toml");
    const ShapeModel& juno = *jupiter_scene.bodies[static_cast<size_t>(jupiter_scene.find("Juno"))].shape;
    glm::dvec3 juno_low(1e9), juno_high(-1e9);
    for (const glm::vec3& p : juno.positions) {
        juno_low = glm::min(juno_low, glm::dvec3(p) * 1000.0);
        juno_high = glm::max(juno_high, glm::dvec3(p) * 1000.0);
    }
    check(juno_high.y - juno_low.y > 18.0 && juno_high.x - juno_low.x > 17.0, "Juno's panels span about 20 m", juno_high.y - juno_low.y);
    const auto juno_panels = lit_share(juno, [](const glm::dvec3& m) { return std::abs(m.y) > 4.0 && m.z > 0.095; });
    check(juno_panels.second > 0 && juno_panels.first > 0.9, "the Sun reaches Juno's panels", juno_panels.first);
    const ShapeModel& galileo = *jupiter_scene.bodies[static_cast<size_t>(jupiter_scene.find("Galileo"))].shape;
    double galileo_low = 1e9, galileo_high = -1e9;
    for (const glm::vec3& p : galileo.positions) {
        galileo_low = std::min(galileo_low, static_cast<double>(p.z) * 1000.0);
        galileo_high = std::max(galileo_high, static_cast<double>(p.z) * 1000.0);
    }
    check(std::abs(galileo_high - galileo_low - 5.3) < 0.02, "Galileo is 5.3 m tall", galileo_high - galileo_low);
    const auto galileo_bus =
        lit_share(galileo, [](const glm::dvec3& m) { return m.z < 0.0 && std::hypot(m.x, m.y) < 1.0; });
    check(galileo_bus.second > 0 && galileo_bus.first < 0.05, "the sunshield shades the Galileo bus", galileo_bus.first);

    // Hubble: 4.3 m wide and about 13 m long [HBN] without its open door; the Sun reaches the
    // front of its solar panels, which stand clear of the tube (the layers behind them are shaded).
    Scene hubble_scene = load_scene_or_die("space_telescopes.toml");
    const ShapeModel& hubble = *hubble_scene.bodies[static_cast<size_t>(hubble_scene.find("Hubble"))].shape;
    double hubble_width = 0.0, hubble_low = 1e9, hubble_high = -1e9;
    for (const glm::vec3& p : hubble.positions) {
        const glm::dvec3 m = glm::dvec3(p) * 1000.0;
        hubble_low = std::min(hubble_low, m.y);
        if (m.y < 7.0) {   // the door is open above that
            hubble_high = std::max(hubble_high, m.y);
        }
        if (m.y < -4.0 && std::hypot(m.x, m.z) < 3.0) {   // the aft shroud, not the panels' ends
            hubble_width = std::max(hubble_width, 2.0 * std::abs(m.x));
        }
    }
    check(std::abs(hubble_width - 4.3) < 0.15, "Hubble is 4.3 m wide", hubble_width);
    check(std::abs(hubble_high - hubble_low - 13.0) < 0.5, "Hubble is 13 m long", hubble_high - hubble_low);
    const auto hubble_panels =
        lit_share(hubble, [](const glm::dvec3& m) { return std::abs(m.x) > 2.6 && std::abs(m.y) < 2.0 && m.z > 0.1 && m.z < 0.3; });
    check(hubble_panels.second > 0 && hubble_panels.first > 0.8, "the Sun reaches Hubble's solar panels", hubble_panels.first);

    // The rest, to scale. The widest span in the plane across the Sun line (the panels, on most),
    // or the length along it, against the published or the model's own figure.
    struct Span {
        const char* file;
        const char* craft;
        double metres;
        double tolerance;
    };
    const Span spans[] = {{"space_telescopes.toml", "Chandra", 19.5, 0.3},   // the model's arrays, in metres
                          {"space_telescopes.toml", "Roman", 12.7, 0.3},     // [RMS]
                          {"space_telescopes.toml", "TESS", 3.9, 0.1},       // Horizons' notes
                          {"space_telescopes.toml", "WMAP", 5.1, 0.1},       // [WMF]
                          {"mars.toml", "MAVEN", 13.1, 0.3},                 // the model's panels, in metres
                          {"mars.toml", "MRO", 13.2, 0.5},                   // 13.6 m in Horizons' notes
                          {"mars.toml", "TGO", 17.5, 0.2},                   // Horizons' notes
                          {"mars.toml", "Europa Clipper", 30.5, 0.2},        // Horizons' notes
                          {"didymos.toml", "DART", 19.0, 0.1},               // the scene's figure
                          {"rosetta.toml", "Rosetta", 32.0, 0.2}};           // [ROS]: the scene's 32 m
    for (const Span& span : spans) {
        Scene s = load_scene_or_die(span.file);
        const ShapeModel& model = *s.bodies[static_cast<size_t>(s.find(span.craft))].shape;
        glm::dvec3 low(1e9), high(-1e9);
        for (const glm::vec3& p : model.positions) {
            low = glm::min(low, glm::dvec3(p) * 1000.0);
            high = glm::max(high, glm::dvec3(p) * 1000.0);
        }
        const double across = std::max(high.x - low.x, high.y - low.y);
        check(std::abs(across - span.metres) < span.tolerance, (std::string(span.craft) + " is to scale").c_str(), across);
    }
    // LADEE: 1.85 m across its body [LDF], which is on the body's y axis here.
    Scene moon_scene = load_scene_or_die("earth_moon.toml");
    const ShapeModel& ladee = *moon_scene.bodies[static_cast<size_t>(moon_scene.find("LADEE"))].shape;
    double ladee_width = 0.0;
    for (const glm::vec3& p : ladee.positions) {
        if (std::abs(p.y) * 1000.0 < 0.8) {   // the body, not what stands on its ends
            ladee_width = std::max(ladee_width, 2.0 * std::hypot(p.x, p.z) * 1000.0 * std::cos(22.5 * kDegToRad));
        }
    }
    check(ladee_width > 1.6 && ladee_width < 2.6, "LADEE's body is about 1.85 m across", ladee_width);

    // What keeps to the shade: WMAP's telescopes behind its array, Spitzer's behind its panels and
    // shield, MESSENGER's bus behind its sunshade, Rosetta's and MRO's panels in the Sun.
    const ShapeModel& wmap = *hubble_scene.bodies[static_cast<size_t>(hubble_scene.find("WMAP"))].shape;
    const auto wmap_scope = lit_share(wmap, [](const glm::dvec3& m) { return m.z < -1.3; });
    check(wmap_scope.second > 0 && wmap_scope.first < 0.05, "the array shades WMAP's telescopes", wmap_scope.first);
    const ShapeModel& spitzer = *hubble_scene.bodies[static_cast<size_t>(hubble_scene.find("Spitzer"))].shape;
    const auto spitzer_scope = lit_share(spitzer, [](const glm::dvec3& m) { return m.z < 0.0 && std::abs(m.x) < 0.5 && std::abs(m.y) < 1.5; });
    check(spitzer_scope.second > 0 && spitzer_scope.first < 0.05, "the shield shades Spitzer's telescope", spitzer_scope.first);
    Scene mercury_scene = load_scene_or_die("mercury.toml");
    const ShapeModel& messenger = *mercury_scene.bodies[static_cast<size_t>(mercury_scene.find("MESSENGER"))].shape;
    const auto messenger_bus = lit_share(messenger, [](const glm::dvec3& m) { return m.z < 0.3 && std::hypot(m.x, m.y) < 0.6; });
    check(messenger_bus.second > 0 && messenger_bus.first < 0.05, "the sunshade shades MESSENGER's bus", messenger_bus.first);

    // TGO and Europa Clipper: the Sun reaches the cell side of the arrays.
    Scene mars_scene = load_scene_or_die("mars.toml");
    const ShapeModel& tgo = *mars_scene.bodies[static_cast<size_t>(mars_scene.find("TGO"))].shape;
    const auto tgo_panels = lit_share(tgo, [](const glm::dvec3& m) { return std::abs(m.x) > 3.0 && m.z > 0.015 && m.z < 0.05; });
    check(tgo_panels.second > 0 && tgo_panels.first > 0.8, "the Sun reaches TGO's solar panels", tgo_panels.first);
    const ShapeModel& clipper = *mars_scene.bodies[static_cast<size_t>(mars_scene.find("Europa Clipper"))].shape;
    const auto clipper_arrays = lit_share(clipper, [](const glm::dvec3& m) { return std::abs(m.x) > 2.5 && m.z > -0.001 && m.z < 0.05; });
    check(clipper_arrays.second > 0 && clipper_arrays.first > 0.8, "the Sun reaches Europa Clipper's arrays", clipper_arrays.first);

    // DART: scaled to 19 m across, its bus comes out at APL's 1.2 m (model's +x side, the body's -x:
    // the other side has a fitting on it). The cells are on the sunward face of both arrays, which
    // sit at different heights.
    Scene didymos_scene = load_scene_or_die("didymos.toml");
    const ShapeModel& dart = *didymos_scene.bodies[static_cast<size_t>(didymos_scene.find("DART"))].shape;
    double dart_bus = 0.0;
    for (const glm::vec3& p : dart.positions) {
        const glm::dvec3 m = glm::dvec3(p) * 1000.0;
        if (m.x < 0.0 && m.x > -1.0 && std::abs(m.y) < 0.4 && std::abs(m.z) < 0.55) {   // not the array's far end
            dart_bus = std::max(dart_bus, -m.x);
        }
    }
    check(std::abs(2.0 * dart_bus - 1.2) < 0.05, "DART's bus is 1.2 m wide", 2.0 * dart_bus);
    const auto dart_cells = lit_share(dart, [](const glm::dvec3& m) {
        return (m.x < -2.0 && m.z > 0.645 && m.z < 0.67) || (m.x > 2.0 && m.z > -0.78 && m.z < -0.75);
    });
    check(dart_cells.second > 0 && dart_cells.first > 0.9, "the Sun reaches DART's solar arrays", dart_cells.first);

    // Huygens: 2.7 m across, its nose (+z) at Titan, and nothing baked: it does not keep to the Sun.
    const Body& huygens_body = saturn_scene.bodies[static_cast<size_t>(saturn_scene.find("Huygens"))];
    const double t_huygens = tdb_from_julian_year(2005.03);
    const glm::dvec3 to_titan = scene_direction(saturn_scene, saturn_scene.find("Huygens"), saturn_scene.find("Titan"), t_huygens);
    check(glm::dot(saturn_scene.body_to_icrf(huygens_body, t_huygens)[2], to_titan) > 1.0 - 1e-12, "Huygens points its shield at Titan");
    double huygens_width = 0.0, huygens_nose = 0.0;
    bool huygens_unbaked = true;
    for (size_t k = 0; k < huygens_body.shape->positions.size(); ++k) {
        const glm::dvec3 m = glm::dvec3(huygens_body.shape->positions[k]) * 1000.0;
        huygens_width = std::max(huygens_width, 2.0 * std::hypot(m.x, m.y));
        if (m.z > 0.3) {
            huygens_nose = std::max(huygens_nose, std::hypot(m.x, m.y));
        }
        huygens_unbaked = huygens_unbaked && huygens_body.shape->albedo[k] == 1.0f;
    }
    check(std::abs(huygens_width - 2.7) < 0.03, "Huygens is 2.7 m across", huygens_width);
    check(huygens_nose < 0.6, "the Huygens shield narrows to its nose on +z", huygens_nose);
    check(huygens_unbaked, "Huygens has no baked sunlight");

    // Gaia Sky's models. Gaia: the sunshield 10.2 m across [GSS], lit underneath, the payload above it
    // in its shade. Euclid: 3.7 m across [EUO] (the model 3.5 m), the telescope behind the sunshield.
    const ShapeModel& gaia = *telescopes_scene_model(hubble_scene, "Gaia");
    glm::dvec3 gaia_low(1e9), gaia_high(-1e9);
    for (const glm::vec3& p : gaia.positions) {
        gaia_low = glm::min(gaia_low, glm::dvec3(p) * 1000.0);
        gaia_high = glm::max(gaia_high, glm::dvec3(p) * 1000.0);
    }
    check(std::abs(gaia_high.x - gaia_low.x - 10.2) < 0.1, "Gaia's sunshield is 10.2 m across", gaia_high.x - gaia_low.x);
    const auto gaia_payload = lit_share(gaia, [](const glm::dvec3& m) { return m.z < -0.5; });
    check(gaia_payload.second > 0 && gaia_payload.first < 0.05, "the sunshield shades Gaia's payload", gaia_payload.first);
    const auto gaia_cells = lit_share(gaia, [](const glm::dvec3& m) { return m.z > 0.2; });   // the lowest layer
    check(gaia_cells.second > 0 && gaia_cells.first > 0.8, "the Sun reaches the underside of Gaia's sunshield", gaia_cells.first);
    const ShapeModel& euclid = *telescopes_scene_model(hubble_scene, "Euclid");
    double euclid_low = 1e9, euclid_high = -1e9;
    for (const glm::vec3& p : euclid.positions) {
        euclid_low = std::min(euclid_low, static_cast<double>(p.x) * 1000.0);
        euclid_high = std::max(euclid_high, static_cast<double>(p.x) * 1000.0);
    }
    check(std::abs(euclid_high - euclid_low - 3.7) < 0.3, "Euclid is about 3.7 m across", euclid_high - euclid_low);
    const auto euclid_scope = lit_share(euclid, [](const glm::dvec3& m) { return m.z < -0.5 && m.y > 0.0; });
    check(euclid_scope.second > 0 && euclid_scope.first < 0.05, "the sunshield shades Euclid's telescope", euclid_scope.first);

    // Orion [MDLO]: the crew module 5 m across on +z, the four arrays about 19 m tip to tip across the
    // Sun line, their cells on the sunward side in the Sun (bake and attitude: the cases above).
    Scene earth_moon = load_scene_or_die("earth_moon.toml");
    const Body& artemis = earth_moon.bodies[static_cast<size_t>(earth_moon.find("Artemis II"))];
    check(artemis.shape != nullptr && !artemis.texture.empty(), "Artemis II has a model");
    if (artemis.shape) {
        double crew_width = 0.0, array_span = 0.0, stack_low = 1e9, stack_high = -1e9;
        for (const glm::vec3& p : artemis.shape->positions) {
            const glm::dvec3 m = glm::dvec3(p) * 1000.0;
            const double r = std::hypot(m.x, m.y);
            if (m.z > 0.5) {
                crew_width = std::max(crew_width, 2.0 * r);
            }
            array_span = std::max(array_span, 2.0 * r);
            stack_low = std::min(stack_low, m.z);
            stack_high = std::max(stack_high, m.z);
        }
        check(std::abs(crew_width - 5.0) < 0.1, "Orion's crew module is 5 m across", crew_width);
        check(std::abs(array_span - 19.0) < 0.3, "Orion's arrays are about 19 m across", array_span);
        check(std::abs(stack_high - stack_low - 7.4) < 0.1, "the Orion stack is about 7.4 m long", stack_high - stack_low);
        const auto cells = lit_share(*artemis.shape, [](const glm::dvec3& m) {
            return std::hypot(m.x, m.y) > 3.5 && m.z > -1.48 && m.z < -1.4;
        });
        check(cells.second > 0 && cells.first > 0.9, "the Sun reaches Orion's solar cells", cells.first);
        const auto rear = lit_share(*artemis.shape, [](const glm::dvec3& m) {
            return std::hypot(m.x, m.y) > 3.5 && m.z < -1.5;
        });
        check(rear.second > 0 && rear.first < 0.1, "the backs of Orion's arrays are in their own shade", rear.first);
        const Body& artemis1 = earth_moon.bodies[static_cast<size_t>(earth_moon.find("Artemis I"))];
        check(artemis1.shape != nullptr && artemis1.shape->positions.size() == artemis.shape->positions.size(),
              "Artemis I has the Orion model too");
    }

    // ISEE-3 [MDLI]: the drum 1.58 m tall (the origin halfway up; thruster pods hang below it) about the spin axis (the body's y, toward the ecliptic north
    // [SPIN]), the four wires 49 m beyond it in the plane across, the axial antennas 7.01 m; a side to
    // the Sun; nothing baked, as it spins.
    Scene isee3_scene = load_scene_or_die("isee3.toml");
    const int isee3_index = isee3_scene.find("ISEE-3");
    const Body& isee3 = isee3_scene.bodies[static_cast<size_t>(isee3_index)];
    check(isee3.shape != nullptr && !isee3.texture.empty(), "ISEE-3 has a model");
    if (isee3.shape) {
        double drum_high = -1e9, drum_radius = 0.0, wire_reach = 0.0, axial_reach = 0.0;
        bool isee3_unbaked = true;
        for (size_t k = 0; k < isee3.shape->positions.size(); ++k) {
            const glm::dvec3 m = glm::dvec3(isee3.shape->positions[k]) * 1000.0;
            const double r = std::hypot(m.x, m.z);
            if (r > 0.5 && r < 0.9) {
                drum_high = std::max(drum_high, m.y);
                drum_radius = std::max(drum_radius, r);
            }
            wire_reach = std::max(wire_reach, r);
            axial_reach = std::max(axial_reach, std::abs(m.y));
            isee3_unbaked = isee3_unbaked && isee3.shape->albedo[k] == 1.0f;
        }
        check(std::abs(2.0 * drum_high - 1.58) < 0.03, "the ISEE-3 drum is 1.58 m tall", 2.0 * drum_high);
        check(std::abs(2.0 * drum_radius - 1.77) < 0.03, "the ISEE-3 drum is 1.77 m across", 2.0 * drum_radius);
        check(std::abs(wire_reach - drum_radius - 49.0) < 0.1, "ISEE-3's wires reach 49 m beyond the drum", wire_reach - drum_radius);
        check(std::abs(axial_reach - drum_high - 7.01) < 0.05, "ISEE-3's axial antennas are 7.01 m", axial_reach - drum_high);
        check(isee3_unbaked, "ISEE-3 has no baked sunlight");
        const double t_isee3 = tdb_from_julian_year(1980.0);
        const glm::dmat3 axes = isee3_scene.body_to_icrf(isee3, t_isee3);
        const glm::dvec3 self = isee3_scene.icrf_state_at(isee3_index, t_isee3).position;
        const glm::dvec3 sun = glm::normalize(isee3_scene.sun_icrf_state_at(t_isee3).position - self);
        check(glm::dot(axes[2], sun) > 1.0 - 1e-12, "ISEE-3 turns a side to the Sun");
        check(glm::dot(axes[1], ecliptic_pole_icrf()) > 0.999, "ISEE-3's spin axis points to the ecliptic north",
              glm::dot(axes[1], ecliptic_pole_icrf()));
    }

    // Webb: the shield about 21 m long, the telescope on the side away from the Sun, in its shade.
    Scene telescopes = load_scene_or_die("space_telescopes.toml");
    const ShapeModel& jwst = *telescopes.bodies[static_cast<size_t>(telescopes.find("JWST"))].shape;
    glm::dvec3 low(1e9), high(-1e9);
    for (const glm::vec3& p : jwst.positions) {
        low = glm::min(low, glm::dvec3(p) * 1000.0);
        high = glm::max(high, glm::dvec3(p) * 1000.0);
    }
    check(std::abs(high.x - low.x - 21.0) < 0.5, "the Webb sunshield is about 21 m long", high.x - low.x);
    check(low.z < -7.0 && high.z < 4.0, "the Webb telescope is on the side away from the Sun", low.z);
    const auto telescope = lit_share(jwst, [](const glm::dvec3& m) { return m.z < -3.5; });
    check(telescope.second > 0 && telescope.first < 0.05, "the sunshield shades the Webb telescope", telescope.first);
    const auto shield = lit_share(jwst, [](const glm::dvec3& m) { return m.z > 2.0; });
    check(shield.second > 0 && shield.first > 0.5, "the Sun reaches Webb's sunward side", shield.first);

    // Parker: the shield's front in the Sun and the bus behind it in its shade, the solar panels'
    // front outside the shadow (the shield is about 2.3 m across; layers behind them are shaded).
    Scene parker = load_scene_or_die("parker.toml");
    const ShapeModel& probe = *parker.bodies[static_cast<size_t>(parker.find("Parker Solar Probe"))].shape;
    const auto front = lit_share(probe, [](const glm::dvec3& m) { return m.z > 1.47; });
    check(front.second > 0 && front.first > 0.95, "the Sun reaches Parker's heat shield", front.first);
    const auto behind = lit_share(probe, [](const glm::dvec3& m) { return m.z < 1.0 && std::hypot(m.x, m.y) < 0.8; });
    check(behind.second > 0 && behind.first < 0.05, "the heat shield shades Parker's bus", behind.first);
    const auto panels = lit_share(probe, [](const glm::dvec3& m) { return std::abs(m.x) > 1.5 && std::abs(m.y) < 0.6 && m.z > 0.02; });
    check(panels.second > 0 && panels.first > 0.8, "Parker's solar panels reach out of the shadow", panels.first);
}

void run_shape_tests()
{
    test_proteus_shape();
    test_arrokoth_shape();
    test_mars_moon_shapes();
    test_vesta_shape();
    test_pallas_shape();
    test_jupiter_saturn_small_moon_shapes();
    test_hyperion_spin_state();
    test_halley_shape();
    test_shape_occlusion();
    test_didymos_shapes();
    test_endurance_shape();
    test_spacecraft_models();
}
