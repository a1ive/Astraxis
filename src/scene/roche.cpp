#include "scene/roche.hpp"

#include "core/math.hpp"
#include "scene/light_curve.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

double potential(double q, const glm::dvec3& p)
{
    const double r1 = glm::length(p);
    const double r2 = glm::length(p - glm::dvec3(1.0, 0.0, 0.0));
    return 1.0 / r1 + q * (1.0 / r2 - p.x) + 0.5 * (1.0 + q) * (p.x * p.x + p.y * p.y);
}

glm::dvec3 potential_gradient(double q, const glm::dvec3& p)
{
    const glm::dvec3 d2 = p - glm::dvec3(1.0, 0.0, 0.0);
    const double r1 = glm::length(p);
    const double r2 = glm::length(d2);
    const glm::dvec3 g = -p / (r1 * r1 * r1) - q * d2 / (r2 * r2 * r2);
    return g + glm::dvec3(-q + (1.0 + q) * p.x, (1.0 + q) * p.y, 0.0);
}

// Surface radius along a unit direction: Omega falls from infinity at the
// center to the critical value on the surface, and stays below it out to the
// distance of L1 (the farthest point of the lobe).
double surface_radius(double q, double l1, double critical, const glm::dvec3& direction)
{
    double lo = 1e-3;
    double hi = l1;
    for (int i = 0; i < 64; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (potential(q, mid * direction) > critical) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

} // namespace

double roche_l1(double q)
{
    // dOmega/dx on the axis rises from -infinity at the star to +infinity at the companion.
    double lo = 1e-6;
    double hi = 1.0 - 1e-6;
    for (int i = 0; i < 100; ++i) {
        const double x = 0.5 * (lo + hi);
        const double slope = -1.0 / (x * x) + q * (1.0 / ((1.0 - x) * (1.0 - x)) - 1.0) + (1.0 + q) * x;
        if (slope < 0.0) {
            lo = x;
        } else {
            hi = x;
        }
    }
    return 0.5 * (lo + hi);
}

double roche_radius(double q, const glm::dvec3& direction)
{
    const double l1 = roche_l1(q);
    return surface_radius(q, l1, potential(q, glm::dvec3(l1, 0.0, 0.0)), direction);
}

RocheSurface make_roche_surface(const RocheStar& star, int rings, int segments)
{
    const double q = star.mass_ratio;
    const double a = star.separation_km;
    const double l1 = roche_l1(q);
    const double critical = potential(q, glm::dvec3(l1, 0.0, 0.0));

    // Vertices: L1 (theta = 0 from the x axis), rings - 1 circles, the far
    // side (theta = pi). Directions (cos theta, sin theta cos phi, sin theta sin phi).
    RocheSurface out;
    ShapeModel& mesh = out.mesh;
    std::vector<glm::dvec3> unit_positions; // units of a
    std::vector<double> gravity;
    unit_positions.push_back(glm::dvec3(l1, 0.0, 0.0));
    for (int k = 1; k < rings; ++k) {
        const double theta = kPi * k / rings;
        for (int m = 0; m < segments; ++m) {
            const double phi = kTwoPi * m / segments;
            const glm::dvec3 d(std::cos(theta), std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi));
            unit_positions.push_back(surface_radius(q, l1, critical, d) * d);
        }
    }
    unit_positions.push_back(surface_radius(q, l1, critical, glm::dvec3(-1.0, 0.0, 0.0)) *
                             glm::dvec3(-1.0, 0.0, 0.0));
    const auto ring_vertex = [&](int k, int m) { return static_cast<uint32_t>(1 + (k - 1) * segments + m % segments); };
    const auto back = static_cast<uint32_t>(unit_positions.size() - 1);

    // Counter-clockwise seen from outside (d/dtheta x d/dphi points outward).
    for (int m = 0; m < segments; ++m) {
        mesh.indices.insert(mesh.indices.end(), {0u, ring_vertex(1, m), ring_vertex(1, m + 1)});
    }
    for (int k = 1; k + 1 < rings; ++k) {
        for (int m = 0; m < segments; ++m) {
            const uint32_t v00 = ring_vertex(k, m);
            const uint32_t v10 = ring_vertex(k + 1, m);
            const uint32_t v11 = ring_vertex(k + 1, m + 1);
            const uint32_t v01 = ring_vertex(k, m + 1);
            mesh.indices.insert(mesh.indices.end(), {v00, v10, v11, v00, v11, v01});
        }
    }
    for (int m = 0; m < segments; ++m) {
        mesh.indices.insert(mesh.indices.end(), {ring_vertex(rings - 1, m), back, ring_vertex(rings - 1, m + 1)});
    }

    // Normals and gravity from the potential (at L1 gravity vanishes: the
    // cusp's normal is taken along +x).
    std::vector<glm::dvec3> normals;
    for (const glm::dvec3& p : unit_positions) {
        const glm::dvec3 g = potential_gradient(q, p);
        const double magnitude = glm::length(g);
        gravity.push_back(magnitude);
        normals.push_back(magnitude > 1e-9 ? -g / magnitude : glm::dvec3(1.0, 0.0, 0.0));
    }

    // Gravity darkening, scaled to the effective temperature's luminosity.
    const size_t n = unit_positions.size();
    std::vector<double> t4(n);
    for (size_t v = 0; v < n; ++v) {
        t4[v] = std::pow(gravity[v], star.gravity_darkening);
    }
    double area = 0.0;
    double weighted = 0.0;
    double volume = 0.0;
    for (size_t i = 0; i < mesh.indices.size(); i += 3) {
        const glm::dvec3& p0 = unit_positions[mesh.indices[i]];
        const glm::dvec3& p1 = unit_positions[mesh.indices[i + 1]];
        const glm::dvec3& p2 = unit_positions[mesh.indices[i + 2]];
        const double facet_area = 0.5 * glm::length(glm::cross(p1 - p0, p2 - p0));
        area += facet_area;
        weighted += facet_area * (t4[mesh.indices[i]] + t4[mesh.indices[i + 1]] + t4[mesh.indices[i + 2]]) / 3.0;
        volume += glm::dot(p0, glm::cross(p1, p2)) / 6.0;
    }
    const double t_eff4 = std::pow(star.temperature_k, 4.0);
    const double scale = t_eff4 * area / weighted;
    const double companion_t4 = std::pow(star.companion_temperature_k, 4.0);
    const double companion_radius = star.companion_radius_km / a;
    for (size_t v = 0; v < n; ++v) {
        t4[v] *= scale;
        const glm::dvec3 s = glm::dvec3(1.0, 0.0, 0.0) - unit_positions[v];
        const double distance = glm::length(s);
        const double cos_gamma = std::max(glm::dot(normals[v], s) / distance, 0.0);
        t4[v] += star.albedo * companion_t4 * companion_radius * companion_radius / (distance * distance) * cos_gamma;
    }

    const double reference = planck_shape(star.temperature_k, kRocheRadianceWavelengthNm);
    mesh.positions.reserve(n);
    mesh.normals.reserve(n);
    mesh.albedo.reserve(n);
    out.vertex_temperature_k.reserve(n);
    for (size_t v = 0; v < n; ++v) {
        const double temperature = std::max(std::pow(t4[v], 0.25), 1.0);
        mesh.positions.push_back(glm::vec3(unit_positions[v] * a));
        mesh.normals.push_back(glm::vec3(normals[v]));
        mesh.albedo.push_back(static_cast<float>(planck_shape(temperature, kRocheRadianceWavelengthNm) / reference));
        out.vertex_temperature_k.push_back(temperature);
    }
    for (size_t i = 0; i < mesh.indices.size(); i += 3) {
        const glm::dvec3 p0 = unit_positions[mesh.indices[i]] * a;
        const glm::dvec3 p1 = unit_positions[mesh.indices[i + 1]] * a;
        const glm::dvec3 p2 = unit_positions[mesh.indices[i + 2]] * a;
        const glm::dvec3 cross = glm::cross(p1 - p0, p2 - p0);
        RocheSurface::Facet f;
        f.center = (p0 + p1 + p2) / 3.0;
        f.normal = glm::normalize(cross);
        f.area_km2 = 0.5 * glm::length(cross);
        f.temperature_k =
            std::max(std::pow((t4[mesh.indices[i]] + t4[mesh.indices[i + 1]] + t4[mesh.indices[i + 2]]) / 3.0, 0.25),
                     1.0);
        f.radiance = planck_shape(f.temperature_k, kRocheRadianceWavelengthNm);
        out.facets.push_back(f);
    }

    out.l1_km = l1 * a;
    out.back_km = glm::length(unit_positions.back()) * a;
    out.side_km = surface_radius(q, l1, critical, glm::dvec3(0.0, 1.0, 0.0)) * a;
    out.pole_km = surface_radius(q, l1, critical, glm::dvec3(0.0, 0.0, 1.0)) * a;
    out.volume_radius_km = std::cbrt(3.0 * volume / (4.0 * kPi)) * a;
    return out;
}

} // namespace astraxis
