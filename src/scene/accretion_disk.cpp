#include "scene/accretion_disk.hpp"

#include "core/math.hpp"
#include "scene/light_curve.hpp"
#include "scene/roche.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>

namespace astraxis {

namespace {

constexpr int kSegments = 128;
constexpr int kFaceRings = 24;
constexpr int kBevelRings = 5;
constexpr int kWallRings = 12;
constexpr int kPhaseSamples = 256;
constexpr int kImageSize = 320;

struct ProjectedVertex {
    glm::vec2 pixel;
    float depth;
    float intensity;
};

struct PhotometryVertex {
    glm::vec3 position;
    glm::vec3 normal;
    float radiance;
};

std::vector<PhotometryVertex> photometry_vertices(const ShapeModel& mesh, const std::vector<double>& temperatures,
                                                  double wavelength_nm, double separation_km, bool donor)
{
    std::vector<PhotometryVertex> vertices;
    vertices.reserve(mesh.positions.size());
    const float scale = static_cast<float>(1.0 / separation_km);
    // Center the image between the stars. The donor's x and y axes are opposite the disk's.
    const glm::vec3 axes = donor ? glm::vec3(-1.0f, -1.0f, 1.0f) : glm::vec3(1.0f);
    const glm::vec3 center(donor ? 0.5f : -0.5f, 0.0f, 0.0f);
    for (size_t i = 0; i < mesh.positions.size(); ++i) {
        vertices.push_back({mesh.positions[i] * axes * scale + center, mesh.normals[i] * axes,
                            static_cast<float>(planck_shape(temperatures[i], wavelength_nm))});
    }
    return vertices;
}

float edge(const glm::vec2& a, const glm::vec2& b)
{
    return a.x * b.y - a.y * b.x;
}

// Orthographic depth buffer for the disk and donor. The closest triangle
// supplies the pixel's radiance, including where one surface eclipses the other.
class PhotometryImage {
public:
    explicit PhotometryImage(float half_span) : m_pixels_per_unit((kImageSize - 2.0f) / (2.0f * half_span)) {}

    void clear()
    {
        std::fill(m_depth.begin(), m_depth.end(), -std::numeric_limits<float>::infinity());
        std::fill(m_intensity.begin(), m_intensity.end(), 0.0f);
    }

    void draw(const std::vector<PhotometryVertex>& vertices, const std::vector<uint32_t>& indices,
              const glm::vec3& toward, float limb_darkening)
    {
        const glm::vec3 helper = std::abs(toward.z) < 0.9f ? glm::vec3(0, 0, 1) : glm::vec3(1, 0, 0);
        const glm::vec3 right = glm::normalize(glm::cross(helper, toward));
        const glm::vec3 up = glm::cross(toward, right);
        m_projected.clear();
        m_projected.reserve(vertices.size());
        for (const PhotometryVertex& v : vertices) {
            const glm::vec2 p(glm::dot(v.position, right), glm::dot(v.position, up));
            const float mu = std::clamp(glm::dot(v.normal, toward), 0.0f, 1.0f);
            m_projected.push_back({p * m_pixels_per_unit + glm::vec2(kImageSize / 2.0f - 0.5f),
                                   glm::dot(v.position, toward), v.radiance * (1.0f - limb_darkening * (1.0f - mu))});
        }
        for (size_t i = 0; i < indices.size(); i += 3) {
            triangle(m_projected[indices[i]], m_projected[indices[i + 1]], m_projected[indices[i + 2]]);
        }
    }

    double flux() const { return std::accumulate(m_intensity.begin(), m_intensity.end(), 0.0); }

private:
    void triangle(const ProjectedVertex& a, const ProjectedVertex& b, const ProjectedVertex& c)
    {
        const float area = edge(b.pixel - a.pixel, c.pixel - a.pixel);
        if (area <= 0.0f) {
            return;
        }
        const glm::vec2 lo = glm::min(a.pixel, glm::min(b.pixel, c.pixel));
        const glm::vec2 hi = glm::max(a.pixel, glm::max(b.pixel, c.pixel));
        const int x0 = std::max(0, static_cast<int>(std::ceil(lo.x)));
        const int y0 = std::max(0, static_cast<int>(std::ceil(lo.y)));
        const int x1 = std::min(kImageSize - 1, static_cast<int>(std::floor(hi.x)));
        const int y1 = std::min(kImageSize - 1, static_cast<int>(std::floor(hi.y)));
        const float inverse_area = 1.0f / area;
        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                const glm::vec2 p(static_cast<float>(x), static_cast<float>(y));
                const float wa = edge(b.pixel - p, c.pixel - p);
                const float wb = edge(c.pixel - p, a.pixel - p);
                const float wc = edge(a.pixel - p, b.pixel - p);
                if (wa < 0.0f || wb < 0.0f || wc < 0.0f) {
                    continue;
                }
                const size_t pixel = static_cast<size_t>(y * kImageSize + x);
                const float depth = (wa * a.depth + wb * b.depth + wc * c.depth) * inverse_area;
                if (depth > m_depth[pixel]) {
                    m_depth[pixel] = depth;
                    m_intensity[pixel] = (wa * a.intensity + wb * b.intensity + wc * c.intensity) * inverse_area;
                }
            }
        }
    }

    float m_pixels_per_unit;
    std::vector<float> m_depth = std::vector<float>(kImageSize * kImageSize);
    std::vector<float> m_intensity = std::vector<float>(kImageSize * kImageSize);
    std::vector<ProjectedVertex> m_projected;
};

} // namespace

AccretionDisk make_accretion_disk(double radius_km, double half_height_km, double core_radius_km,
                                double rim_temperature_k, double inner_temperature_k, double temperature_power,
                                double rim_heating_k, double rim_heating_radius_km)
{
    AccretionDisk disk;
    disk.outer_radius_km = radius_km;
    disk.half_height_km = half_height_km;
    disk.rim_temperature_k = rim_temperature_k;
    disk.inner_temperature_k = inner_temperature_k;

    ShapeModel& mesh = disk.mesh;
    const double bevel = 0.44 * half_height_km;
    const double face_radius = radius_km - bevel;
    const double wall_height = half_height_km - bevel;
    const double reference = planck_shape(rim_temperature_k, kRocheRadianceWavelengthNm);
    const double core_to_rim = std::hypot(radius_km, core_radius_km);
    const auto vertex = [&](double r, double z, double nr, double nz, double azimuth) {
        const double x = std::cos(azimuth);
        const double y = std::sin(azimuth);
        const double temperature = rim_temperature_k * std::pow(core_to_rim / std::hypot(r, core_radius_km),
                                                                temperature_power);
        const glm::dvec3 from_impact(r * x - radius_km, r * y, z);
        const double distance = glm::length(from_impact) / rim_heating_radius_km;
        const double heated = std::min(temperature, inner_temperature_k) + rim_heating_k * std::exp(-0.5 * distance * distance);
        mesh.positions.emplace_back(r * x, r * y, z);
        mesh.normals.emplace_back(nr * x, nr * y, nz);
        mesh.albedo.push_back(static_cast<float>(planck_shape(heated, kRocheRadianceWavelengthNm) / reference));
        disk.vertex_temperature_k.push_back(heated);
    };
    const auto ring = [&](double r, double z, double nr, double nz) {
        for (int i = 0; i < kSegments; ++i) {
            vertex(r, z, nr, nz, kTwoPi * i / kSegments);
        }
    };

    vertex(0.0, half_height_km, 0.0, 1.0, 0.0);
    for (int i = 1; i <= kFaceRings; ++i) {
        ring(face_radius * i / kFaceRings, half_height_km, 0.0, 1.0);
    }
    for (int i = 1; i <= kBevelRings; ++i) {
        const double angle = kPi * 0.5 * i / kBevelRings;
        ring(face_radius + bevel * std::sin(angle), wall_height + bevel * std::cos(angle),
             std::sin(angle), std::cos(angle));
    }
    for (int i = 1; i <= kWallRings; ++i) {
        ring(radius_km, wall_height * (1.0 - 2.0 * i / kWallRings), 1.0, 0.0);
    }
    for (int i = 1; i <= kBevelRings; ++i) {
        const double angle = kPi * 0.5 * i / kBevelRings;
        ring(face_radius + bevel * std::cos(angle), -wall_height - bevel * std::sin(angle),
             std::cos(angle), -std::sin(angle));
    }
    for (int i = kFaceRings - 1; i > 0; --i) {
        ring(face_radius * i / kFaceRings, -half_height_km, 0.0, -1.0);
    }
    const auto bottom = static_cast<uint32_t>(mesh.positions.size());
    const auto rings = (bottom - 1) / kSegments;
    vertex(0.0, -half_height_km, 0.0, -1.0, 0.0);

    const auto index = [](uint32_t j, uint32_t i) { return 1 + j * kSegments + i % kSegments; };
    for (uint32_t i = 0; i < kSegments; ++i) {
        mesh.indices.insert(mesh.indices.end(), {0, index(0, i), index(0, i + 1)});
        for (uint32_t j = 0; j + 1 < rings; ++j) {
            const uint32_t a = index(j, i);
            const uint32_t b = index(j + 1, i);
            const uint32_t c = index(j + 1, i + 1);
            const uint32_t d = index(j, i + 1);
            mesh.indices.insert(mesh.indices.end(), {a, b, c, a, c, d});
        }
        mesh.indices.insert(mesh.indices.end(), {index(rings - 1, i), bottom, index(rings - 1, i + 1)});
    }
    return disk;
}

std::vector<float> accretion_binary_curve(const AccretionDisk& disk, const RocheSurface& donor, double separation_km,
                                         double observer_z, double wavelength_nm, double limb_darkening)
{
    const auto disk_vertices = photometry_vertices(disk.mesh, disk.vertex_temperature_k, wavelength_nm,
                                                   separation_km, false);
    const auto donor_vertices = photometry_vertices(donor.mesh, donor.vertex_temperature_k, wavelength_nm,
                                                     separation_km, true);
    const float xy = static_cast<float>(std::sqrt(std::max(1.0 - observer_z * observer_z, 0.0)));
    const float u = static_cast<float>(limb_darkening);
    const float half_span = 0.5f + static_cast<float>(
        std::max(std::hypot(disk.outer_radius_km, disk.half_height_km), donor.l1_km) / separation_km);
    PhotometryImage image(half_span);
    std::array<double, kPhaseSamples> flux;
    for (int i = 0; i < kPhaseSamples; ++i) {
        const double azimuth = kTwoPi * i / kPhaseSamples;
        const glm::vec3 toward(xy * std::cos(azimuth), xy * std::sin(azimuth), observer_z);
        image.clear();
        image.draw(disk_vertices, disk.mesh.indices, toward, u);
        image.draw(donor_vertices, donor.mesh.indices, toward, u);
        flux[static_cast<size_t>(i)] = image.flux();
    }
    const double maximum = *std::max_element(flux.begin(), flux.end());
    std::vector<float> magnitudes;
    magnitudes.reserve(flux.size());
    for (double sample : flux) {
        magnitudes.push_back(static_cast<float>(-2.5 * std::log10(std::max(sample / maximum, 1e-12))));
    }
    return magnitudes;
}

} // namespace astraxis
