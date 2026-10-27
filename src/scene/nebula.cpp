#include "scene/nebula.hpp"

#include "core/color.hpp"
#include "core/math.hpp"
#include "scene/shape_model.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

// A young shell shows at most this many times brighter than at the
// reference time (nebula_exposure).
constexpr double kNebulaMaxBrightening = 2.0;

// A patch's Gaussian size along the shell in units of the grains' mean
// spacing: they overlap enough to sum smoothly though placed at random-like
// points rather than on a grid (a pinwheel's rings use 0.5; at 0.6 the shell
// seen from inside is blotchy).
constexpr double kNebulaPatch = 0.85;

} // namespace

glm::dmat3 nebula_model_to_icrf(double ra_deg, double dec_deg, double axis_pa_deg, double axis_tilt_deg,
                                double roll_deg)
{
    const double ra = ra_deg * kDegToRad;
    const double dec = dec_deg * kDegToRad;
    const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
    const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
    const glm::dvec3 toward(-std::cos(dec) * std::cos(ra), -std::cos(dec) * std::sin(ra), -std::sin(dec));

    // The model's +z axis on the sky, and the direction perpendicular to it on
    // the observer's side (the line of sight's part across the axis).
    const double pa = axis_pa_deg * kDegToRad;
    const double tilt = axis_tilt_deg * kDegToRad;
    const glm::dvec3 axis = std::sin(tilt) * (north * std::cos(pa) + east * std::sin(pa)) + std::cos(tilt) * toward;
    const glm::dvec3 facing = std::abs(std::sin(tilt)) > 1e-12
                                  ? glm::normalize(toward - glm::dot(toward, axis) * axis)
                                  : north; // looking down the axis: any direction across it
    const glm::dvec3 side = glm::cross(axis, facing);

    // The model's direction at azimuth roll faces the observer: its +x axis is
    // roll behind that.
    const double roll = roll_deg * kDegToRad;
    const glm::dvec3 x = facing * std::cos(roll) - side * std::sin(roll);
    const glm::dvec3 y = facing * std::sin(roll) + side * std::cos(roll);
    return glm::dmat3(x, y, axis);
}

void build_nebula_grains(const ShapeModel& shape, const glm::dvec3& model_center, double km_per_unit,
                         double polar_scale, const glm::dmat3& model_to_icrf, int count, Nebula& nebula)
{
    nebula.grains.clear();
    nebula.extent_km = 0.0;
    nebula.mean_distance_km = 0.0;
    // Model units to km, the z axis stretched by polar_scale; normals by the
    // inverse transpose.
    const glm::dvec3 stretch(km_per_unit, km_per_unit, km_per_unit * polar_scale);
    const glm::dvec3 normal_stretch(1.0, 1.0, 1.0 / polar_scale);
    const size_t triangles = shape.indices.size() / 3;
    std::vector<double> area(triangles);
    double total = 0.0;
    for (size_t k = 0; k < triangles; ++k) {
        const glm::dvec3 a(shape.positions[shape.indices[3 * k]]);
        const glm::dvec3 b(shape.positions[shape.indices[3 * k + 1]]);
        const glm::dvec3 c(shape.positions[shape.indices[3 * k + 2]]);
        area[k] = 0.5 * glm::length(glm::cross((b - a) * stretch, (c - a) * stretch));
        total += area[k];
    }
    if (total <= 0.0 || count <= 0) {
        return;
    }
    nebula.grains.reserve(static_cast<size_t>(count) + 1);
    // Grains per triangle in proportion to its area, carrying the remainder to
    // the next (so the density is even down to single triangles); within
    // each, the R2 sequence (Roberts 2018) folded into the triangle.
    constexpr double kR2x = 0.7548776662466927;
    constexpr double kR2y = 0.5698402909980532;
    const double per_area = count / total;
    double carry = 0.0;
    size_t sequence = 0;
    for (size_t k = 0; k < triangles; ++k) {
        carry += area[k] * per_area;
        const int n = static_cast<int>(std::floor(carry));
        carry -= n;
        const uint32_t ia = shape.indices[3 * k];
        const uint32_t ib = shape.indices[3 * k + 1];
        const uint32_t ic = shape.indices[3 * k + 2];
        for (int j = 0; j < n; ++j, ++sequence) {
            double u = std::fmod(0.5 + kR2x * static_cast<double>(sequence), 1.0);
            double v = std::fmod(0.5 + kR2y * static_cast<double>(sequence), 1.0);
            if (u + v > 1.0) {
                u = 1.0 - u;
                v = 1.0 - v;
            }
            const double w = 1.0 - u - v;
            const glm::dvec3 p = w * glm::dvec3(shape.positions[ia]) + u * glm::dvec3(shape.positions[ib]) +
                                 v * glm::dvec3(shape.positions[ic]);
            glm::dvec3 n_model = w * glm::dvec3(shape.normals[ia]) + u * glm::dvec3(shape.normals[ib]) +
                                 v * glm::dvec3(shape.normals[ic]);
            if (glm::length(n_model) < 1e-12) {
                const glm::dvec3 a(shape.positions[ia]);
                n_model = glm::cross(glm::dvec3(shape.positions[ib]) - a, glm::dvec3(shape.positions[ic]) - a);
            }
            NebulaGrain g;
            g.position = glm::vec3(model_to_icrf * ((p - model_center) * stretch));
            g.normal = glm::vec3(glm::normalize(model_to_icrf * (n_model * normal_stretch)));
            nebula.grains.push_back(g);
            const double d = glm::length(glm::dvec3(g.position));
            nebula.extent_km = std::max(nebula.extent_km, d);
            nebula.mean_distance_km += d;
        }
    }
    if (nebula.grains.empty()) {
        return;
    }
    nebula.mean_distance_km /= static_cast<double>(nebula.grains.size());
    // The area is in km^2 already.
    nebula.patch_sigma_km = kNebulaPatch * std::sqrt(total / static_cast<double>(nebula.grains.size()));
}

double nebula_scale(const Nebula& nebula, double t_tdb)
{
    return std::max(0.0, (t_tdb - nebula.t_ejection_tdb) / (nebula.t_reference_tdb - nebula.t_ejection_tdb));
}

double nebula_extent_km(const Nebula& nebula, double t_tdb)
{
    return nebula.extent_km * nebula_scale(nebula, t_tdb);
}

double henyey_greenstein(double g, double cos_theta)
{
    const double d = 1.0 + g * g - 2.0 * g * cos_theta;
    return (1.0 - g * g) / (d * std::sqrt(d));
}

double nebula_exposure(const Nebula& nebula, double t_tdb)
{
    const double s = nebula_scale(nebula, t_tdb);
    if (s <= 0.0) {
        return 1.0;
    }
    // Emission lacks the additional column dilution of scattered light.
    const double power = nebula.light == NebulaLight::Emission ? 2.0 : 4.0;
    return std::min(1.0, kNebulaMaxBrightening * std::pow(s, power * nebula.display_gamma));
}

double nebula_dust_temperature(const Nebula& nebula, const NebulaLayer& layer, double distance_km)
{
    const double reference = layer.radius_share * nebula.mean_distance_km;
    if (reference <= 0.0 || distance_km <= 0.0) {
        return layer.temperature_k;
    }
    return layer.temperature_k * std::pow(distance_km / reference, -nebula.temperature_index);
}

double nebula_infrared_emission(double kelvin)
{
    const double t = kelvin / 100.0;
    return t * t * t * t * blackbody_band_share(kelvin, kNebulaInfraredLoNm, kNebulaInfraredHiNm);
}

} // namespace astraxis
