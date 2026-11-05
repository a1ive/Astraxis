#include "scene/disk.hpp"

#include "core/math.hpp"
#include "scene/nebula.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr double kInverseSqrtTwoPi = 0.3989422804014327;

// The steps along a ray, the same as disk_march in shaders/disk.hlsli: a step
// takes kVerticalStep of the distance to the midplane (of the Gaussian's sigma
// inside it) and kRadialStep of the distance to the star, so thin disks seen
// face-on are crossed in sigma-sized steps and edge-on ones in steps that shrink
// toward the star, in few samples whichever way the ray goes.
constexpr double kVerticalStep = 0.5;
constexpr double kRadialStep = 0.3;
constexpr int kMaxSteps = 48;

double march(const Disk& disk, const DiskComponent& c, DiskView view, const glm::dvec3& origin,
             const glm::dvec3& direction, double t0, double t1)
{
    const glm::dmat3 to_plane = glm::transpose(c.to_icrf);
    const glm::dvec3 o = to_plane * origin;
    const glm::dvec3 d = to_plane * direction;
    double sum = 0.0;
    double t = t0;
    for (int i = 0; i < kMaxSteps && t < t1; ++i) {
        const glm::dvec3 p = o + t * d;
        const double r = std::hypot(p.x, p.y);
        const double sigma = c.sigma_au * disk_flare(c, r);
        double dt = std::min(kVerticalStep * std::max(std::abs(p.z), sigma) / std::max(std::abs(d.z), 1e-4),
                             kRadialStep * std::max(glm::length(p), 1.0));
        dt = std::min(dt, t1 - t);
        const double middle = t + 0.5 * dt;
        const glm::dvec3 q = o + middle * d;
        const double rq = std::hypot(q.x, q.y);
        const double light = disk_table_light(disk, c, view, rq);
        if (light > 0.0) {
            const double distance = glm::length(q);
            const double along = view == DiskView::Scattered
                                     ? disk_phase(disk, -glm::dot(q, d) / distance) / (distance * distance)
                                     : disk_thermal_light(c, view, distance);
            sum += light * disk_vertical(c, rq, q.z) * dt * along;
        }
        t += dt;
    }
    return sum;
}

} // namespace

double disk_radial_au(const Disk& disk, int i)
{
    return disk.inner_au * std::pow(disk.outer_au / disk.inner_au, static_cast<double>(i) / (kDiskRadialCount - 1));
}

double disk_column(const Disk& disk, const DiskComponent& component, double r_au)
{
    if (!(r_au > disk.inner_au) || r_au >= disk.outer_au) {
        return 0.0;
    }
    const double f = (kDiskRadialCount - 1) * std::log(r_au / disk.inner_au) / std::log(disk.outer_au / disk.inner_au);
    const int i = std::min(static_cast<int>(f), kDiskRadialCount - 2);
    const double a = component.column[static_cast<size_t>(i)];
    const double b = component.column[static_cast<size_t>(i) + 1];
    return a + (b - a) * (f - i);
}

double disk_flare(const DiskComponent& component, double r_au)
{
    const double x = r_au / component.flare_radius_au;
    return std::sqrt(std::sqrt(1.0 + x * x * x * x));
}

double disk_vertical(const DiskComponent& component, double r_au, double z_au)
{
    const double flare = disk_flare(component, r_au);
    const double s1 = component.sigma_au * flare;
    const double a1 = z_au / s1;
    double out = (1.0 - component.share2) * std::exp(-0.5 * a1 * a1) / s1;
    if (component.share2 > 0.0) {
        const double s2 = component.sigma2_au * flare;
        const double a2 = z_au / s2;
        out += component.share2 * std::exp(-0.5 * a2 * a2) / s2;
    }
    return kInverseSqrtTwoPi * out;
}

double disk_temperature(const DiskComponent& component, double d_au)
{
    return component.temperature_k * std::pow(std::max(d_au, 1e-3) / component.temperature_radius_au,
                                              -component.temperature_index);
}

double disk_phase(const Disk& disk, double cos_theta)
{
    double out = 0.0;
    for (size_t i = 0; i < disk.phase_g.size(); ++i) {
        out += disk.phase_w[i] * henyey_greenstein(disk.phase_g[i], cos_theta);
    }
    return out;
}

double disk_table_light(const Disk& disk, const DiskComponent& component, DiskView view, double r_au)
{
    const double gain = component.gain[static_cast<size_t>(view)];
    return gain > 0.0 ? gain * disk_column(disk, component, r_au) : 0.0;
}

double disk_thermal_light(const DiskComponent& component, DiskView view, double d_au)
{
    const double t = disk_temperature(component, d_au);
    return view == DiskView::Infrared ? nebula_infrared_emission(t) : t;
}

double disk_radiance(const Disk& disk, DiskView view, const glm::dvec3& origin, const glm::dvec3& direction,
                     double t_max)
{
    // Only where the ray is inside the disk's sphere.
    const double along = -glm::dot(origin, direction);
    const glm::dvec3 closest = origin + along * direction;
    const double extent2 = disk.outer_au * disk.outer_au - glm::dot(closest, closest);
    if (extent2 <= 0.0) {
        return 0.0;
    }
    const double extent = std::sqrt(extent2);
    const double t0 = std::max(0.0, along - extent);
    const double t1 = std::min(t_max, along + extent);
    if (t1 <= t0) {
        return 0.0;
    }
    double sum = 0.0;
    for (const DiskComponent& c : disk.components) {
        if (c.gain[static_cast<size_t>(view)] > 0.0) {
            sum += march(disk, c, view, origin, direction, t0, t1);
        }
    }
    return sum;
}

void disk_set_references(Disk& disk)
{
    disk.reference.fill(0.0);
    if (disk.components.empty()) {
        return;
    }
    // Along the midplane of the view's first component, past the reference radius.
    for (int v = 0; v < kDiskViewCount; ++v) {
        const DiskView view = static_cast<DiskView>(v);
        for (const DiskComponent& c : disk.components) {
            if (c.gain[static_cast<size_t>(v)] <= 0.0) {
                continue;
            }
            const glm::dvec3 origin = c.to_icrf * glm::dvec3(disk.display_reference_au, -2.0 * disk.outer_au, 0.0);
            disk.reference[static_cast<size_t>(v)] = disk_radiance(disk, view, origin, c.to_icrf[1]);
            break;
        }
    }
}

void build_dust_stream(const DustStream& stream, int count, double patch_sigma_km, Nebula& nebula)
{
    nebula.grains.clear();
    nebula.extent_km = 0.0;
    nebula.mean_distance_km = 0.0;
    nebula.patch_sigma_km = patch_sigma_km;
    if (count <= 0 || stream.age_s <= 0.0) {
        return;
    }
    const double u = stream.argument_deg * kDegToRad;
    const glm::dvec3 x = stream.orbit_to_icrf[0];
    const glm::dvec3 y = stream.orbit_to_icrf[1];
    const glm::dvec3 z = stream.orbit_to_icrf[2];
    const glm::dvec3 start = stream.radius_km * (std::cos(u) * x + std::sin(u) * y);
    const double speed = std::sqrt(stream.gm_km3_s2 / stream.radius_km);
    const glm::dvec3 velocity = speed * (-std::sin(u) * x + std::cos(u) * y);

    // A year per step: the orbit takes centuries, and a grain with beta of 10
    // is still accelerating only over the first few decades.
    constexpr double kStepS = 365.25 * 86400.0;
    const int steps = std::max(1, static_cast<int>(std::ceil(stream.age_s / kStepS)));
    const double dt = stream.age_s / steps;

    // Low-discrepancy kicks: the R3 sequence (Roberts 2018), each coordinate
    // turned into a Gaussian by the logistic approximation of its inverse CDF.
    constexpr double kR3[3] = {0.8191725133961645, 0.6710436067037893, 0.5497004779019703};
    const auto gaussian = [](double v) {
        const double clamped = std::clamp(v, 1e-4, 1.0 - 1e-4);
        return std::log(clamped / (1.0 - clamped)) / 1.702;
    };
    nebula.grains.reserve(static_cast<size_t>(count));
    double weights = 0.0;
    for (int k = 0; k < count; ++k) {
        const double beta = stream.beta_lo + (stream.beta_hi - stream.beta_lo) * (k + 0.5) / count;
        const double gm = stream.gm_km3_s2 * (1.0 - beta);
        const auto acceleration = [gm](const glm::dvec3& p) {
            const double r2 = glm::dot(p, p);
            return -gm * p / (r2 * std::sqrt(r2));
        };
        const double sequence = static_cast<double>(k);
        const glm::dvec3 kick =
            stream.kick_share * speed *
            (gaussian(std::fmod(0.5 + kR3[0] * sequence, 1.0)) * x + gaussian(std::fmod(0.5 + kR3[1] * sequence, 1.0)) * y +
             gaussian(std::fmod(0.5 + kR3[2] * sequence, 1.0)) * z);
        glm::dvec3 p = start;
        glm::dvec3 v = velocity + kick;
        for (int s = 0; s < steps; ++s) {
            const glm::dvec3 a1 = acceleration(p);
            const glm::dvec3 v2 = v + 0.5 * dt * a1;
            const glm::dvec3 a2 = acceleration(p + 0.5 * dt * v);
            const glm::dvec3 v3 = v + 0.5 * dt * a2;
            const glm::dvec3 a3 = acceleration(p + 0.5 * dt * v2);
            const glm::dvec3 v4 = v + dt * a3;
            const glm::dvec3 a4 = acceleration(p + dt * v3);
            p += dt / 6.0 * (v + 2.0 * v2 + 2.0 * v3 + v4);
            v += dt / 6.0 * (a1 + 2.0 * a2 + 2.0 * a3 + a4);
        }
        NebulaGrain grain;
        grain.position = glm::vec3(p);
        grain.normal = glm::vec3(glm::normalize(p));
        grain.weight = static_cast<float>(std::pow(beta / stream.beta_lo, -stream.weight_index));
        weights += grain.weight;
        nebula.grains.push_back(grain);
        const double d = glm::length(p);
        nebula.extent_km = std::max(nebula.extent_km, d);
        nebula.mean_distance_km += d;
    }
    nebula.mean_distance_km /= static_cast<double>(nebula.grains.size());
    const float normalize = static_cast<float>(count / weights);
    for (NebulaGrain& grain : nebula.grains) {
        grain.weight *= normalize;
    }
}

} // namespace astraxis
