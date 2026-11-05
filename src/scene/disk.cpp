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

// The steps in a thick disk, as disk_march_thick in the shader. Where there is dust
// to speak of (the optical depth through a Gaussian's width over kThickDustDepth) a
// step takes a fraction of the height over which the density changes, of the optical
// depth and of the distance to the star, since the light of the disk's surface comes
// from a thin layer and a coarser step shows as moire across it. Elsewhere the ray
// jumps to kThickEnvelope sigmas from the midplane, so the budget lasts.
constexpr double kThickDustDepth = 1e-3;
constexpr double kThickEnvelope = 5.0;
constexpr double kThickDepthStep = 0.15;
constexpr double kThickHeightStep = 0.25;
constexpr double kThickRadialStep = 0.06;
constexpr double kThickEmptyHeightStep = 0.5;
constexpr double kThickEmptyRadialStep = 0.2;
constexpr int kThickMaxSteps = 96;
constexpr double kThickMinTransmission = 0.005;

const DiskComponent* thick_component(const Disk& disk)
{
    for (const DiskComponent& c : disk.components) {
        if (c.tau_per_column > 0.0) {
            return &c;
        }
    }
    return nullptr;
}

// Light scattered toward the ray's origin (phase function / d^2 units), and
// through `transmission` what it lets through.
double march_thick(const Disk& disk, const glm::dvec3& origin, const glm::dvec3& direction, double t0, double t1,
                   double& transmission)
{
    transmission = 1.0;
    const DiskComponent* c = thick_component(disk);
    if (!c) {
        return 0.0;
    }
    const glm::dmat3 to_plane = glm::transpose(c->to_icrf);
    const glm::dvec3 o = to_plane * origin;
    const glm::dvec3 d = to_plane * direction;
    double light = 0.0;
    double t = t0;
    for (int i = 0; i < kThickMaxSteps && t < t1 && transmission > kThickMinTransmission; ++i) {
        const glm::dvec3 p = o + t * d;
        const double r = std::hypot(p.x, p.y);
        // High over the axis the disk's height at r says nothing of the path
        // ahead, and would make the steps vanish.
        const double sigma = c->sigma_au * disk_flare(*c, std::max(r, 0.5 * std::abs(p.z)));
        const double height = std::max(std::abs(p.z), sigma);
        const double climb = std::max(std::abs(d.z), 1e-3);
        const double extinction = disk_extinction(disk, r, p.z);
        const double to_star = std::max(glm::length(p), 1.0);
        double dt = 0.0;
        if (extinction * sigma * 2.5 > kThickDustDepth) {
            // sigma^2 / |z|: where the density changes by a factor e on the disk's surface.
            dt = std::min({kThickHeightStep * sigma * std::min(1.0, sigma / height) / climb, kThickDepthStep / extinction,
                           kThickRadialStep * to_star});
        } else {
            dt = std::min(std::max(kThickEmptyHeightStep * sigma, std::abs(p.z) - kThickEnvelope * sigma) / climb,
                          kThickEmptyRadialStep * to_star);
        }
        dt = std::min(dt, t1 - t);
        const glm::dvec3 q = o + (t + 0.5 * dt) * d;
        const double rq = std::hypot(q.x, q.y);
        const double depth = disk_extinction(disk, rq, q.z) * dt;
        if (depth > 0.0) {
            const double distance = glm::length(q);
            const double through = std::exp(-depth);
            light += transmission * (1.0 - through) * disk_phase(disk, -glm::dot(q, d) / distance) *
                     std::exp(-disk_star_depth(disk, rq, q.z)) / (distance * distance);
            transmission *= through;
        }
        t += dt;
    }
    return light;
}

// Where a ray is inside the disk's sphere.
bool clip_ray(const Disk& disk, const glm::dvec3& origin, const glm::dvec3& direction, double t_max, double& t0,
              double& t1)
{
    const double along = -glm::dot(origin, direction);
    const glm::dvec3 closest = origin + along * direction;
    const double extent2 = disk.outer_au * disk.outer_au - glm::dot(closest, closest);
    if (extent2 <= 0.0) {
        return false;
    }
    const double extent = std::sqrt(extent2);
    t0 = std::max(0.0, along - extent);
    t1 = std::min(t_max, along + extent);
    return t1 > t0;
}

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
    return std::pow(1.0 + x * x * x * x, 0.25 * component.flare_index);
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
    double t0 = 0.0;
    double t1 = 0.0;
    if (!clip_ray(disk, origin, direction, t_max, t0, t1)) {
        return 0.0;
    }
    if (disk.thick && view == DiskView::Scattered) {
        double transmission = 1.0;
        return march_thick(disk, origin, direction, t0, t1, transmission);
    }
    double sum = 0.0;
    for (const DiskComponent& c : disk.components) {
        if (c.gain[static_cast<size_t>(view)] > 0.0) {
            sum += march(disk, c, view, origin, direction, t0, t1);
        }
    }
    return sum;
}

double disk_extinction(const Disk& disk, double r_au, double z_au)
{
    double out = 0.0;
    for (const DiskComponent& c : disk.components) {
        if (c.tau_per_column > 0.0) {
            out += c.tau_per_column * disk_column(disk, c, r_au) * disk_vertical(c, r_au, z_au);
        }
    }
    return out;
}

void disk_build_star_depth(Disk& disk)
{
    disk.star_depth.assign(static_cast<size_t>(kDiskDepthWidth) * kDiskDepthHeight, 0.0f);
    const double ln_inner = std::log(disk.inner_au);
    const double ln_step = std::log(disk.outer_au / disk.inner_au) / (kDiskDepthWidth - 1);
    for (int row = 0; row < kDiskDepthHeight; ++row) {
        const double v = static_cast<double>(row) / (kDiskDepthHeight - 1);
        const double elevation = v * v * 0.5 * kPi;
        const double slope = std::tan(std::min(elevation, 1.55));
        const double per_r = std::sqrt(1.0 + slope * slope); // path per unit of r
        double depth = 0.0;
        for (int column = 0; column < kDiskDepthWidth; ++column) {
            disk.star_depth[static_cast<size_t>(row) * kDiskDepthWidth + static_cast<size_t>(column)] =
                static_cast<float>(std::log1p(depth));
            // To the next entry, in a few sub-steps to catch sharp edges.
            const double r0 = std::exp(ln_inner + column * ln_step);
            const double r1 = std::exp(ln_inner + (column + 1) * ln_step);
            constexpr int kSubSteps = 8;
            for (int s = 0; s < kSubSteps; ++s) {
                const double r = r0 + (r1 - r0) * (s + 0.5) / kSubSteps;
                depth += disk_extinction(disk, r, r * slope) * (r1 - r0) / kSubSteps * per_r;
            }
        }
    }
}

double disk_star_depth(const Disk& disk, double r_au, double z_au)
{
    if (disk.star_depth.empty()) {
        return 0.0;
    }
    const double r = std::max(r_au, 1e-6);
    const double u = std::clamp(std::log(r / disk.inner_au) / std::log(disk.outer_au / disk.inner_au), 0.0, 1.0) *
                     (kDiskDepthWidth - 1);
    const double v = std::sqrt(std::atan2(std::abs(z_au), r) / (0.5 * kPi)) * (kDiskDepthHeight - 1);
    const int iu = std::min(static_cast<int>(u), kDiskDepthWidth - 2);
    const int iv = std::min(static_cast<int>(v), kDiskDepthHeight - 2);
    const double fu = u - iu;
    const double fv = v - iv;
    const auto at = [&](int row, int column) {
        return static_cast<double>(disk.star_depth[static_cast<size_t>(row) * kDiskDepthWidth + static_cast<size_t>(column)]);
    };
    const double low = at(iv, iu) * (1.0 - fu) + at(iv, iu + 1) * fu;
    const double high = at(iv + 1, iu) * (1.0 - fu) + at(iv + 1, iu + 1) * fu;
    return std::expm1(low * (1.0 - fv) + high * fv);
}

double disk_optical_depth(const Disk& disk, const glm::dvec3& origin, const glm::dvec3& direction, double t_max)
{
    double t0 = 0.0;
    double t1 = 0.0;
    if (!disk.thick || !clip_ray(disk, origin, direction, t_max, t0, t1)) {
        return 0.0;
    }
    double transmission = 1.0;
    march_thick(disk, origin, direction, t0, t1, transmission);
    return -std::log(transmission);
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
            // Straight through a thick disk (along its midplane nothing gets out).
            const glm::dvec3 start = disk.thick ? glm::dvec3(disk.display_reference_au, 0.0, 2.0 * disk.outer_au)
                                                : glm::dvec3(disk.display_reference_au, -2.0 * disk.outer_au, 0.0);
            const glm::dvec3 along = disk.thick ? -c.to_icrf[2] : c.to_icrf[1];
            disk.reference[static_cast<size_t>(v)] = disk_radiance(disk, view, c.to_icrf * start, along);
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
