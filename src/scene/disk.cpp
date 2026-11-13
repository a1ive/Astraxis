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
constexpr double kThickDustDepth = 1e-4;
constexpr double kThickEnvelope = 5.0;
constexpr double kThickDepthStep = 0.15;
constexpr double kThickHeightStep = 0.25;
constexpr double kThickRadialStep = 0.06;
constexpr double kThickEmptyHeightStep = 0.5;
constexpr double kThickEmptyRadialStep = 0.1;
constexpr int kThickMaxSteps = 96;
constexpr double kThickMinTransmission = 0.005;

// The grid of the later scatterings' source. Rows and columns as star_depth's, but
// the source changes slowly: the tables only have to follow the extinction.
constexpr int kDiffuseColumns = 80;
constexpr int kDiffuseRows = 40;
constexpr int kDiffuseRays = 40;
constexpr int kDiffuseSteps = 48;
constexpr double kDiffuseLeastDistance = 0.25; // of the point's own distance from the star

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
            light += transmission * (1.0 - through) *
                     (disk.single_scatter_gain * disk_phase(disk, -glm::dot(q, d) / distance) *
                          std::exp(-disk_star_depth(disk, rq, q.z)) +
                      disk_diffuse(disk, rq, q.z)) /
                     (distance * distance);
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
             const glm::dvec3& direction, double t0, double t1, double t_tdb, bool r2_scaled)
{
    const double reference2 = disk.display_reference_au * disk.display_reference_au;
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
                                     ? disk_phase(disk, -glm::dot(q, d) / distance) /
                                           (r2_scaled ? reference2 : distance * distance)
                                 : view == DiskView::Gas ? disk_clump_factor(disk, rq, std::atan2(q.y, q.x), t_tdb)
                                                         : disk_thermal_light(c, view, distance);
            double skew = 0.0;
            if (view == DiskView::Scattered && (c.warp_over_radius != 0.0 || c.contrast != 0.0)) {
                skew = glm::smoothstep(c.asymmetry_start_au, c.asymmetry_full_au, rq) * q.x;
            }
            const double contrast = 1.0 + c.contrast * skew / std::max(rq, disk.inner_au);
            sum += light * disk_vertical(c, rq, q.z, c.warp_over_radius * skew) * contrast * dt * along;
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

double disk_vertical(const DiskComponent& component, double r_au, double z_au, double broad_shift_au)
{
    const double flare = disk_flare(component, r_au);
    const double s1 = component.sigma_au * flare;
    const double a1 = z_au / s1;
    double out = (1.0 - component.share2) * std::exp(-0.5 * a1 * a1) / s1;
    if (component.share2 > 0.0) {
        const double s2 = component.sigma2_au * flare;
        const double a2 = (z_au - broad_shift_au) / s2;
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

double disk_clump_azimuth(const DiskClump& clump, double t_tdb)
{
    return clump.azimuth_deg * kDegToRad + kTwoPi * (t_tdb - clump.epoch_tdb) / clump.period_s;
}

double disk_clump_factor(const Disk& disk, double r_au, double phi, double t_tdb)
{
    double out = 1.0;
    for (const DiskClump& clump : disk.clumps) {
        const double a = (r_au - clump.radius_au) / clump.sigma_au;
        const double from = phi - disk_clump_azimuth(clump, t_tdb);
        const double delta = std::atan2(std::sin(from), std::cos(from));
        const double shape = delta >= 0.0 ? std::exp(-0.5 * (delta / (clump.sigma_deg * kDegToRad)) * (delta / (clump.sigma_deg * kDegToRad)))
                                          : std::exp(delta / (clump.tail_deg * kDegToRad));
        out += clump.gain * std::exp(-0.5 * a * a) * shape;
    }
    return out;
}

double disk_line_of_sight_km_s(const DiskComponent& component, double gm_km3_s2, const glm::dvec3& position_au,
                               const glm::dvec3& direction)
{
    const glm::dmat3 to_plane = glm::transpose(component.to_icrf);
    const glm::dvec3 p = to_plane * position_au;
    const glm::dvec3 d = to_plane * direction;
    const double r = std::max(std::hypot(p.x, p.y), 1e-6);
    const double speed = std::sqrt(gm_km3_s2 / (r * kAuKm));
    return speed * (-p.y * d.x + p.x * d.y) / r;
}

double disk_radiance(const Disk& disk, DiskView view, const glm::dvec3& origin, const glm::dvec3& direction,
                     double t_max, double t_tdb, bool r2_scaled)
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
            sum += march(disk, c, view, origin, direction, t0, t1, t_tdb, r2_scaled && view == DiskView::Scattered);
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

void disk_apply_similarity(Disk& disk)
{
    if (disk.albedo <= 0.0) {
        disk.single_scatter_gain = 1.0;
        return;
    }
    const double denominator = 1.0 - disk.albedo * disk.asymmetry;
    for (DiskComponent& c : disk.components) {
        c.tau_per_column *= denominator;
    }
    disk.single_scatter_gain = disk.albedo / denominator;
}

double disk_scaled_albedo(const Disk& disk)
{
    return disk.albedo * (1.0 - disk.asymmetry) / (1.0 - disk.albedo * disk.asymmetry);
}

namespace {

// Bilinear on a polar table, in whatever the table stores.
struct PolarGrid {
    int columns = 0;
    int rows = 0;
    double ln_inner = 0.0;
    double ln_step = 1.0;

    double radius(int column) const { return std::exp(ln_inner + column * ln_step); }
    static double elevation(int row, int rows)
    {
        const double v = static_cast<double>(row) / (rows - 1);
        return v * v * 0.5 * kPi;
    }
    double sample(const std::vector<float>& table, double r, double z_abs) const
    {
        const double u = std::clamp((std::log(std::max(r, 1e-6)) - ln_inner) / ln_step, 0.0, columns - 1.0);
        const double v = std::sqrt(std::atan2(z_abs, std::max(r, 1e-6)) / (0.5 * kPi)) * (rows - 1);
        const int iu = std::min(static_cast<int>(u), columns - 2);
        const int iv = std::min(static_cast<int>(v), rows - 2);
        const double fu = u - iu;
        const double fv = v - iv;
        const auto at = [&](int row, int column) {
            return static_cast<double>(table[static_cast<size_t>(row) * columns + static_cast<size_t>(column)]);
        };
        return (at(iv, iu) * (1.0 - fu) + at(iv, iu + 1) * fu) * (1.0 - fv) +
               (at(iv + 1, iu) * (1.0 - fu) + at(iv + 1, iu + 1) * fu) * fv;
    }
};

PolarGrid disk_table_grid(const Disk& disk, int columns, int rows)
{
    return {columns, rows, std::log(disk.inner_au), std::log(disk.outer_au / disk.inner_au) / (columns - 1)};
}

} // namespace

void disk_build_diffuse(Disk& disk)
{
    disk.diffuse.clear();
    const double omega = disk_scaled_albedo(disk);
    const DiskComponent* c = thick_component(disk);
    if (disk.albedo <= 0.0 || omega <= 0.0 || !c) {
        return;
    }
    const PolarGrid grid = disk_table_grid(disk, kDiffuseColumns, kDiffuseRows);
    const size_t cells = static_cast<size_t>(kDiffuseColumns) * kDiffuseRows;
    constexpr double kFloorExtinction = 1e-30;
    constexpr double kFloorSource = 1e-30;

    // Fibonacci directions on the sphere.
    std::array<glm::dvec3, kDiffuseRays> directions{};
    for (int i = 0; i < kDiffuseRays; ++i) {
        const double z = 1.0 - (2.0 * i + 1.0) / kDiffuseRays;
        const double s = std::sqrt(1.0 - z * z);
        const double phi = i * 2.399963229728653; // the golden angle
        directions[static_cast<size_t>(i)] = glm::dvec3(s * std::cos(phi), s * std::sin(phi), z);
    }

    // On the grid: ln extinction, the disk's height, and ln of the source times d^2
    // (which has none of the 1 / d^2 to resolve). The first scattering's source is
    // omega * exp(-tau_star) / d^2. It falls off exponentially into the disk, where
    // interpolating the source itself would light the dust just inside the shadow.
    std::vector<float> ln_extinction(cells);
    std::vector<float> first(cells);
    std::vector<double> sigma(kDiffuseColumns);
    for (int column = 0; column < kDiffuseColumns; ++column) {
        sigma[static_cast<size_t>(column)] = c->sigma_au * disk_flare(*c, grid.radius(column));
    }
    for (int row = 0; row < kDiffuseRows; ++row) {
        const double slope = std::tan(std::min(PolarGrid::elevation(row, kDiffuseRows), 1.55));
        for (int column = 0; column < kDiffuseColumns; ++column) {
            const double r = grid.radius(column);
            const size_t k = static_cast<size_t>(row) * kDiffuseColumns + static_cast<size_t>(column);
            ln_extinction[k] = static_cast<float>(std::log(disk_extinction(disk, r, r * slope) + kFloorExtinction));
            first[k] = static_cast<float>(std::log(omega * std::exp(-disk_star_depth(disk, r, r * slope)) + kFloorSource));
        }
    }

    // The mean intensity at each grid point from a source on the grid. A ray stops
    // when little of its light would get through. What lies closer to the star than a
    // fraction of the point's own distance is too small and bright for so few rays to
    // find except by chance, so it counts as one more star at the center: its
    // scattered light, summed, dimmed as the starlight is on the way to the point.
    const auto next_order = [&](const std::vector<float>& source) {
        std::vector<float> out(cells);
        // The light scattered within each distance from the star, on both sides of the
        // disk: extinction * source * volume, in the steps of ln r.
        std::vector<double> within(kDiffuseColumns + 1, 0.0);
        const double dv = 1.0 / (kDiffuseRows - 1);
        for (int row = 0; row < kDiffuseRows; ++row) {
            const double v = static_cast<double>(row) * dv;
            const double elevation = std::min(PolarGrid::elevation(row, kDiffuseRows), 1.55);
            const double slope = std::tan(elevation);
            const double per_elevation = 1.0 + slope * slope;
            for (int column = 0; column < kDiffuseColumns; ++column) {
                const size_t k = static_cast<size_t>(row) * kDiffuseColumns + static_cast<size_t>(column);
                const double r = grid.radius(column);
                const double distance2 = r * r * per_elevation;
                const double dr = r * grid.ln_step;
                const double dz = r * per_elevation * kPi * v * dv;
                const double volume = 2.0 * kTwoPi * r * dr * dz;
                const double light = std::exp(static_cast<double>(ln_extinction[k]) + source[k]) / distance2 * volume;
                const int bin = static_cast<int>((0.5 * std::log(distance2) - grid.ln_inner) / grid.ln_step) + 1;
                within[static_cast<size_t>(std::clamp(bin, 0, kDiffuseColumns))] += light;
            }
        }
        for (int bin = 1; bin <= kDiffuseColumns; ++bin) {
            within[static_cast<size_t>(bin)] += within[static_cast<size_t>(bin) - 1];
        }
        for (int row = 0; row < kDiffuseRows; ++row) {
            const double elevation = std::min(PolarGrid::elevation(row, kDiffuseRows), 1.55);
            const double slope = std::tan(elevation);
            for (int column = 0; column < kDiffuseColumns; ++column) {
                const double r0 = grid.radius(column);
                const glm::dvec3 start(r0, 0.0, r0 * slope);
                const double d0 = glm::length(start);
                const double near_star = kDiffuseLeastDistance * d0;
                const int near_bin = std::clamp(static_cast<int>((std::log(near_star) - grid.ln_inner) / grid.ln_step) + 1, 0,
                                                kDiffuseColumns);
                double sum = 0.0;
                for (const glm::dvec3& d : directions) {
                    double t = 0.0;
                    double tau = 0.0;
                    for (int step = 0; step < kDiffuseSteps && tau < 5.0; ++step) {
                        const glm::dvec3 p = start + t * d;
                        const double r = std::hypot(p.x, p.y);
                        const double u = std::clamp((std::log(std::max(r, 1e-6)) - grid.ln_inner) / grid.ln_step, 0.0,
                                                    kDiffuseColumns - 1.0);
                        const double height = std::max(std::abs(p.z), sigma[static_cast<size_t>(u)]);
                        const double to_star = std::max(glm::length(p), 1e-3);
                        if (to_star > disk.outer_au) {
                            break;
                        }
                        const double dt = std::min(0.25 * height / std::max(std::abs(d.z), 1e-3), 0.3 * to_star);
                        const glm::dvec3 q = p + 0.5 * dt * d;
                        const double rq = std::hypot(q.x, q.y);
                        const double zq = std::abs(q.z);
                        const double depth = std::exp(grid.sample(ln_extinction, rq, zq)) * dt;
                        const double dq = glm::length(q);
                        if (dq > near_star) {
                            sum += std::exp(-tau) * (1.0 - std::exp(-depth)) * std::exp(grid.sample(source, rq, zq)) / (dq * dq);
                        }
                        tau += depth;
                        t += dt;
                    }
                }
                // Mostly light off a flat surface: Lambert's law, with a little that gets
                // out sideways. The mean over the sphere is 1.
                constexpr double kCentralSideways = 0.25;
                const double flat = kCentralSideways + (1.0 - kCentralSideways) * 2.0 * std::sin(elevation);
                const double central = flat * within[static_cast<size_t>(near_bin)] / (2.0 * kTwoPi * d0 * d0) *
                                       std::exp(-disk_star_depth(disk, r0, r0 * slope));
                // Back to the same units: source * d^2 at this grid point.
                out[static_cast<size_t>(row) * kDiffuseColumns + static_cast<size_t>(column)] =
                    static_cast<float>(std::log(omega * (sum / kDiffuseRays + central) * d0 * d0 + kFloorSource));
            }
        }
        return out;
    };
    const std::vector<float> second = next_order(first);
    const std::vector<float> third = next_order(second);

    // Onto the same grid as star_depth. What follows the third order is taken as a
    // geometric series with the ratio omega.
    std::vector<float> ln_source(cells);
    const double rest = 1.0 / (1.0 - omega);
    for (size_t k = 0; k < cells; ++k) {
        ln_source[k] = static_cast<float>(std::log(std::exp(static_cast<double>(second[k])) +
                                                   std::exp(static_cast<double>(third[k])) * rest + kDiskDiffuseFloor));
    }
    const PolarGrid fine = disk_table_grid(disk, kDiskDepthWidth, kDiskDepthHeight);
    disk.diffuse.resize(static_cast<size_t>(kDiskDepthWidth) * kDiskDepthHeight);
    for (int row = 0; row < kDiskDepthHeight; ++row) {
        const double slope = std::tan(std::min(PolarGrid::elevation(row, kDiskDepthHeight), 1.55));
        for (int column = 0; column < kDiskDepthWidth; ++column) {
            const double r = fine.radius(column);
            disk.diffuse[static_cast<size_t>(row) * kDiskDepthWidth + static_cast<size_t>(column)] =
                static_cast<float>(grid.sample(ln_source, r, r * slope));
        }
    }
}

double disk_diffuse(const Disk& disk, double r_au, double z_au)
{
    if (disk.diffuse.empty()) {
        return 0.0;
    }
    const PolarGrid fine = disk_table_grid(disk, kDiskDepthWidth, kDiskDepthHeight);
    return std::exp(fine.sample(disk.diffuse, r_au, std::abs(z_au))) - kDiskDiffuseFloor;
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
    // The reference does not move with the clumps.
    std::vector<DiskClump> clumps;
    clumps.swap(disk.clumps);
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
    clumps.swap(disk.clumps);
}

double dust_stream_drag_per_km(double gm_km3_s2, double luminosity_solar, double radiation_q, double ism_h_per_cm3)
{
    constexpr double kSolarLuminosityKgKm2S3 = 3.828e20; // 3.828e26 W, IAU 2015 B3
    constexpr double kHydrogenKg = 1.6735e-27;
    constexpr double kCm3PerKm3 = 1e15;
    const double density = ism_h_per_cm3 * kHydrogenKg * kCm3PerKm3; // kg per km^3
    return 4.0 * kPi * gm_km3_s2 * kSpeedOfLightKmS * density / (radiation_q * luminosity_solar * kSolarLuminosityKgKm2S3);
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
        const double drag = stream.beta_hi * std::pow(beta / stream.beta_hi, stream.ism_beta_index) * stream.ism_drag_per_km;
        const auto acceleration = [&stream, gm, drag](const glm::dvec3& p, const glm::dvec3& v) {
            const double r2 = glm::dot(p, p);
            const glm::dvec3 wind = stream.ism_velocity_km_s - v;
            return -gm * p / (r2 * std::sqrt(r2)) + drag * glm::length(wind) * wind;
        };
        const double sequence = static_cast<double>(k);
        const glm::dvec3 kick =
            stream.kick_share * speed *
            (gaussian(std::fmod(0.5 + kR3[0] * sequence, 1.0)) * x + gaussian(std::fmod(0.5 + kR3[1] * sequence, 1.0)) * y +
             gaussian(std::fmod(0.5 + kR3[2] * sequence, 1.0)) * z);
        glm::dvec3 p = start;
        glm::dvec3 v = velocity + kick;
        for (int s = 0; s < steps; ++s) {
            const glm::dvec3 a1 = acceleration(p, v);
            const glm::dvec3 v2 = v + 0.5 * dt * a1;
            const glm::dvec3 a2 = acceleration(p + 0.5 * dt * v, v2);
            const glm::dvec3 v3 = v + 0.5 * dt * a2;
            const glm::dvec3 a3 = acceleration(p + 0.5 * dt * v2, v3);
            const glm::dvec3 v4 = v + dt * a3;
            const glm::dvec3 a4 = acceleration(p + dt * v3, v4);
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
