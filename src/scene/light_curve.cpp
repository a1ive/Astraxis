#include "scene/light_curve.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cmath>

namespace astraxis {

namespace {

// h c / k in nm K (CODATA 2018: h, c and k are exact).
constexpr double kHcOverKNmK = 1.438776877e7;

// Each disk is summed over this many rings of equal width: the eclipse of
// Algol A by B is then smooth to well under 0.001 mag.
constexpr int kRings = 128;
// Disks in front of one disk that are taken into account (a ring's covered
// arcs are merged from at most this many).
constexpr int kMaxOccluders = 8;

struct Occluder {
    double distance = 0.0; // between the centers, on the sky (km)
    double angle = 0.0;    // direction of its center seen from the disk's (radians)
    double radius = 0.0;
};

// Fraction of the circle of radius r around a disk's center hidden by the
// occluders: the union of their arcs.
double covered_fraction(double r, std::span<const Occluder> occluders)
{
    std::array<std::array<double, 2>, 2 * kMaxOccluders> arcs; // [start, end) in [0, 2 pi)
    int count = 0;
    for (const Occluder& o : occluders) {
        if (o.distance + r <= o.radius) {
            return 1.0; // the whole circle is inside it
        }
        if (o.distance >= r + o.radius || o.distance + o.radius <= r) {
            continue; // apart, or the occluder lies inside the circle without touching it
        }
        const double cos_half = (r * r + o.distance * o.distance - o.radius * o.radius) / (2.0 * r * o.distance);
        const double half = std::acos(std::clamp(cos_half, -1.0, 1.0));
        double start = std::fmod(o.angle - half, kTwoPi);
        if (start < 0.0) {
            start += kTwoPi;
        }
        const double end = start + 2.0 * half;
        if (end > kTwoPi) {
            arcs[static_cast<size_t>(count++)] = {start, kTwoPi};
            arcs[static_cast<size_t>(count++)] = {0.0, end - kTwoPi};
        } else {
            arcs[static_cast<size_t>(count++)] = {start, end};
        }
    }
    if (count == 0) {
        return 0.0;
    }
    std::sort(arcs.begin(), arcs.begin() + count);
    double covered = 0.0;
    double from = arcs[0][0];
    double to = arcs[0][1];
    for (int i = 1; i < count; ++i) {
        if (arcs[static_cast<size_t>(i)][0] > to) {
            covered += to - from;
            from = arcs[static_cast<size_t>(i)][0];
        }
        to = std::max(to, arcs[static_cast<size_t>(i)][1]);
    }
    covered += to - from;
    return std::min(covered / kTwoPi, 1.0);
}

// Light of disk `index`, hiding what the disks in front cover (none if `hide` is false).
double disk_flux(std::span<const StarDisk> disks, size_t index, const glm::dvec3& toward, const glm::dvec3& e1,
                 const glm::dvec3& e2, double wavelength_nm, double limb_darkening, bool hide)
{
    const StarDisk& disk = disks[index];
    std::array<Occluder, kMaxOccluders> occluders;
    int count = 0;
    if (hide) {
        const double depth = glm::dot(disk.position, toward);
        for (size_t k = 0; k < disks.size() && count < kMaxOccluders; ++k) {
            const glm::dvec3 d = disks[k].position - disk.position;
            const double x = glm::dot(d, e1);
            const double y = glm::dot(d, e2);
            const double sky = std::hypot(x, y);
            if (k == index || glm::dot(disks[k].position, toward) <= depth ||
                sky >= disk.radius_km + disks[k].radius_km) {
                continue;
            }
            occluders[static_cast<size_t>(count++)] = {sky, std::atan2(y, x), disks[k].radius_km};
        }
    }
    const std::span<const Occluder> in_front(occluders.data(), static_cast<size_t>(count));
    const double dr = disk.radius_km / kRings;
    double sum = 0.0;
    for (int i = 0; i < kRings; ++i) {
        const double r = (i + 0.5) * dr;
        const double s = r / disk.radius_km;
        const double mu = std::sqrt(std::max(1.0 - s * s, 0.0));
        const double ring = kTwoPi * r * dr * (1.0 - limb_darkening * (1.0 - mu));
        sum += count > 0 ? ring * (1.0 - covered_fraction(r, in_front)) : ring;
    }
    return sum * planck_shape(disk.temperature_k, wavelength_nm);
}

void sky_axes(const glm::dvec3& toward, glm::dvec3& e1, glm::dvec3& e2)
{
    const glm::dvec3 helper = std::abs(toward.z) < 0.9 ? glm::dvec3(0.0, 0.0, 1.0) : glm::dvec3(1.0, 0.0, 0.0);
    e1 = glm::normalize(glm::cross(helper, toward));
    e2 = glm::cross(toward, e1);
}

} // namespace

double planck_shape(double temperature_k, double wavelength_nm)
{
    return 1.0 / std::expm1(kHcOverKNmK / (wavelength_nm * temperature_k));
}

double disks_flux(std::span<const StarDisk> disks, const glm::dvec3& toward_observer, double wavelength_nm,
                  double limb_darkening)
{
    glm::dvec3 e1, e2;
    sky_axes(toward_observer, e1, e2);
    double sum = 0.0;
    for (size_t i = 0; i < disks.size(); ++i) {
        sum += disk_flux(disks, i, toward_observer, e1, e2, wavelength_nm, limb_darkening, true);
    }
    return sum;
}

double disks_full_flux(std::span<const StarDisk> disks, double wavelength_nm, double limb_darkening)
{
    const glm::dvec3 toward(0.0, 0.0, 1.0);
    glm::dvec3 e1, e2;
    sky_axes(toward, e1, e2);
    double sum = 0.0;
    for (size_t i = 0; i < disks.size(); ++i) {
        sum += disk_flux(disks, i, toward, e1, e2, wavelength_nm, limb_darkening, false);
    }
    return sum;
}

double light_curve_magnitude(const Scene& scene, double t_tdb)
{
    const LightCurve& lc = scene.light_curve;
    if (!lc.enabled()) {
        return 0.0;
    }
    std::array<StarDisk, kMaxOccluders + 1> disks;
    const size_t count = std::min(lc.stars.size(), disks.size());
    for (size_t i = 0; i < count; ++i) {
        const Body& star = scene.bodies[static_cast<size_t>(lc.stars[i])];
        disks[i] = {scene.icrf_state_at(lc.stars[i], t_tdb).position, star.equatorial_radius_km, star.temperature_k};
    }
    const std::span<const StarDisk> span(disks.data(), count);
    const glm::dvec3 toward = glm::normalize(scene.icrf_state_at(lc.observer, t_tdb).position - disks[0].position);
    const double flux = disks_flux(span, toward, lc.wavelength_nm, lc.limb_darkening);
    const double full = disks_full_flux(span, lc.wavelength_nm, lc.limb_darkening);
    return std::max(-2.5 * std::log10(std::max(flux / full, 1e-12)), 0.0);
}

void measure_light_curve_depth(Scene& scene, double t_tdb, double days)
{
    LightCurve& lc = scene.light_curve;
    if (!lc.enabled()) {
        return;
    }
    const double step = lc.span_days * kSecondsPerDay / (LightCurveSamples::kCount - 1);
    const int count = static_cast<int>(std::ceil(days * kSecondsPerDay / step));
    double deepest = 0.0;
    for (int i = 0; i <= count; ++i) {
        deepest = std::max(deepest, light_curve_magnitude(scene, t_tdb + i * step));
    }
    lc.depth_mag = std::max(std::ceil(deepest * 10.0 - 1e-6) / 10.0, 0.1);
}

void LightCurveSamples::update(const Scene& scene, double t_tdb)
{
    const LightCurve& lc = scene.light_curve;
    if (!lc.enabled()) {
        reset();
        return;
    }
    const double step = lc.span_days * kSecondsPerDay / (kCount - 1);
    if (step != m_step_s) {
        reset();
        m_step_s = step;
    }
    const long long first = static_cast<long long>(std::floor(t_tdb / step - kBefore * (kCount - 1)));
    if (!m_magnitudes.empty() && first == m_first) {
        return;
    }
    m_scratch.assign(kCount, 0.0f);
    for (int i = 0; i < kCount; ++i) {
        const long long old = first + i - m_first; // the same sample's index before
        if (!m_magnitudes.empty() && old >= 0 && old < kCount) {
            m_scratch[static_cast<size_t>(i)] = m_magnitudes[static_cast<size_t>(old)];
        } else {
            m_scratch[static_cast<size_t>(i)] =
                static_cast<float>(light_curve_magnitude(scene, static_cast<double>(first + i) * step));
        }
    }
    m_magnitudes.swap(m_scratch);
    m_first = first;
}

} // namespace astraxis
