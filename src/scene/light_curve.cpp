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
// Directions a Roche surface is seen from to find where it is brightest:
// circles of latitude about its x axis (both poles included), points on each.
constexpr int kBrightestLatitudes = 9;
constexpr int kBrightestLongitudes = 32;

struct Occluder {
    double x = 0.0; // its center on the sky, from the hidden disk's (km)
    double y = 0.0;
    double distance = 0.0; // between the centers, on the sky (km)
    double angle = 0.0;    // direction of its center seen from the disk's (radians)
    double radius = 0.0;
};

// A Roche surface's outline on the sky: a convex polygon, counter-clockwise,
// around the hidden disk's center (km).
using Outline = std::vector<glm::dvec2>;
using Arc = std::array<double, 2>; // [start, end) in [0, 2 pi)

double cross(const glm::dvec2& a, const glm::dvec2& b) { return a.x * b.y - a.y * b.x; }

// Andrew's monotone chain over the projected vertices of the surface.
Outline roche_outline(const StarDisk& star, const glm::dvec3& origin, const glm::dvec3& e1, const glm::dvec3& e2)
{
    std::vector<glm::dvec2> points;
    points.reserve(star.roche->mesh.positions.size());
    for (const glm::vec3& v : star.roche->mesh.positions) {
        const glm::dvec3 p = star.position + star.axes * glm::dvec3(v) - origin;
        points.emplace_back(glm::dot(p, e1), glm::dot(p, e2));
    }
    std::sort(points.begin(), points.end(),
              [](const glm::dvec2& a, const glm::dvec2& b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
    Outline hull(2 * points.size());
    size_t k = 0;
    for (size_t i = 0; i < points.size(); ++i) { // lower hull
        while (k >= 2 && cross(hull[k - 1] - hull[k - 2], points[i] - hull[k - 2]) <= 0.0) {
            --k;
        }
        hull[k++] = points[i];
    }
    for (size_t i = points.size() - 1, lower = k + 1; i > 0; --i) { // upper hull
        while (k >= lower && cross(hull[k - 1] - hull[k - 2], points[i - 1] - hull[k - 2]) <= 0.0) {
            --k;
        }
        hull[k++] = points[i - 1];
    }
    hull.resize(k - 1); // the last point is the first
    return hull;
}

bool inside(const Outline& outline, const glm::dvec2& p)
{
    for (size_t i = 0; i < outline.size(); ++i) {
        const glm::dvec2& a = outline[i];
        const glm::dvec2& b = outline[(i + 1) % outline.size()];
        if (cross(b - a, p - a) < 0.0) {
            return false;
        }
    }
    return true;
}

void add_arc(double start, double end, std::vector<Arc>& arcs)
{
    const double length = std::min(end - start, kTwoPi);
    start = std::fmod(start, kTwoPi);
    if (start < 0.0) {
        start += kTwoPi;
    }
    end = start + length;
    if (end > kTwoPi) {
        arcs.push_back({start, kTwoPi});
        arcs.push_back({0.0, end - kTwoPi});
    } else {
        arcs.push_back({start, end});
    }
}

// Adds the arcs of the circle of radius r around the origin that lie inside
// the outline; true if the whole circle does.
bool outline_arcs(double r, const Outline& outline, std::vector<double>& crossings, std::vector<Arc>& arcs)
{
    crossings.clear();
    for (size_t i = 0; i < outline.size(); ++i) {
        const glm::dvec2& p = outline[i];
        const glm::dvec2 d = outline[(i + 1) % outline.size()] - p;
        // |p + s d| = r for s in [0, 1)
        const double a = glm::dot(d, d);
        const double b = glm::dot(p, d);
        const double c = glm::dot(p, p) - r * r;
        const double disc = b * b - a * c;
        if (a <= 0.0 || disc < 0.0) {
            continue;
        }
        const double root = std::sqrt(disc);
        for (const double s : {(-b - root) / a, (-b + root) / a}) {
            if (s >= 0.0 && s < 1.0) {
                const glm::dvec2 q = p + s * d;
                crossings.push_back(std::atan2(q.y, q.x));
            }
        }
    }
    if (crossings.empty()) {
        return inside(outline, glm::dvec2(r, 0.0));
    }
    std::sort(crossings.begin(), crossings.end());
    for (size_t i = 0; i < crossings.size(); ++i) {
        const double start = crossings[i];
        const double end = i + 1 < crossings.size() ? crossings[i + 1] : crossings[0] + kTwoPi;
        const double mid = 0.5 * (start + end);
        if (end > start && inside(outline, glm::dvec2(r * std::cos(mid), r * std::sin(mid)))) {
            add_arc(start, end, arcs);
        }
    }
    return false;
}

// Fraction of the circle of radius r around a disk's center hidden by the
// occluders: the union of their arcs.
double covered_fraction(double r, std::span<const Occluder> occluders, std::span<const Outline> outlines,
                        std::vector<double>& crossings, std::vector<Arc>& arcs)
{
    arcs.clear();
    for (const Occluder& o : occluders) {
        if (o.distance + r <= o.radius) {
            return 1.0; // the whole circle is inside it
        }
        if (o.distance >= r + o.radius || o.distance + o.radius <= r) {
            continue; // apart, or the occluder lies inside the circle without touching it
        }
        const double cos_half = (r * r + o.distance * o.distance - o.radius * o.radius) / (2.0 * r * o.distance);
        const double half = std::acos(std::clamp(cos_half, -1.0, 1.0));
        add_arc(o.angle - half, o.angle + half, arcs);
    }
    for (const Outline& outline : outlines) {
        if (outline_arcs(r, outline, crossings, arcs)) {
            return 1.0;
        }
    }
    if (arcs.empty()) {
        return 0.0;
    }
    std::sort(arcs.begin(), arcs.end());
    double covered = 0.0;
    double from = arcs[0][0];
    double to = arcs[0][1];
    for (size_t i = 1; i < arcs.size(); ++i) {
        if (arcs[i][0] > to) {
            covered += to - from;
            from = arcs[i][0];
        }
        to = std::max(to, arcs[i][1]);
    }
    covered += to - from;
    return std::min(covered / kTwoPi, 1.0);
}

// Light of a Roche surface: its facets facing the observer, each limb
// darkened, less those behind the disks in front.
double roche_flux(const StarDisk& star, const glm::dvec3& toward, const glm::dvec3& e1, const glm::dvec3& e2,
                  double wavelength_nm, double limb_darkening, std::span<const Occluder> in_front)
{
    const bool cached = wavelength_nm == kRocheRadianceWavelengthNm;
    double sum = 0.0;
    for (const RocheSurface::Facet& f : star.roche->facets) {
        const double mu = glm::dot(star.axes * f.normal, toward);
        if (mu <= 0.0) {
            continue;
        }
        if (!in_front.empty()) {
            const glm::dvec3 c = star.axes * f.center;
            const double x = glm::dot(c, e1);
            const double y = glm::dot(c, e2);
            const auto hides = [&](const Occluder& o) {
                return (x - o.x) * (x - o.x) + (y - o.y) * (y - o.y) < o.radius * o.radius;
            };
            if (std::any_of(in_front.begin(), in_front.end(), hides)) {
                continue;
            }
        }
        const double radiance = cached ? f.radiance : planck_shape(f.temperature_k, wavelength_nm);
        sum += f.area_km2 * mu * (1.0 - limb_darkening * (1.0 - mu)) * radiance;
    }
    return sum;
}

void sky_axes(const glm::dvec3& toward, glm::dvec3& e1, glm::dvec3& e2)
{
    const glm::dvec3 helper = std::abs(toward.z) < 0.9 ? glm::dvec3(0.0, 0.0, 1.0) : glm::dvec3(1.0, 0.0, 0.0);
    e1 = glm::normalize(glm::cross(helper, toward));
    e2 = glm::cross(toward, e1);
}

// A Roche surface at its brightest, over directions all around it.
double roche_brightest(const RocheSurface& surface, double wavelength_nm, double limb_darkening)
{
    StarDisk star;
    star.roche = &surface;
    double brightest = 0.0;
    for (int i = 0; i < kBrightestLatitudes; ++i) {
        const double theta = kPi * i / (kBrightestLatitudes - 1);
        const int count = i == 0 || i == kBrightestLatitudes - 1 ? 1 : kBrightestLongitudes;
        for (int j = 0; j < count; ++j) {
            const double phi = kTwoPi * j / count;
            const glm::dvec3 toward(std::cos(theta), std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi));
            glm::dvec3 e1, e2;
            sky_axes(toward, e1, e2);
            brightest = std::max(brightest, roche_flux(star, toward, e1, e2, wavelength_nm, limb_darkening, {}));
        }
    }
    return brightest;
}

// Light of disk `index`, hiding what the disks in front cover (none if `hide` is false).
double disk_flux(std::span<const StarDisk> disks, size_t index, const glm::dvec3& toward, const glm::dvec3& e1,
                 const glm::dvec3& e2, double wavelength_nm, double limb_darkening, bool hide)
{
    const StarDisk& disk = disks[index];
    std::array<Occluder, kMaxOccluders> occluders;
    std::vector<Outline> outlines;
    int count = 0;
    if (hide) {
        const double depth = glm::dot(disk.position, toward);
        for (size_t k = 0; k < disks.size() && count < kMaxOccluders; ++k) {
            const StarDisk& other = disks[k];
            const glm::dvec3 d = other.position - disk.position;
            const double x = glm::dot(d, e1);
            const double y = glm::dot(d, e2);
            const double sky = std::hypot(x, y);
            if (k == index || glm::dot(other.position, toward) <= depth || sky >= disk.radius_km + other.radius_km) {
                continue;
            }
            if (other.roche && !disk.roche) {
                outlines.push_back(roche_outline(other, disk.position, e1, e2));
            } else {
                // A Roche surface hiding another as the sphere of its volume.
                const double radius = other.roche ? other.roche->volume_radius_km : other.radius_km;
                occluders[static_cast<size_t>(count++)] = {x, y, sky, std::atan2(y, x), radius};
            }
        }
    }
    const std::span<const Occluder> in_front(occluders.data(), static_cast<size_t>(count));
    if (disk.roche) {
        return hide ? roche_flux(disk, toward, e1, e2, wavelength_nm, limb_darkening, in_front)
                    : roche_brightest(*disk.roche, wavelength_nm, limb_darkening);
    }
    std::vector<double> crossings;
    std::vector<Arc> arcs;
    const double dr = disk.radius_km / kRings;
    double sum = 0.0;
    for (int i = 0; i < kRings; ++i) {
        const double r = (i + 0.5) * dr;
        const double s = r / disk.radius_km;
        const double mu = std::sqrt(std::max(1.0 - s * s, 0.0));
        const double ring = kTwoPi * r * dr * (1.0 - limb_darkening * (1.0 - mu));
        const bool covered = count > 0 || !outlines.empty();
        sum += covered ? ring * (1.0 - covered_fraction(r, in_front, outlines, crossings, arcs)) : ring;
    }
    return sum * planck_shape(disk.temperature_k, wavelength_nm);
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
        disks[i] = star_disk_at(scene, lc.stars[i], t_tdb);
    }
    const std::span<const StarDisk> span(disks.data(), count);
    const glm::dvec3 toward = glm::normalize(scene.icrf_state_at(lc.observer, t_tdb).position - disks[0].position);
    const double flux = disks_flux(span, toward, lc.wavelength_nm, lc.limb_darkening);
    const double reference =
        lc.reference_flux > 0.0 ? lc.reference_flux : disks_full_flux(span, lc.wavelength_nm, lc.limb_darkening);
    return std::max(-2.5 * std::log10(std::max(flux / reference, 1e-12)), 0.0);
}

StarDisk star_disk_at(const Scene& scene, int star, double t_tdb)
{
    const Body& body = scene.bodies[static_cast<size_t>(star)];
    const glm::dvec3 position = scene.icrf_state_at(star, t_tdb).position;
    if (body.pulsation) {
        const Pulsation::Sample s = body.pulsation->at(t_tdb);
        return {position, s.radius_km, s.temperature_k};
    }
    if (body.roche) {
        return {position, body.equatorial_radius_km, body.temperature_k, body.roche.get(),
                scene.body_to_icrf(body, t_tdb)};
    }
    return {position, body.equatorial_radius_km, body.temperature_k};
}

void set_light_curve_reference(Scene& scene)
{
    LightCurve& lc = scene.light_curve;
    lc.reference_flux = 0.0;
    for (int star : lc.stars) {
        const Body& body = scene.bodies[static_cast<size_t>(star)];
        double brightest = 0.0;
        if (body.pulsation) {
            for (const Pulsation::Sample& s : body.pulsation->samples) {
                const StarDisk disk{glm::dvec3(0.0), s.radius_km, s.temperature_k};
                brightest = std::max(brightest, disks_full_flux(std::span(&disk, 1), lc.wavelength_nm,
                                                                lc.limb_darkening));
            }
        } else if (body.roche) {
            brightest = roche_brightest(*body.roche, lc.wavelength_nm, lc.limb_darkening);
        } else {
            const StarDisk disk{glm::dvec3(0.0), body.equatorial_radius_km, body.temperature_k};
            brightest = disks_full_flux(std::span(&disk, 1), lc.wavelength_nm, lc.limb_darkening);
        }
        lc.reference_flux += brightest;
    }
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

double light_curve_visibility(const LightCurve& curve, double rate)
{
    if (!curve.enabled()) {
        return 0.0;
    }
    constexpr double kFastSpanSeconds = 2.0;
    const double span_s = curve.span_days * kSecondsPerDay;
    const double fast = std::clamp((std::abs(rate) * kFastSpanSeconds / span_s - 0.5) * 2.0, 0.0, 1.0);
    return 1.0 - fast;
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
