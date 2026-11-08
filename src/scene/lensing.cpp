#include "scene/lensing.hpp"

#include "core/math.hpp"
#include "ephem/kerr_null.hpp"
#include "scene/camera.hpp"
#include "scene/scene.hpp"

#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/matrix.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

// A body closer to the camera than this fraction of its distance from the hole is
// not lensed: the light from it barely bends (and the solve would lose precision).
constexpr double kNearCamera = 1e-4;
// The camera must be outside the ergosphere for its static frame (r > 2 M).
constexpr double kMinCameraRadius = 3.0;
constexpr double kToleranceInPixels = 0.05;
// A secondary image this faint (relative to the body unlensed) is left out: its label
// dot would be under 5% opaque (LabelLayout), a drawn body under 5% as wide.
constexpr double kMinSecondaryMagnification = 0.0025;
// A body whose image could move less than this from where it is is not lensed (nor its
// secondary image looked for, which is then fainter than (shift / angle from the hole)^2).
constexpr double kMinShiftPixels = 0.1;
// Frames between searches for a secondary image that was not found.
constexpr int kSearchInterval = 30;
// Depth margins about the hole's: in front of it, or behind it.
constexpr double kFrontOfHole = 0.999;
constexpr double kBehindHole = 1.001;

bool orbits(const Scene& scene, int body, int hole)
{
    for (int p = scene.bodies[static_cast<size_t>(body)].parent; p >= 0;
         p = scene.bodies[static_cast<size_t>(p)].parent) {
        if (p == hole) {
            return true;
        }
    }
    return false;
}

// Rotation taking unit vector a onto unit vector b.
glm::dmat3 rotation_between(const glm::dvec3& a, const glm::dvec3& b)
{
    const glm::dvec3 v = glm::cross(a, b);
    const double s = glm::length(v);
    const double c = glm::dot(a, b);
    if (s < 1e-12) {
        if (c > 0.0) {
            return glm::dmat3(1.0);
        }
        const glm::dvec3 helper = std::abs(a.z) < 0.9 ? glm::dvec3(0.0, 0.0, 1.0) : glm::dvec3(1.0, 0.0, 0.0);
        const glm::dvec3 p = glm::normalize(glm::cross(a, helper));
        return 2.0 * glm::outerProduct(p, p) - glm::dmat3(1.0); // half a turn about p
    }
    return glm::dmat3(glm::rotate(glm::dmat4(1.0), std::atan2(s, c), v / s));
}

} // namespace

glm::dmat3 black_hole_axes(const Scene& scene, const Body& hole)
{
    const glm::dvec3 z = glm::normalize(scene.current_transform().direction_to_display(hole.pole));
    const glm::dvec3 helper = std::abs(z.z) < 0.9 ? glm::dvec3(0.0, 0.0, 1.0) : glm::dvec3(1.0, 0.0, 0.0);
    const glm::dvec3 x = glm::normalize(glm::cross(helper, z));
    return glm::dmat3(x, glm::cross(z, x), z);
}

void lens_bodies(const Scene& scene, const CameraView& view, std::vector<BodyLens>& lens)
{
    lens.resize(scene.bodies.size());
    int hole_index = -1;
    for (size_t i = 0; i < scene.bodies.size(); ++i) {
        const Body& b = scene.bodies[i];
        if (b.kind == BodyKind::BlackHole && b.visible && b.gm_km3_s2 > 0.0) {
            hole_index = static_cast<int>(i);
            break; // one hole per scene, as the ray tracer
        }
    }
    if (hole_index < 0) {
        std::fill(lens.begin(), lens.end(), BodyLens{});
        return;
    }
    const Body& hole = scene.bodies[static_cast<size_t>(hole_index)];
    const double m_km = hole.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
    const double a = hole.spin;
    const glm::dmat3 axes = black_hole_axes(scene, hole);
    const glm::dmat3 to_hole = glm::transpose(axes);
    const glm::dvec3 cam = to_hole * (view.position - hole.world_position) / m_km;
    const double hole_distance = glm::length(view.position - hole.world_position);
    const bool camera_ok = kerr_schild_r(cam, a) > kMinCameraRadius;

    // The camera's static frame, with the view's axes as the ray tracer builds it: a ray's
    // coordinate direction d is seen in local direction n = normalize(to_local d), drawn
    // where the view's Euclidean axes put n.
    const glm::dvec3 right(view.view[0][0], view.view[1][0], view.view[2][0]);
    const glm::dvec3 up(view.view[0][1], view.view[1][1], view.view[2][1]);
    const glm::dvec3 forward(-view.view[0][2], -view.view[1][2], -view.view[2][2]);
    const glm::dmat3 screen(right, up, forward);
    glm::dmat3 to_local(1.0);
    if (camera_ok) {
        const ObserverFrame obs = static_observer_frame(cam, a, to_hole * forward, to_hole * up, to_hole * right);
        to_local = glm::inverse(glm::dmat3(obs.right_dir, obs.up_dir, obs.forward_dir));
    }

    LensSettings settings;
    settings.spin = a;
    settings.disk_r_in = std::max(isco_radius(a), hole.disk_inner_m);
    settings.disk_r_out = hole.disk_outer_m > settings.disk_r_in && hole.disk_brightness > 0.0 ? hole.disk_outer_m : 0.0;
    const double pixel = 2.0 / (static_cast<double>(view.focal_y) * view.viewport.y);
    settings.tolerance_rad = kToleranceInPixels * pixel;
    settings.min_secondary_magnification = kMinSecondaryMagnification;

    for (size_t i = 0; i < scene.bodies.size(); ++i) {
        BodyLens& out = lens[i];
        const Body& body = scene.bodies[i];
        if (!camera_ok || !body.visible || body.kind == BodyKind::Barycenter ||
            !orbits(scene, static_cast<int>(i), hole_index)) {
            out = BodyLens{};
            continue;
        }
        const glm::dvec3 target = to_hole * (body.world_position - hole.world_position) / m_km;
        const glm::dvec3 path = target - cam;
        const double path_length = glm::length(path);
        if (path_length < kNearCamera * glm::length(target)) {
            out = BodyLens{};
            continue;
        }
        // Weak field: light passing the hole at b bends by 4 M / b in all, and the image
        // moves by at most that, less when the body is near the camera. The straight
        // line's distance from the hole stands for b (it is never larger).
        const double b = glm::length(glm::cross(cam, path)) / path_length;
        if (4.0 / b * std::min(1.0, path_length / glm::length(cam)) < kMinShiftPixels * pixel) {
            out = BodyLens{};
            continue;
        }
        const LensImage guesses[2] = {out.images[0].solution, out.images[1].solution};
        LensImage solved[2];
        LensSettings body_settings = settings;
        body_settings.search_secondary = !out.lensed || out.search_wait <= 0;
        const int found = find_lens_images(cam, target, body_settings, guesses, out.lensed ? out.found : 0, solved);
        if (!(found & 1)) {
            out = BodyLens{};
            continue;
        }
        // A secondary image that a search did not find (or was too faint to look for) is
        // looked for again only now and then: the search takes dozens of rays.
        out.search_wait = (found & 2) ? 0 : body_settings.search_secondary ? kSearchInterval : out.search_wait - 1;
        out.lensed = true;
        out.found = found;
        const double distance = glm::length(body.world_position - view.position);
        for (int k = 0; k < 2; ++k) {
            if (!(found & (1 << k))) {
                continue;
            }
            const LensImage& s = solved[k];
            LensedView& image = out.images[k];
            image.solution = s;

            const glm::dvec3 v = to_local * s.direction;
            const double v_len = glm::length(v);
            const glm::dvec3 n = v / v_len;
            const glm::dvec3 dir = glm::normalize(screen * n);
            // The ellipse's axes on screen, as angles per km of body radius.
            glm::dvec3 w[2];
            double angle[2];
            for (int j = 0; j < 2; ++j) {
                const glm::dvec3 t = to_local * s.axis[j];
                w[j] = screen * ((t - glm::dot(t, n) * n) * (s.angle_per_m[j] / v_len));
                angle[j] = glm::length(w[j]) / m_km;
            }
            glm::dvec3 axis0 = glm::normalize(w[0] - glm::dot(w[0], dir) * dir);
            glm::dvec3 axis1 = glm::normalize(glm::cross(dir, axis0));

            double depth = distance;
            if (s.disk_first) {
                depth = std::max(depth, kBehindHole * hole_distance);
            } else {
                depth = std::min(depth, kFrontOfHole * hole_distance);
            }
            const double s0 = depth * angle[0];
            const double s1 = depth * angle[1];
            image.offset = depth * dir;
            image.stretch = std::sqrt(s0 * s1) * glm::outerProduct(dir, dir) + s0 * glm::outerProduct(axis0, axis0) +
                            s1 * glm::outerProduct(axis1, axis1);
            image.turn = rotation_between(glm::normalize(axes * s.emission), dir);
            image.angular_radius_per_km = std::sqrt(angle[0] * angle[1]);
            image.magnification = s.magnification;
            image.behind_disk = s.disk_first;
        }
    }
}

} // namespace astraxis
