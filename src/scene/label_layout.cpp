#include "scene/label_layout.hpp"

#include "scene/scene.hpp"
#include "scene/simulation.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr float kLabelPadding = 2.0f;     // points around each label when checking overlaps
constexpr float kLabelFadeSeconds = 0.7f; // labels ease in/out over this long
// A label that was covered comes back only once its place has stayed free this
// long, so that it does not blink as bodies graze each other.
constexpr float kLabelReshowSeconds = 1.0f;
constexpr float kPickMinAlpha = 0.3f;     // a label must be at least this visible to be picked

// True if the segment from the camera (origin) to `target` passes through a sphere.
bool occluded(const glm::dvec3& target, const glm::dvec3& center, double radius)
{
    const double dist = glm::length(target);
    const glm::dvec3 dir = target / dist;
    const double along = glm::dot(center, dir);
    if (along <= 0.0 || along >= dist) {
        return false;
    }
    const double perp2 = glm::dot(center, center) - along * along;
    return perp2 < radius * radius;
}

} // namespace

void LabelLayout::update(const Scene& scene, const OutputView& view, int focus, glm::vec2 screen_size,
                         float font_size, const MeasureText& measure_text, float dt)
{
    const CameraView& camera = view.camera;
    const glm::dvec3& cam = camera.position;
    m_markers.clear();
    m_bodies.clear();

    auto project = [&](const glm::dvec3& rel, glm::vec2* screen, float* w) {
        const glm::vec4 clip = camera.view_proj * glm::vec4(glm::vec3(rel), 1.0f);
        if (clip.w <= 0.0f) {
            return false;
        }
        *screen = glm::vec2((clip.x / clip.w * 0.5f + 0.5f) * screen_size.x,
                            (0.5f - clip.y / clip.w * 0.5f) * screen_size.y);
        *w = clip.w;
        return true;
    };

    for (size_t k = 0; k < scene.markers.size(); ++k) {
        glm::vec2 s;
        float w = 0.0f;
        if (project(scene.markers[k].world_position - cam, &s, &w)) {
            m_markers.push_back({static_cast<int>(k), s});
        }
    }

    struct Candidate {
        int body;
        int rank;
        glm::vec2 screen;
        float radius_px;
    };
    std::vector<Candidate> candidates;
    m_labels.resize(scene.bodies.size());
    for (size_t i = 0; i < scene.bodies.size(); ++i) {
        const Body& body = scene.bodies[i];
        const float fade = i < view.body_fades.size() ? view.body_fades[i] : 1.0f;
        if (!body.visible || body.kind == BodyKind::Barycenter || fade <= 0.0f ||
            (!body.label && static_cast<int>(i) != focus)) {
            m_labels[i] = {};
            continue;
        }
        const glm::dvec3 rel = body.world_position - cam;

        bool hidden = false;
        for (size_t k = 0; k < scene.bodies.size() && !hidden; ++k) {
            const Body& other = scene.bodies[k];
            if (k != i && other.visible && other.kind != BodyKind::Spacecraft &&
                other.kind != BodyKind::Barycenter && other.kind != BodyKind::Ghost) {
                hidden = occluded(rel, other.world_position - cam, other.equatorial_radius_km);
            }
        }
        glm::vec2 screen;
        float w = 0.0f;
        if (hidden || !project(rel, &screen, &w)) {
            m_labels[i] = {};
            continue;
        }
        const float radius_px =
            static_cast<float>(body.equatorial_radius_km) / w * camera.focal_y * 0.5f * screen_size.y;
        const int host = scene.satellite_host(static_cast<int>(i));
        int rank = 4; // moons
        if (static_cast<int>(i) == focus) {
            rank = 0;
        } else if (body.kind == BodyKind::Star) {
            rank = 1;
        } else if (body.kind == BodyKind::Spacecraft) {
            rank = 3;
        } else if (body.kind == BodyKind::Ghost) {
            rank = 5;
        } else if (host < 0 || scene.bodies[static_cast<size_t>(host)].kind == BodyKind::Star) {
            rank = 2;
        }
        candidates.push_back({static_cast<int>(i), rank, screen, radius_px});
    }
    std::sort(candidates.begin(), candidates.end(), [&](const Candidate& a, const Candidate& b) {
        if (a.rank != b.rank) {
            return a.rank < b.rank;
        }
        return scene.bodies[static_cast<size_t>(a.body)].equatorial_radius_km >
               scene.bodies[static_cast<size_t>(b.body)].equatorial_radius_km;
    });

    std::vector<glm::vec4> placed; // label rectangles: min x, min y, max x, max y
    const float step = dt / kLabelFadeSeconds;
    for (const Candidate& c : candidates) {
        const size_t i = static_cast<size_t>(c.body);
        const Body& body = scene.bodies[i];
        const float fade = i < view.body_fades.size() ? view.body_fades[i] : 1.0f;

        const float offset = std::max(c.radius_px, kMarkerRadius) * 0.7071f + 4.0f;
        const glm::vec2 text_pos(c.screen.x + offset, c.screen.y - offset - font_size * 0.5f);
        const glm::vec2 text_size = measure_text(body.name);
        const glm::vec4 rect(text_pos.x - kLabelPadding, text_pos.y - kLabelPadding,
                             text_pos.x + text_size.x + kLabelPadding, text_pos.y + text_size.y + kLabelPadding);
        bool overlaps = false;
        for (const glm::vec4& r : placed) {
            if (rect.x < r.z && r.x < rect.z && rect.y < r.w && r.y < rect.w) {
                overlaps = true;
                break;
            }
        }
        if (!overlaps && fade > 0.05f) {
            placed.push_back(rect);
        }
        LabelState& s = m_labels[i];
        s.wait = overlaps ? kLabelReshowSeconds : std::max(s.wait - dt, 0.0f);
        const bool show = !overlaps && s.wait <= 0.0f;
        s.alpha = std::clamp(s.alpha + (show ? step : -step), 0.0f, 1.0f);

        BodyMark mark;
        mark.body = c.body;
        mark.position = c.screen;
        mark.radius = c.radius_px;
        mark.fade = fade;
        mark.text_position = text_pos;
        mark.text_alpha = s.alpha * fade;
        // Spacecraft and ghosts are never drawn as bodies, so they always get a dot.
        if (body.kind == BodyKind::Ghost) {
            mark.dot = Dot::Hollow;
        } else if (c.radius_px < kMarkerRadius || body.kind == BodyKind::Spacecraft) {
            mark.dot = Dot::Filled;
        }
        m_bodies.push_back(mark);
    }
}

int LabelLayout::pick(glm::vec2 point, float pick_radius) const
{
    int picked = -1;
    for (const BodyMark& mark : m_bodies) {
        if (mark.text_alpha > kPickMinAlpha && glm::length(point - mark.position) < std::max(mark.radius, pick_radius)) {
            picked = mark.body;
        }
    }
    return picked;
}

} // namespace astraxis
