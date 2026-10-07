#include "scene/label_layout.hpp"

#include "scene/scene.hpp"
#include "scene/simulation.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr float kLabelPadding = 2.0f;     // points around each label when checking overlaps
constexpr float kLabelFadeSeconds = 0.3f; // labels ease in/out over this long
constexpr float kPickMinAlpha = 0.3f;     // a label must be at least this visible to be picked
constexpr float kLabelGap = 4.0f;         // points between a body (or its dot) and its label
// A label moves back to a more preferred anchor only after it has been free
// this long, so that labels do not hop back and forth as bodies pass; it
// slides to a new anchor over kAnchorSlideSeconds.
constexpr float kAnchorSettleSeconds = 0.5f;
constexpr float kAnchorSlideSeconds = 0.2f;

// Label anchors around a body (screen y down), in order of preference: the
// diagonals first (NE, SE, NW, SW; a body at the right edge keeps its label at
// the same height in NW), then E, W, N, S.
constexpr int kAnchorCount = 8;
constexpr glm::ivec2 kAnchors[kAnchorCount] = {{1, -1}, {1, 1}, {-1, -1}, {-1, 1}, {1, 0}, {-1, 0}, {0, -1}, {0, 1}};

// Top left of a label of `size` at `anchor` around a body at `center` with
// on-screen radius `radius`. Labels to the side are centered on the font's
// line height, diagonal ones on the 45 degree point of the body's rim.
glm::vec2 anchor_position(int anchor, glm::vec2 center, float radius, glm::vec2 size, float font_size)
{
    const glm::ivec2 dir = kAnchors[anchor];
    const float r = std::max(radius, LabelLayout::kMarkerRadius);
    const float off = dir.x != 0 && dir.y != 0 ? r * 0.7071f + kLabelGap : r + kLabelGap;
    glm::vec2 pos;
    if (dir.x == 0) {
        pos.x = center.x - 0.5f * size.x;
        pos.y = dir.y < 0 ? center.y - off - size.y : center.y + off;
    } else {
        pos.x = dir.x > 0 ? center.x + off : center.x - off - size.x;
        pos.y = center.y + static_cast<float>(dir.y) * off - 0.5f * font_size;
    }
    return pos;
}

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
            m_labels[i].alpha = 0.0f;
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
            m_labels[i].alpha = 0.0f;
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

        // Each anchor's rectangle: free if it overlaps no placed label, usable
        // if also on screen (or if no free anchor is on screen, e.g. for a
        // body at the edge).
        const glm::vec2 text_size = measure_text(body.name);
        glm::vec2 positions[kAnchorCount];
        glm::vec4 rects[kAnchorCount];
        bool free[kAnchorCount];
        bool on_screen[kAnchorCount];
        bool any_on_screen = false;
        for (int a = 0; a < kAnchorCount; ++a) {
            const glm::vec2 p = anchor_position(a, c.screen, c.radius_px, text_size, font_size);
            positions[a] = p;
            rects[a] = glm::vec4(p.x - kLabelPadding, p.y - kLabelPadding, p.x + text_size.x + kLabelPadding,
                                 p.y + text_size.y + kLabelPadding);
            free[a] = true;
            for (const glm::vec4& r : placed) {
                if (rects[a].x < r.z && r.x < rects[a].z && rects[a].y < r.w && r.y < rects[a].w) {
                    free[a] = false;
                    break;
                }
            }
            on_screen[a] = p.x >= 0.0f && p.y >= 0.0f && p.x + text_size.x <= screen_size.x &&
                           p.y + text_size.y <= screen_size.y;
            any_on_screen = any_on_screen || (free[a] && on_screen[a]);
        }
        auto usable = [&](int a) { return free[a] && (on_screen[a] || !any_on_screen); };
        int best = -1;
        for (int a = 0; a < kAnchorCount && best < 0; ++a) {
            best = usable(a) ? a : -1;
        }

        // Keep the current anchor while it is usable; move at once when it is
        // not, and back to a better one only once that has stayed free.
        LabelState& s = m_labels[i];
        auto move_to = [&](int a) {
            s.from = s.anchor;
            s.anchor = a;
            s.slide = 0.0f;
            s.better_for = 0.0f;
        };
        if (s.anchor < 0 || s.alpha <= 0.0f) {
            s.anchor = best >= 0 ? best : 0; // not shown yet: no need to slide
            s.from = s.anchor;
            s.slide = 1.0f;
            s.better_for = 0.0f;
        } else if (!usable(s.anchor)) {
            if (best >= 0) {
                move_to(best);
            }
        } else if (best >= 0 && best < s.anchor) {
            s.better_for += dt;
            if (s.better_for >= kAnchorSettleSeconds) {
                move_to(best);
            }
        } else {
            s.better_for = 0.0f;
        }
        s.slide = std::min(s.slide + dt / kAnchorSlideSeconds, 1.0f);

        const bool overlaps = !free[s.anchor];
        if (!overlaps && fade > 0.05f) {
            placed.push_back(rects[s.anchor]);
        }
        s.alpha = std::clamp(s.alpha + (overlaps ? -step : step), 0.0f, 1.0f);

        BodyMark mark;
        mark.body = c.body;
        mark.position = c.screen;
        mark.radius = c.radius_px;
        mark.fade = fade;
        const float t = s.slide * s.slide * (3.0f - 2.0f * s.slide);
        mark.text_position = glm::mix(positions[s.from], positions[s.anchor], t);
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
