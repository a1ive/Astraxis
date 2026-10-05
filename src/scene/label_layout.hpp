#pragma once

#include <glm/vec2.hpp>

#include <functional>
#include <string>
#include <vector>

namespace astraxis {

class Scene;
struct OutputView;

// Where the on-screen annotations go: scene markers (e.g. Lagrange points),
// body markers and body labels. Labels are placed by priority (the focus,
// stars, bodies orbiting a star, spacecraft, moons, then ghosts, larger bodies
// first); a label that would overlap one already placed fades out, and labels
// ease in and out so that they do not flicker as bodies pass each other.
//
// Screen coordinates are in "points" from the top left, in whatever unit the
// host draws text with (`screen_size`, `font_size` and `measure_text` agree);
// the layout itself draws nothing.
class LabelLayout {
public:
    static constexpr float kMarkerRadius = 2.5f; // a body smaller than this on screen gets a dot

    struct SceneMarker {
        int marker;         // Scene::markers index
        glm::vec2 position; // on screen
    };

    enum class Dot { None, Filled, Hollow };

    struct BodyMark {
        int body;
        glm::vec2 position;  // body center on screen
        float radius = 0.0f; // body radius on screen
        Dot dot = Dot::None; // spacecraft and bodies under kMarkerRadius: filled; ghosts: hollow
        float fade = 1.0f;   // the body's fade (satellite_fades); the dot's opacity
        glm::vec2 text_position{0.0f}; // top left of the label
        float text_alpha = 0.0f;       // label opacity (eased, times the fade)
    };

    // Width and height of `text` on screen.
    using MeasureText = std::function<glm::vec2(const std::string& text)>;

    // Lays out one frame. `dt`: real seconds since the last update (label easing).
    void update(const Scene& scene, const OutputView& view, int focus, glm::vec2 screen_size, float font_size,
                const MeasureText& measure_text, float dt);
    // Forgets the eased label opacities (e.g. on a scene change).
    void reset() { m_alpha.clear(); }

    const std::vector<SceneMarker>& markers() const { return m_markers; }
    // In placement order (highest priority first).
    const std::vector<BodyMark>& bodies() const { return m_bodies; }

    // The body whose mark is under `point` (within its disk or `pick_radius`)
    // and whose label is shown enough to be picked; -1 if none. The last such
    // mark in placement order wins.
    int pick(glm::vec2 point, float pick_radius) const;

private:
    std::vector<float> m_alpha; // per body: label opacity, eased toward its target
    std::vector<SceneMarker> m_markers;
    std::vector<BodyMark> m_bodies;
};

} // namespace astraxis
