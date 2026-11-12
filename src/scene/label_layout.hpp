#pragma once

#include <glm/vec2.hpp>
#include <glm/vec4.hpp>

#include <functional>
#include <string>
#include <vector>

namespace astraxis {

class Scene;
struct OutputView;

// Where the on-screen annotations go: scene markers (e.g. Lagrange points),
// body markers and body labels. Labels are placed by priority (a spacecraft
// whose burn's delta-v is shown below its name, the focus,
// stars, bodies orbiting a star, spacecraft, moons, then ghosts, larger bodies
// first); a label that would overlap one already placed fades out, and comes
// back only after its place has stayed free for a moment (it keeps that place
// meanwhile). Labels never move relative to their body; they only fade, slowly,
// so that they do not draw the eye as bodies pass each other.
//
// Screen coordinates are in "points" from the top left, in whatever unit the
// host draws text with (`screen_size`, `font_size` and `measure_text` agree);
// the layout itself draws nothing.
class LabelLayout {
public:
    static constexpr float kMarkerRadius = 2.5f; // a body smaller than this on screen gets a dot
    // A second line below a label (Scene::burn_delta_v): its offset from the
    // name's top and the room it takes, in font sizes.
    static constexpr float kReadoutLineSpacing = 1.15f;
    static constexpr float kReadoutWidthEm = 5.5f;

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
        glm::vec4 text_bounds{0.0f};   // padded bounds of the label
        float text_alpha = 0.0f;       // label opacity (eased, times the fade)
        // A body bent around a black hole is marked at its primary image (`position`);
        // a filled dot also goes at its secondary image, at this opacity relative to the
        // first's (0: none).
        glm::vec2 image2{0.0f};
        float image2_alpha = 0.0f;
    };

    // Width and height of `text` on screen.
    using MeasureText = std::function<glm::vec2(const std::string& text)>;

    // Lays out one frame. `dt`: real seconds since the last update (label easing).
    void update(const Scene& scene, const OutputView& view, int focus, glm::vec2 screen_size, float font_size,
                const MeasureText& measure_text, float dt);
    // Forgets the eased label opacities (e.g. on a scene change).
    void reset() { m_labels.clear(); }

    const std::vector<SceneMarker>& markers() const { return m_markers; }
    // In placement order (highest priority first).
    const std::vector<BodyMark>& bodies() const { return m_bodies; }

    // Pick a visible label under `point`, then fall back to the body disk or
    // `pick_radius`. The last eligible mark in placement order wins in each
    // case. Returns -1 when there is no hit.
    int pick(glm::vec2 point, float pick_radius) const;

private:
    struct LabelState {
        float alpha = 0.0f; // opacity, eased toward its target
        float wait = 0.0f;  // seconds before a label that was covered may show again
    };
    std::vector<LabelState> m_labels; // per body
    std::vector<SceneMarker> m_markers;
    std::vector<BodyMark> m_bodies;
};

} // namespace astraxis
