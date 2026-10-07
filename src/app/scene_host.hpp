#pragma once

#include "app/captions.hpp"
#include "app/settings.hpp"
#include "platform/window.hpp"
#include "render/gpu_device.hpp"
#include "render/offscreen_output.hpp"
#include "render/render_output.hpp"
#include "scene/label_layout.hpp"
#include "scene/simulation.hpp"
#include "view/scene_renderer.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace astraxis {

// One Simulation on the tour, drawn on several outputs with the labels and
// the scene's title card as the only overlays: what the screensaver and the wallpaper hosts share. The
// host creates the windows and runs the event loop; SceneHost draws.
class SceneHost {
public:
    struct Output {
        Window window;
        RenderOutput output;       // a window with a swapchain
        OffscreenOutput offscreen; // or: read back and drawn with GDI (window is a native child)
        bool gdi = false;
        bool scene = true;  // false: drawn black
        bool active = true; // false: not drawn at all (e.g. covered); the window keeps its last frame
        OutputView view;
        LabelLayout labels;
        uint64_t frames = 0; // frames presented, for the frame rate log
    };

    bool init_gpu();
    // Takes an output whose window the host has created and claims it for the
    // GPU (or sets up the offscreen output). False if that fails.
    bool add_output(std::unique_ptr<Output> out);
    // Shuts down and drops all outputs (the scene stays loaded).
    void remove_outputs();

    // After the first outputs: loads the scene of `settings` (or `scene`),
    // starts the tour. `labels` draws them, `title` the scene's title card;
    // `error` gets the message on failure.
    bool start(const Settings& settings, const std::optional<std::string>& scene, bool labels, bool title,
               std::string* error);
    void set_fps_limit(int fps) { m_fps_limit = fps; }
    // The scene's resolution relative to each output's (window outputs only).
    void set_render_scale(float scale);

    // Advances the simulation by the real time since the last frame, draws
    // every active output and waits out the frame rate cap.
    void frame();
    bool any_active() const;
    // Restarts the real-time clock (after a pause, so the next frame does not jump).
    void reset_clock();
    void shutdown();

    std::vector<std::unique_ptr<Output>>& outputs() { return m_outputs; }
    Simulation& sim() { return m_sim; }

private:
    void render_output(Output& out);
    void draw_scene(Output& out, const Frame& frame);
    void log_frame_rates();

    GpuDevice m_gpu;
    SceneRenderer m_scene_renderer;
    bool m_renderer_ready = false;
    Simulation m_sim;
    std::vector<std::unique_ptr<Output>> m_outputs;
    ViewOptions m_view_options;
    bool m_show_labels = false;
    bool m_show_title = false;
    bool m_imgui_ready = false;
    CaptionFonts m_fonts;
    int m_fps_limit = 0;
    float m_render_scale = 1.0f;

    uint64_t m_last_counter = 0;
    uint64_t m_frame_start_ns = 0;
    uint64_t m_log_start_ns = 0;
    double m_real_dt = 0.0;
};

} // namespace astraxis
