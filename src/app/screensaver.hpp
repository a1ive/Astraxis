#pragma once

#include "app/settings.hpp"
#include "platform/window.hpp"
#include "render/gpu_device.hpp"
#include "render/offscreen_output.hpp"
#include "render/render_output.hpp"
#include "scene/label_layout.hpp"
#include "scene/simulation.hpp"
#include "view/scene_renderer.hpp"

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace astraxis {

struct ScreensaverOptions {
    bool preview = false;   // /p: draw into `parent` until it goes away
    void* parent = nullptr; // preview: the HWND given by Windows
    std::filesystem::path config; // empty: config.toml next to the executable
    // Command line overrides of [screensaver], for this run only.
    std::optional<std::string> scene;
    std::optional<int> display;
    std::optional<int> fps;
};

// The screensaver host (/s and /p). /s covers every display with a window,
// draws the [screensaver] scene on the chosen display(s) and black on the
// others, runs the auto tour, and quits on any key, click or mouse movement.
// /p draws a small view into the preview window of the Screen Saver Settings
// dialog until that window is destroyed.
class Screensaver {
public:
    bool init(const ScreensaverOptions& options);
    void run();
    void shutdown();

private:
    struct Output {
        Window window;
        RenderOutput output;       // a display's window
        OffscreenOutput offscreen; // the preview, drawn with GDI
        bool preview = false;
        OutputView view;
        LabelLayout labels;
        bool scene = true;   // false: drawn black (not a chosen display)
        uint64_t frames = 0; // frames presented, for the frame rate log
    };

    bool create_outputs(int display);
    void handle_events();
    void quit(const char* reason);
    void render();
    void render_output(Output& out, bool wait);
    // Records the scene and the labels into a frame of `out`.
    void draw_scene(Output& out, const Frame& frame);
    void log_frame_rates();

    ScreensaverOptions m_options;
    GpuDevice m_gpu;
    SceneRenderer m_scene_renderer;
    Simulation m_sim;
    std::vector<std::unique_ptr<Output>> m_outputs;
    ViewOptions m_view_options;
    bool m_show_labels = true;
    int m_fps_limit = 0;

    bool m_running = false;
    bool m_imgui_ready = false;
    uint64_t m_start_ns = 0;
    uint64_t m_last_counter = 0;
    uint64_t m_frame_start_ns = 0;
    uint64_t m_log_start_ns = 0;
    double m_real_dt = 0.0;

    // Mouse position (desktop pixels) at the first motion event.
    bool m_have_mouse = false;
    float m_mouse_x = 0.0f;
    float m_mouse_y = 0.0f;
};

} // namespace astraxis
