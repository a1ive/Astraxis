#pragma once

#include "platform/window.hpp"
#include "render/gpu_device.hpp"
#include "render/render_output.hpp"
#include "scene/label_layout.hpp"
#include "scene/simulation.hpp"
#include "view/scene_renderer.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

union SDL_Event;

namespace astraxis {

// Command line: --scene <file stem> (default solar_system), --event <n> (1-based,
// in the order of the Events list) to start at that event.
struct LaunchOptions {
    std::string scene = "solar_system";
    int event = 0; // 0 = none
};

// The windowed host: an SDL window with mouse/keyboard camera control and the
// ImGui panels around a Simulation drawn by SceneRenderer.
class App {
public:
    bool init(const LaunchOptions& options);
    void run();
    void shutdown();

private:
    void handle_events();
    void handle_key(const SDL_Event& event);
    void update(double real_dt);
    void build_ui();
    void build_control_panel();
    void build_time_bar();   // playback controls at the bottom
    void build_info_panel(); // details of the focused body
    void build_labels();
    void render();

    bool load_scene(size_t index);

    // Auto tour (idle mode).
    void start_tour();
    // Any user input; `takes_camera` stops the tour (camera control returns to the user).
    void note_activity(bool takes_camera);
    bool panel_visible() const;

    Window m_window;
    GpuDevice m_gpu;
    RenderOutput m_output;
    SceneRenderer m_scene_renderer;
    Simulation m_sim;
    OutputView m_view; // the window's
    ViewOptions m_view_options;
    LabelLayout m_labels;

    std::filesystem::path m_asset_dir;
    std::vector<std::filesystem::path> m_scene_files;
    size_t m_scene_index = 0;
    std::string m_scene_error;

    bool m_running = false;
    bool m_imgui_ready = false;
    bool m_dragging = false;
    uint64_t m_last_counter = 0;

    // UI options.
    bool m_show_ui = true;
    bool m_show_labels = true;
    bool m_show_info = true;
    bool m_show_demo = false;
    bool m_auto_tour = true;     // start the tour after kIdleSeconds without input
    double m_idle_time = 0.0;    // seconds since the last input
    double m_pointer_idle = 0.0; // seconds since the mouse last moved
    bool m_cursor_hidden = false;

    // "Go to date" input.
    char m_goto_text[32] = "";
    bool m_goto_error = false;
};

} // namespace astraxis
