#pragma once

#include "app/captions.hpp"
#include "app/settings.hpp"
#include "platform/window.hpp"
#include "render/gpu_device.hpp"
#include "render/render_output.hpp"
#include "scene/label_layout.hpp"
#include "scene/simulation.hpp"
#include "view/scene_renderer.hpp"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

union SDL_Event;

namespace astraxis {

// Command line; each option overrides config.toml for this run only.
struct LaunchOptions {
    std::filesystem::path config; // empty: config.toml next to the executable
    std::optional<std::string> scene; // file stem
    int event = 0; // 1-based, in the order of the Events list; 0 = none
    std::optional<bool> fullscreen;
    std::optional<int> display; // 0 = primary, n = the n-th display
    std::optional<int> fps;     // frame rate cap, 0 = none
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
    void draw_captions();  // the scene title card and the last event's caption
    void toggle_caption(); // dismisses the event caption, or shows it again
    void render();

    bool load_scene(size_t index);
    void save_settings();

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

    // config.toml: the [window] settings as read, and what this run started
    // with after the command line. Only what the user changes in the app is
    // saved, not the command line's overrides.
    std::filesystem::path m_config_path;
    Settings m_settings;
    bool m_config_ok = false; // false: the file could not be parsed; never overwrite it
    std::string m_start_scene;
    uint32_t m_start_display = 0;
    bool m_start_fullscreen = false;
    int m_fps_limit = 0;
    uint64_t m_frame_start_ns = 0;

    std::filesystem::path m_asset_dir;
    std::vector<std::filesystem::path> m_scene_files;
    std::vector<std::string> m_scene_names; // each file's own `name`, for the scene list
    size_t m_scene_index = 0;
    std::string m_scene_error;

    bool m_initialized = false; // init succeeded: save the settings on shutdown
    bool m_running = false;
    bool m_imgui_ready = false;
    CaptionFonts m_fonts; // the scene title card and the event captions
    bool m_dragging = false;
    uint64_t m_last_counter = 0;

    float m_time_bar_top = 0.0f; // where the time bar was last drawn (event captions stand above it)

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
