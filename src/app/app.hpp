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
#include <future>
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
    // False if the first scene could not be loaded.
    bool run();
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
    void draw_light_curve_plot();
    void draw_loading_spinner(); // while a scene is read in the background
    void toggle_caption(); // dismisses the event caption, or shows it again
    void set_msaa(int samples);
    void set_render_scale(float scale);
    void render();

    // Scenes are read on a worker thread (the file, its integrations and its
    // maps) while the window keeps drawing the previous one; the latest request wins.
    struct SceneLoad {
        bool ok = false;
        Scene scene;
        SceneImages images;
        std::string error;
    };
    void request_scene(size_t index);
    void poll_scene_load(); // starts the scene once it has been read
    void start_scene(SceneLoad&& load);
    bool scene_loading() const { return m_loading.valid(); }
    size_t requested_scene() const; // the one loading or queued, else the one shown
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
    int m_msaa = 4;              // as chosen (the device may use fewer samples)
    float m_render_scale = 1.0f;

    std::filesystem::path m_asset_dir;
    std::vector<std::filesystem::path> m_scene_files;
    std::vector<std::string> m_scene_names; // each file's own `name`, for the scene list
    size_t m_scene_index = 0; // the scene shown
    std::string m_scene_error;
    bool m_has_scene = false; // false until the first scene has loaded: only the spinner is drawn
    std::future<SceneLoad> m_loading;
    size_t m_loading_index = 0;
    std::optional<size_t> m_next_scene; // requested while another was loading
    double m_loading_age = 0.0;         // seconds since the load started (the spinner fades in)
    int m_start_event = 0;              // the command line's --event, for the first scene

    bool m_initialized = false; // init succeeded: save the settings on shutdown
    bool m_running = false;
    bool m_imgui_ready = false;
    CaptionFonts m_fonts; // the scene title card and the event captions
    bool m_dragging = false;
    uint64_t m_last_counter = 0;

    float m_time_bar_top = 0.0f;          // where the time bar was last drawn (event captions stand above it)
    std::vector<ScreenRect> m_panel_rects; // the side panels drawn this frame (the title card moves aside)

    // UI options.
    bool m_show_ui = true;
    bool m_show_labels = true;
    bool m_show_captions = true; // event captions (E dismisses one)
    bool m_show_info = true;
    bool m_auto_tour = true;     // start the tour after kIdleSeconds without input
    double m_idle_time = 0.0;    // seconds since the last input
    double m_pointer_idle = 0.0; // seconds since the mouse last moved
    bool m_cursor_hidden = false;

    // "Go to date" input.
    char m_goto_text[32] = "";
    bool m_goto_error = false;
};

} // namespace astraxis
