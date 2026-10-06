#pragma once

#include "app/scene_host.hpp"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>

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
    bool create_outputs(int display);
    void handle_events();
    void quit(const char* reason);

    ScreensaverOptions m_options;
    SceneHost m_host;
    bool m_sdl_ready = false;
    bool m_running = false;
    uint64_t m_start_ns = 0;

    // Mouse position (desktop pixels) at the first motion event.
    bool m_have_mouse = false;
    float m_mouse_x = 0.0f;
    float m_mouse_y = 0.0f;
};

} // namespace astraxis
