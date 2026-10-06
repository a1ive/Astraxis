#pragma once

#include "app/scene_host.hpp"
#include "platform/power_monitor.hpp"

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
//
// Power (/s): nothing is drawn while the display is off; on battery the frame
// rate is capped or drawing stops ([screensaver] battery); [screensaver]
// render_scale lowers the scene's resolution.
class Screensaver {
public:
    bool init(const ScreensaverOptions& options);
    void run();
    void shutdown();

private:
    bool create_outputs(int display);
    void handle_events();
    void quit(const char* reason);
    // From the power state: whether to draw, and at what frame rate.
    void update_activity();

    ScreensaverOptions m_options;
    SceneHost m_host;
    PowerMonitor m_power;
    BatteryPolicy m_battery = BatteryPolicy::Limit;
    int m_fps = 0;       // the cap on AC power
    bool m_idle = false; // nothing is drawn
    uint64_t m_last_check_ns = 0;
    bool m_sdl_ready = false;
    bool m_running = false;
    uint64_t m_start_ns = 0;

    // Mouse position (desktop pixels) at the first motion event.
    bool m_have_mouse = false;
    float m_mouse_x = 0.0f;
    float m_mouse_y = 0.0f;
};

} // namespace astraxis
