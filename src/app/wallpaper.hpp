#pragma once

#include "app/scene_host.hpp"
#include "platform/power_monitor.hpp"
#include "platform/tray_icon.hpp"
#include "platform/wallpaper_layer.hpp"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace astraxis {

struct WallpaperOptions {
    std::filesystem::path config; // empty: config.toml next to the executable
    // Command line overrides of [wallpaper], for this run only.
    std::optional<std::string> scene;
    std::optional<int> display;
    std::optional<int> fps;
};

// The live wallpaper host (Windows): a window per chosen display, placed
// behind the desktop icons (WallpaperLayer), running the auto tour. A tray
// icon pauses it or ends it. Re-attaches when Explorer restarts or the
// displays change; repaints Explorer's wallpaper on exit.
//
// Power: nothing is drawn while the session is locked or the display is
// off; a display covered by a maximized or full-screen window is not drawn
// ([wallpaper] pause_covered); on battery the frame rate is capped or
// drawing stops ([wallpaper] battery); [wallpaper] render_scale lowers the
// scene's resolution.
class Wallpaper {
public:
    bool init(const WallpaperOptions& options);
    void run();
    void shutdown();

private:
    bool build_outputs();
    void rebuild(const char* reason);
    void handle_events();
    // Decides what is drawn, and at what frame rate, from the tray, the
    // power state and the windows covering the displays.
    void update_activity();

    WallpaperOptions m_options;
    Settings m_settings;
    int m_display = 0; // 0 = all, n = the n-th display
    InstanceLock m_lock;
    SceneHost m_host;
    WallpaperLayer m_layer;
    TrayIcon m_tray;
    PowerMonitor m_power;
    std::vector<void*> m_monitors; // HMONITOR per output
    int m_fps = 0;                 // the cap on AC power
    bool m_idle = false;           // nothing is drawn
    std::string m_activity;        // the last logged state
    bool m_sdl_ready = false;
    bool m_running = false;
    bool m_paused = false;
    bool m_rebuild = false;
    uint64_t m_last_check_ns = 0;
};

} // namespace astraxis
