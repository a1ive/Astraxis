#pragma once

#include "view/view_options.hpp"

#include <filesystem>
#include <string>

namespace astraxis {

// The section of config.toml that a host reads on top of the shared [view].
enum class SettingsSection { Window, Wallpaper, Screensaver };

// Settings kept in config.toml next to the executable. [view] holds the
// display options shared by all modes; each mode's section ([window],
// [wallpaper], [screensaver]) holds that mode's own keys and may override
// any [view] key.
struct Settings {
    // Mode keys.
    std::string scene = "solar_system";
    int display = 0;          // 0 = primary, n = the n-th display (1-based)
    std::string display_name; // finds the display again if the order changes
    int fps = 0;              // frame rate cap, 0 = none (the display's refresh rate)
    bool fullscreen = false;  // [window] only
    bool auto_tour = true;    // [window] only: "Start when idle"
    bool info = true;         // [window] only: the focused body's info panel

    // View keys.
    bool labels = true;
    ViewOptions view;
};

// Reads the settings of `section` into `out` (defaults where the file has no
// value; out-of-range values are clamped). A missing file is not an error.
// False if the file cannot be read or parsed; `out` then keeps the defaults.
bool load_settings(const std::filesystem::path& path, SettingsSection section, Settings& out,
                   std::string* error = nullptr);

// Writes the view keys to [view] (dropping the section's overrides of them)
// and the mode keys to `section`, keeping the file's other sections and keys.
// Refuses to overwrite a file that cannot be parsed.
bool save_settings(const std::filesystem::path& path, SettingsSection section, const Settings& settings,
                   std::string* error = nullptr);

} // namespace astraxis
