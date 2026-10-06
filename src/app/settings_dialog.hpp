#pragma once

#include <filesystem>

namespace astraxis {

enum class SettingsDialogKind {
    Wallpaper,   // the launcher: Apply (saves, starts or restarts the wallpaper), Stop, Close
    Screensaver, // /c: OK (saves), Cancel
};

struct SettingsDialogContext {
    SettingsDialogKind kind = SettingsDialogKind::Wallpaper;
    std::filesystem::path config;    // config.toml
    std::filesystem::path asset_dir; // for the scene list
    std::filesystem::path program;   // astraxis.exe, which runs the wallpaper
};

// The settings dialog for the [wallpaper] or [screensaver] section (Windows,
// modal to `parent`, an HWND or null). The view options are saved as
// overrides of [view] where they differ from it. Initializes SDL's video
// subsystem for the display list if it is not up.
void run_settings_dialog(const SettingsDialogContext& context, void* parent);

} // namespace astraxis
