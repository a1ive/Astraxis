// astraxis_wallpaper.exe: the settings dialog of the live wallpaper. Apply
// saves [wallpaper] in config.toml and (re)starts `astraxis --mode wallpaper`
// next to it; the wallpaper's tray menu opens this program.

#include "app/settings_dialog.hpp"
#include "platform/paths.hpp"

#include <SDL3/SDL.h>
#define SDL_MAIN_HANDLED // our own wWinMain
#include <SDL3/SDL_main.h>

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

int WINAPI wWinMain(HINSTANCE, HINSTANCE, PWSTR, int)
{
    // Only one dialog: bring an open one to the front.
    if (const HWND open = FindWindowW(L"#32770", L"Astraxis wallpaper")) {
        SetForegroundWindow(open);
        return 0;
    }

    SDL_SetMainReady(); // SDL is used for the display list, without SDL_main
    astraxis::SettingsDialogContext context;
    context.kind = astraxis::SettingsDialogKind::Wallpaper;
    context.config = astraxis::config_path();
    context.asset_dir = astraxis::asset_directory();
    context.program = context.config.parent_path() / "astraxis.exe";
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        MessageBoxA(nullptr, SDL_GetError(), "Astraxis", MB_ICONERROR);
        return 1;
    }
    astraxis::run_settings_dialog(context, nullptr);
    SDL_Quit();
    return 0;
}
