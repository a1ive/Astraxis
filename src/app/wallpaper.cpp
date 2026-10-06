#include "app/wallpaper.hpp"

#include "platform/paths.hpp"

#include <SDL3/SDL.h>

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <shellapi.h>

namespace astraxis {

namespace {

constexpr uint64_t kCheckPeriodNs = 1'000'000'000; // how often the layer is checked
constexpr int kPausedWaitMs = 200;

} // namespace

bool Wallpaper::init(const WallpaperOptions& options)
{
    m_options = options;

    // The wallpaper must not keep the display on, and its windows disappear
    // when Explorer restarts: that is not a reason to quit.
    SDL_SetHint(SDL_HINT_VIDEO_ALLOW_SCREENSAVER, "1");
    SDL_SetHint(SDL_HINT_VIDEO_MINIMIZE_ON_FOCUS_LOSS, "0");
    SDL_SetHint(SDL_HINT_QUIT_ON_LAST_WINDOW_CLOSE, "0");
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_Init failed: %s", SDL_GetError());
        return false;
    }
    m_sdl_ready = true;

    if (!m_lock.acquire("Local\\AstraxisWallpaper")) {
        SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_INFORMATION, "Astraxis",
                                 "The Astraxis wallpaper is already running (see the notification area).", nullptr);
        return false;
    }

    const std::filesystem::path config = options.config.empty() ? config_path() : options.config;
    std::string error;
    if (!load_settings(config, SettingsSection::Wallpaper, m_settings, &error)) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Settings not read: %s", error.c_str());
    }
    m_host.set_fps_limit(options.fps.value_or(m_settings.fps));
    m_display = options.display.value_or(0);
    if (!options.display && m_settings.display > 0) {
        m_display = display_index(find_display(m_settings.display, m_settings.display_name));
    }

    if (!m_host.init_gpu() || !build_outputs()) {
        SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Astraxis",
                                 "The wallpaper could not be placed on the desktop.", nullptr);
        return false;
    }
    if (!m_host.start(m_settings, options.scene, m_settings.labels, &error)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "%s", error.c_str());
        SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Astraxis", error.c_str(), nullptr);
        return false;
    }
    m_tray.create("Astraxis wallpaper");

    m_last_check_ns = SDL_GetTicksNS();
    m_running = true;
    return true;
}

bool Wallpaper::build_outputs()
{
    if (!m_layer.prepare()) {
        return false;
    }
    int count = 0;
    SDL_DisplayID* displays = SDL_GetDisplays(&count);
    for (int i = 0; displays && i < count; ++i) {
        if (m_display > 0 && m_display <= count && m_display != i + 1) {
            continue; // Explorer's wallpaper stays on the other displays
        }
        void* monitor = native_display(displays[i]);
        int width = 0;
        int height = 0;
        auto out = std::make_unique<SceneHost::Output>();
        if (!WallpaperLayer::monitor_size(monitor, width, height) ||
            !out->window.create_hidden("Astraxis wallpaper", width, height)) {
            continue;
        }
        out->window.set_content_scale(SDL_GetDisplayContentScale(displays[i]));
        if (!m_layer.attach(out->window.native_handle(), monitor)) {
            out->window.destroy();
            continue;
        }
        int sdl_width = 0;
        int sdl_height = 0;
        SDL_GetWindowSizeInPixels(out->window.handle(), &sdl_width, &sdl_height);
        SDL_LogInfo(SDL_LOG_CATEGORY_APPLICATION, "Display %d: %dx%d (SDL: %dx%d), scale %.2f", i + 1, width, height,
                    sdl_width, sdl_height, static_cast<double>(out->window.content_scale()));
        if (sdl_width != width || sdl_height != height) {
            SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "SDL has another size for the window than the monitor");
        }
        m_host.add_output(std::move(out));
    }
    SDL_free(displays);
    if (m_host.outputs().empty()) {
        m_layer.release();
        return false;
    }
    return true;
}

void Wallpaper::rebuild(const char* reason)
{
    SDL_LogInfo(SDL_LOG_CATEGORY_APPLICATION, "Rebuilding the wallpaper: %s", reason);
    m_host.remove_outputs();
    m_layer.release();
    if (!build_outputs()) {
        // Explorer may still be starting: try again at the next check.
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "The wallpaper could not be placed yet");
        m_rebuild = true;
        return;
    }
    m_rebuild = false;
    m_host.reset_clock();
}

void Wallpaper::run()
{
    while (m_running) {
        handle_events();

        switch (m_tray.take_command()) {
        case TrayIcon::Command::Settings: {
            // The launcher, next to us: its Apply restarts this process.
            const std::filesystem::path launcher = config_path().parent_path() / "astraxis_wallpaper.exe";
            if (reinterpret_cast<INT_PTR>(ShellExecuteW(nullptr, L"open", launcher.c_str(), nullptr, nullptr,
                                                        SW_SHOWNORMAL)) <= 32) {
                SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Astraxis",
                                         "astraxis_wallpaper.exe was not found next to astraxis.exe.", nullptr);
            }
            break;
        }
        case TrayIcon::Command::TogglePause:
            m_paused = !m_paused;
            m_tray.set_paused(m_paused);
            m_host.reset_clock();
            break;
        case TrayIcon::Command::Exit:
            SDL_LogInfo(SDL_LOG_CATEGORY_APPLICATION, "Wallpaper ends: exit from the tray");
            m_running = false;
            break;
        case TrayIcon::Command::None:
            break;
        }
        if (!m_running) {
            break;
        }

        const uint64_t now = SDL_GetTicksNS();
        if (now - m_last_check_ns >= kCheckPeriodNs) {
            m_last_check_ns = now;
            if (m_rebuild || !m_layer.check()) {
                rebuild(m_rebuild ? "the displays changed (or retrying)" : "the desktop changed");
            }
        }

        if (m_paused || m_host.outputs().empty()) {
            SDL_WaitEventTimeout(nullptr, kPausedWaitMs);
            m_host.reset_clock();
        } else {
            m_host.frame();
        }
    }
}

void Wallpaper::shutdown()
{
    m_host.shutdown();
    m_layer.release();
    m_tray.destroy();
    if (m_sdl_ready) {
        SDL_Quit();
        m_sdl_ready = false;
    }
}

void Wallpaper::handle_events()
{
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
        case SDL_EVENT_QUIT:
            SDL_LogInfo(SDL_LOG_CATEGORY_APPLICATION, "Wallpaper ends: quit requested");
            m_running = false;
            break;
        case SDL_EVENT_DISPLAY_ADDED:
        case SDL_EVENT_DISPLAY_REMOVED:
        case SDL_EVENT_DISPLAY_MOVED:
        case SDL_EVENT_DISPLAY_CURRENT_MODE_CHANGED:
        case SDL_EVENT_DISPLAY_CONTENT_SCALE_CHANGED:
            m_rebuild = true; // at the next check, once the changes have settled
            break;
        default:
            break;
        }
    }
}

} // namespace astraxis
