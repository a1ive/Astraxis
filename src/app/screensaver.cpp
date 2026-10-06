#include "app/screensaver.hpp"

#include "platform/paths.hpp"

#include <SDL3/SDL.h>

#include <cmath>

namespace astraxis {

namespace {

constexpr uint64_t kGraceNs = 300'000'000; // input right after the start is ignored
constexpr float kMouseSlack = 8.0f;        // pixels the mouse may drift before it counts
constexpr int kPreviewFps = 30;            // the preview is tiny; spare the GPU

} // namespace

bool Screensaver::init(const ScreensaverOptions& options)
{
    m_options = options;

    // Let the display power down while the screensaver runs (SDL otherwise
    // keeps it on), and keep every display's window up when another one
    // takes the focus.
    SDL_SetHint(SDL_HINT_VIDEO_ALLOW_SCREENSAVER, "1");
    SDL_SetHint(SDL_HINT_VIDEO_MINIMIZE_ON_FOCUS_LOSS, "0");
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_Init failed: %s", SDL_GetError());
        return false;
    }
    m_sdl_ready = true;

    const std::filesystem::path config = options.config.empty() ? config_path() : options.config;
    Settings settings;
    std::string error;
    if (!load_settings(config, SettingsSection::Screensaver, settings, &error)) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Settings not read: %s", error.c_str());
    }
    m_host.set_fps_limit(options.preview ? kPreviewFps : options.fps.value_or(settings.fps));

    if (!m_host.init_gpu()) {
        return false;
    }
    int display = options.display.value_or(0);
    if (!options.display && settings.display > 0) {
        display = display_index(find_display(settings.display, settings.display_name));
    }
    if (!create_outputs(display)) {
        return false;
    }
    if (!m_host.start(settings, options.scene, settings.labels && !options.preview, &error)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "%s", error.c_str());
        if (!options.preview) {
            SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Astraxis", error.c_str(), nullptr);
        }
        return false;
    }

    if (!options.preview) {
        SDL_HideCursor();
    }
    m_start_ns = SDL_GetTicksNS();
    m_running = true;
    return true;
}

bool Screensaver::create_outputs(int display)
{
    if (m_options.preview) {
        // No swapchain: the frames are read back and drawn with GDI.
        auto out = std::make_unique<SceneHost::Output>();
        out->gdi = true;
        return out->window.create_child(m_options.parent) && m_host.add_output(std::move(out));
    }

    // One window per display; `display` (1-based) alone shows the scene, or all if 0.
    int count = 0;
    SDL_DisplayID* displays = SDL_GetDisplays(&count);
    for (int i = 0; displays && i < count; ++i) {
        auto out = std::make_unique<SceneHost::Output>();
        out->scene = display <= 0 || display > count || display == i + 1;
        if (out->window.create_cover("Astraxis", displays[i])) {
            m_host.add_output(std::move(out));
        }
    }
    SDL_free(displays);
    if (m_host.outputs().empty()) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "No screensaver window could be created");
        return false;
    }
    // The first window takes the keyboard.
    SDL_RaiseWindow(m_host.outputs().front()->window.handle());
    return true;
}

void Screensaver::run()
{
    while (m_running) {
        handle_events();
        if (m_options.preview && !m_host.outputs().front()->window.parent_alive()) {
            quit("preview window gone");
        }
        if (m_running) {
            m_host.frame();
        }
    }
}

void Screensaver::shutdown()
{
    m_host.shutdown();
    if (m_sdl_ready) {
        SDL_ShowCursor();
        SDL_Quit();
        m_sdl_ready = false;
    }
}

void Screensaver::handle_events()
{
    // Input that arrives with the windows (focus changes, the click that
    // started a preview) does not count.
    const bool armed = !m_options.preview && SDL_GetTicksNS() - m_start_ns > kGraceNs;

    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
        case SDL_EVENT_QUIT:
        case SDL_EVENT_WINDOW_CLOSE_REQUESTED:
        case SDL_EVENT_WINDOW_DESTROYED:
            quit("window closed");
            break;
        case SDL_EVENT_KEY_DOWN:
            if (armed) {
                quit("key pressed");
            }
            break;
        case SDL_EVENT_MOUSE_BUTTON_DOWN:
        case SDL_EVENT_MOUSE_WHEEL:
            if (armed) {
                quit("mouse button or wheel");
            }
            break;
        case SDL_EVENT_MOUSE_MOTION: {
            // In desktop coordinates: motion events are relative to whichever
            // window the mouse is over.
            float x = 0.0f;
            float y = 0.0f;
            SDL_GetGlobalMouseState(&x, &y);
            if (!m_have_mouse) {
                m_have_mouse = true;
                m_mouse_x = x;
                m_mouse_y = y;
            } else if (armed && std::fabs(x - m_mouse_x) + std::fabs(y - m_mouse_y) > kMouseSlack) {
                quit("mouse moved");
            }
            break;
        }
        default:
            break;
        }
    }
}

void Screensaver::quit(const char* reason)
{
    if (m_running) {
        SDL_LogInfo(SDL_LOG_CATEGORY_APPLICATION, "Screensaver ends: %s", reason);
    }
    m_running = false;
}

} // namespace astraxis
