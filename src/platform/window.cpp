#include "platform/window.hpp"

#include <SDL3/SDL.h>

namespace astraxis {

uint32_t find_display(int index, const std::string& name)
{
    int count = 0;
    SDL_DisplayID* displays = SDL_GetDisplays(&count);
    SDL_DisplayID found = 0;
    auto named = [&](int i) {
        const char* n = SDL_GetDisplayName(displays[i]);
        return n && name == n;
    };
    if (displays && index >= 1 && index <= count && (name.empty() || named(index - 1))) {
        found = displays[index - 1];
    }
    for (int i = 0; displays && !found && !name.empty() && i < count; ++i) {
        if (named(i)) {
            found = displays[i];
        }
    }
    if (displays && !found && index >= 1 && index <= count) {
        found = displays[index - 1];
    }
    SDL_free(displays);
    return found ? found : SDL_GetPrimaryDisplay();
}

int display_index(uint32_t display)
{
    int count = 0;
    SDL_DisplayID* displays = SDL_GetDisplays(&count);
    int index = 0;
    for (int i = 0; displays && i < count; ++i) {
        if (displays[i] == display) {
            index = i + 1;
        }
    }
    SDL_free(displays);
    return index;
}

std::string display_name(uint32_t display)
{
    const char* name = SDL_GetDisplayName(display);
    return name ? name : "";
}

bool Window::create(const char* title, int width, int height, uint32_t display, bool fullscreen)
{
    const float scale = SDL_GetDisplayContentScale(display);
    const SDL_WindowFlags flags = SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIDDEN | SDL_WINDOW_HIGH_PIXEL_DENSITY;

    m_window = SDL_CreateWindow(title, static_cast<int>(width * scale), static_cast<int>(height * scale), flags);
    if (!m_window) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_CreateWindow failed: %s", SDL_GetError());
        return false;
    }

    SDL_SetWindowPosition(m_window, SDL_WINDOWPOS_CENTERED_DISPLAY(display), SDL_WINDOWPOS_CENTERED_DISPLAY(display));
    if (fullscreen) {
        set_fullscreen(true);
    }
    SDL_ShowWindow(m_window);
    return true;
}

void Window::destroy()
{
    if (m_window) {
        SDL_DestroyWindow(m_window);
        m_window = nullptr;
    }
}

float Window::content_scale() const
{
    const float scale = SDL_GetWindowDisplayScale(m_window);
    return scale > 0.0f ? scale : 1.0f;
}

bool Window::is_minimized() const
{
    return (SDL_GetWindowFlags(m_window) & SDL_WINDOW_MINIMIZED) != 0;
}

bool Window::is_fullscreen() const
{
    return (SDL_GetWindowFlags(m_window) & SDL_WINDOW_FULLSCREEN) != 0;
}

void Window::set_fullscreen(bool fullscreen)
{
    // No fullscreen mode set: borderless fullscreen at the desktop resolution.
    if (!SDL_SetWindowFullscreen(m_window, fullscreen)) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "SDL_SetWindowFullscreen failed: %s", SDL_GetError());
    }
}

uint32_t Window::display() const
{
    const SDL_DisplayID display = SDL_GetDisplayForWindow(m_window);
    return display ? display : SDL_GetPrimaryDisplay();
}

} // namespace astraxis
