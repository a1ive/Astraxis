#include "platform/window.hpp"

#include <SDL3/SDL.h>

namespace astraxis {

bool Window::create(const char* title, int width, int height)
{
    const float scale = SDL_GetDisplayContentScale(SDL_GetPrimaryDisplay());
    const SDL_WindowFlags flags = SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIDDEN | SDL_WINDOW_HIGH_PIXEL_DENSITY;

    m_window = SDL_CreateWindow(title, static_cast<int>(width * scale), static_cast<int>(height * scale), flags);
    if (!m_window) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_CreateWindow failed: %s", SDL_GetError());
        return false;
    }

    SDL_SetWindowPosition(m_window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
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

} // namespace astraxis
