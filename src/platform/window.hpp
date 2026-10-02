#pragma once

struct SDL_Window;

namespace astraxis {

// Thin owner of the main SDL window.
class Window {
public:
    bool create(const char* title, int width, int height);
    void destroy();

    SDL_Window* handle() const { return m_window; }
    float content_scale() const;
    bool is_minimized() const;

private:
    SDL_Window* m_window = nullptr;
};

} // namespace astraxis
