#pragma once

#include <cstdint>
#include <string>

struct SDL_Window;

namespace astraxis {

// A display as settings name it: `index` 0 = primary, n = the n-th display
// (1-based); `name` picks the display again if the order has changed.
// Falls back to the primary display. Returns an SDL_DisplayID.
uint32_t find_display(int index, const std::string& name);

// The inverse of find_display: 1-based index and name of a display.
int display_index(uint32_t display);
std::string display_name(uint32_t display);

// Thin owner of the main SDL window.
class Window {
public:
    // Centered on `display`; `fullscreen` is borderless desktop fullscreen.
    bool create(const char* title, int width, int height, uint32_t display, bool fullscreen);
    void destroy();

    SDL_Window* handle() const { return m_window; }
    float content_scale() const;
    bool is_minimized() const;

    bool is_fullscreen() const;
    void set_fullscreen(bool fullscreen);
    uint32_t display() const; // the display the window is on (SDL_DisplayID)

private:
    SDL_Window* m_window = nullptr;
};

} // namespace astraxis
