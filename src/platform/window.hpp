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

// Thin owner of an SDL window.
class Window {
public:
    // Centered on `display`; `fullscreen` is borderless desktop fullscreen.
    bool create(const char* title, int width, int height, uint32_t display, bool fullscreen);
    // Covers `display`, above all other windows (screensaver).
    bool create_cover(const char* title, uint32_t display);
    // A native child window filling the client area of a foreign native
    // window (Win32 HWND; the screensaver preview), drawn with draw_child.
    // Not an SDL window: SDL would move and resize it by its own DPI logic,
    // meant for top-level windows. Windows only.
    bool create_child(void* parent);
    void destroy();

    SDL_Window* handle() const { return m_window; }
    float content_scale() const; // 1 for a child window
    bool is_minimized() const;

    bool is_fullscreen() const;
    void set_fullscreen(bool fullscreen);
    uint32_t display() const; // the display the window is on (SDL_DisplayID)

    // The foreign parent of a child window still exists.
    bool parent_alive() const;
    // Size of a child window, kept filling its parent's client area.
    void child_size(int& width, int& height);
    // Draws pixels (B, G, R, x; rows top to bottom) into a child window with
    // GDI: a swapchain does not show in a child of another process's window.
    void draw_child(const uint8_t* pixels, int width, int height);

private:
    SDL_Window* m_window = nullptr;
    void* m_parent = nullptr; // create_child: the foreign parent
    void* m_child = nullptr;  // create_child: our native child window (no SDL window)
};

} // namespace astraxis
