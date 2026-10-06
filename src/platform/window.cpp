#include "platform/window.hpp"

#include <SDL3/SDL.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif

namespace astraxis {

#ifdef _WIN32
namespace {

constexpr const wchar_t* kChildClass = L"AstraxisChild";

bool register_child_class()
{
    static const bool registered = [] {
        WNDCLASSEXW wc = {};
        wc.cbSize = sizeof(wc);
        wc.lpfnWndProc = DefWindowProcW; // SDL_PumpEvents dispatches its messages
        wc.hInstance = GetModuleHandleW(nullptr);
        wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
        wc.lpszClassName = kChildClass;
        return RegisterClassExW(&wc) != 0;
    }();
    return registered;
}

} // namespace
#endif

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

void* native_display(uint32_t display)
{
#ifdef _WIN32
    return SDL_GetPointerProperty(SDL_GetDisplayProperties(display), SDL_PROP_DISPLAY_WINDOWS_HMONITOR_POINTER, nullptr);
#else
    (void)display;
    return nullptr;
#endif
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

bool Window::create_cover(const char* title, uint32_t display)
{
    SDL_Rect bounds = {};
    if (!SDL_GetDisplayBounds(display, &bounds)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_GetDisplayBounds failed: %s", SDL_GetError());
        return false;
    }
    const SDL_WindowFlags flags =
        SDL_WINDOW_HIDDEN | SDL_WINDOW_BORDERLESS | SDL_WINDOW_ALWAYS_ON_TOP | SDL_WINDOW_HIGH_PIXEL_DENSITY;
    m_window = SDL_CreateWindow(title, bounds.w, bounds.h, flags);
    if (!m_window) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_CreateWindow failed: %s", SDL_GetError());
        return false;
    }
    SDL_SetWindowPosition(m_window, bounds.x, bounds.y);
    set_fullscreen(true);
    SDL_ShowWindow(m_window);
    return true;
}

bool Window::create_hidden(const char* title, int width, int height)
{
    const SDL_WindowFlags flags = SDL_WINDOW_HIDDEN | SDL_WINDOW_BORDERLESS | SDL_WINDOW_HIGH_PIXEL_DENSITY;
    m_window = SDL_CreateWindow(title, width, height, flags);
    if (!m_window) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_CreateWindow failed: %s", SDL_GetError());
        return false;
    }
    return true;
}

void* Window::native_handle() const
{
#ifdef _WIN32
    if (m_child) {
        return m_child;
    }
    return m_window ? SDL_GetPointerProperty(SDL_GetWindowProperties(m_window), SDL_PROP_WINDOW_WIN32_HWND_POINTER,
                                             nullptr)
                    : nullptr;
#else
    return nullptr;
#endif
}

bool Window::create_child(void* parent)
{
#ifdef _WIN32
    const HWND parent_hwnd = static_cast<HWND>(parent);
    RECT rect = {};
    if (!IsWindow(parent_hwnd) || !GetClientRect(parent_hwnd, &rect)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Parent window %p does not exist", parent);
        return false;
    }
    if (!register_child_class()) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "RegisterClassExW failed: %lu", GetLastError());
        return false;
    }
    const HWND child = CreateWindowExW(0, kChildClass, L"Astraxis", WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS, 0, 0,
                                       rect.right, rect.bottom, parent_hwnd, nullptr, GetModuleHandleW(nullptr),
                                       nullptr);
    if (!child) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "CreateWindowExW failed: %lu", GetLastError());
        return false;
    }

    m_parent = parent;
    m_child = child;
    return true;
#else
    (void)parent;
    SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Child windows of foreign windows are only supported on Windows");
    return false;
#endif
}

void Window::destroy()
{
    if (m_window) {
        SDL_DestroyWindow(m_window);
        m_window = nullptr;
    }
#ifdef _WIN32
    if (m_child && IsWindow(static_cast<HWND>(m_child))) {
        DestroyWindow(static_cast<HWND>(m_child));
    }
#endif
    m_child = nullptr;
    m_parent = nullptr;
}

void Window::child_size(int& width, int& height)
{
    width = 0;
    height = 0;
#ifdef _WIN32
    RECT parent = {};
    RECT child = {};
    if (!m_child || !GetClientRect(static_cast<HWND>(m_parent), &parent) ||
        !GetClientRect(static_cast<HWND>(m_child), &child)) {
        return;
    }
    if (child.right != parent.right || child.bottom != parent.bottom) {
        MoveWindow(static_cast<HWND>(m_child), 0, 0, parent.right, parent.bottom, FALSE);
    }
    width = parent.right;
    height = parent.bottom;
#endif
}

void Window::draw_child(const uint8_t* pixels, int width, int height)
{
#ifdef _WIN32
    if (!m_child || !pixels) {
        return;
    }
    const HWND hwnd = static_cast<HWND>(m_child);
    BITMAPINFO info = {};
    info.bmiHeader.biSize = sizeof(info.bmiHeader);
    info.bmiHeader.biWidth = width;
    info.bmiHeader.biHeight = -height; // top-down rows
    info.bmiHeader.biPlanes = 1;
    info.bmiHeader.biBitCount = 32;
    info.bmiHeader.biCompression = BI_RGB;
    if (const HDC dc = GetDC(hwnd)) {
        SetDIBitsToDevice(dc, 0, 0, static_cast<DWORD>(width), static_cast<DWORD>(height), 0, 0, 0,
                          static_cast<UINT>(height), pixels, &info, DIB_RGB_COLORS);
        ReleaseDC(hwnd, dc);
    }
#else
    (void)pixels;
    (void)width;
    (void)height;
#endif
}

bool Window::parent_alive() const
{
#ifdef _WIN32
    return !m_parent || IsWindow(static_cast<HWND>(m_parent));
#else
    return true;
#endif
}

float Window::content_scale() const
{
    if (m_scale > 0.0f) {
        return m_scale;
    }
    if (!m_window) {
        return 1.0f;
    }
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
