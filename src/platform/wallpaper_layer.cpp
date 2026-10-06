#include "platform/wallpaper_layer.hpp"

#include <SDL3/SDL_log.h>
#include <SDL3/SDL_timer.h>

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <dwmapi.h>

#include <cwchar>
#include <initializer_list>

namespace astraxis {

namespace {

constexpr const wchar_t* kHolderClass = L"AstraxisWallpaperHolder";
constexpr UINT kSpawnWorkerW = 0x052C; // undocumented; Explorer splits off the wallpaper layer
constexpr int kSpawnWaitMs = 2000;

HWND hwnd(void* p)
{
    return static_cast<HWND>(p);
}

bool register_holder_class()
{
    static const bool registered = [] {
        WNDCLASSEXW wc = {};
        wc.cbSize = sizeof(wc);
        wc.lpfnWndProc = DefWindowProcW;
        wc.hInstance = GetModuleHandleW(nullptr);
        wc.hbrBackground = static_cast<HBRUSH>(GetStockObject(BLACK_BRUSH));
        wc.lpszClassName = kHolderClass;
        return RegisterClassExW(&wc) != 0;
    }();
    return registered;
}

// Classic layout: the top-level WorkerW right behind the one holding the icons.
struct ClassicSearch {
    HWND icons = nullptr;
    HWND workerw = nullptr;
};

BOOL CALLBACK find_classic(HWND top, LPARAM param)
{
    auto* search = reinterpret_cast<ClassicSearch*>(param);
    if (HWND icons = FindWindowExW(top, nullptr, L"SHELLDLL_DefView", nullptr)) {
        search->icons = icons;
        search->workerw = FindWindowExW(nullptr, top, L"WorkerW", nullptr);
        return FALSE;
    }
    return TRUE;
}

bool find_layer(HWND progman, bool raised, HWND& icons, HWND& workerw)
{
    if (raised) {
        icons = FindWindowExW(progman, nullptr, L"SHELLDLL_DefView", nullptr);
        workerw = FindWindowExW(progman, nullptr, L"WorkerW", nullptr);
    } else {
        ClassicSearch search;
        EnumWindows(find_classic, reinterpret_cast<LPARAM>(&search));
        icons = search.icons;
        workerw = search.workerw;
    }
    return icons && workerw;
}

struct CoverSearch {
    RECT work;
    DWORD own_process;
    bool covered = false;
};

BOOL CALLBACK find_cover(HWND top, LPARAM param)
{
    auto* search = reinterpret_cast<CoverSearch*>(param);
    if (!IsWindowVisible(top) || IsIconic(top)) {
        return TRUE;
    }
    DWORD process = 0;
    GetWindowThreadProcessId(top, &process);
    if (process == search->own_process) {
        return TRUE;
    }
    // Click-through overlays and tool windows do not hide the desktop; nor do
    // windows on another virtual desktop (cloaked).
    const LONG_PTR ex = GetWindowLongPtrW(top, GWL_EXSTYLE);
    if (ex & (WS_EX_TRANSPARENT | WS_EX_TOOLWINDOW | WS_EX_NOACTIVATE)) {
        return TRUE;
    }
    DWORD cloaked = 0;
    if (SUCCEEDED(DwmGetWindowAttribute(top, DWMWA_CLOAKED, &cloaked, sizeof(cloaked))) && cloaked) {
        return TRUE;
    }
    wchar_t cls[64] = {};
    GetClassNameW(top, cls, 64);
    for (const wchar_t* shell : {L"Progman", L"WorkerW", L"Shell_TrayWnd", L"Shell_SecondaryTrayWnd"}) {
        if (wcscmp(cls, shell) == 0) {
            return TRUE;
        }
    }
    RECT r = {};
    GetWindowRect(top, &r);
    const RECT& w = search->work;
    if (r.left <= w.left && r.top <= w.top && r.right >= w.right && r.bottom >= w.bottom) {
        search->covered = true;
        return FALSE;
    }
    return TRUE;
}

} // namespace

bool WallpaperLayer::monitor_covered(void* monitor)
{
    MONITORINFO info = {};
    info.cbSize = sizeof(info);
    if (!GetMonitorInfoW(static_cast<HMONITOR>(monitor), &info)) {
        return false;
    }
    CoverSearch search;
    search.work = info.rcWork;
    search.own_process = GetCurrentProcessId();
    EnumWindows(find_cover, reinterpret_cast<LPARAM>(&search));
    return search.covered;
}

bool WallpaperLayer::monitor_size(void* monitor, int& width, int& height)
{
    MONITORINFO info = {};
    info.cbSize = sizeof(info);
    if (!GetMonitorInfoW(static_cast<HMONITOR>(monitor), &info)) {
        return false;
    }
    width = info.rcMonitor.right - info.rcMonitor.left;
    height = info.rcMonitor.bottom - info.rcMonitor.top;
    return true;
}

bool WallpaperLayer::prepare()
{
    const HWND progman = FindWindowW(L"Progman", nullptr);
    if (!progman) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "No desktop (Progman) to put the wallpaper on");
        return false;
    }
    const bool raised = (GetWindowLongPtrW(progman, GWL_EXSTYLE) & WS_EX_NOREDIRECTIONBITMAP) != 0;

    HWND icons = nullptr;
    HWND workerw = nullptr;
    if (!find_layer(progman, raised, icons, workerw)) {
        // 0xD, 1: what Lively and others send; with 0, 0 Windows 11 24H2 no longer splits.
        SendMessageTimeoutW(progman, kSpawnWorkerW, 0xD, 1, SMTO_NORMAL, 1000, nullptr);
        for (int waited = 0; waited < kSpawnWaitMs && !find_layer(progman, raised, icons, workerw); waited += 50) {
            SDL_Delay(50);
        }
    }
    if (!icons || !workerw) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Explorer did not create the wallpaper layer (%s desktop)",
                     raised ? "raised" : "classic");
        return false;
    }
    m_progman = progman;
    m_icons = icons;
    m_workerw = workerw;
    m_raised = raised;
    SDL_LogInfo(SDL_LOG_CATEGORY_APPLICATION, "Wallpaper layer: %s desktop", raised ? "raised (24H2+)" : "classic");
    return true;
}

bool WallpaperLayer::attach(void* window, void* monitor)
{
    MONITORINFO info = {};
    info.cbSize = sizeof(info);
    if (!m_progman || !GetMonitorInfoW(static_cast<HMONITOR>(monitor), &info)) {
        return false;
    }
    const RECT& r = info.rcMonitor;
    const int width = r.right - r.left;
    const int height = r.bottom - r.top;

    HWND parent = nullptr;
    POINT origin = {r.left, r.top};
    if (m_raised) {
        if (!register_holder_class()) {
            return false;
        }
        ScreenToClient(hwnd(m_progman), &origin);
        const HWND holder =
            CreateWindowExW(WS_EX_LAYERED | WS_EX_NOACTIVATE, kHolderClass, L"", WS_CHILD | WS_VISIBLE | WS_CLIPCHILDREN,
                            origin.x, origin.y, width, height, hwnd(m_progman), nullptr, GetModuleHandleW(nullptr),
                            nullptr);
        if (!holder) {
            SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Creating the wallpaper holder failed: %lu", GetLastError());
            return false;
        }
        SetLayeredWindowAttributes(holder, 0, 255, LWA_ALPHA);
        SetWindowPos(holder, hwnd(m_icons), 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE | SWP_NOACTIVATE);
        m_holders.push_back(holder);
        parent = holder;
        origin = {0, 0};
    } else {
        parent = hwnd(m_workerw);
        ScreenToClient(parent, &origin);
    }

    const HWND child = hwnd(window);
    LONG_PTR style = GetWindowLongPtrW(child, GWL_STYLE);
    style &= ~static_cast<LONG_PTR>(WS_POPUP | WS_CAPTION | WS_THICKFRAME | WS_SYSMENU | WS_MINIMIZEBOX |
                                    WS_MAXIMIZEBOX);
    style |= WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS;
    SetWindowLongPtrW(child, GWL_STYLE, style);
    if (!SetParent(child, parent)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SetParent into the wallpaper layer failed: %lu", GetLastError());
        return false;
    }
    SetWindowPos(child, nullptr, origin.x, origin.y, width, height,
                 SWP_NOZORDER | SWP_NOACTIVATE | SWP_FRAMECHANGED | SWP_SHOWWINDOW);
    m_windows.push_back(window);
    return true;
}

bool WallpaperLayer::check()
{
    if (!m_progman || !IsWindow(hwnd(m_progman)) || FindWindowW(L"Progman", nullptr) != hwnd(m_progman) ||
        !IsWindow(hwnd(m_icons)) || !IsWindow(hwnd(m_workerw))) {
        return false;
    }
    for (size_t i = 0; i < m_windows.size(); ++i) {
        const HWND parent = m_raised ? hwnd(m_holders[i]) : hwnd(m_workerw);
        if (!IsWindow(hwnd(m_windows[i])) || GetParent(hwnd(m_windows[i])) != parent) {
            return false;
        }
    }
    if (m_raised) {
        if (GetParent(hwnd(m_icons)) != hwnd(m_progman)) {
            return false;
        }
        // Each holder between the icons and Explorer's wallpaper; put back any
        // that Explorer moved.
        for (void* holder : m_holders) {
            if (!IsWindow(hwnd(holder))) {
                return false;
            }
            bool placed = false;
            for (HWND h = GetWindow(hwnd(m_icons), GW_HWNDNEXT); h && h != hwnd(m_workerw);
                 h = GetWindow(h, GW_HWNDNEXT)) {
                if (h == hwnd(holder)) {
                    placed = true;
                    break;
                }
            }
            if (!placed) {
                SetWindowPos(hwnd(holder), hwnd(m_icons), 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE | SWP_NOACTIVATE);
            }
        }
    }
    return true;
}

void WallpaperLayer::release()
{
    for (void* holder : m_holders) {
        if (IsWindow(hwnd(holder))) {
            DestroyWindow(hwnd(holder));
        }
    }
    const bool had_windows = !m_windows.empty();
    m_holders.clear();
    m_windows.clear();

    if (had_windows && !m_raised) {
        // The classic WorkerW keeps showing our last frame until the
        // wallpaper is set again; set the current one.
        wchar_t path[MAX_PATH] = {};
        SystemParametersInfoW(SPI_GETDESKWALLPAPER, MAX_PATH, path, 0);
        SystemParametersInfoW(SPI_SETDESKWALLPAPER, 0, path, 0);
    } else if (m_workerw && IsWindow(hwnd(m_workerw))) {
        RedrawWindow(hwnd(m_workerw), nullptr, nullptr, RDW_INVALIDATE | RDW_ERASE | RDW_ALLCHILDREN);
    }
    m_progman = nullptr;
    m_icons = nullptr;
    m_workerw = nullptr;
}

} // namespace astraxis
