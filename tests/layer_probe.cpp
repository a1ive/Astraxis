// astraxis_layer_probe.exe: checks WallpaperLayer on a desktop without a GPU
// (a virtual machine). It prints the desktop's window tree, puts a plain GDI
// window on every monitor into the wallpaper layer, keeps it there for a few
// seconds (checking the layer every second, as the wallpaper does), releases
// it and prints the tree again. A working layer shows a numbered color panel
// per monitor behind the desktop icons, and the old wallpaper afterwards.
//
// usage: astraxis_layer_probe [seconds]   (default 15)

#include "platform/wallpaper_layer.hpp"

#include <SDL3/SDL_log.h>

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cwchar>
#include <vector>

namespace {

constexpr const wchar_t* kProbeClass = L"AstraxisLayerProbe";

struct Monitor {
    HMONITOR handle;
    RECT rect;
};

BOOL CALLBACK add_monitor(HMONITOR monitor, HDC, LPRECT, LPARAM param)
{
    MONITORINFO info = {};
    info.cbSize = sizeof(info);
    GetMonitorInfoW(monitor, &info);
    reinterpret_cast<std::vector<Monitor>*>(param)->push_back({monitor, info.rcMonitor});
    return TRUE;
}

LRESULT CALLBACK probe_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam)
{
    if (msg == WM_PAINT) {
        PAINTSTRUCT ps;
        const HDC dc = BeginPaint(hwnd, &ps);
        RECT r;
        GetClientRect(hwnd, &r);
        // A color per monitor, and a frame to show that the whole monitor is covered.
        const int index = static_cast<int>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
        static const COLORREF kColors[] = {RGB(20, 60, 140), RGB(120, 30, 90), RGB(20, 110, 60), RGB(140, 90, 20)};
        const HBRUSH fill = CreateSolidBrush(kColors[index % 4]);
        FillRect(dc, &r, fill);
        DeleteObject(fill);
        const HBRUSH frame = CreateSolidBrush(RGB(255, 220, 0));
        for (int i = 0; i < 8; ++i) {
            RECT f = {r.left + i, r.top + i, r.right - i, r.bottom - i};
            FrameRect(dc, &f, frame);
        }
        DeleteObject(frame);
        wchar_t text[128];
        std::swprintf(text, 128, L"Astraxis layer probe: monitor %d (%ldx%ld)", index + 1, r.right, r.bottom);
        SetBkMode(dc, TRANSPARENT);
        SetTextColor(dc, RGB(255, 255, 255));
        DrawTextW(dc, text, -1, &r, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
        EndPaint(hwnd, &ps);
        return 0;
    }
    return DefWindowProcW(hwnd, msg, wparam, lparam);
}

void print_window(HWND h, int depth)
{
    wchar_t cls[64] = {};
    GetClassNameW(h, cls, 64);
    RECT r = {};
    GetWindowRect(h, &r);
    DWORD pid = 0;
    GetWindowThreadProcessId(h, &pid);
    std::wprintf(L"%*s%ls hwnd=%p pid=%lu rect=%ld,%ld-%ld,%ld visible=%d ex=0x%llX\n", depth * 2, L"", cls,
                 static_cast<void*>(h), pid, r.left, r.top, r.right, r.bottom, IsWindowVisible(h) ? 1 : 0,
                 static_cast<unsigned long long>(GetWindowLongPtrW(h, GWL_EXSTYLE)));
}

void print_children(HWND parent, int depth)
{
    for (HWND c = GetWindow(parent, GW_CHILD); c; c = GetWindow(c, GW_HWNDNEXT)) {
        print_window(c, depth);
        if (depth < 2) {
            print_children(c, depth + 1);
        }
    }
}

// Progman and the visible top-level WorkerW windows, with their children.
void print_desktop(const wchar_t* title)
{
    std::wprintf(L"--- %ls\n", title);
    for (HWND top = GetTopWindow(nullptr); top; top = GetWindow(top, GW_HWNDNEXT)) {
        wchar_t cls[64] = {};
        GetClassNameW(top, cls, 64);
        const bool progman = wcscmp(cls, L"Progman") == 0;
        if (progman || (wcscmp(cls, L"WorkerW") == 0 && IsWindowVisible(top))) {
            print_window(top, 0);
            print_children(top, 1);
        }
    }
}

void pump_messages(DWORD ms)
{
    const ULONGLONG end = GetTickCount64() + ms;
    while (GetTickCount64() < end) {
        MSG msg;
        while (PeekMessageW(&msg, nullptr, 0, 0, PM_REMOVE)) {
            TranslateMessage(&msg);
            DispatchMessageW(&msg);
        }
        Sleep(15);
    }
}

} // namespace

int main(int argc, char* argv[])
{
    const int seconds = argc > 1 ? std::max(1, std::atoi(argv[1])) : 15;
    // As SDL makes the wallpaper: per-monitor DPI aware (physical pixels).
    SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);
    SDL_SetLogPriorities(SDL_LOG_PRIORITY_INFO);

    OSVERSIONINFOW version = {};
    version.dwOSVersionInfoSize = sizeof(version);
    using RtlGetVersion = LONG(WINAPI*)(OSVERSIONINFOW*);
    if (auto get = reinterpret_cast<RtlGetVersion>(GetProcAddress(GetModuleHandleW(L"ntdll.dll"), "RtlGetVersion"))) {
        get(&version);
    }
    const HWND progman = FindWindowW(L"Progman", nullptr);
    std::wprintf(L"Windows %lu.%lu build %lu; Progman %p, WS_EX_NOREDIRECTIONBITMAP %ls\n", version.dwMajorVersion,
                 version.dwMinorVersion, version.dwBuildNumber, static_cast<void*>(progman),
                 progman && (GetWindowLongPtrW(progman, GWL_EXSTYLE) & WS_EX_NOREDIRECTIONBITMAP) ? L"yes (raised)"
                                                                                                   : L"no (classic)");
    print_desktop(L"desktop before");

    astraxis::WallpaperLayer layer;
    if (!layer.prepare()) {
        std::wprintf(L"FAILED: no wallpaper layer\n");
        print_desktop(L"desktop after asking Explorer for the layer");
        return 1;
    }
    std::wprintf(L"layer: %ls\n", layer.raised() ? L"raised" : L"classic");

    WNDCLASSEXW wc = {};
    wc.cbSize = sizeof(wc);
    wc.lpfnWndProc = probe_proc;
    wc.hInstance = GetModuleHandleW(nullptr);
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.lpszClassName = kProbeClass;
    RegisterClassExW(&wc);

    std::vector<Monitor> monitors;
    EnumDisplayMonitors(nullptr, nullptr, add_monitor, reinterpret_cast<LPARAM>(&monitors));
    std::vector<HWND> windows;
    for (size_t i = 0; i < monitors.size(); ++i) {
        const RECT& r = monitors[i].rect;
        std::wprintf(L"monitor %zu: %ld,%ld-%ld,%ld\n", i + 1, r.left, r.top, r.right, r.bottom);
        // A hidden top-level window first, as the wallpaper's SDL windows are.
        const HWND hwnd = CreateWindowExW(0, kProbeClass, L"Astraxis layer probe", WS_POPUP, 0, 0,
                                          r.right - r.left, r.bottom - r.top, nullptr, nullptr,
                                          GetModuleHandleW(nullptr), nullptr);
        SetWindowLongPtrW(hwnd, GWLP_USERDATA, static_cast<LONG_PTR>(i));
        if (!layer.attach(hwnd, monitors[i].handle)) {
            std::wprintf(L"FAILED: attaching to monitor %zu\n", i + 1);
            DestroyWindow(hwnd);
            continue;
        }
        windows.push_back(hwnd);
    }
    print_desktop(L"desktop with the probe attached");

    std::wprintf(L"keeping it for %d s (look behind the icons)...\n", seconds);
    bool intact = true;
    for (int s = 0; s < seconds; ++s) {
        pump_messages(1000);
        const bool now = layer.check();
        if (now != intact) {
            std::wprintf(L"after %d s: the layer %ls\n", s + 1, now ? L"is back" : L"is no longer intact");
            intact = now;
        }
    }

    for (HWND hwnd : windows) {
        DestroyWindow(hwnd);
    }
    layer.release();
    pump_messages(500);
    print_desktop(L"desktop after release");
    std::wprintf(L"%ls\n", windows.size() == monitors.size() ? L"done" : L"done, with failures");
    return windows.size() == monitors.size() ? 0 : 1;
}
