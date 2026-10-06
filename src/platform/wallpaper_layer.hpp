#pragma once

#include <vector>

namespace astraxis {

// Where a live wallpaper's windows go: behind the desktop icons, in front of
// Explorer's own wallpaper. Windows only.
//
// Explorer builds the desktop in one of two ways:
// - Windows 11 24H2 and later ("raised desktop", Progman has
//   WS_EX_NOREDIRECTIONBITMAP): Progman holds the icons (SHELLDLL_DefView)
//   and, once asked with 0x052C, a WorkerW below them that paints the
//   wallpaper. Each of our windows goes into a layered child of Progman
//   placed between the two: a layered child has a surface of its own, which
//   Progman, without a redirection bitmap, cannot give its children.
// - Before: 0x052C splits off a top-level WorkerW behind the one that holds
//   the icons; our windows become its children.
class WallpaperLayer {
public:
    // Size of a monitor (HMONITOR) in physical pixels.
    static bool monitor_size(void* monitor, int& width, int& height);
    // Another program's visible window covers the monitor's work area (a
    // maximized or full-screen application): the wallpaper is out of sight.
    static bool monitor_covered(void* monitor);

    // Finds the layer, asking Explorer to create it if needed (waits up to 2 s).
    bool prepare();
    // Reparents a top-level window (HWND) into the layer, covering `monitor`
    // (HMONITOR). Coordinates are physical pixels (the process is per-monitor
    // DPI aware).
    bool attach(void* window, void* monitor);
    // Keeps our windows below the icons; false once the layer is gone or our
    // windows are no longer in it (Explorer restarted or rebuilt the desktop):
    // then release() and attach everything again.
    bool check();
    // After our windows are destroyed: removes the holders and repaints
    // Explorer's wallpaper.
    void release();

    bool raised() const { return m_raised; }

private:
    void* m_progman = nullptr;
    void* m_icons = nullptr;   // SHELLDLL_DefView
    void* m_workerw = nullptr; // raised: Explorer's wallpaper; classic: our parent
    bool m_raised = false;
    std::vector<void*> m_holders; // raised: one layered child of Progman per window
    std::vector<void*> m_windows;
};

} // namespace astraxis
