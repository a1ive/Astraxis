#pragma once

namespace astraxis {

// A notification area icon with a context menu (Settings, Pause / Resume, Exit).
// Windows only. Its messages are dispatched by the thread's message loop
// (SDL_PumpEvents).
class TrayIcon {
public:
    enum class Command { None, Settings, TogglePause, Exit };

    bool create(const char* tooltip);
    void destroy();

    void set_paused(bool paused) { m_paused = paused; }
    // The command chosen since the last call.
    Command take_command();

    // For the window procedure: a message to the icon's hidden window; true
    // if handled (WM_CLOSE asks to exit instead of destroying the window).
    bool handle_message(unsigned msg, long long lparam);

private:
    bool add_icon();
    void show_menu();

    void* m_hwnd = nullptr;
    unsigned m_taskbar_created = 0; // registered message: Explorer restarted
    char m_tooltip[128] = {};
    bool m_paused = false;
    Command m_command = Command::None;
};

// Holds a named lock while alive; false from acquire() if another process holds it.
class InstanceLock {
public:
    bool acquire(const char* name);
    ~InstanceLock();

private:
    void* m_handle = nullptr;
};

} // namespace astraxis
