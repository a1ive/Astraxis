#pragma once

namespace astraxis {

// The states that make drawing pointless or costly: the session is locked,
// the display is off, the computer runs on battery (or the battery saver is
// on). Windows; elsewhere it reports an unlocked session, a display that is
// on and AC power.
class PowerMonitor {
public:
    // Creates the hidden window that receives the lock and display
    // notifications; they are dispatched by the thread's message loop
    // (SDL_PumpEvents).
    bool start();
    void stop();
    // Reads the battery state (call every second or so).
    void poll();

    bool locked() const { return m_locked; }
    bool display_off() const { return m_display_off; }
    bool on_battery() const { return m_on_battery; }

    // For the window procedure.
    void handle_message(unsigned msg, unsigned long long wparam, long long lparam);

private:
#ifdef _WIN32
    void* m_hwnd = nullptr;
    void* m_display_notify = nullptr;
#endif
    bool m_locked = false;
    bool m_display_off = false;
    bool m_on_battery = false;
};

} // namespace astraxis
