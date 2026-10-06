#include "platform/power_monitor.hpp"

#include <SDL3/SDL_log.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <wtsapi32.h>
#endif

namespace astraxis {

#ifdef _WIN32

namespace {

constexpr const wchar_t* kPowerClass = L"AstraxisPower";

// GUID_CONSOLE_DISPLAY_STATE (winnt.h), spelled out to avoid initguid.h.
constexpr GUID kConsoleDisplayState = {0x6fe69556, 0x704a, 0x47a0, {0x8f, 0x24, 0xc2, 0x8d, 0x93, 0x6f, 0xda, 0x47}};

LRESULT CALLBACK power_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam)
{
    if (auto* monitor = reinterpret_cast<PowerMonitor*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA))) {
        monitor->handle_message(msg, wparam, lparam);
    }
    return DefWindowProcW(hwnd, msg, wparam, lparam);
}

} // namespace

bool PowerMonitor::start()
{
    WNDCLASSEXW wc = {};
    wc.cbSize = sizeof(wc);
    wc.lpfnWndProc = power_proc;
    wc.hInstance = GetModuleHandleW(nullptr);
    wc.lpszClassName = kPowerClass;
    RegisterClassExW(&wc); // fails harmlessly if already registered

    const HWND hwnd = CreateWindowExW(0, kPowerClass, L"", 0, 0, 0, 0, 0, HWND_MESSAGE, nullptr,
                                      GetModuleHandleW(nullptr), nullptr);
    if (!hwnd) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "No power notifications: %lu", GetLastError());
        return false;
    }
    SetWindowLongPtrW(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    m_hwnd = hwnd;
    WTSRegisterSessionNotification(hwnd, NOTIFY_FOR_THIS_SESSION);
    // Windows sends the current display state right away.
    m_display_notify = RegisterPowerSettingNotification(hwnd, &kConsoleDisplayState, DEVICE_NOTIFY_WINDOW_HANDLE);
    poll();
    return true;
}

void PowerMonitor::stop()
{
    if (m_display_notify) {
        UnregisterPowerSettingNotification(static_cast<HPOWERNOTIFY>(m_display_notify));
        m_display_notify = nullptr;
    }
    if (m_hwnd) {
        WTSUnRegisterSessionNotification(static_cast<HWND>(m_hwnd));
        DestroyWindow(static_cast<HWND>(m_hwnd));
        m_hwnd = nullptr;
    }
}

void PowerMonitor::poll()
{
    SYSTEM_POWER_STATUS status = {};
    if (GetSystemPowerStatus(&status)) {
        // ACLineStatus 0: on battery; SystemStatusFlag bit 0: the battery saver is on.
        m_on_battery = status.ACLineStatus == 0 || (status.SystemStatusFlag & 1) != 0;
    }
}

void PowerMonitor::handle_message(unsigned msg, unsigned long long wparam, long long lparam)
{
    if (msg == WM_WTSSESSION_CHANGE) {
        if (wparam == WTS_SESSION_LOCK) {
            m_locked = true;
        } else if (wparam == WTS_SESSION_UNLOCK) {
            m_locked = false;
        }
    } else if (msg == WM_POWERBROADCAST && wparam == PBT_POWERSETTINGCHANGE) {
        const auto* setting = reinterpret_cast<const POWERBROADCAST_SETTING*>(lparam);
        if (setting && IsEqualGUID(setting->PowerSetting, kConsoleDisplayState) && setting->DataLength >= 1) {
            m_display_off = setting->Data[0] == 0; // 0 off, 1 on, 2 dimmed
        }
    }
}

#else

bool PowerMonitor::start()
{
    return true;
}

void PowerMonitor::stop() {}

void PowerMonitor::poll() {}

void PowerMonitor::handle_message(unsigned, unsigned long long, long long) {}

#endif

} // namespace astraxis
