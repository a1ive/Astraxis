#include "platform/tray_icon.hpp"

#include <SDL3/SDL_log.h>

#include <cstdio>
#include <string_view>

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <shellapi.h>

namespace astraxis {

namespace {

constexpr const wchar_t* kTrayClass = L"AstraxisTray";
constexpr UINT kCallbackMessage = WM_APP + 1;
constexpr UINT kIconId = 1;
constexpr UINT kMenuPause = 1;
constexpr UINT kMenuExit = 2;
constexpr UINT kMenuSettings = 3;

LRESULT CALLBACK tray_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam)
{
    auto* tray = reinterpret_cast<TrayIcon*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (tray && tray->handle_message(msg, lparam)) {
        return 0;
    }
    return DefWindowProcW(hwnd, msg, wparam, lparam);
}

std::wstring_view widen(const char* utf8, wchar_t* buffer, int size)
{
    const int n = MultiByteToWideChar(CP_UTF8, 0, utf8, -1, buffer, size);
    return n > 0 ? std::wstring_view(buffer, static_cast<size_t>(n - 1)) : std::wstring_view();
}

} // namespace

bool TrayIcon::create(const char* tooltip)
{
    std::snprintf(m_tooltip, sizeof(m_tooltip), "%s", tooltip);
    m_taskbar_created = RegisterWindowMessageW(L"TaskbarCreated");

    WNDCLASSEXW wc = {};
    wc.cbSize = sizeof(wc);
    wc.lpfnWndProc = tray_proc;
    wc.hInstance = GetModuleHandleW(nullptr);
    wc.lpszClassName = kTrayClass;
    RegisterClassExW(&wc); // fails harmlessly if already registered

    // A hidden top-level window (not message-only: TaskbarCreated is broadcast
    // to top-level windows).
    const HWND hwnd = CreateWindowExW(WS_EX_TOOLWINDOW, kTrayClass, L"Astraxis", WS_POPUP, 0, 0, 0, 0, nullptr,
                                      nullptr, GetModuleHandleW(nullptr), nullptr);
    if (!hwnd) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Creating the tray window failed: %lu", GetLastError());
        return false;
    }
    SetWindowLongPtrW(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    m_hwnd = hwnd;
    return add_icon();
}

bool TrayIcon::add_icon()
{
    NOTIFYICONDATAW data = {};
    data.cbSize = sizeof(data);
    data.hWnd = static_cast<HWND>(m_hwnd);
    data.uID = kIconId;
    data.uFlags = NIF_MESSAGE | NIF_ICON | NIF_TIP;
    data.uCallbackMessage = kCallbackMessage;
    // Icon group 1 of the resources (res/astraxis.rc), at the small icon size.
    data.hIcon = static_cast<HICON>(LoadImageW(GetModuleHandleW(nullptr), MAKEINTRESOURCEW(1), IMAGE_ICON,
                                               GetSystemMetrics(SM_CXSMICON), GetSystemMetrics(SM_CYSMICON), 0));
    widen(m_tooltip, data.szTip, static_cast<int>(sizeof(data.szTip) / sizeof(data.szTip[0])));
    const bool ok = Shell_NotifyIconW(NIM_ADD, &data) != FALSE;
    if (data.hIcon) {
        DestroyIcon(data.hIcon);
    }
    if (!ok) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Adding the tray icon failed");
    }
    return ok;
}

void TrayIcon::destroy()
{
    if (!m_hwnd) {
        return;
    }
    NOTIFYICONDATAW data = {};
    data.cbSize = sizeof(data);
    data.hWnd = static_cast<HWND>(m_hwnd);
    data.uID = kIconId;
    Shell_NotifyIconW(NIM_DELETE, &data);
    DestroyWindow(static_cast<HWND>(m_hwnd));
    m_hwnd = nullptr;
}

TrayIcon::Command TrayIcon::take_command()
{
    const Command command = m_command;
    m_command = Command::None;
    return command;
}

bool TrayIcon::handle_message(unsigned msg, long long lparam)
{
    if (msg == WM_CLOSE) {
        m_command = Command::Exit;
        return true;
    }
    if (msg == m_taskbar_created && m_taskbar_created != 0) {
        add_icon(); // Explorer restarted: the icon is gone
    } else if (msg == kCallbackMessage) {
        const UINT event = static_cast<UINT>(lparam);
        if (event == WM_RBUTTONUP || event == WM_CONTEXTMENU) {
            show_menu();
        } else if (event == WM_LBUTTONDBLCLK) {
            m_command = Command::TogglePause;
        }
    }
    return false;
}

void TrayIcon::show_menu()
{
    const HMENU menu = CreatePopupMenu();
    AppendMenuW(menu, MF_STRING, kMenuSettings, L"Settings...");
    SetMenuDefaultItem(menu, kMenuSettings, FALSE);
    AppendMenuW(menu, MF_STRING, kMenuPause, m_paused ? L"Resume" : L"Pause");
    AppendMenuW(menu, MF_SEPARATOR, 0, nullptr);
    AppendMenuW(menu, MF_STRING, kMenuExit, L"Exit");

    POINT cursor = {};
    GetCursorPos(&cursor);
    // The menu closes when it loses the focus only if our window is in front.
    const HWND hwnd = static_cast<HWND>(m_hwnd);
    SetForegroundWindow(hwnd);
    const UINT chosen = static_cast<UINT>(
        TrackPopupMenu(menu, TPM_RETURNCMD | TPM_NONOTIFY | TPM_RIGHTBUTTON, cursor.x, cursor.y, 0, hwnd, nullptr));
    PostMessageW(hwnd, WM_NULL, 0, 0);
    DestroyMenu(menu);

    if (chosen == kMenuSettings) {
        m_command = Command::Settings;
    } else if (chosen == kMenuPause) {
        m_command = Command::TogglePause;
    } else if (chosen == kMenuExit) {
        m_command = Command::Exit;
    }
}

bool InstanceLock::acquire(const char* name)
{
    wchar_t wide[128];
    widen(name, wide, 128);
    m_handle = CreateMutexW(nullptr, TRUE, wide);
    if (m_handle && GetLastError() == ERROR_ALREADY_EXISTS) {
        CloseHandle(static_cast<HANDLE>(m_handle));
        m_handle = nullptr;
        return false;
    }
    return m_handle != nullptr;
}

InstanceLock::~InstanceLock()
{
    if (m_handle) {
        ReleaseMutex(static_cast<HANDLE>(m_handle));
        CloseHandle(static_cast<HANDLE>(m_handle));
    }
}

} // namespace astraxis
