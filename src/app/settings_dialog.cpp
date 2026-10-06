#include "app/settings_dialog.hpp"

#include "app/settings.hpp"
#include "platform/window.hpp"

#include <SDL3/SDL.h>

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <commctrl.h>

#include "settings_dialog.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <string>
#include <vector>

namespace astraxis {

namespace {

constexpr const wchar_t* kRunKey = L"Software\\Microsoft\\Windows\\CurrentVersion\\Run";
constexpr const wchar_t* kRunValue = L"Astraxis Wallpaper";
constexpr const wchar_t* kTrayClass = L"AstraxisTray"; // the running wallpaper's tray window
constexpr UINT_PTR kStatusTimer = 1;
constexpr int kFpsChoices[] = {0, 60, 30, 24, 15};

std::wstring widen(const std::string& utf8)
{
    const int n = MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(), -1, nullptr, 0);
    std::wstring out(n > 0 ? static_cast<size_t>(n - 1) : 0, L'\0');
    if (n > 1) {
        MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(), -1, out.data(), n);
    }
    return out;
}

// A slider: a float range mapped onto integer trackbar positions, linearly
// or logarithmically.
struct Slider {
    int id;
    int value_id;
    float lo;
    float hi;
    int steps;
    bool log;
    const wchar_t* format;

    int position(float v) const
    {
        v = std::clamp(v, lo, hi);
        const float t = log ? std::log(v / lo) / std::log(hi / lo) : (v - lo) / (hi - lo);
        return static_cast<int>(std::lround(t * steps));
    }
    float value(int pos) const
    {
        const float t = static_cast<float>(pos) / steps;
        return log ? lo * std::pow(hi / lo, t) : lo + t * (hi - lo);
    }
};

// Ranges as on the window's control panel.
constexpr Slider kSliders[] = {
    {IDC_STARS, IDC_STARS_VALUE, 0.0f, 2.0f, 40, false, L"%.2f"},
    {IDC_LINE_WIDTH, IDC_LINE_WIDTH_VALUE, 0.5f, 4.0f, 35, false, L"%.1f px"},
    {IDC_EXPOSURE, IDC_EXPOSURE_VALUE, 0.1f, 8.0f, 100, true, L"%.2f"},
    {IDC_BLOOM, IDC_BLOOM_VALUE, 0.0f, 0.2f, 40, false, L"%.3f"},
};

float* slider_field(Settings& s, int id)
{
    switch (id) {
    case IDC_STARS:
        return &s.view.star_brightness;
    case IDC_LINE_WIDTH:
        return &s.view.line_width;
    case IDC_EXPOSURE:
        return &s.view.post.exposure;
    default:
        return &s.view.post.bloom_strength;
    }
}

struct Check {
    int id;
    bool Settings::*field;
    bool ViewOptions::*view_field;
};

constexpr Check kChecks[] = {
    {IDC_LABELS, &Settings::labels, nullptr},
    {IDC_ORBITS, nullptr, &ViewOptions::orbits},
    {IDC_BELTS, nullptr, &ViewOptions::belts},
    {IDC_ATMOSPHERES, nullptr, &ViewOptions::atmospheres},
    {IDC_PLUMES, nullptr, &ViewOptions::plumes},
    {IDC_COMETS, nullptr, &ViewOptions::comets},
};

bool& check_field(Settings& s, const Check& c)
{
    return c.field ? s.*(c.field) : s.view.*(c.view_field);
}

struct DisplayEntry {
    int index; // 1-based, as in the settings
    std::string name;
};

struct DialogState {
    const SettingsDialogContext* context = nullptr;
    SettingsSection section = SettingsSection::Wallpaper;
    Settings settings;
    bool config_ok = true;
    std::vector<SceneEntry> scenes;
    std::vector<DisplayEntry> displays;
    std::vector<int> fps_values; // per FPS combo item
    // Slider positions as set up: an untouched slider keeps its exact value
    // (positions are rounded: exposure 1.0 would come back as 1.02).
    int slider_start[std::size(kSliders)] = {};
};

std::vector<DisplayEntry> list_displays(std::vector<std::wstring>& labels)
{
    std::vector<DisplayEntry> displays;
    const bool init = !SDL_WasInit(SDL_INIT_VIDEO);
    if (init && !SDL_InitSubSystem(SDL_INIT_VIDEO)) {
        return displays;
    }
    int count = 0;
    SDL_DisplayID* ids = SDL_GetDisplays(&count);
    for (int i = 0; ids && i < count; ++i) {
        DisplayEntry entry{i + 1, display_name(ids[i])};
        wchar_t label[256];
        const SDL_DisplayMode* mode = SDL_GetDesktopDisplayMode(ids[i]);
        std::swprintf(label, 256, L"%d: %ls (%dx%d)", i + 1, widen(entry.name).c_str(), mode ? mode->w : 0,
                      mode ? mode->h : 0);
        labels.emplace_back(label);
        displays.push_back(std::move(entry));
    }
    SDL_free(ids);
    if (init) {
        SDL_QuitSubSystem(SDL_INIT_VIDEO);
    }
    return displays;
}

// --- Wallpaper process and autostart ---

HWND running_wallpaper()
{
    return FindWindowW(kTrayClass, nullptr);
}

std::wstring autostart_command(const std::filesystem::path& program)
{
    return L"\"" + program.wstring() + L"\" --mode wallpaper";
}

bool autostart_enabled()
{
    HKEY key = nullptr;
    if (RegOpenKeyExW(HKEY_CURRENT_USER, kRunKey, 0, KEY_QUERY_VALUE, &key) != ERROR_SUCCESS) {
        return false;
    }
    const bool found = RegQueryValueExW(key, kRunValue, nullptr, nullptr, nullptr, nullptr) == ERROR_SUCCESS;
    RegCloseKey(key);
    return found;
}

void set_autostart(bool enabled, const std::filesystem::path& program)
{
    HKEY key = nullptr;
    if (RegOpenKeyExW(HKEY_CURRENT_USER, kRunKey, 0, KEY_SET_VALUE, &key) != ERROR_SUCCESS) {
        return;
    }
    if (enabled) {
        const std::wstring command = autostart_command(program);
        RegSetValueExW(key, kRunValue, 0, REG_SZ, reinterpret_cast<const BYTE*>(command.c_str()),
                       static_cast<DWORD>((command.size() + 1) * sizeof(wchar_t)));
    } else {
        RegDeleteValueW(key, kRunValue);
    }
    RegCloseKey(key);
}

// Asks a running wallpaper to exit (through its tray window); waits up to 5 s.
void stop_wallpaper()
{
    const HWND tray = running_wallpaper();
    if (!tray) {
        return;
    }
    PostMessageW(tray, WM_CLOSE, 0, 0);
    for (int waited = 0; waited < 5000 && running_wallpaper(); waited += 50) {
        Sleep(50);
    }
}

bool start_wallpaper(const std::filesystem::path& program)
{
    std::wstring command = autostart_command(program);
    const std::wstring dir = program.parent_path().wstring();
    STARTUPINFOW startup = {};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process = {};
    if (!CreateProcessW(program.c_str(), command.data(), nullptr, nullptr, FALSE, 0, nullptr, dir.c_str(), &startup,
                        &process)) {
        return false;
    }
    CloseHandle(process.hThread);
    CloseHandle(process.hProcess);
    return true;
}

// --- The dialog ---

void set_slider_label(HWND dialog, const Slider& slider, float value)
{
    wchar_t text[32];
    std::swprintf(text, 32, slider.format, static_cast<double>(value));
    SetDlgItemTextW(dialog, slider.value_id, text);
}

void update_status(HWND dialog)
{
    SetDlgItemTextW(dialog, IDC_STATUS,
                    running_wallpaper() ? L"The wallpaper is running." : L"The wallpaper is not running.");
    EnableWindow(GetDlgItem(dialog, IDC_STOP), running_wallpaper() != nullptr);
}

void init_dialog(HWND dialog, DialogState& state)
{
    const bool wallpaper = state.context->kind == SettingsDialogKind::Wallpaper;
    SetWindowTextW(dialog, wallpaper ? L"Astraxis wallpaper" : L"Astraxis screensaver");
    const Settings& s = state.settings;

    // Scenes: by their own names; one missing from the list is kept as it is.
    const HWND scene = GetDlgItem(dialog, IDC_SCENE);
    int selected = -1;
    for (size_t i = 0; i < state.scenes.size(); ++i) {
        SendMessageW(scene, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(widen(state.scenes[i].name).c_str()));
        if (state.scenes[i].stem == s.scene) {
            selected = static_cast<int>(i);
        }
    }
    if (selected < 0) {
        state.scenes.push_back({s.scene, s.scene});
        SendMessageW(scene, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(widen(s.scene).c_str()));
        selected = static_cast<int>(state.scenes.size() - 1);
    }
    SendMessageW(scene, CB_SETCURSEL, static_cast<WPARAM>(selected), 0);

    // Displays: "All displays", then each; the saved one found again by name.
    const HWND display = GetDlgItem(dialog, IDC_DISPLAY);
    SendMessageW(display, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(L"All displays"));
    std::vector<std::wstring> labels;
    state.displays = list_displays(labels);
    int display_item = 0;
    for (size_t i = 0; i < labels.size(); ++i) {
        SendMessageW(display, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(labels[i].c_str()));
    }
    if (s.display > 0) {
        for (size_t i = 0; i < state.displays.size(); ++i) {
            if (state.displays[i].name == s.display_name && !s.display_name.empty()) {
                display_item = static_cast<int>(i + 1);
                break;
            }
        }
        if (display_item == 0 && s.display <= static_cast<int>(state.displays.size())) {
            display_item = s.display;
        }
    }
    SendMessageW(display, CB_SETCURSEL, static_cast<WPARAM>(display_item), 0);

    // Frame rate.
    const HWND fps = GetDlgItem(dialog, IDC_FPS);
    state.fps_values.assign(std::begin(kFpsChoices), std::end(kFpsChoices));
    if (std::find(state.fps_values.begin(), state.fps_values.end(), s.fps) == state.fps_values.end()) {
        state.fps_values.push_back(s.fps);
    }
    for (int value : state.fps_values) {
        wchar_t text[64];
        if (value == 0) {
            std::swprintf(text, 64, L"The display's refresh rate");
        } else {
            std::swprintf(text, 64, L"At most %d FPS", value);
        }
        SendMessageW(fps, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(text));
    }
    const auto fps_item = std::find(state.fps_values.begin(), state.fps_values.end(), s.fps);
    SendMessageW(fps, CB_SETCURSEL, static_cast<WPARAM>(fps_item - state.fps_values.begin()), 0);

    Settings& editable = state.settings;
    for (const Check& c : kChecks) {
        CheckDlgButton(dialog, c.id, check_field(editable, c) ? BST_CHECKED : BST_UNCHECKED);
    }
    for (size_t i = 0; i < std::size(kSliders); ++i) {
        const Slider& slider = kSliders[i];
        state.slider_start[i] = slider.position(*slider_field(editable, slider.id));
        SendDlgItemMessageW(dialog, slider.id, TBM_SETRANGE, FALSE, MAKELPARAM(0, slider.steps));
        SendDlgItemMessageW(dialog, slider.id, TBM_SETPOS, TRUE, state.slider_start[i]);
        set_slider_label(dialog, slider, *slider_field(editable, slider.id));
    }

    if (wallpaper) {
        CheckDlgButton(dialog, IDC_AUTOSTART, autostart_enabled() ? BST_CHECKED : BST_UNCHECKED);
        SetDlgItemTextW(dialog, IDOK, L"Apply");
        SetDlgItemTextW(dialog, IDCANCEL, L"Close");
        update_status(dialog);
        SetTimer(dialog, kStatusTimer, 1000, nullptr);
    } else {
        ShowWindow(GetDlgItem(dialog, IDC_AUTOSTART), SW_HIDE);
        ShowWindow(GetDlgItem(dialog, IDC_STATUS), SW_HIDE);
        ShowWindow(GetDlgItem(dialog, IDC_STOP), SW_HIDE);
    }
}

// Reads the controls into state.settings.
void read_dialog(HWND dialog, DialogState& state)
{
    Settings& s = state.settings;
    const LRESULT scene = SendDlgItemMessageW(dialog, IDC_SCENE, CB_GETCURSEL, 0, 0);
    if (scene >= 0 && static_cast<size_t>(scene) < state.scenes.size()) {
        s.scene = state.scenes[static_cast<size_t>(scene)].stem;
    }
    const LRESULT display = SendDlgItemMessageW(dialog, IDC_DISPLAY, CB_GETCURSEL, 0, 0);
    if (display > 0 && static_cast<size_t>(display) <= state.displays.size()) {
        const DisplayEntry& entry = state.displays[static_cast<size_t>(display - 1)];
        s.display = entry.index;
        s.display_name = entry.name;
    } else {
        s.display = 0;
        s.display_name.clear();
    }
    const LRESULT fps = SendDlgItemMessageW(dialog, IDC_FPS, CB_GETCURSEL, 0, 0);
    if (fps >= 0 && static_cast<size_t>(fps) < state.fps_values.size()) {
        s.fps = state.fps_values[static_cast<size_t>(fps)];
    }
    for (const Check& c : kChecks) {
        check_field(s, c) = IsDlgButtonChecked(dialog, c.id) == BST_CHECKED;
    }
    for (size_t i = 0; i < std::size(kSliders); ++i) {
        const Slider& slider = kSliders[i];
        const int pos = static_cast<int>(SendDlgItemMessageW(dialog, slider.id, TBM_GETPOS, 0, 0));
        if (pos != state.slider_start[i]) {
            *slider_field(s, slider.id) = slider.value(pos);
            state.slider_start[i] = pos; // applied: the new value is exact from now on
        }
    }
}

bool save(HWND dialog, DialogState& state)
{
    if (!state.config_ok) {
        MessageBoxW(dialog, L"config.toml could not be read, so it is not overwritten. Fix or remove it first.",
                    L"Astraxis", MB_ICONWARNING);
        return false;
    }
    read_dialog(dialog, state);
    std::string error;
    if (!save_settings(state.context->config, state.section, state.settings, &error, ViewScope::Overrides)) {
        MessageBoxW(dialog, widen("The settings could not be saved:\n" + error).c_str(), L"Astraxis", MB_ICONERROR);
        return false;
    }
    return true;
}

INT_PTR CALLBACK dialog_proc(HWND dialog, UINT msg, WPARAM wparam, LPARAM lparam)
{
    auto* state = reinterpret_cast<DialogState*>(GetWindowLongPtrW(dialog, DWLP_USER));
    switch (msg) {
    case WM_INITDIALOG:
        SetWindowLongPtrW(dialog, DWLP_USER, lparam);
        init_dialog(dialog, *reinterpret_cast<DialogState*>(lparam));
        return TRUE;
    case WM_HSCROLL:
        for (const Slider& slider : kSliders) {
            if (reinterpret_cast<HWND>(lparam) == GetDlgItem(dialog, slider.id)) {
                const int pos = static_cast<int>(SendDlgItemMessageW(dialog, slider.id, TBM_GETPOS, 0, 0));
                set_slider_label(dialog, slider, slider.value(pos));
            }
        }
        return TRUE;
    case WM_TIMER:
        if (wparam == kStatusTimer) {
            update_status(dialog);
        }
        return TRUE;
    case WM_COMMAND:
        switch (LOWORD(wparam)) {
        case IDOK:
            if (!state || !save(dialog, *state)) {
                return TRUE;
            }
            if (state->context->kind == SettingsDialogKind::Screensaver) {
                EndDialog(dialog, IDOK);
                return TRUE;
            }
            // Apply: autostart, and (re)start the wallpaper with the new settings.
            set_autostart(IsDlgButtonChecked(dialog, IDC_AUTOSTART) == BST_CHECKED, state->context->program);
            stop_wallpaper();
            if (!start_wallpaper(state->context->program)) {
                MessageBoxW(dialog, L"The wallpaper could not be started.", L"Astraxis", MB_ICONERROR);
            }
            update_status(dialog);
            return TRUE;
        case IDC_STOP:
            stop_wallpaper();
            update_status(dialog);
            return TRUE;
        case IDCANCEL:
            KillTimer(dialog, kStatusTimer);
            EndDialog(dialog, IDCANCEL);
            return TRUE;
        default:
            break;
        }
        break;
    default:
        break;
    }
    return FALSE;
}

} // namespace

void run_settings_dialog(const SettingsDialogContext& context, void* parent)
{
    INITCOMMONCONTROLSEX controls = {};
    controls.dwSize = sizeof(controls);
    controls.dwICC = ICC_BAR_CLASSES | ICC_STANDARD_CLASSES;
    InitCommonControlsEx(&controls);

    DialogState state;
    state.context = &context;
    state.section =
        context.kind == SettingsDialogKind::Wallpaper ? SettingsSection::Wallpaper : SettingsSection::Screensaver;
    std::string error;
    state.config_ok = load_settings(context.config, state.section, state.settings, &error);
    state.scenes = list_scenes(context.asset_dir);
    if (!state.config_ok) {
        MessageBoxW(static_cast<HWND>(parent),
                    widen("config.toml could not be read; the defaults are shown and nothing will be saved:\n" +
                          error)
                        .c_str(),
                    L"Astraxis", MB_ICONWARNING);
    }
    DialogBoxParamW(GetModuleHandleW(nullptr), MAKEINTRESOURCEW(IDD_SETTINGS), static_cast<HWND>(parent), dialog_proc,
                    reinterpret_cast<LPARAM>(&state));
}

} // namespace astraxis
