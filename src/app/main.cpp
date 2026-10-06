#include "app/app.hpp"
#include "app/screensaver.hpp"
#include "platform/paths.hpp"

#include <SDL3/SDL_log.h>
#include <SDL3/SDL_main.h>
#include <SDL3/SDL_messagebox.h>

#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <string_view>

namespace {

constexpr const char* kUsage =
    "usage: astraxis [--scene <name>] [--event <n>] [--mode window|fullscreen|screensaver]\n"
    "                [--display <n>] [--fps <n>] [--config <path>]\n"
    "       astraxis /s | /p <hwnd> | /c   (the Windows screensaver protocol)";

bool parse_int(const char* text, int& out)
{
    char* end = nullptr;
    errno = 0;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || errno != 0 || value < 0 || value > 1000000) {
        return false;
    }
    out = static_cast<int>(value);
    return true;
}

// The screensaver protocol: Windows runs a .scr with /s (run), /p <hwnd>
// (preview in that window), /c[:<hwnd>] (settings) or nothing (settings), in
// any case, with '-' for '/' and ':' or a space before the window handle.
struct ScreensaverCommand {
    char verb = 0; // 's', 'p', 'c' or 'a' (change password: Windows 9x only)
    uint64_t hwnd = 0;
};

bool parse_screensaver_command(int argc, char* argv[], ScreensaverCommand& out)
{
    if (argc < 2) {
        return false;
    }
    const std::string_view arg = argv[1];
    if (arg.size() < 2 || (arg[0] != '/' && arg[0] != '-')) {
        return false;
    }
    const char verb = static_cast<char>(std::tolower(static_cast<unsigned char>(arg[1])));
    if ((verb != 's' && verb != 'p' && verb != 'c' && verb != 'a') || (arg.size() > 2 && arg[2] != ':')) {
        return false;
    }
    out.verb = verb;
    const char* handle = arg.size() > 3 ? argv[1] + 3 : (argc > 2 ? argv[2] : nullptr);
    if (handle) {
        out.hwnd = std::strtoull(handle, nullptr, 10);
    }
    return true;
}

bool running_as_scr(const char* argv0)
{
    std::string ext = astraxis::path_from_utf8(argv0).extension().string();
    for (char& c : ext) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return ext == ".scr";
}

int run_screensaver(const astraxis::ScreensaverOptions& options)
{
    astraxis::Screensaver saver;
    const bool ok = saver.init(options);
    if (ok) {
        saver.run();
    }
    saver.shutdown();
    return ok ? 0 : 1;
}

void show_screensaver_settings()
{
    const std::string text = "Astraxis has no screensaver settings dialog yet.\n\n"
                             "The screensaver reads the [view] and [screensaver] sections of\n" +
                             astraxis::config_path().string();
    SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_INFORMATION, "Astraxis", text.c_str(), nullptr);
}

// ASTRAXIS_LOG=<file> appends the log to that file: the Release build and a
// screensaver started by Windows have no console.
void log_to_file(void* userdata, int /*category*/, SDL_LogPriority priority, const char* message)
{
    static constexpr const char* kNames[] = {"", "TRACE", "VERBOSE", "DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"};
    const int p = static_cast<int>(priority);
    std::FILE* file = static_cast<std::FILE*>(userdata);
    std::fprintf(file, "%s: %s\n", p >= 0 && p < 8 ? kNames[p] : "", message);
    std::fflush(file);
}

void init_log_file()
{
    const char* path = SDL_getenv("ASTRAXIS_LOG");
    if (!path || !*path) {
        return;
    }
    if (std::FILE* file = std::fopen(path, "a")) {
        SDL_SetLogOutputFunction(log_to_file, file); // the file stays open until exit
    }
}

} // namespace

int main(int argc, char* argv[])
{
    init_log_file();
    ScreensaverCommand command;
    if (parse_screensaver_command(argc, argv, command)) {
        switch (command.verb) {
        case 's':
            return run_screensaver({});
        case 'p': {
            if (command.hwnd == 0) {
                SDL_Log("/p needs a window handle");
                return 1;
            }
            astraxis::ScreensaverOptions options;
            options.preview = true;
            options.parent = reinterpret_cast<void*>(static_cast<uintptr_t>(command.hwnd));
            return run_screensaver(options);
        }
        case 'c':
            show_screensaver_settings();
            return 0;
        default:
            return 0; // 'a': no passwords
        }
    }
    if (argc < 2 && running_as_scr(argv[0])) {
        show_screensaver_settings();
        return 0;
    }

    astraxis::LaunchOptions options;
    bool screensaver = false;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        const char* value = i + 1 < argc ? argv[i + 1] : nullptr;
        int number = 0;
        bool ok = value != nullptr;
        if (arg == "--scene" && ok) {
            options.scene = value;
        } else if (arg == "--event" && ok) {
            ok = parse_int(value, options.event);
        } else if (arg == "--mode" && ok) {
            const std::string_view mode = value;
            ok = mode == "window" || mode == "fullscreen" || mode == "screensaver";
            if (ok) {
                screensaver = mode == "screensaver";
                options.fullscreen = mode == "fullscreen";
            }
        } else if (arg == "--display" && ok) {
            ok = parse_int(value, number);
            if (ok) {
                options.display = number;
            }
        } else if (arg == "--fps" && ok) {
            ok = parse_int(value, number);
            if (ok) {
                options.fps = number;
            }
        } else if (arg == "--config" && ok) {
            options.config = astraxis::path_from_utf8(value);
        } else {
            SDL_Log("Ignoring argument '%s'\n%s", argv[i], kUsage);
            continue;
        }
        if (!ok) {
            SDL_Log("Invalid value for %s: '%s'\n%s", argv[i], value ? value : "", kUsage);
        }
        ++i; // the value
    }

    if (screensaver) {
        astraxis::ScreensaverOptions saver;
        saver.config = options.config;
        saver.scene = options.scene;
        saver.display = options.display;
        saver.fps = options.fps;
        return run_screensaver(saver);
    }

    astraxis::App app;
    const bool ok = app.init(options);
    if (ok) {
        app.run();
    }
    app.shutdown();
    return ok ? 0 : 1;
}
