#include "app/app.hpp"
#include "platform/paths.hpp"

#include <SDL3/SDL_log.h>
#include <SDL3/SDL_main.h>

#include <cerrno>
#include <cstdlib>
#include <string_view>

namespace {

constexpr const char* kUsage = "usage: astraxis [--scene <name>] [--event <n>] [--mode window|fullscreen]\n"
                               "                [--display <n>] [--fps <n>] [--config <path>]";

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

} // namespace

int main(int argc, char* argv[])
{
    astraxis::LaunchOptions options;
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
            ok = mode == "window" || mode == "fullscreen";
            if (ok) {
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

    astraxis::App app;
    const bool ok = app.init(options);
    if (ok) {
        app.run();
    }
    app.shutdown();
    return ok ? 0 : 1;
}
