#include "app/app.hpp"

#include <SDL3/SDL_log.h>
#include <SDL3/SDL_main.h>

#include <cstdlib>
#include <string_view>

int main(int argc, char* argv[])
{
    astraxis::LaunchOptions options;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        if (arg == "--scene" && i + 1 < argc) {
            options.scene = argv[++i];
        } else if (arg == "--event" && i + 1 < argc) {
            options.event = std::atoi(argv[++i]);
        } else {
            SDL_Log("Ignoring argument '%s' (usage: astraxis [--scene <name>] [--event <n>])", argv[i]);
        }
    }

    astraxis::App app;
    const bool ok = app.init(options);
    if (ok) {
        app.run();
    }
    app.shutdown();
    return ok ? 0 : 1;
}
