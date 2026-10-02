#include "app/app.hpp"

#include <SDL3/SDL_main.h>

int main(int /*argc*/, char* /*argv*/[])
{
    astraxis::App app;
    const bool ok = app.init();
    if (ok) {
        app.run();
    }
    app.shutdown();
    return ok ? 0 : 1;
}
