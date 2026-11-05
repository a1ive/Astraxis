#pragma once

#include "render/post_settings.hpp"

namespace astraxis {

// What to draw and how, chosen by the host (SceneRenderer).
struct ViewOptions {
    bool orbits = true;
    bool belts = true;
    bool atmospheres = true; // otherwise the surface under Venus' and Titan's clouds; also Earth's aurora
    int aurora_view = 0; // 0: event setting (off by default), 1: visible light, 2: UV and IR overlay + visible
    int magnetic_view = 0; // field lines. 0: event setting (off by default), 1: always, 2: never
    int magnetic_map_view = 0; // radial field maps. 0: event setting (off by default), 1: surface, 2: deep, 3: never
    bool plumes = true;
    bool comets = true;      // comae, ion and dust tails
    float star_brightness = 1.0f;
    float line_width = 1.6f;    // pixels
    PostSettings post;
};

} // namespace astraxis
