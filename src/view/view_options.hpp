#pragma once

#include "render/post_settings.hpp"

namespace astraxis {

// What to draw and how, chosen by the host (SceneRenderer).
struct ViewOptions {
    bool orbits = true;
    bool belts = true;
    bool atmospheres = true; // otherwise the surface under Venus' and Titan's clouds
    bool plumes = true;
    bool comets = true;      // comae, ion and dust tails
    float star_brightness = 1.0f;
    float line_width = 1.6f;    // pixels
    PostSettings post;
};

} // namespace astraxis
