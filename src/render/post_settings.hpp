#pragma once

namespace astraxis {

// Settings of PostProcess. Kept free of SDL so that settings files can be
// read without the renderer.
struct PostSettings {
    float exposure = 1.0f;
    float bloom_strength = 0.04f;
    float bloom_threshold = 1.5f; // linear HDR; only brighter light blooms
};

} // namespace astraxis
