#pragma once

#include "render/black_hole_uniforms.hpp"
#include "render/corona_uniforms.hpp"
#include "render/gpu_device.hpp"

#include <SDL3/SDL_gpu.h>

#include <cstdint>
#include <vector>

namespace astraxis {

// PostProcess's bloom chain (halving sizes), sized with its output.
struct BloomChain {
    struct Level {
        SDL_GPUTexture* texture = nullptr;
        uint32_t width = 0;
        uint32_t height = 0;
    };
    std::vector<Level> levels;
    uint32_t width = 0; // the scene size the chain was made for
    uint32_t height = 0;
};

// BlackHolePass's off-screen trace result (reduced resolution).
struct BlackHoleTargets {
    SDL_GPUTexture* color = nullptr; // RGBA16F: disk light * fade, sky transmission * fade
    SDL_GPUTexture* sky = nullptr;   // RGBA32F: sky direction, fade (negative where opaque)
    uint32_t width = 0;
    uint32_t height = 0;
    bool traced = false;         // a trace waits to be composited
    BlackHoleUniforms uniforms;  // of that trace (the composite traces again at edges)
};

// A SplatPass's off-screen nebula (reduced resolution, one sample per pixel).
struct SplatTargets {
    SDL_GPUTexture* color = nullptr; // RGBA16F
    uint32_t width = 0;
    uint32_t height = 0;
    bool drawn = false; // waits to be composited
};

struct CoronaTargets {
    SDL_GPUTexture* color = nullptr;
    uint32_t width = 0;
    uint32_t height = 0;
    bool drawn = false;
    CoronaUniforms uniforms;
};

// The render targets of one output, all of its size: MSAA HDR color resolved
// into a sampleable HDR texture, reversed-Z depth, and the size-dependent
// textures of the passes (bloom chain, black-hole trace). Each output keeps
// its own, so that outputs of different sizes drawn in turn do not make the
// passes reallocate every frame.
class SceneTargets {
public:
    // (Re)creates the scene targets for `width` x `height` pixels (and the
    // format's sample count) if needed.
    bool ensure(SDL_GPUDevice* device, const SceneTargetFormat& format, uint32_t width, uint32_t height);
    // Releases every texture, the passes' included.
    void release();

    // Begins the scene pass: clears color and depth (to 0, reversed-Z) and
    // resolves MSAA into hdr() at the end.
    SDL_GPURenderPass* begin_scene_pass(SDL_GPUCommandBuffer* cmd, SDL_FColor clear_color) const;
    SDL_GPUTexture* hdr() const { return m_hdr; }

    CoronaTargets corona;
    BloomChain bloom;              // managed by PostProcess
    BlackHoleTargets black_hole;   // managed by BlackHolePass
    SplatTargets pinwheel;         // managed by PinwheelPass
    SplatTargets nebula;           // managed by NebulaPass

private:
    void release_scene();

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUTexture* m_msaa_color = nullptr;
    SDL_GPUTexture* m_hdr = nullptr;
    SDL_GPUTexture* m_depth = nullptr;
    uint32_t m_width = 0;
    uint32_t m_height = 0;
    SDL_GPUSampleCount m_samples = SDL_GPU_SAMPLECOUNT_1;
};

} // namespace astraxis
