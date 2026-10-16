#pragma once

#include <SDL3/SDL_gpu.h>

namespace astraxis {

// Formats and sample count shared by all scene pipelines. The scene is
// rendered in linear HDR and tonemapped to the swapchain by PostProcess.
struct SceneTargetFormat {
    SDL_GPUTextureFormat color = SDL_GPU_TEXTUREFORMAT_INVALID;
    SDL_GPUTextureFormat depth = SDL_GPU_TEXTUREFORMAT_INVALID;
    SDL_GPUSampleCount samples = SDL_GPU_SAMPLECOUNT_1;
};

// Owns the SDL_GPU device, shared by every output (window) it draws to.
class GpuDevice {
public:
    bool init(bool debug);
    void shutdown();

    SDL_GPUDevice* device() const { return m_device; }
    const SceneTargetFormat& scene_format() const { return m_scene_format; }
    const char* driver_name() const;

    // The scene's MSAA sample count: 1 (off), 2 or 4, lowered to what the
    // device supports for the scene formats (4 after init). True if it
    // changed; the scene pipelines must then be recreated
    // (SceneRenderer::refresh_pipelines).
    bool set_msaa(int samples);
    int msaa() const;
    int max_msaa() const { return m_max_msaa; }

private:
    bool supports_msaa(int samples) const;

    SDL_GPUDevice* m_device = nullptr;
    SceneTargetFormat m_scene_format;
    int m_max_msaa = 1;
};

} // namespace astraxis
