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

private:
    SDL_GPUDevice* m_device = nullptr;
    SceneTargetFormat m_scene_format;
};

} // namespace astraxis
