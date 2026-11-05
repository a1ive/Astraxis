#pragma once

#include "render/disk_uniforms.hpp"
#include "render/gpu_device.hpp"
#include "render/scene_targets.hpp"

namespace astraxis {

// Integrates an optically thin disk's light along each ray off-screen at half
// resolution (one sample per pixel), then adds it to the scene, integrating
// again at full resolution around what blocks the rays (stars and planets).
class DiskPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    bool set_format(const SceneTargetFormat& format);
    void shutdown();
    // Nothing is drawn when uniforms.params.x is 0.
    void render(SDL_GPUCommandBuffer* cmd, DiskTargets& targets, const DiskUniforms& uniforms, uint32_t width,
                uint32_t height);
    void composite(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, DiskTargets& targets);

private:
    bool create_pipelines(const SceneTargetFormat& format);

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_trace = nullptr;
    SDL_GPUGraphicsPipeline* m_composite = nullptr;
    SDL_GPUSampler* m_linear = nullptr;
};

} // namespace astraxis
