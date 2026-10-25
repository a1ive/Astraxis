#pragma once

#include "render/corona_uniforms.hpp"
#include "render/gpu_device.hpp"
#include "render/scene_targets.hpp"

namespace astraxis {

// Integrates the optically thin corona off-screen, then adds it to the scene.
class CoronaPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    bool set_format(const SceneTargetFormat& format);
    void shutdown();
    void render(SDL_GPUCommandBuffer* cmd, CoronaTargets& targets, const CoronaUniforms& uniforms,
                uint32_t width, uint32_t height);
    void composite(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, CoronaTargets& targets);

private:
    bool create_pipelines(const SceneTargetFormat& format);
    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_trace = nullptr;
    SDL_GPUGraphicsPipeline* m_composite = nullptr;
    SDL_GPUSampler* m_linear = nullptr;
};

} // namespace astraxis
