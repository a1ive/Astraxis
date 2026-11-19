#pragma once

#include "render/jet_uniforms.hpp"
#include "render/scene_targets.hpp"

#include <SDL3/SDL_gpu.h>

namespace astraxis {

// A black hole's jet in visible light (shaders/jet.frag.hlsl): a full-screen pass in
// the scene pass, adding the light integrated along each view ray. It draws after the
// hole's composite (the jet is not lensed; rays into the shadow end at the hole).
class JetPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();
    // Recreates the pipeline for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format);

    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const JetUniforms& u) const;

private:
    bool create_pipeline(const SceneTargetFormat& format);

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
};

} // namespace astraxis
