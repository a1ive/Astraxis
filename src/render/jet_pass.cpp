#include "render/jet_pass.hpp"

#include "render/shader.hpp"

#include <shaders/fullscreen.vert.h>
#include <shaders/jet.frag.h>

#include <SDL3/SDL_log.h>

namespace astraxis {

bool JetPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;
    return create_pipeline(format);
}

bool JetPass::set_format(const SceneTargetFormat& format)
{
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_pipeline);
    m_pipeline = nullptr;
    return create_pipeline(format);
}

bool JetPass::create_pipeline(const SceneTargetFormat& format)
{
    SDL_GPUShader* vs = create_shader(
        m_device, {.dxil = kFullscreenVertDxil, .spirv = kFullscreenVertSpirv, .stage = SDL_GPU_SHADERSTAGE_VERTEX});
    SDL_GPUShader* fs = create_shader(m_device, {.dxil = kJetFragDxil,
                                                 .spirv = kJetFragSpirv,
                                                 .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                                 .num_uniform_buffers = 1});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(m_device, vs);
        SDL_ReleaseGPUShader(m_device, fs);
        return false;
    }

    // Emission only: added, no depth.
    SDL_GPUColorTargetDescription color = {};
    color.format = format.color;
    color.blend_state.enable_blend = true;
    color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ZERO;
    color.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    SDL_GPUGraphicsPipelineCreateInfo info = {};
    info.vertex_shader = vs;
    info.fragment_shader = fs;
    info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
    info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
    info.multisample_state.sample_count = format.samples;
    info.depth_stencil_state.enable_depth_test = false;
    info.depth_stencil_state.enable_depth_write = false;
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = true;
    info.target_info.depth_stencil_format = format.depth;
    m_pipeline = SDL_CreateGPUGraphicsPipeline(m_device, &info);

    SDL_ReleaseGPUShader(m_device, vs);
    SDL_ReleaseGPUShader(m_device, fs);
    if (!m_pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Jet pipeline creation failed: %s", SDL_GetError());
        return false;
    }
    return true;
}

void JetPass::shutdown()
{
    if (!m_device) {
        return;
    }
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_pipeline);
    m_pipeline = nullptr;
    m_device = nullptr;
}

void JetPass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const JetUniforms& u) const
{
    SDL_BindGPUGraphicsPipeline(pass, m_pipeline);
    SDL_PushGPUFragmentUniformData(cmd, 0, &u, sizeof(u));
    SDL_DrawGPUPrimitives(pass, 3, 1, 0, 0);
}

} // namespace astraxis
