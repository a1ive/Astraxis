#include "render/beam_pass.hpp"

#include "render/shader.hpp"

#include <shaders/beam.frag.h>
#include <shaders/beam.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

#include <cmath>

namespace astraxis {

namespace {

struct VertexUniforms {
    glm::mat4 view_proj;
    glm::vec4 apex;
    glm::vec4 axis;
};

struct FragmentUniforms {
    glm::vec4 color;
};

} // namespace

bool BeamPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;

    SDL_GPUShader* vs = create_shader(device, {.dxil = kBeamVertDxil,
                                               .spirv = kBeamVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kBeamFragDxil,
                                               .spirv = kBeamFragSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                               .num_uniform_buffers = 1});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, fs);
        return false;
    }

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
    info.depth_stencil_state.enable_depth_test = true;
    info.depth_stencil_state.enable_depth_write = false;
    info.depth_stencil_state.compare_op = SDL_GPU_COMPAREOP_GREATER_OR_EQUAL; // reversed-Z
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = true;
    info.target_info.depth_stencil_format = format.depth;

    m_pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, fs);
    if (!m_pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Beam pipeline creation failed: %s", SDL_GetError());
        return false;
    }
    return true;
}

void BeamPass::shutdown()
{
    if (m_device) {
        SDL_ReleaseGPUGraphicsPipeline(m_device, m_pipeline);
    }
    m_pipeline = nullptr;
    m_device = nullptr;
}

void BeamPass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                    std::span<const BeamDrawItem> items) const
{
    if (items.empty()) {
        return;
    }
    SDL_BindGPUGraphicsPipeline(pass, m_pipeline);
    for (const BeamDrawItem& item : items) {
        VertexUniforms vu;
        vu.view_proj = view.view_proj;
        vu.apex = glm::vec4(item.apex, item.length_km);
        vu.axis = glm::vec4(item.axis, std::tan(item.half_angle_rad));
        FragmentUniforms fu;
        fu.color = glm::vec4(item.color, item.intensity);
        SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));
        SDL_DrawGPUPrimitives(pass, 6, 1, 0, 0);
    }
}

} // namespace astraxis
