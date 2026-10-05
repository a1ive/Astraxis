#include "render/dust_pass.hpp"

#include "render/shader.hpp"

#include <shaders/dust.frag.h>
#include <shaders/dust.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

namespace astraxis {

namespace {

struct VertexUniforms {
    glm::mat4 view_proj;
    glm::vec4 params;
};

} // namespace

bool DustPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;
    m_buffer.init(device, SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ);

    SDL_GPUShader* vs = create_shader(device, {.dxil = kDustVertDxil,
                                               .spirv = kDustVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_storage_buffers = 1,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kDustFragDxil,
                                               .spirv = kDustFragSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_FRAGMENT});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, fs);
        return false;
    }

    // Emission, added to what lies behind.
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
    info.depth_stencil_state.compare_op = SDL_GPU_COMPAREOP_GREATER; // reversed-Z
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = true;
    info.target_info.depth_stencil_format = format.depth;

    m_pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, fs);
    if (!m_pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Dust pipeline creation failed: %s", SDL_GetError());
        return false;
    }
    return true;
}

void DustPass::shutdown()
{
    if (!m_device) {
        return;
    }
    m_buffer.shutdown();
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_pipeline);
    m_pipeline = nullptr;
    m_device = nullptr;
}

void DustPass::upload(SDL_GPUCommandBuffer* cmd)
{
    m_uploaded = 0;
    if (m_splats.empty()) {
        return;
    }
    if (m_buffer.upload(cmd, m_splats.data(), static_cast<uint32_t>(m_splats.size() * sizeof(DustSplat)))) {
        m_uploaded = static_cast<uint32_t>(m_splats.size());
    }
}

void DustPass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                    float min_sigma_px) const
{
    if (m_uploaded == 0 || !m_buffer.buffer()) {
        return;
    }
    SDL_BindGPUGraphicsPipeline(pass, m_pipeline);
    SDL_GPUBuffer* storage = m_buffer.buffer();
    SDL_BindGPUVertexStorageBuffers(pass, 0, &storage, 1);

    VertexUniforms vu;
    vu.view_proj = view.view_proj;
    vu.params = glm::vec4(0.5f * view.viewport.y * view.focal_y, min_sigma_px, view.viewport);
    SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));
    SDL_DrawGPUPrimitives(pass, m_uploaded * 6, 1, 0, 0);
}

} // namespace astraxis
