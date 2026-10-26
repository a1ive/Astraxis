#include "render/splat_pass.hpp"

#include "render/gpu_buffer.hpp"

#include <shaders/dust.frag.h>
#include <shaders/fullscreen.vert.h>
#include <shaders/splat_composite.frag.h>

#include <SDL3/SDL_log.h>

#include <algorithm>

namespace astraxis {

namespace {

constexpr uint32_t kDivisor = 2; // off-screen resolution: 1 / kDivisor per axis
constexpr SDL_GPUTextureFormat kSplatFormat = SDL_GPU_TEXTUREFORMAT_R16G16B16A16_FLOAT;

// Additive (emission of optically thin dust).
SDL_GPUColorTargetDescription additive(SDL_GPUTextureFormat format)
{
    SDL_GPUColorTargetDescription color = {};
    color.format = format;
    color.blend_state.enable_blend = true;
    color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ZERO;
    color.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;
    return color;
}

} // namespace

bool SplatPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format, const ShaderDesc& vertex_shader)
{
    m_device = device;
    m_vertex_shader = vertex_shader;
    SDL_GPUSamplerCreateInfo info = {};
    info.min_filter = SDL_GPU_FILTER_LINEAR;
    info.mag_filter = SDL_GPU_FILTER_LINEAR;
    info.mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_NEAREST;
    info.address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    info.address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    info.address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    m_linear = SDL_CreateGPUSampler(device, &info);
    return m_linear && create_pipelines(format);
}

bool SplatPass::set_format(const SceneTargetFormat& format)
{
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_splat_pipeline);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_composite_pipeline);
    m_splat_pipeline = nullptr;
    m_composite_pipeline = nullptr;
    return create_pipelines(format);
}

bool SplatPass::create_pipelines(const SceneTargetFormat& format)
{
    SDL_GPUDevice* device = m_device;
    SDL_GPUShader* splat_vs = create_shader(device, m_vertex_shader);
    SDL_GPUShader* splat_fs = create_shader(device, {.dxil = kDustFragDxil,
                                                     .spirv = kDustFragSpirv,
                                                     .stage = SDL_GPU_SHADERSTAGE_FRAGMENT});
    SDL_GPUShader* full_vs = create_shader(device, {.dxil = kFullscreenVertDxil,
                                                    .spirv = kFullscreenVertSpirv,
                                                    .stage = SDL_GPU_SHADERSTAGE_VERTEX});
    SDL_GPUShader* composite_fs = create_shader(device, {.dxil = kSplatCompositeFragDxil,
                                                         .spirv = kSplatCompositeFragSpirv,
                                                         .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                                         .num_samplers = 1});
    if (splat_vs && splat_fs && full_vs && composite_fs) {
        // Splats: into the off-screen target, one sample, no depth (the dust
        // is thin; the stars are drawn over it).
        const SDL_GPUColorTargetDescription dust = additive(kSplatFormat);
        SDL_GPUGraphicsPipelineCreateInfo info = {};
        info.vertex_shader = splat_vs;
        info.fragment_shader = splat_fs;
        info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
        info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
        info.target_info.color_target_descriptions = &dust;
        info.target_info.num_color_targets = 1;
        m_splat_pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);

        // Composite: a fullscreen triangle in the scene pass, per pixel.
        const SDL_GPUColorTargetDescription scene = additive(format.color);
        info = {};
        info.vertex_shader = full_vs;
        info.fragment_shader = composite_fs;
        info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
        info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
        info.multisample_state.sample_count = format.samples;
        info.target_info.color_target_descriptions = &scene;
        info.target_info.num_color_targets = 1;
        info.target_info.has_depth_stencil_target = true;
        info.target_info.depth_stencil_format = format.depth;
        m_composite_pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    }
    SDL_ReleaseGPUShader(device, splat_vs);
    SDL_ReleaseGPUShader(device, splat_fs);
    SDL_ReleaseGPUShader(device, full_vs);
    SDL_ReleaseGPUShader(device, composite_fs);
    if (!m_splat_pipeline || !m_composite_pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Splat pipeline creation failed: %s", SDL_GetError());
        return false;
    }
    return true;
}

void SplatPass::shutdown()
{
    if (!m_device) {
        return;
    }
    clear();
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_splat_pipeline);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_composite_pipeline);
    SDL_ReleaseGPUSampler(m_device, m_linear);
    m_splat_pipeline = nullptr;
    m_composite_pipeline = nullptr;
    m_linear = nullptr;
    m_device = nullptr;
}

void SplatPass::clear()
{
    SDL_ReleaseGPUBuffer(m_device, m_buffer);
    m_buffer = nullptr;
    m_count = 0;
}

bool SplatPass::set_points(const void* data, uint32_t bytes, uint32_t count)
{
    clear();
    if (count == 0) {
        return true;
    }
    m_buffer = create_static_buffer(m_device, SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ, data, bytes);
    if (!m_buffer) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Splat buffer creation failed: %s", SDL_GetError());
        return false;
    }
    m_count = count;
    return true;
}

uint32_t SplatPass::target_size(uint32_t scene_size)
{
    return std::max(1u, (scene_size + kDivisor - 1) / kDivisor);
}

bool SplatPass::ensure_targets(SplatTargets& targets, uint32_t width, uint32_t height)
{
    if (targets.color && targets.width == width && targets.height == height) {
        return true;
    }
    SDL_ReleaseGPUTexture(m_device, targets.color);
    targets = {};
    SDL_GPUTextureCreateInfo info = {};
    info.type = SDL_GPU_TEXTURETYPE_2D;
    info.format = kSplatFormat;
    info.width = width;
    info.height = height;
    info.layer_count_or_depth = 1;
    info.num_levels = 1;
    info.usage = SDL_GPU_TEXTUREUSAGE_COLOR_TARGET | SDL_GPU_TEXTUREUSAGE_SAMPLER;
    targets.color = SDL_CreateGPUTexture(m_device, &info);
    if (!targets.color) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Splat target: %s", SDL_GetError());
        return false;
    }
    targets.width = width;
    targets.height = height;
    return true;
}

void SplatPass::render(SDL_GPUCommandBuffer* cmd, SplatTargets& targets, const void* uniforms, uint32_t size,
                       uint32_t width, uint32_t height)
{
    targets.drawn = false;
    if (m_count == 0 || !m_buffer || !ensure_targets(targets, target_size(width), target_size(height))) {
        return;
    }
    SDL_GPUColorTargetInfo color = {};
    color.texture = targets.color;
    color.load_op = SDL_GPU_LOADOP_CLEAR;
    color.store_op = SDL_GPU_STOREOP_STORE;
    color.cycle = true;
    SDL_GPURenderPass* pass = SDL_BeginGPURenderPass(cmd, &color, 1, nullptr);
    if (!pass) {
        return;
    }
    SDL_BindGPUGraphicsPipeline(pass, m_splat_pipeline);
    SDL_GPUBuffer* storage = m_buffer;
    SDL_BindGPUVertexStorageBuffers(pass, 0, &storage, 1);
    SDL_PushGPUVertexUniformData(cmd, 0, uniforms, size);
    SDL_DrawGPUPrimitives(pass, m_count * 6u, 1, 0, 0);
    SDL_EndGPURenderPass(pass);
    targets.drawn = true;
}

void SplatPass::composite(SDL_GPURenderPass* pass, SplatTargets& targets)
{
    if (!targets.drawn) {
        return;
    }
    SDL_BindGPUGraphicsPipeline(pass, m_composite_pipeline);
    const SDL_GPUTextureSamplerBinding binding = {.texture = targets.color, .sampler = m_linear};
    SDL_BindGPUFragmentSamplers(pass, 0, &binding, 1);
    SDL_DrawGPUPrimitives(pass, 3, 1, 0, 0);
    targets.drawn = false;
}

} // namespace astraxis
