#include "render/pinwheel_pass.hpp"

#include "render/gpu_buffer.hpp"
#include "render/shader.hpp"

#include <shaders/dust.frag.h>
#include <shaders/fullscreen.vert.h>
#include <shaders/pinwheel.vert.h>
#include <shaders/pinwheel_composite.frag.h>

#include <SDL3/SDL_log.h>
#include <glm/mat4x4.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr uint32_t kDivisor = 2; // off-screen resolution: 1 / kDivisor per axis
constexpr SDL_GPUTextureFormat kDustFormat = SDL_GPU_TEXTUREFORMAT_R16G16B16A16_FLOAT;

struct VertexUniforms {
    glm::mat4 view_proj;
    glm::mat4 view;
    glm::vec4 center;
    glm::vec4 display_x;
    glm::vec4 display_y;
    glm::vec4 display_z;
    glm::vec4 cycle;
    glm::vec4 shape;
    glm::vec4 light;
    glm::vec4 screen;
    std::array<glm::vec4, kPinwheelColorCount> colors;
    std::array<float, kPinwheelRadiusCount> radius;
};

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

bool PinwheelPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;
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

bool PinwheelPass::set_format(const SceneTargetFormat& format)
{
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_splat_pipeline);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_composite_pipeline);
    m_splat_pipeline = nullptr;
    m_composite_pipeline = nullptr;
    return create_pipelines(format);
}

bool PinwheelPass::create_pipelines(const SceneTargetFormat& format)
{
    SDL_GPUDevice* device = m_device;
    SDL_GPUShader* splat_vs = create_shader(device, {.dxil = kPinwheelVertDxil,
                                                     .spirv = kPinwheelVertSpirv,
                                                     .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                                     .num_storage_buffers = 1,
                                                     .num_uniform_buffers = 1});
    SDL_GPUShader* splat_fs = create_shader(device, {.dxil = kDustFragDxil,
                                                     .spirv = kDustFragSpirv,
                                                     .stage = SDL_GPU_SHADERSTAGE_FRAGMENT});
    SDL_GPUShader* full_vs = create_shader(device, {.dxil = kFullscreenVertDxil,
                                                    .spirv = kFullscreenVertSpirv,
                                                    .stage = SDL_GPU_SHADERSTAGE_VERTEX});
    SDL_GPUShader* composite_fs = create_shader(device, {.dxil = kPinwheelCompositeFragDxil,
                                                         .spirv = kPinwheelCompositeFragSpirv,
                                                         .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                                         .num_samplers = 1});
    if (splat_vs && splat_fs && full_vs && composite_fs) {
        // Splats: into the off-screen target, one sample, no depth (the dust
        // is thin; the stars are drawn over it).
        const SDL_GPUColorTargetDescription dust = additive(kDustFormat);
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
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Pinwheel pipeline creation failed: %s", SDL_GetError());
        return false;
    }
    return true;
}

void PinwheelPass::shutdown()
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

void PinwheelPass::clear()
{
    SDL_ReleaseGPUBuffer(m_device, m_buffer);
    m_buffer = nullptr;
    m_count = 0;
}

bool PinwheelPass::set_grains(const PinwheelGpuGrain* grains, uint32_t count)
{
    clear();
    if (count == 0) {
        return true;
    }
    m_buffer = create_static_buffer(m_device, SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ, grains, count * static_cast<uint32_t>(sizeof(PinwheelGpuGrain)));
    if (!m_buffer) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Pinwheel buffer creation failed: %s", SDL_GetError());
        return false;
    }
    m_count = count;
    return true;
}

bool PinwheelPass::ensure_targets(PinwheelTargets& targets, uint32_t width, uint32_t height)
{
    if (targets.color && targets.width == width && targets.height == height) {
        return true;
    }
    SDL_ReleaseGPUTexture(m_device, targets.color);
    targets = {};
    SDL_GPUTextureCreateInfo info = {};
    info.type = SDL_GPU_TEXTURETYPE_2D;
    info.format = kDustFormat;
    info.width = width;
    info.height = height;
    info.layer_count_or_depth = 1;
    info.num_levels = 1;
    info.usage = SDL_GPU_TEXTUREUSAGE_COLOR_TARGET | SDL_GPU_TEXTUREUSAGE_SAMPLER;
    targets.color = SDL_CreateGPUTexture(m_device, &info);
    if (!targets.color) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Pinwheel target: %s", SDL_GetError());
        return false;
    }
    targets.width = width;
    targets.height = height;
    return true;
}

void PinwheelPass::render(SDL_GPUCommandBuffer* cmd, PinwheelTargets& targets, const CameraView& view,
                          const PinwheelDrawItem& item, uint32_t width, uint32_t height)
{
    targets.drawn = false;
    const uint32_t w = std::max(1u, (width + kDivisor - 1) / kDivisor);
    const uint32_t h = std::max(1u, (height + kDivisor - 1) / kDivisor);
    if (m_count == 0 || !m_buffer || !ensure_targets(targets, w, h)) {
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

    // glm is column-major: row r of the rotation is (m[0][r], m[1][r], m[2][r]).
    const glm::mat3& m = item.display_from_icrf;
    VertexUniforms vu;
    vu.view_proj = view.view_proj;
    vu.view = view.view;
    vu.center = glm::vec4(item.center, 0.0f);
    vu.display_x = glm::vec4(m[0][0], m[1][0], m[2][0], 0.0f);
    vu.display_y = glm::vec4(m[0][1], m[1][1], m[2][1], 0.0f);
    vu.display_z = glm::vec4(m[0][2], m[1][2], m[2][2], 0.0f);
    vu.cycle = glm::vec4(item.cycle, item.turns, item.reference_km, item.nucleation_km);
    vu.shape = glm::vec4(item.sigma_along, item.sigma_across, item.wall_sigma_km, item.min_sigma_px);
    const float log_span = std::log(item.radius.back() / item.nucleation_km);
    vu.light = glm::vec4(item.radiance, item.gamma, log_span > 0.0f ? 1.0f / log_span : 0.0f, item.cull_radiance);
    // Splat sizes in the off-screen target's pixels (the projection itself is the scene's).
    vu.screen = glm::vec4(static_cast<float>(w), static_cast<float>(h), 0.5f * static_cast<float>(h) * view.focal_y,
                          0.0f);
    vu.colors = item.colors;
    vu.radius = item.radius;
    SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));
    SDL_DrawGPUPrimitives(pass, m_count * 6u, 1, 0, 0);
    SDL_EndGPURenderPass(pass);
    targets.drawn = true;
}

void PinwheelPass::composite(SDL_GPURenderPass* pass, PinwheelTargets& targets)
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
