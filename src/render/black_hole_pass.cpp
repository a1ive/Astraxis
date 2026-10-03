#include "render/black_hole_pass.hpp"

#include "core/color.hpp"
#include "render/shader.hpp"
#include "render/texture.hpp"

#include <shaders/black_hole.frag.h>
#include <shaders/black_hole_composite.frag.h>
#include <shaders/fullscreen.vert.h>

#include <SDL3/SDL_log.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace astraxis {

namespace {

constexpr int kLutSize = 256;
constexpr double kLutMinK = 1000.0;
constexpr double kLutMaxK = 40000.0; // the shader maps log(T / 1000 K) / log(40)

// Rays are traced at 1/kTraceDivisor of the scene resolution per axis; the
// disk and shadow are smooth enough, and the lensed stars are bilinearly upsampled.
constexpr uint32_t kTraceDivisor = 2;
constexpr SDL_GPUTextureFormat kTraceColorFormat = SDL_GPU_TEXTUREFORMAT_R16G16B16A16_FLOAT;
constexpr SDL_GPUTextureFormat kTraceDepthFormat = SDL_GPU_TEXTUREFORMAT_R32_FLOAT;

SDL_GPUSampler* create_sampler(SDL_GPUDevice* device, SDL_GPUFilter filter)
{
    SDL_GPUSamplerCreateInfo info = {};
    info.min_filter = filter;
    info.mag_filter = filter;
    info.mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_NEAREST;
    info.address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    info.address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    info.address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    return SDL_CreateGPUSampler(device, &info);
}

} // namespace

bool BlackHolePass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;

    SDL_GPUShader* vs = create_shader(
        device, {.dxil = kFullscreenVertDxil, .spirv = kFullscreenVertSpirv, .stage = SDL_GPU_SHADERSTAGE_VERTEX});
    SDL_GPUShader* trace_fs = create_shader(device, {.dxil = kBlackHoleFragDxil,
                                                     .spirv = kBlackHoleFragSpirv,
                                                     .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                                     .num_samplers = 2,
                                                     .num_uniform_buffers = 1});
    SDL_GPUShader* composite_fs = create_shader(
        device, {.dxil = kBlackHoleCompositeFragDxil,
                 .spirv = kBlackHoleCompositeFragSpirv,
                 .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                 .num_samplers = 2});
    if (!vs || !trace_fs || !composite_fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, trace_fs);
        SDL_ReleaseGPUShader(device, composite_fs);
        return false;
    }

    // Trace: two plain off-screen targets, no depth, no blending.
    {
        SDL_GPUColorTargetDescription targets[2] = {};
        targets[0].format = kTraceColorFormat;
        targets[1].format = kTraceDepthFormat;

        SDL_GPUGraphicsPipelineCreateInfo info = {};
        info.vertex_shader = vs;
        info.fragment_shader = trace_fs;
        info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
        info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
        info.target_info.color_target_descriptions = targets;
        info.target_info.num_color_targets = 2;
        m_trace_pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    }

    // Composite: premultiplied color over the scene; SV_Depth from the trace, so the
    // shadow and disk occlude what lies behind the hole while the lensed sky (depth 0)
    // only replaces empty background.
    {
        SDL_GPUColorTargetDescription color = {};
        color.format = format.color;
        color.blend_state.enable_blend = true;
        color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
        color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
        color.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
        color.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
        color.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
        color.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

        SDL_GPUGraphicsPipelineCreateInfo info = {};
        info.vertex_shader = vs;
        info.fragment_shader = composite_fs;
        info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
        info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
        info.multisample_state.sample_count = format.samples;
        info.depth_stencil_state.enable_depth_test = true;
        info.depth_stencil_state.enable_depth_write = true;
        info.depth_stencil_state.compare_op = SDL_GPU_COMPAREOP_GREATER_OR_EQUAL;
        info.target_info.color_target_descriptions = &color;
        info.target_info.num_color_targets = 1;
        info.target_info.has_depth_stencil_target = true;
        info.target_info.depth_stencil_format = format.depth;
        m_composite_pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    }

    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, trace_fs);
    SDL_ReleaseGPUShader(device, composite_fs);
    if (!m_trace_pipeline || !m_composite_pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Black hole pipeline creation failed: %s", SDL_GetError());
        return false;
    }

    m_linear = create_sampler(device, SDL_GPU_FILTER_LINEAR);
    m_point = create_sampler(device, SDL_GPU_FILTER_NEAREST);

    // Blackbody chromaticity (linear sRGB, max component 1) on a log temperature axis.
    std::vector<uint8_t> lut(kLutSize * 4);
    for (int i = 0; i < kLutSize; ++i) {
        const double t = kLutMinK * std::pow(kLutMaxK / kLutMinK, static_cast<double>(i) / (kLutSize - 1));
        const glm::dvec3 c = blackbody_linear_srgb(t);
        for (int k = 0; k < 3; ++k) {
            lut[i * 4 + k] = static_cast<uint8_t>(std::lround(std::clamp(c[k], 0.0, 1.0) * 255.0));
        }
        lut[i * 4 + 3] = 255;
    }
    m_blackbody = create_texture_rgba8(device, lut.data(), kLutSize, 1, false);
    return m_linear && m_point && m_blackbody;
}

void BlackHolePass::shutdown()
{
    if (!m_device) {
        return;
    }
    release_targets();
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_trace_pipeline);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_composite_pipeline);
    SDL_ReleaseGPUSampler(m_device, m_linear);
    SDL_ReleaseGPUSampler(m_device, m_point);
    SDL_ReleaseGPUTexture(m_device, m_blackbody);
    m_trace_pipeline = nullptr;
    m_composite_pipeline = nullptr;
    m_linear = nullptr;
    m_point = nullptr;
    m_blackbody = nullptr;
    m_device = nullptr;
}

bool BlackHolePass::ensure_targets(uint32_t width, uint32_t height)
{
    if (m_color && m_width == width && m_height == height) {
        return true;
    }
    release_targets();

    SDL_GPUTextureCreateInfo info = {};
    info.type = SDL_GPU_TEXTURETYPE_2D;
    info.width = width;
    info.height = height;
    info.layer_count_or_depth = 1;
    info.num_levels = 1;
    info.usage = SDL_GPU_TEXTUREUSAGE_COLOR_TARGET | SDL_GPU_TEXTUREUSAGE_SAMPLER;
    info.format = kTraceColorFormat;
    m_color = SDL_CreateGPUTexture(m_device, &info);
    info.format = kTraceDepthFormat;
    m_depth = SDL_CreateGPUTexture(m_device, &info);
    if (!m_color || !m_depth) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Black hole trace targets: %s", SDL_GetError());
        release_targets();
        return false;
    }
    m_width = width;
    m_height = height;
    return true;
}

void BlackHolePass::release_targets()
{
    SDL_ReleaseGPUTexture(m_device, m_color);
    SDL_ReleaseGPUTexture(m_device, m_depth);
    m_color = nullptr;
    m_depth = nullptr;
    m_width = 0;
    m_height = 0;
}

void BlackHolePass::trace(SDL_GPUCommandBuffer* cmd, uint32_t width, uint32_t height, const BlackHoleUniforms& u)
{
    m_traced = false;
    const uint32_t w = std::max(1u, (width + kTraceDivisor - 1) / kTraceDivisor);
    const uint32_t h = std::max(1u, (height + kTraceDivisor - 1) / kTraceDivisor);
    if (!m_sky || !ensure_targets(w, h)) {
        return;
    }

    SDL_GPUColorTargetInfo targets[2] = {};
    targets[0].texture = m_color;
    targets[0].load_op = SDL_GPU_LOADOP_CLEAR; // transparent where the shader discards
    targets[0].store_op = SDL_GPU_STOREOP_STORE;
    targets[1].texture = m_depth;
    targets[1].load_op = SDL_GPU_LOADOP_CLEAR;
    targets[1].store_op = SDL_GPU_STOREOP_STORE;
    SDL_GPURenderPass* pass = SDL_BeginGPURenderPass(cmd, targets, 2, nullptr);
    if (!pass) {
        return;
    }
    SDL_BindGPUGraphicsPipeline(pass, m_trace_pipeline);
    const SDL_GPUTextureSamplerBinding bindings[2] = {
        {.texture = m_sky, .sampler = m_linear},
        {.texture = m_blackbody, .sampler = m_linear},
    };
    SDL_BindGPUFragmentSamplers(pass, 0, bindings, 2);
    SDL_PushGPUFragmentUniformData(cmd, 0, &u, sizeof(u));
    SDL_DrawGPUPrimitives(pass, 3, 1, 0, 0);
    SDL_EndGPURenderPass(pass);
    m_traced = true;
}

void BlackHolePass::composite(SDL_GPURenderPass* pass)
{
    if (!m_traced) {
        return;
    }
    SDL_BindGPUGraphicsPipeline(pass, m_composite_pipeline);
    const SDL_GPUTextureSamplerBinding bindings[2] = {
        {.texture = m_color, .sampler = m_linear},
        {.texture = m_depth, .sampler = m_point},
    };
    SDL_BindGPUFragmentSamplers(pass, 0, bindings, 2);
    SDL_DrawGPUPrimitives(pass, 3, 1, 0, 0);
    m_traced = false; // one composite per trace
}

} // namespace astraxis
