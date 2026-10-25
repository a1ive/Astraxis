#include "render/corona_pass.hpp"

#include "render/shader.hpp"
#include <shaders/corona.frag.h>
#include <shaders/corona_composite.frag.h>
#include <shaders/fullscreen.vert.h>
#include <SDL3/SDL_log.h>

#include <algorithm>
#include <cmath>

namespace astraxis {

bool CoronaPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;
    SDL_GPUSamplerCreateInfo sampler = {};
    sampler.min_filter = sampler.mag_filter = SDL_GPU_FILTER_LINEAR;
    sampler.address_mode_u = sampler.address_mode_v = sampler.address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    m_linear = SDL_CreateGPUSampler(device, &sampler);
    return m_linear && create_pipelines(format);
}

bool CoronaPass::set_format(const SceneTargetFormat& format)
{
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_trace);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_composite);
    m_trace = m_composite = nullptr;
    return create_pipelines(format);
}

bool CoronaPass::create_pipelines(const SceneTargetFormat& format)
{
    SDL_GPUShader* vs = create_shader(m_device, {.dxil = kFullscreenVertDxil, .spirv = kFullscreenVertSpirv,
                                                 .stage = SDL_GPU_SHADERSTAGE_VERTEX});
    SDL_GPUShader* trace = create_shader(m_device, {.dxil = kCoronaFragDxil, .spirv = kCoronaFragSpirv,
                                                    .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                                    .num_uniform_buffers = 1});
    SDL_GPUShader* composite = create_shader(m_device, {.dxil = kCoronaCompositeFragDxil,
                                                        .spirv = kCoronaCompositeFragSpirv,
                                                        .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                                        .num_samplers = 1, .num_uniform_buffers = 1});
    if (vs && trace && composite) {
        SDL_GPUColorTargetDescription color = {};
        color.format = SDL_GPU_TEXTUREFORMAT_R16G16B16A16_FLOAT;
        SDL_GPUGraphicsPipelineCreateInfo info = {};
        info.vertex_shader = vs;
        info.fragment_shader = trace;
        info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
        info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
        info.target_info.color_target_descriptions = &color;
        info.target_info.num_color_targets = 1;
        m_trace = SDL_CreateGPUGraphicsPipeline(m_device, &info);

        color.format = format.color;
        color.blend_state.enable_blend = true;
        color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
        color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
        color.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
        color.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ZERO;
        color.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
        color.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;
        info.fragment_shader = composite;
        info.multisample_state.sample_count = format.samples;
        info.target_info.has_depth_stencil_target = true;
        info.target_info.depth_stencil_format = format.depth;
        m_composite = SDL_CreateGPUGraphicsPipeline(m_device, &info);
    }
    SDL_ReleaseGPUShader(m_device, vs);
    SDL_ReleaseGPUShader(m_device, trace);
    SDL_ReleaseGPUShader(m_device, composite);
    if (!m_trace || !m_composite) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Corona pipeline: %s", SDL_GetError());
        return false;
    }
    return true;
}

void CoronaPass::shutdown()
{
    if (m_device) {
        SDL_ReleaseGPUGraphicsPipeline(m_device, m_trace);
        SDL_ReleaseGPUGraphicsPipeline(m_device, m_composite);
        SDL_ReleaseGPUSampler(m_device, m_linear);
    }
    m_trace = m_composite = nullptr;
    m_linear = nullptr;
    m_device = nullptr;
}

void CoronaPass::render(SDL_GPUCommandBuffer* cmd, CoronaTargets& targets, const CoronaUniforms& uniforms,
                        uint32_t width, uint32_t height)
{
    targets.drawn = false;
    if (uniforms.params.x <= 0.0f) return;
    const uint32_t w = std::max(1u, (width + 1) / 2);
    const uint32_t h = std::max(1u, (height + 1) / 2);
    if (!targets.color || targets.width != w || targets.height != h) {
        SDL_ReleaseGPUTexture(m_device, targets.color);
        targets = {};
        SDL_GPUTextureCreateInfo info = {};
        info.type = SDL_GPU_TEXTURETYPE_2D;
        info.format = SDL_GPU_TEXTUREFORMAT_R16G16B16A16_FLOAT;
        info.width = w;
        info.height = h;
        info.layer_count_or_depth = 1;
        info.num_levels = 1;
        info.usage = SDL_GPU_TEXTUREUSAGE_COLOR_TARGET | SDL_GPU_TEXTUREUSAGE_SAMPLER;
        targets.color = SDL_CreateGPUTexture(m_device, &info);
        if (!targets.color) {
            SDL_LogError(SDL_LOG_CATEGORY_GPU, "Corona target: %s", SDL_GetError());
            return;
        }
        targets.width = w;
        targets.height = h;
    }
    static const auto radial = [] {
        std::array<glm::vec4, kCoronaRadialCount> result{};
        for (int i = 0; i < kCoronaRadialCount; ++i) {
            const double t = static_cast<double>(i) / (kCoronaRadialCount - 1);
            const glm::dvec2 c = corona_scattering(std::exp(t * t * std::log(kCoronaOuterRadius)));
            result[static_cast<size_t>(i)] = glm::vec4(c.x, c.y, 0.0, 0.0);
        }
        return result;
    }();
    targets.uniforms = uniforms;
    targets.uniforms.radial = radial;
    SDL_GPUColorTargetInfo color = {};
    color.texture = targets.color;
    color.load_op = SDL_GPU_LOADOP_CLEAR;
    color.store_op = SDL_GPU_STOREOP_STORE;
    color.cycle = true;
    SDL_GPURenderPass* pass = SDL_BeginGPURenderPass(cmd, &color, 1, nullptr);
    if (!pass) return;
    SDL_BindGPUGraphicsPipeline(pass, m_trace);
    SDL_PushGPUFragmentUniformData(cmd, 0, &targets.uniforms, sizeof(CoronaUniforms));
    SDL_DrawGPUPrimitives(pass, 3, 1, 0, 0);
    SDL_EndGPURenderPass(pass);
    targets.drawn = true;
}

void CoronaPass::composite(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, CoronaTargets& targets)
{
    if (!targets.drawn) return;
    SDL_BindGPUGraphicsPipeline(pass, m_composite);
    SDL_PushGPUFragmentUniformData(cmd, 0, &targets.uniforms, sizeof(CoronaUniforms));
    const SDL_GPUTextureSamplerBinding binding = {.texture = targets.color, .sampler = m_linear};
    SDL_BindGPUFragmentSamplers(pass, 0, &binding, 1);
    SDL_DrawGPUPrimitives(pass, 3, 1, 0, 0);
    targets.drawn = false;
}

} // namespace astraxis
