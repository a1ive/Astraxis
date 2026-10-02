#include "render/post_process.hpp"

#include "render/shader.hpp"

#include <shaders/bloom_downsample.frag.h>
#include <shaders/bloom_upsample.frag.h>
#include <shaders/composite.frag.h>
#include <shaders/fullscreen.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/vec4.hpp>

#include <algorithm>

namespace astraxis {

namespace {

constexpr int kMaxBloomLevels = 7;
// The Karis average suppresses small very bright features, which here is exactly
// the sun glare we want. The scene has no specular fireflies, so it stays off.
constexpr bool kKarisAverage = false;
constexpr uint32_t kMinBloomSize = 8;

SDL_GPUGraphicsPipeline* create_fullscreen_pipeline(SDL_GPUDevice* device, std::span<const unsigned char> frag,
                                                    uint32_t num_samplers, SDL_GPUTextureFormat format,
                                                    bool additive)
{
    SDL_GPUShader* vs = create_shader(device, {.dxil = kFullscreenVertDxil, .stage = SDL_GPU_SHADERSTAGE_VERTEX});
    SDL_GPUShader* fs = create_shader(device, {.dxil = frag,
                                               .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                               .num_samplers = num_samplers,
                                               .num_uniform_buffers = 1});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, fs);
        return nullptr;
    }

    SDL_GPUColorTargetDescription color = {};
    color.format = format;
    if (additive) {
        color.blend_state.enable_blend = true;
        color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
        color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
        color.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
        color.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ZERO;
        color.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
        color.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;
    }

    SDL_GPUGraphicsPipelineCreateInfo info = {};
    info.vertex_shader = vs;
    info.fragment_shader = fs;
    info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
    info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;

    SDL_GPUGraphicsPipeline* pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, fs);
    if (!pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Post pipeline creation failed: %s", SDL_GetError());
    }
    return pipeline;
}

} // namespace

bool PostProcess::init(SDL_GPUDevice* device, SDL_GPUTextureFormat hdr_format, SDL_GPUTextureFormat output_format)
{
    m_device = device;
    m_hdr_format = hdr_format;

    m_downsample = create_fullscreen_pipeline(device, kBloomDownsampleFragDxil, 1, hdr_format, false);
    m_upsample = create_fullscreen_pipeline(device, kBloomUpsampleFragDxil, 1, hdr_format, true);
    m_composite = create_fullscreen_pipeline(device, kCompositeFragDxil, 2, output_format, false);

    SDL_GPUSamplerCreateInfo sampler = {};
    sampler.min_filter = SDL_GPU_FILTER_LINEAR;
    sampler.mag_filter = SDL_GPU_FILTER_LINEAR;
    sampler.mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_NEAREST;
    sampler.address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    sampler.address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    sampler.address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    m_linear_clamp = SDL_CreateGPUSampler(device, &sampler);

    return m_downsample && m_upsample && m_composite && m_linear_clamp;
}

void PostProcess::shutdown()
{
    if (!m_device) {
        return;
    }
    release_chain();
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_downsample);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_upsample);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_composite);
    SDL_ReleaseGPUSampler(m_device, m_linear_clamp);
    m_downsample = m_upsample = m_composite = nullptr;
    m_linear_clamp = nullptr;
    m_device = nullptr;
}

bool PostProcess::ensure_chain(uint32_t width, uint32_t height)
{
    if (!m_levels.empty() && width == m_chain_width && height == m_chain_height) {
        return true;
    }
    release_chain();

    uint32_t w = std::max(width / 2, 1u);
    uint32_t h = std::max(height / 2, 1u);
    for (int i = 0; i < kMaxBloomLevels && w >= kMinBloomSize && h >= kMinBloomSize; ++i) {
        SDL_GPUTextureCreateInfo info = {};
        info.type = SDL_GPU_TEXTURETYPE_2D;
        info.format = m_hdr_format;
        info.usage = SDL_GPU_TEXTUREUSAGE_COLOR_TARGET | SDL_GPU_TEXTUREUSAGE_SAMPLER;
        info.width = w;
        info.height = h;
        info.layer_count_or_depth = 1;
        info.num_levels = 1;
        info.sample_count = SDL_GPU_SAMPLECOUNT_1;
        SDL_GPUTexture* texture = SDL_CreateGPUTexture(m_device, &info);
        if (!texture) {
            SDL_LogError(SDL_LOG_CATEGORY_GPU, "Bloom texture creation failed: %s", SDL_GetError());
            release_chain();
            return false;
        }
        m_levels.push_back({texture, w, h});
        w = std::max(w / 2, 1u);
        h = std::max(h / 2, 1u);
    }

    m_chain_width = width;
    m_chain_height = height;
    return !m_levels.empty();
}

void PostProcess::release_chain()
{
    for (Level& level : m_levels) {
        SDL_ReleaseGPUTexture(m_device, level.texture);
    }
    m_levels.clear();
    m_chain_width = m_chain_height = 0;
}

void PostProcess::fullscreen_pass(SDL_GPUCommandBuffer* cmd, SDL_GPUTexture* target, SDL_GPULoadOp load,
                                  SDL_GPUGraphicsPipeline* pipeline, SDL_GPUTexture* source, const void* uniforms,
                                  uint32_t uniforms_size)
{
    SDL_GPUColorTargetInfo color = {};
    color.texture = target;
    color.load_op = load;
    color.store_op = SDL_GPU_STOREOP_STORE;

    SDL_GPURenderPass* pass = SDL_BeginGPURenderPass(cmd, &color, 1, nullptr);
    SDL_BindGPUGraphicsPipeline(pass, pipeline);
    SDL_GPUTextureSamplerBinding binding = {.texture = source, .sampler = m_linear_clamp};
    SDL_BindGPUFragmentSamplers(pass, 0, &binding, 1);
    SDL_PushGPUFragmentUniformData(cmd, 0, uniforms, uniforms_size);
    SDL_DrawGPUPrimitives(pass, 3, 1, 0, 0);
    SDL_EndGPURenderPass(pass);
}

void PostProcess::run(SDL_GPUCommandBuffer* cmd, SDL_GPUTexture* hdr, uint32_t width, uint32_t height,
                      SDL_GPUTexture* output, const PostSettings& settings)
{
    const bool bloom = settings.bloom_strength > 0.0f && ensure_chain(width, height);

    if (bloom) {
        // Downsample chain: hdr -> level 0 -> level 1 -> ...
        uint32_t src_w = width;
        uint32_t src_h = height;
        SDL_GPUTexture* src = hdr;
        for (size_t i = 0; i < m_levels.size(); ++i) {
            const glm::vec4 params(1.0f / static_cast<float>(src_w), 1.0f / static_cast<float>(src_h),
                                   (i == 0 && kKarisAverage) ? 1.0f : 0.0f,
                                   i == 0 ? settings.bloom_threshold : 0.0f);
            fullscreen_pass(cmd, m_levels[i].texture, SDL_GPU_LOADOP_DONT_CARE, m_downsample, src, &params,
                            sizeof(params));
            src = m_levels[i].texture;
            src_w = m_levels[i].width;
            src_h = m_levels[i].height;
        }
        // Upsample chain: accumulate each level into the next larger one.
        for (size_t i = m_levels.size() - 1; i > 0; --i) {
            const glm::vec4 params(1.0f / static_cast<float>(m_levels[i].width),
                                   1.0f / static_cast<float>(m_levels[i].height), 0.0f, 0.0f);
            fullscreen_pass(cmd, m_levels[i - 1].texture, SDL_GPU_LOADOP_LOAD, m_upsample, m_levels[i].texture,
                            &params, sizeof(params));
        }
    }

    // Composite to the output.
    SDL_GPUColorTargetInfo color = {};
    color.texture = output;
    color.load_op = SDL_GPU_LOADOP_DONT_CARE;
    color.store_op = SDL_GPU_STOREOP_STORE;

    SDL_GPURenderPass* pass = SDL_BeginGPURenderPass(cmd, &color, 1, nullptr);
    SDL_BindGPUGraphicsPipeline(pass, m_composite);
    SDL_GPUTextureSamplerBinding bindings[2] = {
        {.texture = hdr, .sampler = m_linear_clamp},
        {.texture = bloom ? m_levels[0].texture : hdr, .sampler = m_linear_clamp},
    };
    SDL_BindGPUFragmentSamplers(pass, 0, bindings, 2);
    const glm::vec4 params(settings.exposure, bloom ? settings.bloom_strength : 0.0f, 0.0f, 0.0f);
    SDL_PushGPUFragmentUniformData(cmd, 0, &params, sizeof(params));
    SDL_DrawGPUPrimitives(pass, 3, 1, 0, 0);
    SDL_EndGPURenderPass(pass);
}

} // namespace astraxis
