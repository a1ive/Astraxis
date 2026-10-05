#include "render/scene_targets.hpp"

#include <SDL3/SDL_log.h>

namespace astraxis {

bool SceneTargets::ensure(SDL_GPUDevice* device, const SceneTargetFormat& format, uint32_t width, uint32_t height)
{
    if (m_depth && device == m_device && width == m_width && height == m_height) {
        return true;
    }
    release_scene();
    m_device = device;

    SDL_GPUTextureCreateInfo info = {};
    info.type = SDL_GPU_TEXTURETYPE_2D;
    info.width = width;
    info.height = height;
    info.layer_count_or_depth = 1;
    info.num_levels = 1;
    info.sample_count = format.samples;

    if (format.samples != SDL_GPU_SAMPLECOUNT_1) {
        info.format = format.color;
        info.usage = SDL_GPU_TEXTUREUSAGE_COLOR_TARGET;
        m_msaa_color = SDL_CreateGPUTexture(device, &info);
        if (!m_msaa_color) {
            SDL_LogError(SDL_LOG_CATEGORY_GPU, "MSAA color target creation failed: %s", SDL_GetError());
            return false;
        }
    }

    info.sample_count = SDL_GPU_SAMPLECOUNT_1;
    info.format = format.color;
    info.usage = SDL_GPU_TEXTUREUSAGE_COLOR_TARGET | SDL_GPU_TEXTUREUSAGE_SAMPLER;
    m_hdr = SDL_CreateGPUTexture(device, &info);
    if (!m_hdr) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "HDR target creation failed: %s", SDL_GetError());
        release_scene();
        return false;
    }

    info.sample_count = format.samples;
    info.format = format.depth;
    info.usage = SDL_GPU_TEXTUREUSAGE_DEPTH_STENCIL_TARGET;
    m_depth = SDL_CreateGPUTexture(device, &info);
    if (!m_depth) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Depth target creation failed: %s", SDL_GetError());
        release_scene();
        return false;
    }

    m_width = width;
    m_height = height;
    return true;
}

void SceneTargets::release()
{
    release_scene();
    if (!m_device) {
        return;
    }
    for (BloomChain::Level& level : bloom.levels) {
        SDL_ReleaseGPUTexture(m_device, level.texture);
    }
    bloom = {};
    SDL_ReleaseGPUTexture(m_device, black_hole.color);
    SDL_ReleaseGPUTexture(m_device, black_hole.depth);
    black_hole = {};
    m_device = nullptr;
}

SDL_GPURenderPass* SceneTargets::begin_scene_pass(SDL_GPUCommandBuffer* cmd, SDL_FColor clear_color) const
{
    SDL_GPUColorTargetInfo color = {};
    color.clear_color = clear_color;
    color.load_op = SDL_GPU_LOADOP_CLEAR;
    color.cycle = true;
    if (m_msaa_color) {
        color.texture = m_msaa_color;
        color.resolve_texture = m_hdr;
        color.store_op = SDL_GPU_STOREOP_RESOLVE;
        color.cycle_resolve_texture = true;
    } else {
        color.texture = m_hdr;
        color.store_op = SDL_GPU_STOREOP_STORE;
    }

    SDL_GPUDepthStencilTargetInfo depth = {};
    depth.texture = m_depth;
    depth.clear_depth = 0.0f; // reversed-Z: far = 0
    depth.load_op = SDL_GPU_LOADOP_CLEAR;
    depth.store_op = SDL_GPU_STOREOP_DONT_CARE;
    depth.stencil_load_op = SDL_GPU_LOADOP_DONT_CARE;
    depth.stencil_store_op = SDL_GPU_STOREOP_DONT_CARE;
    depth.cycle = true;

    return SDL_BeginGPURenderPass(cmd, &color, 1, &depth);
}

void SceneTargets::release_scene()
{
    if (m_device) {
        SDL_ReleaseGPUTexture(m_device, m_msaa_color);
        SDL_ReleaseGPUTexture(m_device, m_hdr);
        SDL_ReleaseGPUTexture(m_device, m_depth);
    }
    m_msaa_color = nullptr;
    m_hdr = nullptr;
    m_depth = nullptr;
    m_width = m_height = 0;
}

} // namespace astraxis
