#include "render/renderer.hpp"

#include <SDL3/SDL_log.h>
#include <SDL3/SDL_video.h>

namespace astraxis {

bool Renderer::init(SDL_Window* window, bool debug)
{
    m_window = window;

    // DXIL for our shaders; DXBC because the ImGui backend ships DXBC for D3D12.
    m_device = SDL_CreateGPUDevice(SDL_GPU_SHADERFORMAT_DXIL | SDL_GPU_SHADERFORMAT_DXBC, debug, nullptr);
    if (!m_device) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_CreateGPUDevice failed: %s", SDL_GetError());
        return false;
    }

    if (!SDL_ClaimWindowForGPUDevice(m_device, m_window)) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_ClaimWindowForGPUDevice failed: %s", SDL_GetError());
        return false;
    }

    SDL_SetGPUSwapchainParameters(m_device, m_window, SDL_GPU_SWAPCHAINCOMPOSITION_SDR, SDL_GPU_PRESENTMODE_VSYNC);

    m_scene_format.color = SDL_GPU_TEXTUREFORMAT_R16G16B16A16_FLOAT;
    m_scene_format.depth = SDL_GPUTextureSupportsFormat(m_device, SDL_GPU_TEXTUREFORMAT_D32_FLOAT,
                                                        SDL_GPU_TEXTURETYPE_2D,
                                                        SDL_GPU_TEXTUREUSAGE_DEPTH_STENCIL_TARGET)
                               ? SDL_GPU_TEXTUREFORMAT_D32_FLOAT
                               : SDL_GPU_TEXTUREFORMAT_D24_UNORM;
    const bool msaa4 = SDL_GPUTextureSupportsSampleCount(m_device, m_scene_format.color, SDL_GPU_SAMPLECOUNT_4) &&
                       SDL_GPUTextureSupportsSampleCount(m_device, m_scene_format.depth, SDL_GPU_SAMPLECOUNT_4);
    m_scene_format.samples = msaa4 ? SDL_GPU_SAMPLECOUNT_4 : SDL_GPU_SAMPLECOUNT_1;

    SDL_LogInfo(SDL_LOG_CATEGORY_GPU, "GPU driver: %s, MSAA %s", driver_name(), msaa4 ? "4x" : "off");
    return true;
}

void Renderer::shutdown()
{
    if (m_device) {
        SDL_WaitForGPUIdle(m_device);
        release_targets();
        if (m_window) {
            SDL_ReleaseWindowFromGPUDevice(m_device, m_window);
        }
        SDL_DestroyGPUDevice(m_device);
        m_device = nullptr;
    }
    m_window = nullptr;
}

SDL_GPUTextureFormat Renderer::swapchain_format() const
{
    return SDL_GetGPUSwapchainTextureFormat(m_device, m_window);
}

const char* Renderer::driver_name() const
{
    return m_device ? SDL_GetGPUDeviceDriver(m_device) : "none";
}

bool Renderer::begin_frame(Frame& frame)
{
    frame = {};
    frame.cmd = SDL_AcquireGPUCommandBuffer(m_device);
    if (!frame.cmd) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_AcquireGPUCommandBuffer failed: %s", SDL_GetError());
        return false;
    }

    if (!SDL_WaitAndAcquireGPUSwapchainTexture(frame.cmd, m_window, &frame.swapchain, &frame.width, &frame.height)) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_WaitAndAcquireGPUSwapchainTexture failed: %s", SDL_GetError());
        frame.swapchain = nullptr;
    }
    if (frame.swapchain && !ensure_targets(frame.width, frame.height)) {
        frame.swapchain = nullptr;
    }
    return true;
}

void Renderer::end_frame(Frame& frame)
{
    if (frame.cmd) {
        SDL_SubmitGPUCommandBuffer(frame.cmd);
    }
    frame = {};
}

SDL_GPURenderPass* Renderer::begin_scene_pass(const Frame& frame, SDL_FColor clear_color)
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

    return SDL_BeginGPURenderPass(frame.cmd, &color, 1, &depth);
}

SDL_GPURenderPass* Renderer::begin_overlay_pass(const Frame& frame)
{
    SDL_GPUColorTargetInfo color = {};
    color.texture = frame.swapchain;
    color.load_op = SDL_GPU_LOADOP_LOAD;
    color.store_op = SDL_GPU_STOREOP_STORE;
    return SDL_BeginGPURenderPass(frame.cmd, &color, 1, nullptr);
}

bool Renderer::ensure_targets(uint32_t width, uint32_t height)
{
    if (m_depth && width == m_target_width && height == m_target_height) {
        return true;
    }
    release_targets();

    SDL_GPUTextureCreateInfo info = {};
    info.type = SDL_GPU_TEXTURETYPE_2D;
    info.width = width;
    info.height = height;
    info.layer_count_or_depth = 1;
    info.num_levels = 1;
    info.sample_count = m_scene_format.samples;

    if (m_scene_format.samples != SDL_GPU_SAMPLECOUNT_1) {
        info.format = m_scene_format.color;
        info.usage = SDL_GPU_TEXTUREUSAGE_COLOR_TARGET;
        m_msaa_color = SDL_CreateGPUTexture(m_device, &info);
        if (!m_msaa_color) {
            SDL_LogError(SDL_LOG_CATEGORY_GPU, "MSAA color target creation failed: %s", SDL_GetError());
            return false;
        }
    }

    info.sample_count = SDL_GPU_SAMPLECOUNT_1;
    info.format = m_scene_format.color;
    info.usage = SDL_GPU_TEXTUREUSAGE_COLOR_TARGET | SDL_GPU_TEXTUREUSAGE_SAMPLER;
    m_hdr = SDL_CreateGPUTexture(m_device, &info);
    if (!m_hdr) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "HDR target creation failed: %s", SDL_GetError());
        release_targets();
        return false;
    }

    info.sample_count = m_scene_format.samples;

    info.format = m_scene_format.depth;
    info.usage = SDL_GPU_TEXTUREUSAGE_DEPTH_STENCIL_TARGET;
    m_depth = SDL_CreateGPUTexture(m_device, &info);
    if (!m_depth) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Depth target creation failed: %s", SDL_GetError());
        release_targets();
        return false;
    }

    m_target_width = width;
    m_target_height = height;
    return true;
}

void Renderer::release_targets()
{
    if (m_msaa_color) {
        SDL_ReleaseGPUTexture(m_device, m_msaa_color);
        m_msaa_color = nullptr;
    }
    if (m_hdr) {
        SDL_ReleaseGPUTexture(m_device, m_hdr);
        m_hdr = nullptr;
    }
    if (m_depth) {
        SDL_ReleaseGPUTexture(m_device, m_depth);
        m_depth = nullptr;
    }
    m_target_width = m_target_height = 0;
}

} // namespace astraxis
