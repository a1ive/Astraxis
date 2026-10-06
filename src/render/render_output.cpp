#include "render/render_output.hpp"

#include <SDL3/SDL_log.h>
#include <SDL3/SDL_video.h>

#include <algorithm>
#include <cmath>

namespace astraxis {

void scale_scene_size(uint32_t output_width, uint32_t output_height, float scale, uint32_t& width, uint32_t& height)
{
    width = std::max(1u, static_cast<uint32_t>(std::lround(static_cast<float>(output_width) * scale)));
    height = std::max(1u, static_cast<uint32_t>(std::lround(static_cast<float>(output_height) * scale)));
}

bool RenderOutput::init(const GpuDevice& device, SDL_Window* window)
{
    if (!SDL_ClaimWindowForGPUDevice(device.device(), window)) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_ClaimWindowForGPUDevice failed: %s", SDL_GetError());
        return false;
    }
    m_device = &device;
    m_window = window;
    SDL_SetGPUSwapchainParameters(device.device(), window, SDL_GPU_SWAPCHAINCOMPOSITION_SDR,
                                  SDL_GPU_PRESENTMODE_VSYNC);
    return true;
}

void RenderOutput::shutdown()
{
    if (m_device && m_window) {
        SDL_WaitForGPUIdle(m_device->device());
        m_targets.release();
        SDL_ReleaseWindowFromGPUDevice(m_device->device(), m_window);
    }
    m_device = nullptr;
    m_window = nullptr;
}

SDL_GPUTextureFormat RenderOutput::swapchain_format() const
{
    return SDL_GetGPUSwapchainTextureFormat(m_device->device(), m_window);
}

bool RenderOutput::begin_frame(Frame& frame)
{
    frame = {};
    frame.cmd = SDL_AcquireGPUCommandBuffer(m_device->device());
    if (!frame.cmd) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_AcquireGPUCommandBuffer failed: %s", SDL_GetError());
        return false;
    }

    if (!SDL_WaitAndAcquireGPUSwapchainTexture(frame.cmd, m_window, &frame.swapchain, &frame.width, &frame.height)) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_WaitAndAcquireGPUSwapchainTexture failed: %s", SDL_GetError());
        frame.swapchain = nullptr;
    }
    if (frame.swapchain) {
        frame.output_width = frame.width;
        frame.output_height = frame.height;
        scale_scene_size(frame.output_width, frame.output_height, m_render_scale, frame.width, frame.height);
        if (!m_targets.ensure(m_device->device(), m_device->scene_format(), frame.width, frame.height)) {
            frame.swapchain = nullptr;
        }
    }
    frame.targets = &m_targets;
    return true;
}

void RenderOutput::end_frame(Frame& frame)
{
    if (frame.cmd) {
        SDL_SubmitGPUCommandBuffer(frame.cmd);
    }
    frame = {};
}

SDL_GPURenderPass* RenderOutput::begin_overlay_pass(const Frame& frame)
{
    SDL_GPUColorTargetInfo color = {};
    color.texture = frame.swapchain;
    color.load_op = SDL_GPU_LOADOP_LOAD;
    color.store_op = SDL_GPU_STOREOP_STORE;
    return SDL_BeginGPURenderPass(frame.cmd, &color, 1, nullptr);
}

} // namespace astraxis
