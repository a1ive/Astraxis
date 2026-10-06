#include "render/offscreen_output.hpp"

#include <SDL3/SDL_log.h>

#include <cstring>

namespace astraxis {

bool OffscreenOutput::init(const GpuDevice& device, SDL_GPUTextureFormat format)
{
    m_device = &device;
    m_format = format;
    return true;
}

void OffscreenOutput::shutdown()
{
    if (!m_device) {
        return;
    }
    SDL_GPUDevice* device = m_device->device();
    SDL_WaitForGPUIdle(device);
    m_targets.release();
    if (m_texture) {
        SDL_ReleaseGPUTexture(device, m_texture);
        m_texture = nullptr;
    }
    if (m_download) {
        SDL_ReleaseGPUTransferBuffer(device, m_download);
        m_download = nullptr;
    }
    m_width = 0;
    m_height = 0;
    m_device = nullptr;
}

bool OffscreenOutput::ensure_size(uint32_t width, uint32_t height)
{
    if (m_texture && width == m_width && height == m_height) {
        return true;
    }
    SDL_GPUDevice* device = m_device->device();
    SDL_WaitForGPUIdle(device);
    if (m_texture) {
        SDL_ReleaseGPUTexture(device, m_texture);
        m_texture = nullptr;
    }
    if (m_download) {
        SDL_ReleaseGPUTransferBuffer(device, m_download);
        m_download = nullptr;
    }

    SDL_GPUTextureCreateInfo tex = {};
    tex.type = SDL_GPU_TEXTURETYPE_2D;
    tex.format = m_format;
    tex.usage = SDL_GPU_TEXTUREUSAGE_COLOR_TARGET;
    tex.width = width;
    tex.height = height;
    tex.layer_count_or_depth = 1;
    tex.num_levels = 1;
    m_texture = SDL_CreateGPUTexture(device, &tex);

    SDL_GPUTransferBufferCreateInfo transfer = {};
    transfer.usage = SDL_GPU_TRANSFERBUFFERUSAGE_DOWNLOAD;
    transfer.size = width * height * 4;
    m_download = SDL_CreateGPUTransferBuffer(device, &transfer);

    if (!m_texture || !m_download) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Offscreen output %ux%u: %s", width, height, SDL_GetError());
        return false;
    }
    m_width = width;
    m_height = height;
    m_pixels.resize(static_cast<size_t>(width) * height * 4);
    return true;
}

bool OffscreenOutput::begin_frame(Frame& frame, uint32_t width, uint32_t height)
{
    frame = {};
    frame.cmd = SDL_AcquireGPUCommandBuffer(m_device->device());
    if (!frame.cmd) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_AcquireGPUCommandBuffer failed: %s", SDL_GetError());
        return false;
    }
    if (width > 0 && height > 0 && ensure_size(width, height) &&
        m_targets.ensure(m_device->device(), m_device->scene_format(), width, height)) {
        frame.swapchain = m_texture;
        frame.width = width;
        frame.height = height;
        frame.output_width = width;
        frame.output_height = height;
    }
    frame.targets = &m_targets;
    return true;
}

const uint8_t* OffscreenOutput::end_frame(Frame& frame)
{
    if (!frame.cmd) {
        return nullptr;
    }
    if (!frame.swapchain) {
        SDL_SubmitGPUCommandBuffer(frame.cmd);
        frame = {};
        return nullptr;
    }

    SDL_GPUCopyPass* copy = SDL_BeginGPUCopyPass(frame.cmd);
    SDL_GPUTextureRegion src = {};
    src.texture = m_texture;
    src.w = m_width;
    src.h = m_height;
    src.d = 1;
    SDL_GPUTextureTransferInfo dst = {};
    dst.transfer_buffer = m_download;
    dst.pixels_per_row = m_width;
    dst.rows_per_layer = m_height;
    SDL_DownloadFromGPUTexture(copy, &src, &dst);
    SDL_EndGPUCopyPass(copy);

    SDL_GPUDevice* device = m_device->device();
    SDL_GPUFence* fence = SDL_SubmitGPUCommandBufferAndAcquireFence(frame.cmd);
    frame = {};
    if (!fence) {
        return nullptr;
    }
    SDL_WaitForGPUFences(device, true, &fence, 1);
    SDL_ReleaseGPUFence(device, fence);

    const void* mapped = SDL_MapGPUTransferBuffer(device, m_download, false);
    if (!mapped) {
        return nullptr;
    }
    std::memcpy(m_pixels.data(), mapped, m_pixels.size());
    SDL_UnmapGPUTransferBuffer(device, m_download);
    return m_pixels.data();
}

} // namespace astraxis
