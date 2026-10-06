#pragma once

#include "render/gpu_device.hpp"
#include "render/render_output.hpp"
#include "render/scene_targets.hpp"

#include <SDL3/SDL_gpu.h>

#include <cstdint>
#include <vector>

namespace astraxis {

// An output without a swapchain: frames are drawn into a texture and read
// back to the CPU, for a window the host draws into itself. Used for the
// screensaver preview, a child of another process's window: flip-model
// swapchains (the only kind D3D12 has) do not show up there.
class OffscreenOutput {
public:
    // `format`: an SDR swapchain format; the pixels are what a swapchain of
    // that format would show (sRGB encoded by PostProcess).
    bool init(const GpuDevice& device, SDL_GPUTextureFormat format);
    void shutdown();

    SDL_GPUTextureFormat format() const { return m_format; }

    // Like RenderOutput::begin_frame, at the given size.
    bool begin_frame(Frame& frame, uint32_t width, uint32_t height);
    // Submits the frame, waits for it and returns its pixels (4 bytes each,
    // rows top to bottom, no padding); null if the frame had no texture.
    const uint8_t* end_frame(Frame& frame);

private:
    bool ensure_size(uint32_t width, uint32_t height);

    const GpuDevice* m_device = nullptr;
    SDL_GPUTextureFormat m_format = SDL_GPU_TEXTUREFORMAT_INVALID;
    SDL_GPUTexture* m_texture = nullptr;
    SDL_GPUTransferBuffer* m_download = nullptr;
    uint32_t m_width = 0;
    uint32_t m_height = 0;
    SceneTargets m_targets;
    std::vector<uint8_t> m_pixels;
};

} // namespace astraxis
