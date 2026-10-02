#pragma once

#include <SDL3/SDL_gpu.h>

#include <cstdint>

struct SDL_Window;

namespace astraxis {

// Per-frame GPU state. `swapchain` may be null (e.g. window occluded);
// the command buffer must still be submitted via Renderer::end_frame.
struct Frame {
    SDL_GPUCommandBuffer* cmd = nullptr;
    SDL_GPUTexture* swapchain = nullptr;
    uint32_t width = 0;
    uint32_t height = 0;
};

// Formats and sample count shared by all scene pipelines. The scene is
// rendered in linear HDR and tonemapped to the swapchain by PostProcess.
struct SceneTargetFormat {
    SDL_GPUTextureFormat color = SDL_GPU_TEXTUREFORMAT_INVALID;
    SDL_GPUTextureFormat depth = SDL_GPU_TEXTUREFORMAT_INVALID;
    SDL_GPUSampleCount samples = SDL_GPU_SAMPLECOUNT_1;
};

// Owns the SDL_GPU device, the window's swapchain and the scene render
// targets (MSAA HDR color resolved into a sampleable HDR texture, plus
// reversed-Z depth).
class Renderer {
public:
    bool init(SDL_Window* window, bool debug);
    void shutdown();

    SDL_GPUDevice* device() const { return m_device; }
    SDL_GPUTextureFormat swapchain_format() const;
    const SceneTargetFormat& scene_format() const { return m_scene_format; }
    const char* driver_name() const;

    bool begin_frame(Frame& frame);
    void end_frame(Frame& frame);

    // Begins the scene pass: clears color and depth (to 0, reversed-Z) and
    // resolves MSAA into hdr_texture() at the end.
    SDL_GPURenderPass* begin_scene_pass(const Frame& frame, SDL_FColor clear_color);
    SDL_GPUTexture* hdr_texture() const { return m_hdr; }
    // Begins an overlay pass that draws on top of the swapchain (e.g. ImGui).
    SDL_GPURenderPass* begin_overlay_pass(const Frame& frame);

private:
    bool ensure_targets(uint32_t width, uint32_t height);
    void release_targets();

    SDL_Window* m_window = nullptr;
    SDL_GPUDevice* m_device = nullptr;
    SceneTargetFormat m_scene_format;

    SDL_GPUTexture* m_msaa_color = nullptr;
    SDL_GPUTexture* m_hdr = nullptr;
    SDL_GPUTexture* m_depth = nullptr;
    uint32_t m_target_width = 0;
    uint32_t m_target_height = 0;
};

} // namespace astraxis
