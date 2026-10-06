#pragma once

#include "render/gpu_device.hpp"
#include "render/scene_targets.hpp"

#include <SDL3/SDL_gpu.h>

#include <cstdint>

struct SDL_Window;

namespace astraxis {

// Per-frame GPU state of one output. `swapchain` may be null (e.g. window
// occluded); the command buffer must still be submitted via RenderOutput::end_frame.
struct Frame {
    SDL_GPUCommandBuffer* cmd = nullptr;
    SDL_GPUTexture* swapchain = nullptr;
    uint32_t width = 0;  // the scene's size: the swapchain texture's times the render scale
    uint32_t height = 0;
    uint32_t output_width = 0; // the swapchain texture's size
    uint32_t output_height = 0;
    SceneTargets* targets = nullptr; // the output's, sized to the scene
};

// The scene size for an output size and a render scale (at least 1 x 1).
void scale_scene_size(uint32_t output_width, uint32_t output_height, float scale, uint32_t& width, uint32_t& height);

// One window drawn by a GpuDevice: its swapchain and its scene targets. A
// device may drive several outputs (e.g. one per monitor); their swapchains
// all use SDR composition, so they share one format (and one PostProcess).
class RenderOutput {
public:
    bool init(const GpuDevice& device, SDL_Window* window);
    void shutdown();

    SDL_Window* window() const { return m_window; }
    SDL_GPUTextureFormat swapchain_format() const;

    // Renders the scene at this fraction of the swapchain's resolution
    // (PostProcess scales it up); 1 by default.
    void set_render_scale(float scale) { m_render_scale = scale; }

    // Acquires a command buffer and the swapchain texture, and sizes the scene
    // targets to it. False only if no command buffer could be acquired.
    bool begin_frame(Frame& frame);
    void end_frame(Frame& frame);

    // Begins a pass that draws on top of the swapchain texture (e.g. ImGui).
    static SDL_GPURenderPass* begin_overlay_pass(const Frame& frame);

private:
    const GpuDevice* m_device = nullptr;
    SDL_Window* m_window = nullptr;
    SceneTargets m_targets;
    float m_render_scale = 1.0f;
};

} // namespace astraxis
