#pragma once

#include "render/gpu_device.hpp"
#include "render/scene_targets.hpp"
#include "render/shader.hpp"

#include <SDL3/SDL_gpu.h>

#include <cstdint>

namespace astraxis {

// Additive Gaussian splats of a fixed set of points, drawn off-screen at half
// the scene's resolution with one sample per pixel (render), then added to the
// scene pass (composite): smooth nebulae (a pinwheel's dust, the Homunculus)
// made of thousands of overlapping splats, which blended into every MSAA
// sample at full resolution cost ~16 times the fill rate (8 frames per second
// on an Arc A380 at 1280 x 720, 4x MSAA).
//
// The points are a storage buffer (set_points, once per scene) read by the
// vertex shader given to init: six vertices per point, uniform buffer 0, and
// the outputs dust.frag.hlsl takes (the Gaussian's offset in sigmas and its
// radiance at the center).
class SplatPass {
public:
    // One draw of `count` points (a range the vertex shader picks out with
    // its uniforms) with `uniforms` (`size` bytes) as uniform buffer 0.
    struct Draw {
        const void* uniforms = nullptr;
        uint32_t size = 0;
        uint32_t count = 0;
    };

    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format, const ShaderDesc& vertex_shader);
    void shutdown();
    // Recreates the pipelines for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format);

    void clear();
    bool set_points(const void* data, uint32_t bytes, uint32_t count);
    bool has_points() const { return m_count > 0; }

    // Size of the off-screen target for a scene of width x height pixels.
    static uint32_t target_size(uint32_t scene_size);

    // Draws every point into targets for a scene of width x height pixels,
    // with `uniforms` (`size` bytes) as the vertex shader's uniform buffer 0;
    // outside any render pass.
    void render(SDL_GPUCommandBuffer* cmd, SplatTargets& targets, const void* uniforms, uint32_t size,
                uint32_t width, uint32_t height);
    // The same with several draws into the one target (several nebulae, or
    // one in layers).
    void render(SDL_GPUCommandBuffer* cmd, SplatTargets& targets, const Draw* draws, size_t draw_count,
                uint32_t width, uint32_t height);
    // Adds what render() drew to the scene pass.
    void composite(SDL_GPURenderPass* pass, SplatTargets& targets);

private:
    bool create_pipelines(const SceneTargetFormat& format);
    bool ensure_targets(SplatTargets& targets, uint32_t width, uint32_t height);

    SDL_GPUDevice* m_device = nullptr;
    ShaderDesc m_vertex_shader;
    SDL_GPUGraphicsPipeline* m_splat_pipeline = nullptr;
    SDL_GPUGraphicsPipeline* m_composite_pipeline = nullptr;
    SDL_GPUSampler* m_linear = nullptr;
    SDL_GPUBuffer* m_buffer = nullptr;
    uint32_t m_count = 0;
};

} // namespace astraxis
