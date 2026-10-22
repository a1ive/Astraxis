#pragma once

#include "render/gpu_device.hpp"
#include "render/scene_targets.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <array>
#include <cstdint>

namespace astraxis {

inline constexpr int kPinwheelColorCount = 64; // LUT_SIZE in pinwheel.vert.hlsl

// What a pinwheel looks like this frame (see shaders/pinwheel.vert.hlsl).
struct PinwheelDrawItem {
    glm::vec3 center{0.0f};                 // the pair's barycenter relative to the camera (km, display axes)
    glm::mat3 display_from_icrf{1.0f};
    float cycle = 0.0f;                     // frac((t - t_ref) / period / turns)
    float turns = 1.0f;
    float step_km = 1.0f;                   // speed * period
    float nucleation_km = 1.0f;
    float sigma_along = 0.0f;               // patch sigma along the ring, per km of r
    float sigma_across = 0.0f;              // patch sigma across the rings, per km of r
    float wall_sigma_km = 0.0f;
    float min_sigma_px = 0.7f;              // in the off-screen target's pixels
    float radiance = 1.0f;                  // face-on surface radiance at the nucleation radius (linear HDR)
    float gamma = 1.0f;                     // display compression of the falloff
    float cull_radiance = 0.0f;             // splats fainter than this at their center are skipped
    // At ln(r / r_nuc) / ln(r_end / r_nuc) = i / (kPinwheelColorCount - 1): linear
    // RGB (largest component 1), and in alpha a factor on the brightness
    // (already raised to `gamma`).
    std::array<glm::vec4, kPinwheelColorCount> colors{};
};

// A colliding-wind binary's dust spiral (scene/pinwheel.hpp): its grains are
// uploaded once (set_grains, when a scene is loaded) and placed on the GPU
// each frame from the time alone. They are splatted off-screen at half the
// scene's resolution with one sample per pixel (render), then added to the
// scene pass (composite): the nebula is smooth, and blending its overlapping
// splats into every MSAA sample at full resolution costs ~16 times the fill
// rate (8 frames per second on an Arc A380 at 1280 x 720, 4x MSAA).
class PinwheelPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();
    // Recreates the pipelines for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format);

    void clear();
    // `grains`: PinwheelGrain records (32 bytes each).
    bool set_grains(const void* grains, uint32_t count);
    bool has_grains() const { return m_count > 0; }

    // Draws the dust into targets.pinwheel for a scene of width x height
    // pixels; outside any render pass.
    void render(SDL_GPUCommandBuffer* cmd, PinwheelTargets& targets, const CameraView& view,
                const PinwheelDrawItem& item, uint32_t width, uint32_t height);
    // Adds the dust drawn by render() to the scene pass.
    void composite(SDL_GPURenderPass* pass, PinwheelTargets& targets);

private:
    bool create_pipelines(const SceneTargetFormat& format);
    bool ensure_targets(PinwheelTargets& targets, uint32_t width, uint32_t height);

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_splat_pipeline = nullptr;
    SDL_GPUGraphicsPipeline* m_composite_pipeline = nullptr;
    SDL_GPUSampler* m_linear = nullptr;
    SDL_GPUBuffer* m_buffer = nullptr;
    uint32_t m_count = 0;
};

} // namespace astraxis
