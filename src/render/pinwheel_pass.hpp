#pragma once

#include "render/gpu_device.hpp"
#include "render/scene_targets.hpp"
#include "render/splat_pass.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <array>
#include <cstdint>

namespace astraxis {

inline constexpr int kPinwheelColorCount = 64;   // LUT_SIZE in pinwheel.vert.hlsl
inline constexpr int kPinwheelRadiusCount = 128; // RADIUS_SIZE in pinwheel.vert.hlsl

// One grain as the shader reads it (32 bytes): see PinwheelGrain.
struct PinwheelGpuGrain {
    glm::vec3 direction{1.0f, 0.0f, 0.0f};
    float launch = 0.0f; // turn + phase
    glm::vec3 tangent{0.0f, 1.0f, 0.0f};
    float weight = 1.0f;
};
static_assert(sizeof(PinwheelGpuGrain) == 32, "PinwheelGpuGrain must match the shader's raw buffer loads");

// What a pinwheel looks like this frame (see shaders/pinwheel.vert.hlsl).
struct PinwheelDrawItem {
    glm::vec3 center{0.0f};                 // the pair's barycenter relative to the camera (km, display axes)
    glm::mat3 display_from_icrf{1.0f};
    float cycle = 0.0f;                     // frac((t - t_ref) / period / turns)
    float turns = 1.0f;
    float nucleation_km = 1.0f;
    float sigma_along = 0.0f;               // patch sigma along the ring, per km of r
    float sigma_across = 0.0f;              // patch sigma across the rings, per km of r
    float wall_sigma_km = 0.0f;
    float min_sigma_px = 0.7f;              // in the off-screen target's pixels
    float reference_km = 1.0f;              // where a face-on wall has `radiance`
    float radiance = 1.0f;                  // face-on surface radiance at reference_km (linear HDR)
    float gamma = 1.0f;                     // display compression of the falloff
    float cull_radiance = 0.0f;             // splats fainter than this at their center are skipped
    float max_age = 1e30f;                  // older grains (periods) left before the start: not drawn
    // At ln(r / r_nuc) / ln(r_end / r_nuc) = i / (kPinwheelColorCount - 1): linear
    // RGB (largest component 1), and in alpha a factor on the brightness
    // (already raised to `gamma`).
    std::array<glm::vec4, kPinwheelColorCount> colors{};
    // Distance from the center (km) at age = turns * (i / (kPinwheelRadiusCount - 1))^2 periods.
    std::array<float, kPinwheelRadiusCount> radius{};
};

// A colliding-wind binary's dust spiral (scene/pinwheel.hpp): its grains are
// uploaded once (set_grains, when a scene is loaded) and placed on the GPU
// each frame from the time alone, splatted off-screen (SplatPass).
class PinwheelPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown() { m_splats.shutdown(); }
    // Recreates the pipelines for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format) { return m_splats.set_format(format); }

    void clear() { m_splats.clear(); }
    bool set_grains(const PinwheelGpuGrain* grains, uint32_t count);
    bool has_grains() const { return m_splats.has_points(); }

    // Draws the dust into targets.pinwheel for a scene of width x height
    // pixels; outside any render pass.
    void render(SDL_GPUCommandBuffer* cmd, SplatTargets& targets, const CameraView& view,
                const PinwheelDrawItem& item, uint32_t width, uint32_t height);
    // Adds the dust drawn by render() to the scene pass.
    void composite(SDL_GPURenderPass* pass, SplatTargets& targets) { m_splats.composite(pass, targets); }

private:
    SplatPass m_splats;
};

} // namespace astraxis
