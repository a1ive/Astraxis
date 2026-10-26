#pragma once

#include "render/gpu_device.hpp"
#include "render/scene_targets.hpp"
#include "render/splat_pass.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <cstdint>

namespace astraxis {

// One patch as the shader reads it (32 bytes): see NebulaGrain.
struct NebulaGpuGrain {
    glm::vec3 position{0.0f};
    float weight = 1.0f;
    glm::vec3 normal{0.0f, 0.0f, 1.0f};
    float padding = 0.0f;
};
static_assert(sizeof(NebulaGpuGrain) == 32, "NebulaGpuGrain must match the shader's raw buffer loads");

// What an expanding nebula looks like this frame (see shaders/nebula.vert.hlsl).
struct NebulaDrawItem {
    glm::vec3 center{0.0f};         // the star relative to the camera (km, display axes)
    glm::mat3 display_from_icrf{1.0f};
    float scale = 1.0f;             // size relative to the reference time
    float patch_sigma_km = 0.0f;    // at the reference time
    float wall_sigma_km = 0.0f;     // at the reference time
    float min_sigma_px = 0.7f;      // in the off-screen target's pixels
    float asymmetry_g = 0.0f;
    float reference_km = 1.0f;      // where a face-on patch has `radiance` (phase function 1, reference time)
    float radiance = 1.0f;          // linear HDR
    float gamma = 1.0f;             // display compression
    float cull_radiance = 0.0f;     // splats fainter than this at their center are skipped
    glm::vec3 color{1.0f};          // of the scattered light (sRGB)
};

// An expanding dusty shell (scene/nebula.hpp): its patches are uploaded once
// (set_grains, when a scene is loaded), scaled from the star on the GPU each
// frame and splatted off-screen (SplatPass).
class NebulaPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown() { m_splats.shutdown(); }
    // Recreates the pipelines for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format) { return m_splats.set_format(format); }

    void clear() { m_splats.clear(); }
    bool set_grains(const NebulaGpuGrain* grains, uint32_t count);
    bool has_grains() const { return m_splats.has_points(); }

    // Draws the shell into targets.nebula for a scene of width x height
    // pixels; outside any render pass.
    void render(SDL_GPUCommandBuffer* cmd, SplatTargets& targets, const CameraView& view, const NebulaDrawItem& item,
                uint32_t width, uint32_t height);
    // Adds the shell drawn by render() to the scene pass.
    void composite(SDL_GPURenderPass* pass, SplatTargets& targets) { m_splats.composite(pass, targets); }

private:
    SplatPass m_splats;
};

} // namespace astraxis
