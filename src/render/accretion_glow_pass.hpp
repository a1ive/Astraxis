#pragma once

#include "render/dust_pass.hpp"
#include "render/gpu_device.hpp"
#include "render/scene_targets.hpp"
#include "render/splat_pass.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <cstdint>
#include <vector>

namespace astraxis {

inline constexpr int kMaxAccretionGlowOccluders = 2; // MAX_OCCLUDERS in accretion_glow.vert.hlsl

// One accretion system's share of the frame's splats, and the stars that hide them.
struct AccretionGlowItem {
    uint32_t first = 0; // in AccretionGlowPass::splats()
    uint32_t count = 0;
    glm::mat3 to_local{1.0f}; // world -> disk frame
    glm::vec3 center{0.0f};   // the gainer relative to the camera (km)
    int occluder_count = 0;
    glm::vec3 occluder_centers[kMaxAccretionGlowOccluders] = {}; // disk frame (km)
    glm::vec3 occluder_radii[kMaxAccretionGlowOccluders] = {};   // semi-axes along the disk's axes (km)
};

// The glowing gas around an accretion disk (the stream, the hot spot, the jets):
// hundreds of large overlapping splats rebuilt every frame, drawn at half
// resolution (SplatPass) and added to the scene. Usage per frame: the splats
// and items -> render() (outside render passes) -> composite().
class AccretionGlowPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown() { m_splats.shutdown(); }
    // Recreates the pipelines for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format) { return m_splats.set_format(format); }

    std::vector<DustSplat>& splats() { return m_points; }
    std::vector<AccretionGlowItem>& items() { return m_items; }

    // Uploads this frame's splats and draws them into `targets` for a scene of
    // width x height pixels; outside any render pass.
    void render(SDL_GPUCommandBuffer* cmd, SplatTargets& targets, const CameraView& view, uint32_t width,
                uint32_t height);
    // Adds the glow drawn by render() to the scene pass.
    void composite(SDL_GPURenderPass* pass, SplatTargets& targets) { m_splats.composite(pass, targets); }

private:
    SplatPass m_splats;
    std::vector<DustSplat> m_points;
    std::vector<AccretionGlowItem> m_items;
};

} // namespace astraxis
