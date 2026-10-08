#pragma once

#include "render/gpu_device.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <array>
#include <span>

namespace astraxis {

inline constexpr int kMaxRingBands = 16; // drawn one by one (eccentric or inclined rings)

struct RingDrawItem {
    glm::mat4 model{1.0f};          // camera-relative translation * body-fixed axes (km)
    glm::vec3 sun_direction{1.0f};  // body frame, unit
    float sun_angular_radius = 0.0f; // radians
    glm::vec3 camera{0.0f};         // camera position in the body frame (km)
    glm::vec3 color{1.0f};          // sRGB tint
    float gain = 1.0f;
    float phase_g = 0.0f;
    float equatorial_radius = 1.0f; // of the planet, for its shadow (km)
    float polar_radius = 1.0f;
    float inner_km = 0.0f; // radial range of the profile
    float outer_km = 0.0f;
    float mesh_inner_km = 0.0f; // radial range of the drawn annulus
    float mesh_outer_km = 0.0f;
    SDL_GPUTexture* profile = nullptr; // create_profile_texture: optical depth, mu floor

    // Bands drawn in their true shapes instead of the profile (0: the profile),
    // in the body frame (see shaders/ring.frag.hlsl).
    int band_count = 0;
    std::array<glm::vec4, kMaxRingBands> band_edges{};  // a and e of the inner and outer edge
    std::array<glm::vec4, kMaxRingBands> band_shape{};  // xy = unit vector to periapsis, zw = node (a sin i)
    std::array<glm::vec4, kMaxRingBands> band_optics{}; // x = mean tau x mean width (km), y = mu floor
};

// Draws planetary rings as single-scattering annuli in the planet's
// equatorial plane: scattered light is added, and what lies behind is dimmed by
// the ring's transmission (depth-tested, not written). Eccentric and inclined
// narrow rings are found along each pixel's line of sight from that plane.
class RingPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
              std::span<const RingDrawItem> items) const;

private:
    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
    SDL_GPUBuffer* m_vertices = nullptr;
    uint32_t m_vertex_count = 0;
    SDL_GPUSampler* m_sampler = nullptr;
};

} // namespace astraxis
