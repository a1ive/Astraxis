#pragma once

#include "render/renderer.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <span>

namespace astraxis {

inline constexpr int kMaxRingBands = 8;

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
    int band_count = 0;
    glm::vec4 bands[kMaxRingBands] = {}; // inner, outer (km), optical depth, thickness (km)
};

// Draws faint planetary rings as additive, single-scattering annuli in the
// planet's equatorial plane (depth-tested, not written).
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
};

} // namespace astraxis
