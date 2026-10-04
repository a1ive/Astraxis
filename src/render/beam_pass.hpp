#pragma once

#include "render/renderer.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/vec3.hpp>

#include <span>

namespace astraxis {

struct BeamDrawItem {
    glm::vec3 apex{0.0f};      // camera-relative star position (km)
    glm::vec3 axis{0.0f, 0.0f, 1.0f}; // unit beam direction
    float length_km = 1.0f;
    float half_angle_rad = 0.1f;
    glm::vec3 color{1.0f};     // sRGB
    float intensity = 1.0f;    // linear HDR scale
};

// Pulsar beams: additive, depth-tested (bodies in front hide them), no depth writes.
class BeamPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();
    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
              std::span<const BeamDrawItem> items) const;

private:
    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
};

} // namespace astraxis
