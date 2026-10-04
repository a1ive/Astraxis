#pragma once

#include "render/renderer.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat3x3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>

#include <span>

namespace astraxis {

// Optical properties of an atmosphere in units of the planet radius.
struct AtmosphereOptics {
    glm::vec3 rayleigh_depth{0.0f}; // vertical optical depth per channel
    float rayleigh_scale_height = 1.0f;
    // Haze attenuates with tau (1 - w g): what it scatters into its forward
    // peak mostly keeps going (similarity relation). It scatters with tau w.
    glm::vec3 haze_attenuation{0.0f};
    float haze_scale_height = 1.0f;
    glm::vec3 haze_scattering{0.0f};
    float haze_g = 0.0f;
};

struct AtmosphereDrawItem {
    // Camera-relative translation * body-fixed axes * radii of the shell's base
    // (the planet, or its opaque deck): the planet's unit-sphere frame.
    glm::mat4 model{1.0f};
    glm::mat3 to_unit{1.0f};        // inverse of model's linear part
    glm::vec3 camera_unit{0.0f};    // camera in the unit-sphere frame
    float top = 1.0f;               // shell top radius (unit-sphere frame)
    glm::vec3 sun_direction{1.0f, 0.0f, 0.0f}; // world, unit
    float sun_angular_radius = 0.0f;
    AtmosphereOptics optics;
    float gain = 1.0f;
};

// Draws atmospheres as shells around their planets: single-scattered sunlight
// is added and what lies behind is dimmed by the view transmission (depth-
// tested, not written, so nearer bodies hide the shell).
class AtmospherePass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
              std::span<const AtmosphereDrawItem> items) const;

private:
    SDL_GPUDevice* m_device = nullptr;
    // Camera outside the shell: front faces, depth-tested. Inside (or near
    // enough for the front faces to be clipped): back faces, not depth-tested,
    // with the rays starting at the camera.
    SDL_GPUGraphicsPipeline* m_outside = nullptr;
    SDL_GPUGraphicsPipeline* m_inside = nullptr;
    SDL_GPUBuffer* m_vertices = nullptr;
    SDL_GPUBuffer* m_indices = nullptr;
    uint32_t m_index_count = 0;
};

} // namespace astraxis
