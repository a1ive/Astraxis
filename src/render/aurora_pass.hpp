#pragma once

#include "render/gpu_device.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat3x3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <array>
#include <span>

namespace astraxis {

inline constexpr int kAuroraOvalSamples = 48; // AuroraEmission::kOvalSamples

struct AuroraDrawItem {
    // Camera-relative translation * body-fixed axes * radii: the planet's
    // unit-sphere frame (as for AtmosphereDrawItem).
    glm::mat4 model{1.0f};
    glm::mat3 to_unit{1.0f};
    glm::vec3 camera_unit{0.0f};
    float km_per_unit = 1.0f;  // the planet's equatorial radius
    float bottom = 1.0f;       // shell radii (unit-sphere frame)
    float top = 1.0f;
    // Unit-sphere frame: the north geomagnetic pole, magnetic noon (the sun's
    // direction across the dipole axis) and the sun.
    glm::vec3 pole{0.0f, 0.0f, 1.0f};
    glm::vec3 noon{1.0f, 0.0f, 0.0f};
    glm::vec3 sun{1.0f, 0.0f, 0.0f};
    // A vertical column through the oval: the linear RGB of its low layer
    // (green, blue), exposure included, and its 427.8 nm brightness (kR); the
    // red line's linear RGB per kR, exposure included.
    glm::vec3 low_rgb{0.0f};
    float blue_kr = 0.0f;
    glm::vec3 red_per_kr{0.0f};
    // Emission layers (km): low peak, low scale height, high peak, high scale height.
    glm::vec4 layers{115.0f, 10.0f, 221.0f, 40.0f};
    float diffuse_flux = 0.0f; // the diffuse aurora's flux relative to the oval's
    float time = 0.0f;         // animation clock (s)
    float pixel_angle = 0.0f;  // radians per pixel, for the detail of the arcs
    // Boundaries around the oval (AuroraEmission::oval): colatitudes (rad) of
    // the diffuse aurora, the oval's equatorward and poleward edges.
    std::array<glm::vec4, kAuroraOvalSamples> oval{};
};

// Draws Earth's aurora as an emitting shell (Chapman layers in height, the
// auroral ovals around the geomagnetic poles), added to what lies behind.
class AuroraPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();
    // Recreates the pipelines for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format);

    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
              std::span<const AuroraDrawItem> items) const;

private:
    bool create_pipelines(const SceneTargetFormat& format);

    SDL_GPUDevice* m_device = nullptr;
    // As for the atmosphere: outside the shell its front faces, depth-tested;
    // inside, its back faces with the rays starting at the camera.
    SDL_GPUGraphicsPipeline* m_outside = nullptr;
    SDL_GPUGraphicsPipeline* m_inside = nullptr;
    SDL_GPUBuffer* m_vertices = nullptr;
    SDL_GPUBuffer* m_indices = nullptr;
    uint32_t m_index_count = 0;
    SDL_GPUTexture* m_noise = nullptr; // tiling value noise for the arcs and rays
    SDL_GPUSampler* m_sampler = nullptr;
};

} // namespace astraxis
