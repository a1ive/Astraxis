#pragma once

#include "render/renderer.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <span>

namespace astraxis {

inline constexpr int kMaxOccluders = 8;

struct BodyDrawItem {
    glm::mat4 model{1.0f};    // camera-relative translation * rotation * scale
    glm::mat4 rotation{1.0f}; // body-fixed axes (including the texture longitude offset)
    glm::vec3 inv_scale{1.0f};
    glm::vec3 color{1.0f};
    int style = 0;
    glm::vec3 sun_direction{1.0f, 0.0f, 0.0f}; // unit vector from the body toward the sun
    float sun_angular_radius = 0.0f;            // as seen from the body (radians)
    SDL_GPUTexture* texture = nullptr; // optional equirectangular albedo (sRGB)
    bool flip_u = false;               // map longitudes increase westward
    int occluder_count = 0;
    glm::vec4 occluders[kMaxOccluders] = {}; // camera-relative center, radius
};

struct SunLight {
    float ambient = 0.0f;
};

// Draws planets and moons as lit ellipsoids.
class BodyPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, const SunLight& sun,
              std::span<const BodyDrawItem> items) const;

private:
    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
    SDL_GPUBuffer* m_vertices = nullptr;
    SDL_GPUBuffer* m_indices = nullptr;
    uint32_t m_index_count = 0;
    SDL_GPUSampler* m_sampler = nullptr;
    SDL_GPUTexture* m_white = nullptr; // bound when a body has no texture
};

} // namespace astraxis
