#pragma once

#include "render/renderer.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <cstdint>
#include <span>
#include <vector>

namespace astraxis {

inline constexpr int kMaxOccluders = 8;

struct BodyDrawItem {
    int mesh = -1;            // index returned by BodyPass::add_mesh; -1 = the ellipsoid
    glm::mat4 model{1.0f};    // camera-relative translation * rotation * scale
    glm::mat4 rotation{1.0f}; // body-fixed axes (x toward the prime meridian, z along the pole)
    glm::vec3 inv_scale{1.0f};
    glm::vec3 color{1.0f};
    int style = 0;
    glm::vec3 sun_direction{1.0f, 0.0f, 0.0f}; // unit vector from the body toward the sun
    float sun_angular_radius = 0.0f;            // as seen from the body (radians)
    SDL_GPUTexture* texture = nullptr; // optional equirectangular albedo (sRGB), ellipsoids only
    float texture_left_lon_deg = 0.0f; // east longitude of the map's left edge
    bool flip_u = false;               // map longitudes increase westward
    int occluder_count = 0;
    glm::vec4 occluders[kMaxOccluders] = {}; // camera-relative center, radius
    // The body's own rings, which shade it (null: none).
    SDL_GPUTexture* ring_profile = nullptr; // create_profile_texture
    glm::vec3 ring_center{0.0f};            // camera-relative
    glm::vec3 ring_normal{0.0f, 0.0f, 1.0f};
    float ring_inner_km = 0.0f;
    float ring_outer_km = 0.0f;
    float ring_samples = 1.0f;
};

struct SunLight {
    float ambient = 0.0f;
};

// Draws planets and moons as lit ellipsoids, and irregular bodies as meshes.
class BodyPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    // Meshes (shape models) are uploaded once, when a scene is loaded: positions in
    // km in the body-fixed frame, unit normals, albedo relative to the body color,
    // and counter-clockwise triangles seen from outside.
    void clear_meshes();
    int add_mesh(std::span<const glm::vec3> positions, std::span<const glm::vec3> normals,
                 std::span<const float> albedo, std::span<const uint32_t> indices);

    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, const SunLight& sun,
              std::span<const BodyDrawItem> items) const;

private:
    struct Mesh {
        SDL_GPUBuffer* vertices = nullptr;
        SDL_GPUBuffer* indices = nullptr; // 32-bit
        uint32_t index_count = 0;
    };

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
    Mesh m_sphere;
    std::vector<Mesh> m_meshes;
    SDL_GPUSampler* m_sampler = nullptr;
    SDL_GPUSampler* m_profile_sampler = nullptr;
    SDL_GPUTexture* m_white = nullptr;    // bound when a body has no texture
    SDL_GPUTexture* m_no_rings = nullptr; // bound when a body has no rings
};

} // namespace astraxis
