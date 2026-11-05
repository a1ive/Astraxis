#pragma once

#include "render/gpu_device.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat3x3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <array>
#include <span>

namespace astraxis {

inline constexpr int kAuroraOvalSamples = 48; // AuroraEmission::kOvalSamples
inline constexpr int kJovianOvalSamples = 256; // JovianAurora::kSamples
inline constexpr int kJovianOvalPairs = kJovianOvalSamples / 2;
inline constexpr int kAuroraTableVectors = 4; // AuroraProfile::table

struct AuroraDrawItem {
    // Camera-relative translation * body-fixed axes * radii: the planet's
    // unit-sphere frame (as for AtmosphereDrawItem).
    glm::mat4 model{1.0f};
    glm::mat3 to_unit{1.0f};
    glm::vec3 camera_unit{0.0f};
    float km_per_unit = 1.0f;  // the body's equatorial radius
    float bottom = 1.0f;       // shell radii (unit-sphere frame)
    float top = 1.0f;
    float polar_ratio = 1.0f;  // polar over equatorial radius
    // What the shader does: AuroraPattern, AuroraVertical, AuroraLighting.
    int pattern = 0;
    int vertical = 0;
    int lighting = 0;
    // Earth: the north geomagnetic pole and magnetic noon (the sun's direction
    // across the dipole axis). Saturn: the spin axis and solar noon. Jupiter:
    // the north oval's center and east, with south_center, south_east.
    glm::vec3 pole{0.0f, 0.0f, 1.0f};
    glm::vec3 noon{1.0f, 0.0f, 0.0f};
    glm::vec3 south_center{0.0f};
    glm::vec3 south_east{0.0f};
    glm::vec3 sun{1.0f, 0.0f, 0.0f};
    float sun_angular_radius = 0.0f;
    glm::vec4 parent{0.0f}; // xyz = center (unit-sphere frame), w = radius; casts the shadow of NightOrEclipsed
    // The two layers' linear RGB, exposure included, and their parameters
    // (AuroraProfile::low_rgb, high_rgb, layers).
    glm::vec3 low_rgb{0.0f};
    glm::vec3 high_rgb{0.0f};
    glm::vec4 layers{0.0f};
    float blue_kr = 0.0f;      // AuroraProfile::blue_kr
    float diffuse_flux = 0.0f; // AuroraProfile::diffuse_flux
    float map_gain = 0.0f;     // AuroraProfile::map_gain
    glm::vec4 map_sun{0.0f};   // AuroraProfile::map_sun
    float sigma_km = 0.0f;     // AuroraProfile::sigma_km
    float time = 0.0f;         // animation clock (s)
    float pixel_angle = 0.0f;  // radians per pixel, for filtering the pattern
    std::array<glm::vec4, kAuroraOvalSamples> oval{}; // AuroraProfile::oval
    std::array<glm::vec4, kJovianOvalPairs> jovian_oval{}; // two samples each: xy, zw = north, south
    std::array<glm::vec4, 6> footprints{}; // physical unit direction, display gain
    std::array<glm::vec4, kAuroraTableVectors> table{}; // AuroraProfile::table, four per vector
    SDL_GPUTexture* map = nullptr; // BodyMap's brightness
    // Saturn is composited after the rings; each ray sample behind their
    // equatorial plane uses the same radial optical-depth profile.
    bool after_rings = false;
    SDL_GPUTexture* ring_profile = nullptr;
    glm::vec2 ring_radii{0.0f}; // inner, outer (km)
};

// Emitting shells above planets and moons (scene/aurora_model.hpp).
class AuroraPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();
    // Recreates the pipelines for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format);

    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
              std::span<const AuroraDrawItem> items, bool after_rings = false) const;

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
    SDL_GPUSampler* m_ring_sampler = nullptr;
    SDL_GPUSampler* m_map_sampler = nullptr; // longitude wraps, latitude clamps
};

} // namespace astraxis
