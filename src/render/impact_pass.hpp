#pragma once

#include "render/gpu_device.hpp"
#include "scene/camera.hpp"
#include "scene/impact.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat3x3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <span>

namespace astraxis {

inline constexpr int kMaxImpactSites = kImpactSiteLimit;

// The debris of a planet's impacts, drawn over its surface (shaders/impact_scar.frag.hlsl).
struct ImpactScarDrawItem {
    glm::mat4 model{1.0f};    // body frame -> camera-relative world (translation * rotation)
    glm::mat3 to_local{1.0f}; // world vector -> body frame
    glm::vec3 camera_local{0.0f};
    glm::vec3 sun_direction{0.0f, 0.0f, 1.0f}; // body frame, unit
    float equatorial_km = 1.0f;
    float polar_km = 1.0f;
    float km_per_radian = 1.0f;
    float band_depth = 0.0f;
    float band_south_rad = 0.0f;
    float band_north_rad = 0.0f;
    int site_count = 0;
    glm::vec4 sites[kMaxImpactSites] = {};  // east longitude, latitude (rad), age (s), depth at the darkest
    glm::vec4 shapes[kMaxImpactSites] = {}; // ejecta km, core km, expansion km/s, seed + 16 * rings
};

// One impact's flash and plume, in the site's local frame (km): x east, y north, z up
// (shaders/impact_plume.frag.hlsl).
struct ImpactPlumeDrawItem {
    glm::mat4 model{1.0f};
    glm::mat3 to_local{1.0f};
    glm::vec3 camera_local{0.0f};
    glm::vec3 box_min{0.0f};
    glm::vec3 box_max{0.0f};
    glm::vec3 planet_center{0.0f};
    float planet_radius = 1.0f;
    glm::vec3 sun_direction{0.0f, 0.0f, 1.0f};
    float sun_angular_radius = 0.0f;
    float age_s = 0.0f;
    float extinction_km2 = 0.0f; // per unit of the model's density; 0: only the flash
    glm::vec3 glow{0.0f};
    float flash_height_km = 0.0f;
    float flash_sigma_km = 1.0f;
    float flash_radiance = 0.0f;
    float flash_unhidden = 1.0f;
    glm::vec3 flash_color{1.0f};
};

class ImpactPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();
    // Recreates the pipelines for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format);

    // Right after the bodies: multiplies them.
    void draw_scars(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                    std::span<const ImpactScarDrawItem> items) const;
    void draw_plumes(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                     std::span<const ImpactPlumeDrawItem> items) const;

private:
    bool create_pipelines(const SceneTargetFormat& format);

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_scar_outside = nullptr;
    SDL_GPUGraphicsPipeline* m_scar_inside = nullptr;
    SDL_GPUGraphicsPipeline* m_plume_outside = nullptr;
    SDL_GPUGraphicsPipeline* m_plume_inside = nullptr;
    SDL_GPUBuffer* m_vertices = nullptr;
    uint32_t m_vertex_count = 0;
};

} // namespace astraxis
