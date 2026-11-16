#pragma once

#include "render/gpu_device.hpp"
#include "scene/camera.hpp"
#include "scene/habitable_zone.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <array>
#include <span>

namespace astraxis {

struct HabitableZoneDrawItem {
    glm::mat4 model{1.0f}; // camera-relative translation * plane axes (km per au)
    float inner_au = 0.0f; // the annulus drawn
    float outer_au = 0.0f;
    int star_count = 0;
    std::array<glm::vec4, kMaxHabitableZoneStars> stars{};   // xyz = in the plane frame (au)
    std::array<glm::vec4, kMaxHabitableZoneStars> weights{}; // HabitableZone::weights_au2
    float near_start_km = 0.0f; // camera distance over which the zone fades in
    float near_end_km = 0.0f;
};

// Draws habitable zones (scene/habitable_zone.hpp) as annuli in their planes,
// lit where each limit's sum of the stars' insolation puts them inside the zone:
// light added, depth-tested and not written.
class HabitableZonePass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();
    // Recreates the pipeline for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format);

    // Colors: linear radiance of the conservative and the optimistic zone.
    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
              std::span<const HabitableZoneDrawItem> items, const glm::vec3& conservative,
              const glm::vec3& optimistic) const;

private:
    bool create_pipeline(const SceneTargetFormat& format);

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
    SDL_GPUBuffer* m_vertices = nullptr;
    uint32_t m_vertex_count = 0;
};

} // namespace astraxis
