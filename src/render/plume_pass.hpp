#pragma once

#include "render/renderer.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat3x3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <span>

namespace astraxis {

inline constexpr int kMaxPlumeOccluders = 4;
inline constexpr int kMaxPlumeDrawJets = 12;

// One plume component in its local frame (km): origin at the vent on the
// surface, z up. Values match the TYPE_* constants of shaders/plume.frag.hlsl.
struct PlumeDrawItem {
    enum Type : int { Umbrella = 0, Jets = 1, Column = 2, Tail = 3 };
    Type type = Umbrella;
    glm::mat4 model{1.0f};    // local frame -> camera-relative world (translation * rotation)
    glm::mat3 to_local{1.0f}; // world vector -> local frame (the rotation's transpose)
    glm::vec3 camera_local{0.0f};
    glm::vec3 box_min{0.0f};
    glm::vec3 box_max{0.0f};
    glm::vec3 planet_center{0.0f}; // local
    float planet_radius = 1.0f;    // km, at the vent (curvature and shadow)
    // Umbrella: canopy height, canopy radius, shell thickness; jets: scale
    // height, tan(spread), base width; column: height, radius; tail: height,
    // radius, length, widening per km.
    glm::vec4 shape{0.0f};
    float density = 0.0f; // extinction coefficient scale (per km)
    glm::vec3 albedo{1.0f};
    float g = 0.0f;
    glm::vec3 sun_direction{0.0f, 0.0f, 1.0f}; // local, unit
    float sun_angular_radius = 0.0f;
    int occluder_count = 0;
    glm::vec4 occluders[kMaxPlumeOccluders] = {}; // local center, radius
    int jet_count = 0;
    glm::vec4 jet_positions[kMaxPlumeDrawJets] = {}; // local
    glm::vec4 jet_axes[kMaxPlumeDrawJets] = {};
};

// Draws volcanic and cryovolcanic plumes by single scattering inside their
// bounding boxes: scattered sunlight is added and what lies behind is dimmed by
// the view transmission (depth-tested, not written).
class PlumePass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    // `time`: seconds for the plumes' slow animation.
    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, float time,
              std::span<const PlumeDrawItem> items) const;

private:
    SDL_GPUDevice* m_device = nullptr;
    // Camera outside the box: front faces, depth-tested. Inside (or close):
    // back faces, not depth-tested, with the rays starting at the camera.
    SDL_GPUGraphicsPipeline* m_outside = nullptr;
    SDL_GPUGraphicsPipeline* m_inside = nullptr;
    SDL_GPUBuffer* m_vertices = nullptr;
    uint32_t m_vertex_count = 0;
};

} // namespace astraxis
