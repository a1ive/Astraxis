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

inline constexpr int kMaxCometOccluders = 4;

// One glowing component of a comet in its local frame (km): the nucleus at the
// origin, +x along the ion tail. Values match the TYPE_* constants of
// shaders/comet.frag.hlsl.
struct CometDrawItem {
    enum Type : int { Coma = 0, IonTail = 1 };
    Type type = Coma;
    glm::mat4 model{1.0f};    // local frame -> camera-relative world (translation * rotation)
    glm::mat3 to_local{1.0f}; // world vector -> local frame (the rotation's transpose)
    glm::vec3 camera_local{0.0f};
    glm::vec3 box_min{0.0f};
    glm::vec3 box_max{0.0f};
    // Coma: parent and daughter scale lengths, nucleus radius; ion tail: 1/e
    // radius at the nucleus, tan(opening angle), e-folding length (km).
    glm::vec3 shape{0.0f};
    glm::vec3 color{1.0f};  // sRGB
    float brightness = 0.0f; // linear HDR radiance per unit column of the normalized density
    int occluder_count = 0;
    glm::vec4 occluders[kMaxCometOccluders] = {}; // local center, radius: the rays stop there
};

// Draws comae and ion tails as optically thin emission added to what lies
// behind (depth-tested, not written; from inside a box, rays start at the camera).
class CometPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    // `time`: seconds for the tails' slow animation.
    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, float time,
              std::span<const CometDrawItem> items) const;

private:
    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_outside = nullptr;
    SDL_GPUGraphicsPipeline* m_inside = nullptr;
    SDL_GPUBuffer* m_vertices = nullptr;
    uint32_t m_vertex_count = 0;
};

} // namespace astraxis
