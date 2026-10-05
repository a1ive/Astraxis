#pragma once

#include "render/gpu_device.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <cstdint>
#include <vector>

namespace astraxis {

// One belt to draw this frame (see shaders/belt.vert.hlsl).
struct BeltDrawItem {
    int belt = 0;                            // index returned by BeltPass::add_belt
    glm::vec3 sun_relative{0.0f};            // sun relative to the camera (km)
    glm::mat3 display_from_ecliptic{1.0f};   // J2000 ecliptic -> display frame
    float days_since_epoch = 0.0f;
    float point_size_px = 1.5f;
    glm::vec3 color{1.0f};                   // sRGB
    float brightness = 1.0f;
    float reference_h = 15.0f;               // absolute magnitude drawn at weight 1
};

// Point clouds of small bodies on two-body orbits, propagated on the GPU.
// add_belt() uploads a belt's elements once (when a scene is loaded).
class BeltPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    void clear();
    // `elements`: 7 floats per object (a au, e, i, node, arg_peri, mean anomaly rad, H).
    int add_belt(const float* elements, uint32_t count);

    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
              const std::vector<BeltDrawItem>& items) const;

private:
    struct Belt {
        SDL_GPUBuffer* buffer = nullptr;
        uint32_t count = 0;
    };

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
    std::vector<Belt> m_belts;
};

} // namespace astraxis
