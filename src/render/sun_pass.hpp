#pragma once

#include "render/renderer.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/vec3.hpp>

namespace astraxis {

// Stars as HDR sprites: a disk (when too small to be drawn as a sphere) plus
// glare for bright ones. Brightness follows the apparent bolometric magnitude
// on the same compressed scale as the background starfield, so the Sun seen
// from alpha Centauri is a bright star while seen from Jupiter it is blinding.
// Draw after bodies so that they occlude it through the depth test.
class SunPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    struct Star {
        glm::vec3 position{0.0f};  // camera-relative (km) when finite, else unit direction
        bool finite = true;
        float angular_radius = 0.0f; // radians
        glm::vec3 color{1.0f};       // linear, normalized
        float magnitude = 0.0f;      // apparent (bolometric-ish) magnitude
        bool draw_disk = true;       // false when the body pass draws the star as a sphere
    };
    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, const Star& star) const;

private:
    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
};

} // namespace astraxis
