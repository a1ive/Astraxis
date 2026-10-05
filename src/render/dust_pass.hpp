#pragma once

#include "render/gpu_buffer.hpp"
#include "render/renderer.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/vec3.hpp>

#include <vector>

namespace astraxis {

// One dust grain (32 bytes, as read by dust.vert.hlsl).
struct DustSplat {
    glm::vec3 position{0.0f}; // camera-relative (km)
    float sigma_km = 0.0f;    // Gaussian standard deviation at the grain
    glm::vec3 color{1.0f};    // sRGB
    float luminosity = 0.0f;  // total flux x distance^2 (linear HDR, as CometDrawItem::brightness)
};
static_assert(sizeof(DustSplat) == 32, "DustSplat must match the shader's raw buffer loads");

// Draws comet dust grains as additive Gaussian splats (depth-tested, not written).
// Usage per frame: the splats -> upload() (outside render passes) -> draw().
class DustPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    std::vector<DustSplat>& splats() { return m_splats; }
    void upload(SDL_GPUCommandBuffer* cmd);
    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, float min_sigma_px) const;

private:
    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
    StreamBuffer m_buffer;
    std::vector<DustSplat> m_splats;
    uint32_t m_uploaded = 0;
};

} // namespace astraxis
