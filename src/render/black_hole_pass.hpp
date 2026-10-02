#pragma once

#include "render/renderer.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

#include <cstdint>

namespace astraxis {

// Per-frame inputs of shaders/black_hole.frag.hlsl (see there for the meaning).
// Everything is in the hole's frame (spin and disk normal along +z), in units of M.
struct BlackHoleUniforms {
    glm::vec4 camera{0.0f};  // xyz = camera position, w = spin a
    glm::vec4 right{0.0f};   // camera right * tan(fov/2) * aspect
    glm::vec4 up{0.0f};      // camera up * tan(fov/2)
    glm::vec4 forward{0.0f}; // camera forward
    glm::mat4 to_sky{1.0f};  // hole frame -> ICRF
    glm::vec4 radii{0.0f};   // r+, disk inner, disk outer, integration radius
    glm::vec4 params{0.0f};  // b_max, depth, disk brightness, disk time
    glm::vec4 disk{0.0f};    // peak temperature, 1 / profile max, max steps, edge fade start
};

// Ray-traced black hole: shadow, lensed sky and thin accretion disk.
// trace() runs before the scene pass, into off-screen targets at a reduced
// resolution (rays are expensive); composite() then draws the result into the
// scene pass with the depth the trace asked for.
class BlackHolePass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    // The lensed sky (cube texture in ICRF directions); not owned.
    void set_sky(SDL_GPUTexture* cube) { m_sky = cube; }
    bool ready() const { return m_sky != nullptr; }

    // Outside any render pass. width/height: scene size in pixels.
    void trace(SDL_GPUCommandBuffer* cmd, uint32_t width, uint32_t height, const BlackHoleUniforms& u);
    // Inside the scene pass; draws the last trace (no-op if none this frame).
    void composite(SDL_GPURenderPass* pass);

private:
    bool ensure_targets(uint32_t width, uint32_t height);
    void release_targets();

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_trace_pipeline = nullptr;
    SDL_GPUGraphicsPipeline* m_composite_pipeline = nullptr;
    SDL_GPUSampler* m_linear = nullptr;
    SDL_GPUSampler* m_point = nullptr;
    SDL_GPUTexture* m_blackbody = nullptr; // temperature -> color lookup
    SDL_GPUTexture* m_sky = nullptr;

    SDL_GPUTexture* m_color = nullptr; // premultiplied color + coverage
    SDL_GPUTexture* m_depth = nullptr; // R32F scene depth to write
    uint32_t m_width = 0;
    uint32_t m_height = 0;
    bool m_traced = false;
};

} // namespace astraxis
