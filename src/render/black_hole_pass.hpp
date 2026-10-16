#pragma once

#include "render/black_hole_uniforms.hpp"
#include "render/scene_targets.hpp"

#include <SDL3/SDL_gpu.h>

#include <cstdint>

namespace astraxis {

// Ray-traced black hole: shadow, lensed sky and thin accretion disk.
// trace() runs before the scene pass, into the output's off-screen targets
// (SceneTargets::black_hole) at a reduced resolution (rays are expensive);
// composite() then draws the result into the scene pass at full resolution:
// the sky sampled per pixel along interpolated directions, rays traced again
// where the prepass cannot be interpolated, and the depth the trace asked for.
class BlackHolePass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();
    // Recreates the pipelines for another scene format (e.g. MSAA sample count).
    bool set_format(const SceneTargetFormat& format);

    // The lensed sky, none of it owned: the cube texture (ICRF directions, with
    // mipmaps) used where lensing squeezes the sky; the stars and their grid
    // (StarfieldPass::sky_stars(), sky_cells()) drawn as points elsewhere; the
    // Milky Way map (equirectangular, mipmapped; null for none).
    void set_sky(SDL_GPUTexture* cube) { m_sky = cube; }
    void set_stars(SDL_GPUBuffer* stars, SDL_GPUBuffer* cells, SDL_GPUTexture* milky_way)
    {
        m_stars = stars;
        m_cells = cells;
        m_milky_way = milky_way;
    }
    bool ready() const { return m_sky && m_stars && m_cells; }

    // Outside any render pass. width/height: scene size in pixels.
    void trace(SDL_GPUCommandBuffer* cmd, uint32_t width, uint32_t height, const BlackHoleUniforms& u,
               BlackHoleTargets& targets);
    // Inside the scene pass; draws the last trace (no-op if none this frame).
    void composite(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, BlackHoleTargets& targets);

private:
    bool create_pipelines(const SceneTargetFormat& format);

    bool ensure_targets(BlackHoleTargets& targets, uint32_t width, uint32_t height);

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_trace_pipeline = nullptr;
    SDL_GPUGraphicsPipeline* m_composite_pipeline = nullptr;
    SDL_GPUSampler* m_linear = nullptr;
    SDL_GPUSampler* m_point = nullptr;
    SDL_GPUSampler* m_sky_sampler = nullptr;       // trilinear, anisotropic
    SDL_GPUSampler* m_milky_way_sampler = nullptr; // the same, wrapping in right ascension
    SDL_GPUTexture* m_blackbody = nullptr;         // temperature -> color lookup
    SDL_GPUTexture* m_sky = nullptr;
    SDL_GPUTexture* m_milky_way = nullptr;
    SDL_GPUBuffer* m_stars = nullptr;
    SDL_GPUBuffer* m_cells = nullptr;
};

} // namespace astraxis
