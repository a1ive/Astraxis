#pragma once

#include "render/renderer.hpp"
#include "scene/camera.hpp"
#include "scene/star_catalog.hpp"

#include <SDL3/SDL_gpu.h>

#include <span>

namespace astraxis {

// Background stars from a catalog (BSC5), drawn as soft dots at infinity.
// Without a catalog, falls back to procedural stars concentrated toward the
// galactic plane. Optionally a Milky Way map (equirectangular, ICRF) is drawn
// behind them.
class StarfieldPass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format, std::span<const CatalogStar> catalog);
    void shutdown();

    // Equirectangular ICRF map with mipmaps (not owned), or null for none.
    void set_milky_way(SDL_GPUTexture* map, uint32_t map_width, float brightness);

    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, float brightness) const;

    // Renders the stars into a cube texture (ICRF directions, Direct3D face
    // layout) for effects that look up the sky in arbitrary directions, such as
    // gravitational lensing. Submits its own command buffer.
    bool render_cubemap(SDL_GPUTexture* cube, uint32_t size, SDL_GPUTextureFormat format, float brightness,
                        float size_scale) const;

private:
    // pixel_angle: angular size of a target pixel (radians), selects the map's mip level.
    void draw_milky_way(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, SDL_GPUGraphicsPipeline* pipeline,
                        const glm::mat4& view_proj, float pixel_angle) const;

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
    SDL_GPUGraphicsPipeline* m_milky_way_pipeline = nullptr;
    SDL_GPUSampler* m_sampler = nullptr;
    SDL_GPUTexture* m_milky_way = nullptr;
    uint32_t m_milky_way_width = 0;
    float m_milky_way_brightness = 0.0f;
    SDL_GPUBuffer* m_instances = nullptr;
    uint32_t m_count = 0;
};

} // namespace astraxis
