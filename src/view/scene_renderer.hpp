#pragma once

#include "render/atmosphere_pass.hpp"
#include "render/beam_pass.hpp"
#include "render/belt_pass.hpp"
#include "render/black_hole_pass.hpp"
#include "render/body_pass.hpp"
#include "render/comet_pass.hpp"
#include "render/dust_pass.hpp"
#include "render/orbit_line_pass.hpp"
#include "render/plume_pass.hpp"
#include "render/post_process.hpp"
#include "render/ring_pass.hpp"
#include "render/starfield_pass.hpp"
#include "render/sun_pass.hpp"
#include "scene/comet.hpp"
#include "scene/star_catalog.hpp"
#include "view/view_options.hpp"

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <cstdint>
#include <filesystem>
#include <vector>

namespace astraxis {

class GpuDevice;
class Scene;
class Simulation;
struct Frame;
struct OutputView;

// Draws a Simulation: turns its scene into the passes' draw items, renders
// them into the frame's scene targets and tonemaps the result into its
// swapchain texture. Owns the passes and each scene's GPU resources; knows
// nothing of windows, input or UI. One SceneRenderer can draw to several
// outputs of the same swapchain format, each with its own OutputView.
class SceneRenderer {
public:
    // `output_format`: the outputs' swapchain format.
    bool init(const GpuDevice& gpu, SDL_GPUTextureFormat output_format, const std::filesystem::path& asset_dir);
    void shutdown();

    // Replaces the GPU resources of the previous scene (textures, meshes, belts,
    // sky). `output_height` (pixels) sizes the stars of the black-hole lensing sky.
    void load_scene(const Scene& scene, const ViewOptions& options, uint32_t output_height);

    // Records the frame (its swapchain texture must be present). The host may
    // draw an overlay on the swapchain texture afterwards.
    void render(const Frame& frame, const Simulation& sim, const OutputView& view, const ViewOptions& options);

private:
    // Renders the background (stars, Milky Way) into the lensing cube map,
    // creating it on first use.
    void update_sky_cube(uint32_t output_height, float star_brightness);
    void create_sky_cube();
    void release_body_textures();
    // Traces the (first visible) black hole off-screen; composited in the scene pass.
    void trace_black_hole(SDL_GPUCommandBuffer* cmd, uint32_t width, uint32_t height, const Simulation& sim,
                          const OutputView& view, BlackHoleTargets& targets);

    void build_body_items(const Simulation& sim, const OutputView& view, const ViewOptions& options);
    void build_orbit_lines(const Simulation& sim, const OutputView& view, const ViewOptions& options);
    void build_belt_items(const Simulation& sim, const OutputView& view, const ViewOptions& options);
    void build_comet_items(const Simulation& sim, const OutputView& view, const ViewOptions& options);
    // Star sprites and pulsar beams (fills m_beam_items).
    void draw_stars(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const Simulation& sim, const OutputView& view);

    const GpuDevice* m_device = nullptr;
    std::filesystem::path m_asset_dir;

    StarfieldPass m_starfield;
    BodyPass m_bodies;
    AtmospherePass m_atmospheres;
    PlumePass m_plumes;
    CometPass m_comets;
    DustPass m_dust;
    RingPass m_rings;
    OrbitLinePass m_orbits;
    BeltPass m_belts;
    SunPass m_sun;
    BeamPass m_beams;
    BlackHolePass m_black_hole;
    PostProcess m_post;

    SDL_GPUTexture* m_sky_cube = nullptr; // starfield cube map for lensing (created on demand)
    SDL_GPUTextureFormat m_sky_cube_format = SDL_GPU_TEXTUREFORMAT_INVALID;
    uint32_t m_sky_cube_size = 0;
    SDL_GPUTexture* m_milky_way = nullptr; // the scene's Milky Way map, if any
    std::string m_milky_way_path;          // and its asset path (kept across scenes that share it)
    uint32_t m_milky_way_width = 0;
    std::vector<CatalogStar> m_catalog;   // the catalog sky (seen from Earth)
    bool m_scene_sky_stars = false;       // the starfield shows the scene's own stars

    // Per scene.
    std::vector<SDL_GPUTexture*> m_body_textures; // per body, may be null
    std::vector<SDL_GPUTexture*> m_ring_textures; // per body: ring profile, may be null
    std::vector<int> m_body_meshes;               // per body: BodyPass mesh index (-1: ellipsoid)
    std::vector<DustTail> m_dust_tails;           // per body (comets only are updated)
    std::vector<int> m_belt_ids;                  // per scene belt: BeltPass index (-1: not uploaded)
    std::vector<float> m_belt_reference_h;        // per scene belt: median absolute magnitude
    std::vector<double> m_belt_radius_km;         // per scene belt: median semi-major axis

    // Per-frame scratch.
    std::vector<BodyDrawItem> m_body_items;
    std::vector<RingDrawItem> m_ring_items;
    std::vector<AtmosphereDrawItem> m_atmosphere_items;
    std::vector<PlumeDrawItem> m_plume_items;
    std::vector<CometDrawItem> m_comet_items;
    std::vector<AtmosphereOptics> m_atmosphere_optics; // per body (BodyDrawItem::atmosphere points here)
    std::vector<BeltDrawItem> m_belt_items;
    std::vector<BeamDrawItem> m_beam_items;
    std::vector<glm::dvec3> m_trail_points;
    std::vector<float> m_trail_fades;
    std::vector<glm::vec4> m_line_points;
};

} // namespace astraxis
