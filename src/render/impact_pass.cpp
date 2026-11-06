#include "render/impact_pass.hpp"

#include "render/gpu_buffer.hpp"
#include "render/shader.hpp"

#include <shaders/impact_plume.frag.h>
#include <shaders/impact_scar.frag.h>
#include <shaders/plume.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/common.hpp>
#include <glm/vector_relational.hpp>

#include <cstdint>
#include <vector>

namespace astraxis {

namespace {

// As PlumePass: inside when within this fraction of the box size from it.
constexpr float kInsideMargin = 0.02f;

struct VertexUniforms {
    glm::mat4 model;
    glm::mat4 view_proj;
    glm::vec4 box_min;
    glm::vec4 box_max;
};

struct ScarUniforms {
    glm::mat4 to_local;
    glm::mat4 view_proj;
    glm::vec4 camera;
    glm::vec4 sun;
    glm::vec4 radii;
    glm::vec4 band;
    glm::vec4 debris;
    glm::vec4 sites[kMaxImpactSites];
    glm::vec4 shapes[kMaxImpactSites];
};
static_assert(sizeof(ScarUniforms) <= 4096, "the Vulkan backend binds 4 KiB of a uniform buffer");

struct PlumeUniforms {
    glm::mat4 to_local;
    glm::vec4 camera;
    glm::vec4 sun;
    glm::vec4 box_min;
    glm::vec4 box_max;
    glm::vec4 planet;
    glm::vec4 params;
    glm::vec4 flash;
    glm::vec4 flash_color;
    glm::vec4 glow;
};

enum class Kind { Scar, Plume };

SDL_GPUGraphicsPipeline* create_pipeline(SDL_GPUDevice* device, const SceneTargetFormat& format, SDL_GPUShader* vs,
                                         SDL_GPUShader* fs, Kind kind, bool inside)
{
    SDL_GPUVertexBufferDescription vb = {};
    vb.slot = 0;
    vb.pitch = sizeof(glm::vec3);
    vb.input_rate = SDL_GPU_VERTEXINPUTRATE_VERTEX;

    SDL_GPUVertexAttribute attr = {.location = 0, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT3,
                                   .offset = 0};

    SDL_GPUColorTargetDescription color = {};
    color.format = format.color;
    color.blend_state.enable_blend = true;
    if (kind == Kind::Scar) {
        // result = what lies behind * factor.
        color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ZERO;
        color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_SRC_COLOR;
    } else {
        // result = light + transmission (alpha) * what lies behind.
        color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
        color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    }
    color.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ZERO;
    color.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    SDL_GPUGraphicsPipelineCreateInfo info = {};
    info.vertex_shader = vs;
    info.fragment_shader = fs;
    info.vertex_input_state.vertex_buffer_descriptions = &vb;
    info.vertex_input_state.num_vertex_buffers = 1;
    info.vertex_input_state.vertex_attributes = &attr;
    info.vertex_input_state.num_vertex_attributes = 1;
    info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
    info.rasterizer_state.fill_mode = SDL_GPU_FILLMODE_FILL;
    info.rasterizer_state.cull_mode = inside ? SDL_GPU_CULLMODE_FRONT : SDL_GPU_CULLMODE_BACK;
    info.rasterizer_state.front_face = SDL_GPU_FRONTFACE_COUNTER_CLOCKWISE;
    info.multisample_state.sample_count = format.samples;
    // The scars write the depth of the surface they hit: always tested. A plume's box is tested
    // only from outside.
    info.depth_stencil_state.enable_depth_test = kind == Kind::Scar || !inside;
    info.depth_stencil_state.enable_depth_write = false;
    info.depth_stencil_state.compare_op = kind == Kind::Scar ? SDL_GPU_COMPAREOP_GREATER_OR_EQUAL
                                                             : SDL_GPU_COMPAREOP_GREATER; // reversed-Z
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = true;
    info.target_info.depth_stencil_format = format.depth;
    return SDL_CreateGPUGraphicsPipeline(device, &info);
}

} // namespace

bool ImpactPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;
    if (!create_pipelines(format)) {
        return false;
    }
    // Unit cube, counter-clockwise when seen from outside.
    const glm::vec3 c[8] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    const int faces[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {2, 3, 7, 6}, {1, 2, 6, 5}, {0, 4, 7, 3}};
    std::vector<glm::vec3> vertices;
    for (const auto& f : faces) {
        for (int k : {f[0], f[1], f[2], f[0], f[2], f[3]}) {
            vertices.push_back(c[k]);
        }
    }
    m_vertex_count = static_cast<uint32_t>(vertices.size());
    m_vertices = create_static_buffer(device, SDL_GPU_BUFFERUSAGE_VERTEX, vertices.data(),
                                      static_cast<uint32_t>(vertices.size() * sizeof(glm::vec3)));
    return m_vertices != nullptr;
}

bool ImpactPass::set_format(const SceneTargetFormat& format)
{
    for (SDL_GPUGraphicsPipeline** p : {&m_scar_outside, &m_scar_inside, &m_plume_outside, &m_plume_inside}) {
        SDL_ReleaseGPUGraphicsPipeline(m_device, *p);
        *p = nullptr;
    }
    return create_pipelines(format);
}

bool ImpactPass::create_pipelines(const SceneTargetFormat& format)
{
    SDL_GPUDevice* device = m_device;
    SDL_GPUShader* vs = create_shader(device, {.dxil = kPlumeVertDxil,
                                               .spirv = kPlumeVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* scar = create_shader(device, {.dxil = kImpactScarFragDxil,
                                                 .spirv = kImpactScarFragSpirv,
                                                 .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                                 .num_uniform_buffers = 1});
    SDL_GPUShader* plume = create_shader(device, {.dxil = kImpactPlumeFragDxil,
                                                  .spirv = kImpactPlumeFragSpirv,
                                                  .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                                  .num_uniform_buffers = 1});
    bool ok = vs && scar && plume;
    if (ok) {
        m_scar_outside = create_pipeline(device, format, vs, scar, Kind::Scar, false);
        m_scar_inside = create_pipeline(device, format, vs, scar, Kind::Scar, true);
        m_plume_outside = create_pipeline(device, format, vs, plume, Kind::Plume, false);
        m_plume_inside = create_pipeline(device, format, vs, plume, Kind::Plume, true);
        ok = m_scar_outside && m_scar_inside && m_plume_outside && m_plume_inside;
    }
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, scar);
    SDL_ReleaseGPUShader(device, plume);
    if (!ok) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Impact pipeline creation failed: %s", SDL_GetError());
    }
    return ok;
}

void ImpactPass::shutdown()
{
    if (!m_device) {
        return;
    }
    SDL_ReleaseGPUBuffer(m_device, m_vertices);
    for (SDL_GPUGraphicsPipeline** p : {&m_scar_outside, &m_scar_inside, &m_plume_outside, &m_plume_inside}) {
        SDL_ReleaseGPUGraphicsPipeline(m_device, *p);
        *p = nullptr;
    }
    m_vertices = nullptr;
    m_device = nullptr;
}

void ImpactPass::draw_scars(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                            std::span<const ImpactScarDrawItem> items) const
{
    if (items.empty()) {
        return;
    }
    SDL_GPUBufferBinding vb = {.buffer = m_vertices, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vb, 1);
    for (const ImpactScarDrawItem& item : items) {
        const glm::vec3 half(item.equatorial_km, item.equatorial_km, item.polar_km);
        const glm::vec3 margin = 2.0f * half * kInsideMargin;
        const bool inside = glm::all(glm::greaterThan(item.camera_local, -half - margin)) &&
                            glm::all(glm::lessThan(item.camera_local, half + margin));
        SDL_BindGPUGraphicsPipeline(pass, inside ? m_scar_inside : m_scar_outside);

        VertexUniforms vu;
        vu.model = item.model;
        vu.view_proj = view.view_proj;
        vu.box_min = glm::vec4(-half, 0.0f);
        vu.box_max = glm::vec4(half, 0.0f);
        SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));

        ScarUniforms fu = {};
        fu.to_local = glm::mat4(item.to_local);
        fu.view_proj = view.view_proj;
        fu.camera = glm::vec4(item.camera_local, 0.0f);
        fu.sun = glm::vec4(item.sun_direction, 0.0f);
        fu.radii = glm::vec4(item.equatorial_km, item.polar_km, static_cast<float>(item.site_count), item.km_per_radian);
        fu.band = glm::vec4(item.band_depth, item.band_south_rad, item.band_north_rad, 0.0f);
        fu.debris = glm::vec4(static_cast<float>(kDebrisRatioRed), static_cast<float>(kDebrisRatioGreen),
                              static_cast<float>(kDebrisRatioBlue), static_cast<float>(kDebrisDimming));
        for (int k = 0; k < item.site_count && k < kMaxImpactSites; ++k) {
            fu.sites[k] = item.sites[k];
            fu.shapes[k] = item.shapes[k];
        }
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));
        SDL_DrawGPUPrimitives(pass, m_vertex_count, 1, 0, 0);
    }
}

void ImpactPass::draw_plumes(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                             std::span<const ImpactPlumeDrawItem> items) const
{
    if (items.empty()) {
        return;
    }
    SDL_GPUBufferBinding vb = {.buffer = m_vertices, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vb, 1);
    for (const ImpactPlumeDrawItem& item : items) {
        const glm::vec3 margin = (item.box_max - item.box_min) * kInsideMargin;
        // A flash is hidden or not as a whole, and its box goes below the surface: as from inside, no depth test.
        const bool inside = item.flash_radiance > 0.0f ||
                            (glm::all(glm::greaterThan(item.camera_local, item.box_min - margin)) &&
                             glm::all(glm::lessThan(item.camera_local, item.box_max + margin)));
        SDL_BindGPUGraphicsPipeline(pass, inside ? m_plume_inside : m_plume_outside);

        VertexUniforms vu;
        vu.model = item.model;
        vu.view_proj = view.view_proj;
        vu.box_min = glm::vec4(item.box_min, 0.0f);
        vu.box_max = glm::vec4(item.box_max, 0.0f);
        SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));

        PlumeUniforms fu = {};
        fu.to_local = glm::mat4(item.to_local);
        fu.camera = glm::vec4(item.camera_local, inside ? 1.0f : 0.0f);
        fu.sun = glm::vec4(item.sun_direction, item.sun_angular_radius);
        fu.box_min = glm::vec4(item.box_min, 0.0f);
        fu.box_max = glm::vec4(item.box_max, 0.0f);
        fu.planet = glm::vec4(item.planet_center, item.planet_radius);
        fu.params = glm::vec4(item.age_s, item.extinction_km2, 0.0f, 0.0f);
        fu.flash_color = glm::vec4(item.flash_color, 0.0f);
        fu.glow = glm::vec4(item.glow, 0.0f);
        fu.flash = glm::vec4(item.flash_height_km, item.flash_sigma_km, item.flash_radiance, item.flash_unhidden);
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));
        SDL_DrawGPUPrimitives(pass, m_vertex_count, 1, 0, 0);
    }
}

} // namespace astraxis
