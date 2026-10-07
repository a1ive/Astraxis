#include "render/starfield_pass.hpp"

#include "core/color.hpp"
#include "core/math.hpp"
#include "render/gpu_buffer.hpp"
#include "render/shader.hpp"

#include <shaders/fullscreen.vert.h>
#include <shaders/milky_way.frag.h>
#include <shaders/starfield.frag.h>
#include <shaders/starfield.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace astraxis {

namespace {

constexpr uint32_t kStarCount = 9000;

// North galactic pole in ICRS (ESA, The Hipparcos and Tycho Catalogues, SP-1200, vol. 1, sec. 1.5.3).
constexpr double kNgpRaDeg = 192.85948;
constexpr double kNgpDecDeg = 27.12825;

struct StarInstance {
    float dir[3];
    float size_px;
    float color[3];
    float intensity;
};

struct Uniforms {
    glm::mat4 view_proj;
    glm::vec4 viewport;
    glm::vec4 params;
};

// Catalog stars: brightness follows the magnitude scale with a compressed
// dynamic range (exponent kMagnitudeGamma on the flux ratio), so that both
// Sirius (-1.46) and naked-eye-limit stars (kNakedEyeMag) are visible on a display.
constexpr double kMagnitudeGamma = 0.56;
constexpr double kFaintMag = kNakedEyeMag;
constexpr double kFaintIntensity = 0.07;
constexpr double kColorSaturation = 0.65; // stars look less saturated than blackbody color

std::vector<StarInstance> stars_from_catalog(std::span<const CatalogStar> catalog)
{
    std::vector<StarInstance> stars;
    stars.reserve(catalog.size());
    for (const CatalogStar& c : catalog) {
        const glm::dvec3 dir = unit_from_ra_dec(c.ra_deg * kDegToRad, c.dec_deg * kDegToRad);
        const double flux_ratio = std::pow(10.0, -0.4 * (c.vmag - kFaintMag));
        const double intensity = kFaintIntensity * std::pow(flux_ratio, kMagnitudeGamma);
        const glm::dvec3 bb = c.has_bv ? blackbody_linear_srgb(temperature_from_bv(c.bv)) : glm::dvec3(1.0);
        const glm::dvec3 color = glm::mix(glm::dvec3(1.0), bb, kColorSaturation);

        StarInstance s;
        s.dir[0] = static_cast<float>(dir.x);
        s.dir[1] = static_cast<float>(dir.y);
        s.dir[2] = static_cast<float>(dir.z);
        s.size_px = static_cast<float>(2.4 + 1.8 * std::clamp((4.5 - c.vmag) / 6.0, 0.0, 1.0));
        s.color[0] = static_cast<float>(color.x);
        s.color[1] = static_cast<float>(color.y);
        s.color[2] = static_cast<float>(color.z);
        s.intensity = static_cast<float>(intensity);
        stars.push_back(s);
    }
    return stars;
}

std::vector<StarInstance> generate_stars()
{
    std::mt19937 rng(20261002u);
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    const glm::dmat3 galactic = frame_from_pole(unit_from_ra_dec(kNgpRaDeg * kDegToRad, kNgpDecDeg * kDegToRad));

    // Magnitudes from a power law N(<m) ~ 10^(0.45 m) between kMinMag and kMaxMag.
    constexpr double kMinMag = -1.0;
    constexpr double kMaxMag = 7.5;
    constexpr double kSlope = 0.45;
    const double span = std::pow(10.0, kSlope * (kMaxMag - kMinMag)) - 1.0;

    std::vector<StarInstance> stars;
    stars.reserve(kStarCount);
    for (uint32_t i = 0; i < kStarCount; ++i) {
        glm::dvec3 dir;
        if (uni(rng) < 0.45) {
            // Disk population: galactic latitude bunched around 0.
            const double l = uni(rng) * kTwoPi;
            const double b = std::asin(2.0 * uni(rng) - 1.0) * 0.22;
            dir = galactic * glm::dvec3(std::cos(b) * std::cos(l), std::cos(b) * std::sin(l), std::sin(b));
        } else {
            const double z = 2.0 * uni(rng) - 1.0;
            const double phi = uni(rng) * kTwoPi;
            const double r = std::sqrt(1.0 - z * z);
            dir = {r * std::cos(phi), r * std::sin(phi), z};
        }

        const double mag = kMinMag + std::log10(1.0 + uni(rng) * span) / kSlope;
        const double flux = std::pow(10.0, -0.4 * (mag - 1.0));

        // Rough color-temperature spread: bluish white to orange.
        const double t = std::clamp(uni(rng) * 1.2 - 0.1, 0.0, 1.0);
        const glm::dvec3 blue(0.70, 0.80, 1.00);
        const glm::dvec3 white(1.00, 0.98, 0.95);
        const glm::dvec3 orange(1.00, 0.80, 0.58);
        const glm::dvec3 color = t < 0.5 ? glm::mix(blue, white, t * 2.0) : glm::mix(white, orange, t * 2.0 - 1.0);

        StarInstance s;
        s.dir[0] = static_cast<float>(dir.x);
        s.dir[1] = static_cast<float>(dir.y);
        s.dir[2] = static_cast<float>(dir.z);
        // Compressed dynamic range so that faint stars stay visible on typical displays.
        s.size_px = static_cast<float>(2.6 + 2.6 * std::min(1.0, std::sqrt(flux)));
        s.color[0] = static_cast<float>(color.x);
        s.color[1] = static_cast<float>(color.y);
        s.color[2] = static_cast<float>(color.z);
        s.intensity = static_cast<float>(std::clamp(std::pow(flux, 0.45) * 0.9, 0.16, 1.6));
        stars.push_back(s);
    }
    return stars;
}

// Additive point-sprite pipeline for a color target (optionally with a depth
// target it ignores, as in the scene pass).
SDL_GPUGraphicsPipeline* create_starfield_pipeline(SDL_GPUDevice* device, SDL_GPUTextureFormat color_format,
                                                   SDL_GPUSampleCount samples, SDL_GPUTextureFormat depth_format)
{
    SDL_GPUShader* vs = create_shader(device, {.dxil = kStarfieldVertDxil,
                                               .spirv = kStarfieldVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(
        device, {.dxil = kStarfieldFragDxil, .spirv = kStarfieldFragSpirv, .stage = SDL_GPU_SHADERSTAGE_FRAGMENT});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, fs);
        return nullptr;
    }

    SDL_GPUVertexBufferDescription vb = {};
    vb.slot = 0;
    vb.pitch = sizeof(StarInstance);
    vb.input_rate = SDL_GPU_VERTEXINPUTRATE_INSTANCE;

    SDL_GPUVertexAttribute attrs[2] = {};
    attrs[0] = {.location = 0, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT4, .offset = 0};
    attrs[1] = {.location = 1, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT4,
                .offset = offsetof(StarInstance, color)};

    SDL_GPUColorTargetDescription color = {};
    color.format = color_format;
    color.blend_state.enable_blend = true;
    color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ZERO;
    color.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    SDL_GPUGraphicsPipelineCreateInfo info = {};
    info.vertex_shader = vs;
    info.fragment_shader = fs;
    info.vertex_input_state.vertex_buffer_descriptions = &vb;
    info.vertex_input_state.num_vertex_buffers = 1;
    info.vertex_input_state.vertex_attributes = attrs;
    info.vertex_input_state.num_vertex_attributes = 2;
    info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
    info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
    info.multisample_state.sample_count = samples;
    info.depth_stencil_state.enable_depth_test = false;
    info.depth_stencil_state.enable_depth_write = false;
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = depth_format != SDL_GPU_TEXTUREFORMAT_INVALID;
    info.target_info.depth_stencil_format = depth_format;

    SDL_GPUGraphicsPipeline* pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, fs);
    if (!pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Starfield pipeline creation failed: %s", SDL_GetError());
    }
    return pipeline;
}

struct MilkyWayUniforms {
    glm::mat4 inv_view_proj;
    glm::vec4 params;
};

// Fullscreen Milky Way map, added to the target (depth ignored, as for the stars).
SDL_GPUGraphicsPipeline* create_milky_way_pipeline(SDL_GPUDevice* device, SDL_GPUTextureFormat color_format,
                                                   SDL_GPUSampleCount samples, SDL_GPUTextureFormat depth_format)
{
    SDL_GPUShader* vs = create_shader(
        device, {.dxil = kFullscreenVertDxil, .spirv = kFullscreenVertSpirv, .stage = SDL_GPU_SHADERSTAGE_VERTEX});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kMilkyWayFragDxil,
                                               .spirv = kMilkyWayFragSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                               .num_samplers = 1,
                                               .num_uniform_buffers = 1});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, fs);
        return nullptr;
    }

    SDL_GPUColorTargetDescription color = {};
    color.format = color_format;
    color.blend_state.enable_blend = true;
    color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ZERO;
    color.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    SDL_GPUGraphicsPipelineCreateInfo info = {};
    info.vertex_shader = vs;
    info.fragment_shader = fs;
    info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
    info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
    info.multisample_state.sample_count = samples;
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = depth_format != SDL_GPU_TEXTUREFORMAT_INVALID;
    info.target_info.depth_stencil_format = depth_format;

    SDL_GPUGraphicsPipeline* pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, fs);
    if (!pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Milky Way pipeline creation failed: %s", SDL_GetError());
    }
    return pipeline;
}

} // namespace

bool StarfieldPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format,
                         std::span<const CatalogStar> catalog)
{
    m_device = device;
    m_pipeline = create_starfield_pipeline(device, format.color, format.samples, format.depth);
    m_milky_way_pipeline = create_milky_way_pipeline(device, format.color, format.samples, format.depth);
    if (!m_pipeline || !m_milky_way_pipeline) {
        return false;
    }

    SDL_GPUSamplerCreateInfo sampler = {};
    sampler.min_filter = SDL_GPU_FILTER_LINEAR;
    sampler.mag_filter = SDL_GPU_FILTER_LINEAR;
    sampler.mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_LINEAR;
    sampler.address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_REPEAT; // right ascension wraps
    sampler.address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    sampler.address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    sampler.max_lod = 1000.0f;
    m_sampler = SDL_CreateGPUSampler(device, &sampler);
    if (!m_sampler) {
        return false;
    }

    return set_stars(catalog);
}

bool StarfieldPass::set_stars(std::span<const CatalogStar> catalog)
{
    if (m_instances) {
        SDL_ReleaseGPUBuffer(m_device, m_instances);
        m_instances = nullptr;
    }
    const std::vector<StarInstance> stars = catalog.empty() ? generate_stars() : stars_from_catalog(catalog);
    m_count = static_cast<uint32_t>(stars.size());
    m_instances = create_static_buffer(m_device, SDL_GPU_BUFFERUSAGE_VERTEX, stars.data(),
                                       static_cast<uint32_t>(stars.size() * sizeof(StarInstance)));
    return m_instances != nullptr;
}

void StarfieldPass::shutdown()
{
    if (!m_device) {
        return;
    }
    SDL_ReleaseGPUBuffer(m_device, m_instances);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_pipeline);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_milky_way_pipeline);
    SDL_ReleaseGPUSampler(m_device, m_sampler);
    m_instances = nullptr;
    m_pipeline = nullptr;
    m_milky_way_pipeline = nullptr;
    m_sampler = nullptr;
    m_milky_way = nullptr;
    m_device = nullptr;
}

bool StarfieldPass::render_cubemap(SDL_GPUTexture* cube, uint32_t size, SDL_GPUTextureFormat format,
                                   float brightness, float size_scale) const
{
    SDL_GPUGraphicsPipeline* pipeline =
        create_starfield_pipeline(m_device, format, SDL_GPU_SAMPLECOUNT_1, SDL_GPU_TEXTUREFORMAT_INVALID);
    SDL_GPUGraphicsPipeline* milky_way_pipeline =
        m_milky_way ? create_milky_way_pipeline(m_device, format, SDL_GPU_SAMPLECOUNT_1, SDL_GPU_TEXTUREFORMAT_INVALID)
                    : nullptr;
    if (!pipeline) {
        SDL_ReleaseGPUGraphicsPipeline(m_device, milky_way_pipeline);
        return false;
    }

    // Direct3D cube face order and orientation (+X, -X, +Y, -Y, +Z, -Z): face
    // "right" is cross(up, forward), i.e. the mirror of a right-handed camera,
    // hence the negated x scale of the projection.
    const glm::vec3 forward[6] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    const glm::vec3 up[6] = {{0, 1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}, {0, 1, 0}, {0, 1, 0}};
    glm::mat4 proj(0.0f);
    proj[0][0] = -1.0f; // 90 degree field of view, mirrored
    proj[1][1] = 1.0f;
    proj[2][3] = -1.0f;
    proj[3][2] = 1.0f;

    SDL_GPUCommandBuffer* cmd = SDL_AcquireGPUCommandBuffer(m_device);
    for (uint32_t face = 0; face < 6; ++face) {
        SDL_GPUColorTargetInfo target = {};
        target.texture = cube;
        target.layer_or_depth_plane = face;
        target.clear_color = {0.0f, 0.0f, 0.0f, 1.0f};
        target.load_op = SDL_GPU_LOADOP_CLEAR;
        target.store_op = SDL_GPU_STOREOP_STORE;
        SDL_GPURenderPass* pass = SDL_BeginGPURenderPass(cmd, &target, 1, nullptr);

        Uniforms u;
        u.view_proj = proj * glm::lookAt(glm::vec3(0.0f), forward[face], up[face]);
        if (milky_way_pipeline) {
            draw_milky_way(cmd, pass, milky_way_pipeline, u.view_proj, 2.0f / static_cast<float>(size));
        }
        u.viewport = glm::vec4(static_cast<float>(size), static_cast<float>(size), 0.0f, 0.0f);
        u.params = glm::vec4(brightness, size_scale, 0.0f, 0.0f);
        SDL_BindGPUGraphicsPipeline(pass, pipeline);
        SDL_PushGPUVertexUniformData(cmd, 0, &u, sizeof(u));
        SDL_GPUBufferBinding binding = {.buffer = m_instances, .offset = 0};
        SDL_BindGPUVertexBuffers(pass, 0, &binding, 1);
        SDL_DrawGPUPrimitives(pass, 6, m_count, 0, 0);
        SDL_EndGPURenderPass(pass);
    }
    SDL_SubmitGPUCommandBuffer(cmd);
    SDL_ReleaseGPUGraphicsPipeline(m_device, pipeline);
    SDL_ReleaseGPUGraphicsPipeline(m_device, milky_way_pipeline);
    return true;
}

void StarfieldPass::set_milky_way(SDL_GPUTexture* map, uint32_t map_width, float brightness)
{
    m_milky_way = map;
    m_milky_way_width = map_width;
    m_milky_way_brightness = brightness;
}

void StarfieldPass::draw_milky_way(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass,
                                   SDL_GPUGraphicsPipeline* pipeline, const glm::mat4& view_proj,
                                   float pixel_angle) const
{
    // Mip level at which a texel (2 pi / width at the equator) matches a pixel.
    const float texel_angle = static_cast<float>(kTwoPi) / static_cast<float>(std::max(m_milky_way_width, 1u));
    MilkyWayUniforms u;
    u.inv_view_proj = glm::inverse(view_proj);
    u.params = glm::vec4(m_milky_way_brightness, std::max(0.0f, std::log2(pixel_angle / texel_angle)), 0.0f, 0.0f);

    SDL_BindGPUGraphicsPipeline(pass, pipeline);
    const SDL_GPUTextureSamplerBinding binding = {.texture = m_milky_way, .sampler = m_sampler};
    SDL_BindGPUFragmentSamplers(pass, 0, &binding, 1);
    SDL_PushGPUFragmentUniformData(cmd, 0, &u, sizeof(u));
    SDL_DrawGPUPrimitives(pass, 3, 1, 0, 0);
}

void StarfieldPass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                         float brightness) const
{
    if (m_milky_way && m_milky_way_brightness > 0.0f) {
        draw_milky_way(cmd, pass, m_milky_way_pipeline, view.sky_view_proj, 2.0f / (view.focal_y * view.viewport.y));
    }
    if (brightness <= 0.0f) {
        return;
    }
    Uniforms u;
    u.view_proj = view.sky_view_proj;
    u.viewport = glm::vec4(view.viewport, 0.0f, 0.0f);
    u.params = glm::vec4(brightness, 1.0f, 0.0f, 0.0f);

    SDL_BindGPUGraphicsPipeline(pass, m_pipeline);
    SDL_PushGPUVertexUniformData(cmd, 0, &u, sizeof(u));
    SDL_GPUBufferBinding binding = {.buffer = m_instances, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &binding, 1);
    SDL_DrawGPUPrimitives(pass, 6, m_count, 0, 0);
}

} // namespace astraxis
