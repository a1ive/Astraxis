#include "render/aurora_pass.hpp"

#include "core/math.hpp"
#include "render/gpu_buffer.hpp"
#include "render/shader.hpp"
#include "render/texture.hpp"

#include <shaders/atmosphere.vert.h>
#include <shaders/aurora.frag.h>

#include <SDL3/SDL_log.h>
#include <glm/geometric.hpp>

#include <cmath>
#include <cstdint>
#include <vector>

namespace astraxis {

namespace {

// The shell's facets: the sagitta at their edges, 1 - cos(pi / 128) = 3e-4
// radii (2 km), is far below the shell's top.
constexpr int kSlices = 128;
constexpr int kStacks = 64;
constexpr float kInsideMargin = 1.02f;

struct VertexUniforms {
    glm::mat4 model;
    glm::mat4 view_proj;
    glm::vec4 top;
};

struct FragmentUniforms {
    glm::mat4 to_unit;
    glm::vec4 camera; // w = 1 inside the shell
    glm::vec4 pole;   // w = km per unit
    glm::vec4 noon;
    glm::vec4 sun;    // w = angular radius
    glm::vec4 low;    // w = time
    glm::vec4 high;
    glm::vec4 layers;
    glm::vec4 shell;  // bottom, top, polar ratio, pixel angle
    glm::vec4 params; // blue kR, diffuse flux, map gain, sigma km
    glm::ivec4 mode;  // pattern, vertical, lighting
    glm::vec4 parent;
    glm::vec4 oval[kAuroraOvalSamples];
    glm::vec4 south_center;
    glm::vec4 south_east;
    glm::vec4 jovian_oval[kJovianOvalPairs];
    glm::vec4 footprints[6];
    glm::vec4 table[kAuroraTableVectors];
    glm::vec4 ring_radii;
    glm::vec4 map_sun;
};
// SDL's Vulkan backend binds each uniform buffer with a 4 KiB range
// (MAX_UBO_SECTION_SIZE); reads past it are undefined.
static_assert(sizeof(FragmentUniforms) <= 4096);

// Four independent tiling value noises, one per channel: random values on a
// kNoiseCells x kNoiseCells lattice (wrapping), smoothstep-interpolated over
// the kNoiseSize x kNoiseSize tile (shaders/aurora.frag.hlsl).
constexpr int kNoiseSize = 256;
constexpr int kNoiseCells = 32;

std::vector<uint8_t> make_noise()
{
    std::vector<float> lattice(kNoiseCells * kNoiseCells * 4);
    uint32_t state = 0x9e3779b9u;
    for (float& v : lattice) {
        // xorshift32
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        v = static_cast<float>(state >> 8) / static_cast<float>(1u << 24);
    }
    auto at = [&](int i, int j, int c) {
        return lattice[static_cast<size_t>(((j % kNoiseCells) * kNoiseCells + (i % kNoiseCells)) * 4 + c)];
    };
    std::vector<uint8_t> rgba(static_cast<size_t>(kNoiseSize * kNoiseSize * 4));
    constexpr float kPixelsPerCell = static_cast<float>(kNoiseSize) / kNoiseCells;
    for (int y = 0; y < kNoiseSize; ++y) {
        for (int x = 0; x < kNoiseSize; ++x) {
            const float fx = (static_cast<float>(x) + 0.5f) / kPixelsPerCell;
            const float fy = (static_cast<float>(y) + 0.5f) / kPixelsPerCell;
            const int i = static_cast<int>(fx);
            const int j = static_cast<int>(fy);
            float sx = fx - static_cast<float>(i);
            float sy = fy - static_cast<float>(j);
            sx = sx * sx * (3.0f - 2.0f * sx);
            sy = sy * sy * (3.0f - 2.0f * sy);
            for (int c = 0; c < 4; ++c) {
                const float a = at(i, j, c) + sx * (at(i + 1, j, c) - at(i, j, c));
                const float b = at(i, j + 1, c) + sx * (at(i + 1, j + 1, c) - at(i, j + 1, c));
                const float v = a + sy * (b - a);
                rgba[static_cast<size_t>((y * kNoiseSize + x) * 4 + c)] =
                    static_cast<uint8_t>(std::lround(v * 255.0f));
            }
        }
    }
    return rgba;
}

SDL_GPUGraphicsPipeline* create_pipeline(SDL_GPUDevice* device, const SceneTargetFormat& format,
                                         SDL_GPUShader* vs, SDL_GPUShader* fs, bool inside)
{
    SDL_GPUVertexBufferDescription vb = {};
    vb.slot = 0;
    vb.pitch = sizeof(glm::vec3);
    vb.input_rate = SDL_GPU_VERTEXINPUTRATE_VERTEX;

    SDL_GPUVertexAttribute attr = {.location = 0, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT3,
                                   .offset = 0};

    // Emission adds to what lies behind.
    SDL_GPUColorTargetDescription color = {};
    color.format = format.color;
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
    info.vertex_input_state.vertex_attributes = &attr;
    info.vertex_input_state.num_vertex_attributes = 1;
    info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
    info.rasterizer_state.fill_mode = SDL_GPU_FILLMODE_FILL;
    info.rasterizer_state.cull_mode = inside ? SDL_GPU_CULLMODE_FRONT : SDL_GPU_CULLMODE_BACK;
    info.rasterizer_state.front_face = SDL_GPU_FRONTFACE_COUNTER_CLOCKWISE;
    info.multisample_state.sample_count = format.samples;
    info.depth_stencil_state.enable_depth_test = !inside;
    info.depth_stencil_state.enable_depth_write = false;
    info.depth_stencil_state.compare_op = SDL_GPU_COMPAREOP_GREATER; // reversed-Z
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = true;
    info.target_info.depth_stencil_format = format.depth;
    return SDL_CreateGPUGraphicsPipeline(device, &info);
}

} // namespace

bool AuroraPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;
    if (!create_pipelines(format)) {
        return false;
    }

    std::vector<glm::vec3> vertices;
    std::vector<uint32_t> indices;
    for (int j = 0; j <= kStacks; ++j) {
        const double lat = -kPi / 2.0 + kPi * j / kStacks;
        for (int i = 0; i < kSlices; ++i) {
            const double lon = kTwoPi * i / kSlices;
            vertices.emplace_back(static_cast<float>(std::cos(lat) * std::cos(lon)),
                                  static_cast<float>(std::cos(lat) * std::sin(lon)),
                                  static_cast<float>(std::sin(lat)));
        }
    }
    // Counter-clockwise when seen from outside.
    for (int j = 0; j < kStacks; ++j) {
        for (int i = 0; i < kSlices; ++i) {
            const uint32_t a = static_cast<uint32_t>(j * kSlices + i);
            const uint32_t b = static_cast<uint32_t>(j * kSlices + (i + 1) % kSlices);
            const uint32_t c = a + kSlices;
            const uint32_t d = b + kSlices;
            indices.insert(indices.end(), {a, b, d, a, d, c});
        }
    }
    m_index_count = static_cast<uint32_t>(indices.size());
    m_vertices = create_static_buffer(device, SDL_GPU_BUFFERUSAGE_VERTEX, vertices.data(),
                                      static_cast<uint32_t>(vertices.size() * sizeof(glm::vec3)));
    m_indices = create_static_buffer(device, SDL_GPU_BUFFERUSAGE_INDEX, indices.data(),
                                     static_cast<uint32_t>(indices.size() * sizeof(uint32_t)));

    const std::vector<uint8_t> noise = make_noise();
    m_noise = create_texture_rgba8(device, noise.data(), kNoiseSize, kNoiseSize, false);
    SDL_GPUSamplerCreateInfo sampler = {};
    sampler.min_filter = SDL_GPU_FILTER_LINEAR;
    sampler.mag_filter = SDL_GPU_FILTER_LINEAR;
    sampler.mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_NEAREST;
    sampler.address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_REPEAT;
    sampler.address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_REPEAT;
    sampler.address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_REPEAT;
    m_sampler = SDL_CreateGPUSampler(device, &sampler);
    sampler.mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_LINEAR;
    sampler.max_lod = 1000.0f; // left at zero, sampling stops at level 0
    sampler.address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    sampler.address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    sampler.address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    m_ring_sampler = SDL_CreateGPUSampler(device, &sampler);
    sampler.address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_REPEAT;
    m_map_sampler = SDL_CreateGPUSampler(device, &sampler);
    return m_vertices && m_indices && m_noise && m_sampler && m_ring_sampler && m_map_sampler;
}

bool AuroraPass::set_format(const SceneTargetFormat& format)
{
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_outside);
    m_outside = nullptr;
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_inside);
    m_inside = nullptr;
    return create_pipelines(format);
}

bool AuroraPass::create_pipelines(const SceneTargetFormat& format)
{
    SDL_GPUDevice* device = m_device;
    // The shell's vertex shader is the atmosphere's.
    SDL_GPUShader* vs = create_shader(device, {.dxil = kAtmosphereVertDxil,
                                               .spirv = kAtmosphereVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kAuroraFragDxil,
                                               .spirv = kAuroraFragSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                               .num_samplers = 3,
                                               .num_uniform_buffers = 1});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, fs);
        return false;
    }
    m_outside = create_pipeline(device, format, vs, fs, false);
    m_inside = create_pipeline(device, format, vs, fs, true);
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, fs);
    if (!m_outside || !m_inside) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Aurora pipeline creation failed: %s", SDL_GetError());
        return false;
    }
    return true;
}

void AuroraPass::shutdown()
{
    if (!m_device) {
        return;
    }
    SDL_ReleaseGPUBuffer(m_device, m_vertices);
    SDL_ReleaseGPUBuffer(m_device, m_indices);
    SDL_ReleaseGPUTexture(m_device, m_noise);
    SDL_ReleaseGPUSampler(m_device, m_sampler);
    SDL_ReleaseGPUSampler(m_device, m_ring_sampler);
    SDL_ReleaseGPUSampler(m_device, m_map_sampler);
    m_map_sampler = nullptr;
    m_ring_sampler = nullptr;
    m_noise = nullptr;
    m_sampler = nullptr;
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_outside);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_inside);
    m_vertices = nullptr;
    m_indices = nullptr;
    m_outside = nullptr;
    m_inside = nullptr;
    m_device = nullptr;
}

void AuroraPass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                      std::span<const AuroraDrawItem> items, bool after_rings) const
{
    if (items.empty()) {
        return;
    }
    SDL_GPUBufferBinding vb = {.buffer = m_vertices, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vb, 1);
    SDL_GPUBufferBinding ib = {.buffer = m_indices, .offset = 0};
    SDL_BindGPUIndexBuffer(pass, &ib, SDL_GPU_INDEXELEMENTSIZE_32BIT);

    const SDL_GPUGraphicsPipeline* bound = nullptr;
    for (const AuroraDrawItem& item : items) {
        if (item.after_rings != after_rings) continue;
        const bool inside = glm::length(item.camera_unit) < item.top * kInsideMargin;
        SDL_GPUGraphicsPipeline* pipeline = inside ? m_inside : m_outside;
        if (pipeline != bound) {
            SDL_BindGPUGraphicsPipeline(pass, pipeline);
            bound = pipeline;
        }

        VertexUniforms vu;
        vu.model = item.model;
        vu.view_proj = view.view_proj;
        vu.top = glm::vec4(item.top, 0.0f, 0.0f, 0.0f);
        SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));

        FragmentUniforms fu;
        fu.to_unit = glm::mat4(item.to_unit);
        fu.camera = glm::vec4(item.camera_unit, inside ? 1.0f : 0.0f);
        fu.pole = glm::vec4(item.pole, item.km_per_unit);
        fu.noon = glm::vec4(item.noon, 0.0f);
        fu.sun = glm::vec4(item.sun, item.sun_angular_radius);
        fu.low = glm::vec4(item.low_rgb, item.time);
        fu.high = glm::vec4(item.high_rgb, 0.0f);
        fu.layers = item.layers;
        fu.shell = glm::vec4(item.bottom, item.top, item.polar_ratio, item.pixel_angle);
        fu.params = glm::vec4(item.blue_kr, item.diffuse_flux, item.map_gain, item.sigma_km);
        fu.mode = glm::ivec4(item.pattern, item.vertical, item.lighting, 0);
        fu.parent = item.parent;
        for (int k = 0; k < kAuroraOvalSamples; ++k) {
            fu.oval[k] = item.oval[static_cast<size_t>(k)];
        }
        fu.south_center = glm::vec4(item.south_center, 0.0f);
        fu.south_east = glm::vec4(item.south_east, 0.0f);
        for (int k = 0; k < kJovianOvalPairs; ++k) {
            fu.jovian_oval[k] = item.jovian_oval[static_cast<size_t>(k)];
        }
        for (size_t k = 0; k < item.footprints.size(); ++k) fu.footprints[k] = item.footprints[k];
        for (size_t k = 0; k < item.table.size(); ++k) fu.table[k] = item.table[k];
        fu.ring_radii = glm::vec4(item.ring_radii, 0.0f, 0.0f);
        fu.map_sun = item.map_sun;
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));
        const SDL_GPUTextureSamplerBinding textures[] = {
            {.texture = m_noise, .sampler = m_sampler},
            {.texture = item.ring_profile ? item.ring_profile : m_noise, .sampler = m_ring_sampler},
            {.texture = item.map ? item.map : m_noise, .sampler = m_map_sampler}};
        SDL_BindGPUFragmentSamplers(pass, 0, textures, 3);
        SDL_DrawGPUIndexedPrimitives(pass, m_index_count, 1, 0, 0, 0);
    }
}

} // namespace astraxis
