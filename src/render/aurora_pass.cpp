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
    glm::vec4 noon;   // w = shell bottom
    glm::vec4 sun;    // w = shell top
    glm::vec4 low;    // w = time
    glm::vec4 red;    // w = pixel angle
    glm::vec4 blue;   // x = kR
    glm::vec4 layers;
    glm::vec4 misc;   // x = diffuse flux
    glm::vec4 oval[kAuroraOvalSamples];
};

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
    return m_vertices && m_indices && m_noise && m_sampler;
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
                                               .num_samplers = 1,
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
                      std::span<const AuroraDrawItem> items) const
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
        fu.noon = glm::vec4(item.noon, item.bottom);
        fu.sun = glm::vec4(item.sun, item.top);
        fu.low = glm::vec4(item.low_rgb, item.time);
        fu.red = glm::vec4(item.red_per_kr, item.pixel_angle);
        fu.blue = glm::vec4(item.blue_kr, 0.0f, 0.0f, 0.0f);
        fu.layers = item.layers;
        fu.misc = glm::vec4(item.diffuse_flux, 0.0f, 0.0f, 0.0f);
        for (int k = 0; k < kAuroraOvalSamples; ++k) {
            fu.oval[k] = item.oval[static_cast<size_t>(k)];
        }
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));
        SDL_GPUTextureSamplerBinding noise = {.texture = m_noise, .sampler = m_sampler};
        SDL_BindGPUFragmentSamplers(pass, 0, &noise, 1);
        SDL_DrawGPUIndexedPrimitives(pass, m_index_count, 1, 0, 0, 0);
    }
}

} // namespace astraxis
