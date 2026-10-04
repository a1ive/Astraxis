#include "render/atmosphere_pass.hpp"

#include "core/math.hpp"
#include "render/gpu_buffer.hpp"
#include "render/shader.hpp"

#include <shaders/atmosphere.frag.h>
#include <shaders/atmosphere.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/geometric.hpp>
#include <glm/vec4.hpp>

#include <cmath>
#include <cstdint>
#include <vector>

namespace astraxis {

namespace {

// Sagitta at the shell's facet edges: 1 - cos(pi / 128) = 3e-4 radii (2 km for
// the Earth), well above the base of any atmosphere's dense part.
constexpr int kSlices = 128;
constexpr int kStacks = 64;

// Beyond this (relative to the shell top) the camera is outside: its front
// faces are not clipped by the near plane.
constexpr float kInsideMargin = 1.02f;

struct VertexUniforms {
    glm::mat4 model;
    glm::mat4 view_proj;
    glm::vec4 top;
};

struct FragmentUniforms {
    glm::mat4 to_unit;
    glm::vec4 camera;
    glm::vec4 sun;
    glm::vec4 sun_unit;
    glm::vec4 rayleigh;
    glm::vec4 haze;
    glm::vec4 haze_scatter;
    glm::vec4 params;
};

SDL_GPUGraphicsPipeline* create_pipeline(SDL_GPUDevice* device, const SceneTargetFormat& format,
                                         SDL_GPUShader* vs, SDL_GPUShader* fs, bool inside)
{
    SDL_GPUVertexBufferDescription vb = {};
    vb.slot = 0;
    vb.pitch = sizeof(glm::vec3);
    vb.input_rate = SDL_GPU_VERTEXINPUTRATE_VERTEX;

    SDL_GPUVertexAttribute attr = {.location = 0, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT3,
                                   .offset = 0};

    // result = scattered light + transmission (alpha) * what lies behind.
    SDL_GPUColorTargetDescription color = {};
    color.format = format.color;
    color.blend_state.enable_blend = true;
    color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
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

bool AtmospherePass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;

    SDL_GPUShader* vs = create_shader(device, {.dxil = kAtmosphereVertDxil,
                                               .spirv = kAtmosphereVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kAtmosphereFragDxil,
                                               .spirv = kAtmosphereFragSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
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
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Atmosphere pipeline creation failed: %s", SDL_GetError());
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
    return m_vertices && m_indices;
}

void AtmospherePass::shutdown()
{
    if (!m_device) {
        return;
    }
    SDL_ReleaseGPUBuffer(m_device, m_vertices);
    SDL_ReleaseGPUBuffer(m_device, m_indices);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_outside);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_inside);
    m_vertices = nullptr;
    m_indices = nullptr;
    m_outside = nullptr;
    m_inside = nullptr;
    m_device = nullptr;
}

void AtmospherePass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                          std::span<const AtmosphereDrawItem> items) const
{
    if (items.empty()) {
        return;
    }
    SDL_GPUBufferBinding vb = {.buffer = m_vertices, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vb, 1);
    SDL_GPUBufferBinding ib = {.buffer = m_indices, .offset = 0};
    SDL_BindGPUIndexBuffer(pass, &ib, SDL_GPU_INDEXELEMENTSIZE_32BIT);

    const SDL_GPUGraphicsPipeline* bound = nullptr;
    for (const AtmosphereDrawItem& item : items) {
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

        const AtmosphereOptics& o = item.optics;
        FragmentUniforms fu;
        fu.to_unit = glm::mat4(item.to_unit);
        fu.camera = glm::vec4(item.camera_unit, inside ? 1.0f : 0.0f);
        fu.sun = glm::vec4(item.sun_direction, item.sun_angular_radius);
        fu.sun_unit = glm::vec4(glm::normalize(item.to_unit * item.sun_direction), item.top);
        fu.rayleigh = glm::vec4(o.rayleigh_depth, o.rayleigh_scale_height);
        fu.haze = glm::vec4(o.haze_attenuation, o.haze_scale_height);
        fu.haze_scatter = glm::vec4(o.haze_scattering, o.haze_g);
        fu.params = glm::vec4(item.gain, 0.0f, 0.0f, 0.0f);
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));

        SDL_DrawGPUIndexedPrimitives(pass, m_index_count, 1, 0, 0, 0);
    }
}

} // namespace astraxis
