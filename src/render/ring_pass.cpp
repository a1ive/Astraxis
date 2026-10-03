#include "render/ring_pass.hpp"

#include "core/math.hpp"
#include "render/gpu_buffer.hpp"
#include "render/shader.hpp"

#include <shaders/ring.frag.h>
#include <shaders/ring.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/vec2.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace astraxis {

namespace {

// Outer-edge sagitta at 280,000 km: 280000 * (1 - cos(pi / 512)) ~ 5 km.
constexpr int kSegments = 512;

struct RingVertex {
    glm::vec2 direction;
    float edge;
};

struct VertexUniforms {
    glm::mat4 model;
    glm::mat4 view_proj;
    glm::vec4 radii;
};

struct FragmentUniforms {
    glm::vec4 sun;
    glm::vec4 camera;
    glm::vec4 color;
    glm::vec4 params;
    glm::vec4 bands[kMaxRingBands];
};

} // namespace

bool RingPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;

    SDL_GPUShader* vs = create_shader(device, {.dxil = kRingVertDxil,
                                               .spirv = kRingVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kRingFragDxil,
                                               .spirv = kRingFragSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                               .num_uniform_buffers = 1});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, fs);
        return false;
    }

    SDL_GPUVertexBufferDescription vb = {};
    vb.slot = 0;
    vb.pitch = sizeof(RingVertex);
    vb.input_rate = SDL_GPU_VERTEXINPUTRATE_VERTEX;

    SDL_GPUVertexAttribute attrs[2] = {};
    attrs[0] = {.location = 0, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2,
                .offset = offsetof(RingVertex, direction)};
    attrs[1] = {.location = 1, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT,
                .offset = offsetof(RingVertex, edge)};

    // Scattered light only adds (tau << 1: no visible extinction).
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
    info.vertex_input_state.vertex_attributes = attrs;
    info.vertex_input_state.num_vertex_attributes = 2;
    info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
    info.rasterizer_state.fill_mode = SDL_GPU_FILLMODE_FILL;
    info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE; // seen from both sides
    info.multisample_state.sample_count = format.samples;
    info.depth_stencil_state.enable_depth_test = true;
    info.depth_stencil_state.enable_depth_write = false;
    info.depth_stencil_state.compare_op = SDL_GPU_COMPAREOP_GREATER; // reversed-Z
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = true;
    info.target_info.depth_stencil_format = format.depth;

    m_pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, fs);
    if (!m_pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Ring pipeline creation failed: %s", SDL_GetError());
        return false;
    }

    std::vector<RingVertex> vertices;
    vertices.reserve(kSegments * 6);
    for (int i = 0; i < kSegments; ++i) {
        const double a0 = kTwoPi * i / kSegments;
        const double a1 = kTwoPi * (i + 1) / kSegments;
        const glm::vec2 d0(static_cast<float>(std::cos(a0)), static_cast<float>(std::sin(a0)));
        const glm::vec2 d1(static_cast<float>(std::cos(a1)), static_cast<float>(std::sin(a1)));
        vertices.insert(vertices.end(), {{d0, 0.0f}, {d1, 0.0f}, {d1, 1.0f}, {d0, 0.0f}, {d1, 1.0f}, {d0, 1.0f}});
    }
    m_vertex_count = static_cast<uint32_t>(vertices.size());
    m_vertices = create_static_buffer(device, SDL_GPU_BUFFERUSAGE_VERTEX, vertices.data(),
                                      static_cast<uint32_t>(vertices.size() * sizeof(RingVertex)));
    return m_vertices != nullptr;
}

void RingPass::shutdown()
{
    if (!m_device) {
        return;
    }
    SDL_ReleaseGPUBuffer(m_device, m_vertices);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_pipeline);
    m_vertices = nullptr;
    m_pipeline = nullptr;
    m_device = nullptr;
}

void RingPass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                    std::span<const RingDrawItem> items) const
{
    if (items.empty()) {
        return;
    }

    SDL_BindGPUGraphicsPipeline(pass, m_pipeline);
    SDL_GPUBufferBinding vb = {.buffer = m_vertices, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vb, 1);

    for (const RingDrawItem& item : items) {
        const int count = std::min(item.band_count, kMaxRingBands);
        if (count <= 0) {
            continue;
        }
        float inner = item.bands[0].x;
        float outer = item.bands[0].y;
        for (int k = 1; k < count; ++k) {
            inner = std::min(inner, item.bands[k].x);
            outer = std::max(outer, item.bands[k].y);
        }
        // Room for the soft edges (2% of a band's width) on both sides.
        const float pad = 0.05f * (outer - inner);

        VertexUniforms vu;
        vu.model = item.model;
        vu.view_proj = view.view_proj;
        vu.radii = glm::vec4(std::max(0.0f, inner - pad), outer + pad, 0.0f, 0.0f);
        SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));

        FragmentUniforms fu = {};
        fu.sun = glm::vec4(item.sun_direction, item.sun_angular_radius);
        fu.camera = glm::vec4(item.camera, 0.0f);
        fu.color = glm::vec4(item.color, item.gain);
        fu.params = glm::vec4(static_cast<float>(count), item.phase_g, item.equatorial_radius, item.polar_radius);
        for (int k = 0; k < count; ++k) {
            fu.bands[k] = item.bands[k];
        }
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));

        SDL_DrawGPUPrimitives(pass, m_vertex_count, 1, 0, 0);
    }
}

} // namespace astraxis
