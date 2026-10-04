#include "render/plume_pass.hpp"

#include "render/gpu_buffer.hpp"
#include "render/shader.hpp"

#include <shaders/plume.frag.h>
#include <shaders/plume.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/common.hpp>
#include <glm/vector_relational.hpp>

#include <cstdint>
#include <vector>

namespace astraxis {

namespace {

// The camera counts as inside when within this fraction of the box size from
// it (the near plane would clip the front faces).
constexpr float kInsideMargin = 0.02f;

struct VertexUniforms {
    glm::mat4 model;
    glm::mat4 view_proj;
    glm::vec4 box_min;
    glm::vec4 box_max;
};

struct FragmentUniforms {
    glm::mat4 to_local;
    glm::vec4 camera;
    glm::vec4 sun;
    glm::vec4 box_min;
    glm::vec4 box_max;
    glm::vec4 planet;
    glm::vec4 shape;
    glm::vec4 albedo;
    glm::vec4 params;
    glm::vec4 occluders[kMaxPlumeOccluders];
    glm::vec4 jet_positions[kMaxPlumeDrawJets];
    glm::vec4 jet_axes[kMaxPlumeDrawJets];
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

bool PlumePass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;

    SDL_GPUShader* vs = create_shader(device, {.dxil = kPlumeVertDxil,
                                               .spirv = kPlumeVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kPlumeFragDxil,
                                               .spirv = kPlumeFragSpirv,
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
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Plume pipeline creation failed: %s", SDL_GetError());
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

void PlumePass::shutdown()
{
    if (!m_device) {
        return;
    }
    SDL_ReleaseGPUBuffer(m_device, m_vertices);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_outside);
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_inside);
    m_vertices = nullptr;
    m_outside = nullptr;
    m_inside = nullptr;
    m_device = nullptr;
}

void PlumePass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, float time,
                     std::span<const PlumeDrawItem> items) const
{
    if (items.empty()) {
        return;
    }
    SDL_GPUBufferBinding vb = {.buffer = m_vertices, .offset = 0};
    SDL_BindGPUVertexBuffers(pass, 0, &vb, 1);

    const SDL_GPUGraphicsPipeline* bound = nullptr;
    for (const PlumeDrawItem& item : items) {
        const glm::vec3 margin = (item.box_max - item.box_min) * kInsideMargin;
        const bool inside = glm::all(glm::greaterThan(item.camera_local, item.box_min - margin)) &&
                            glm::all(glm::lessThan(item.camera_local, item.box_max + margin));
        SDL_GPUGraphicsPipeline* pipeline = inside ? m_inside : m_outside;
        if (pipeline != bound) {
            SDL_BindGPUGraphicsPipeline(pass, pipeline);
            bound = pipeline;
        }

        VertexUniforms vu;
        vu.model = item.model;
        vu.view_proj = view.view_proj;
        vu.box_min = glm::vec4(item.box_min, 0.0f);
        vu.box_max = glm::vec4(item.box_max, 0.0f);
        SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));

        FragmentUniforms fu = {};
        fu.to_local = glm::mat4(item.to_local);
        fu.camera = glm::vec4(item.camera_local, inside ? 1.0f : 0.0f);
        fu.sun = glm::vec4(item.sun_direction, item.sun_angular_radius);
        fu.box_min = glm::vec4(item.box_min, static_cast<float>(item.type));
        fu.box_max = glm::vec4(item.box_max, static_cast<float>(item.jet_count));
        fu.planet = glm::vec4(item.planet_center, item.planet_radius);
        fu.shape = item.shape;
        fu.albedo = glm::vec4(item.albedo, item.g);
        fu.params = glm::vec4(time, static_cast<float>(item.occluder_count), item.density, 0.0f);
        for (int k = 0; k < item.occluder_count && k < kMaxPlumeOccluders; ++k) {
            fu.occluders[k] = item.occluders[k];
        }
        for (int k = 0; k < item.jet_count && k < kMaxPlumeDrawJets; ++k) {
            fu.jet_positions[k] = item.jet_positions[k];
            fu.jet_axes[k] = item.jet_axes[k];
        }
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));

        SDL_DrawGPUPrimitives(pass, m_vertex_count, 1, 0, 0);
    }
}

} // namespace astraxis
