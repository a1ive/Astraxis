#include "render/orbit_line_pass.hpp"

#include "render/shader.hpp"

#include <shaders/orbit_line.frag.h>
#include <shaders/orbit_line.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/mat4x4.hpp>

namespace astraxis {

namespace {

constexpr float kTailOpacity = 0.12f;

struct VertexUniforms {
    glm::mat4 view_proj;
    glm::vec4 viewport;
    glm::vec4 params;
};

struct FragmentUniforms {
    glm::vec4 color;
    glm::vec4 params;
};

} // namespace

bool OrbitLinePass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;
    m_points_buffer.init(device, SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ);

    SDL_GPUShader* vs = create_shader(device, {.dxil = kOrbitLineVertDxil,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_storage_buffers = 1,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kOrbitLineFragDxil,
                                               .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                               .num_uniform_buffers = 1});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, fs);
        return false;
    }

    SDL_GPUColorTargetDescription color = {};
    color.format = format.color;
    color.blend_state.enable_blend = true;
    color.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_SRC_ALPHA;
    color.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
    color.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
    color.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
    color.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;

    SDL_GPUGraphicsPipelineCreateInfo info = {};
    info.vertex_shader = vs;
    info.fragment_shader = fs;
    info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
    info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
    info.multisample_state.sample_count = format.samples;
    info.depth_stencil_state.enable_depth_test = true;
    info.depth_stencil_state.enable_depth_write = false;
    info.depth_stencil_state.compare_op = SDL_GPU_COMPAREOP_GREATER_OR_EQUAL; // reversed-Z
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = true;
    info.target_info.depth_stencil_format = format.depth;

    m_pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, fs);
    if (!m_pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Orbit line pipeline creation failed: %s", SDL_GetError());
        return false;
    }
    return true;
}

void OrbitLinePass::shutdown()
{
    if (!m_device) {
        return;
    }
    m_points_buffer.shutdown();
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_pipeline);
    m_pipeline = nullptr;
    m_device = nullptr;
}

void OrbitLinePass::begin()
{
    m_points.clear();
    m_lines.clear();
}

void OrbitLinePass::add_line(const std::vector<glm::vec4>& points, const glm::vec4& color)
{
    if (points.size() < 2) {
        return;
    }
    m_lines.push_back({static_cast<uint32_t>(m_points.size()), static_cast<uint32_t>(points.size()), color});
    m_points.insert(m_points.end(), points.begin(), points.end());
}

void OrbitLinePass::upload(SDL_GPUCommandBuffer* cmd)
{
    m_points_buffer.upload(cmd, m_points.data(), static_cast<uint32_t>(m_points.size() * sizeof(glm::vec4)));
}

void OrbitLinePass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                         float width_px) const
{
    if (m_lines.empty() || !m_points_buffer.buffer()) {
        return;
    }

    SDL_BindGPUGraphicsPipeline(pass, m_pipeline);
    SDL_GPUBuffer* storage = m_points_buffer.buffer();
    SDL_BindGPUVertexStorageBuffers(pass, 0, &storage, 1);

    for (const Line& line : m_lines) {
        VertexUniforms vu;
        vu.view_proj = view.view_proj;
        vu.viewport = glm::vec4(view.viewport, 0.0f, 0.0f);
        vu.params = glm::vec4(static_cast<float>(line.first), width_px, 0.0f, 0.0f);
        SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));

        FragmentUniforms fu;
        fu.color = line.color;
        fu.params = glm::vec4(width_px, kTailOpacity, 0.0f, 0.0f);
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));

        SDL_DrawGPUPrimitives(pass, (line.count - 1) * 6, 1, 0, 0);
    }
}

} // namespace astraxis
