#include "render/belt_pass.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "render/gpu_buffer.hpp"
#include "render/shader.hpp"

#include <shaders/belt.frag.h>
#include <shaders/belt.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

#include <cmath>

namespace astraxis {

namespace {

struct VertexUniforms {
    glm::mat4 view_proj;
    glm::vec4 display_x;
    glm::vec4 display_y;
    glm::vec4 display_z;
    glm::vec4 sun;
    glm::vec4 params;
    glm::vec4 viewport;
};

struct FragmentUniforms {
    glm::vec4 color;
};

} // namespace

bool BeltPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;
    SDL_GPUShader* vs = create_shader(device, {.dxil = kBeltVertDxil,
                                               .spirv = kBeltVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_storage_buffers = 1,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kBeltFragDxil,
                                               .spirv = kBeltFragSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                               .num_uniform_buffers = 1});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, fs);
        return false;
    }

    // Additive: overlapping points add up into the belt's glow.
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
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Belt pipeline creation failed: %s", SDL_GetError());
        return false;
    }
    return true;
}

void BeltPass::shutdown()
{
    if (!m_device) {
        return;
    }
    clear();
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_pipeline);
    m_pipeline = nullptr;
    m_device = nullptr;
}

void BeltPass::clear()
{
    for (Belt& b : m_belts) {
        SDL_ReleaseGPUBuffer(m_device, b.buffer);
    }
    m_belts.clear();
}

int BeltPass::add_belt(const float* elements, uint32_t count)
{
    Belt belt;
    belt.count = count;
    belt.buffer = create_static_buffer(m_device, SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ, elements,
                                       count * 7u * static_cast<uint32_t>(sizeof(float)));
    if (!belt.buffer) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Belt buffer creation failed: %s", SDL_GetError());
        return -1;
    }
    m_belts.push_back(belt);
    return static_cast<int>(m_belts.size()) - 1;
}

void BeltPass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view,
                    const std::vector<BeltDrawItem>& items) const
{
    if (items.empty()) {
        return;
    }
    SDL_BindGPUGraphicsPipeline(pass, m_pipeline);
    const float sqrt_gm = static_cast<float>(std::sqrt(kSunGmKm3S2 / (kAuKm * kAuKm * kAuKm)) * kSecondsPerDay); // au^1.5/day
    for (const BeltDrawItem& item : items) {
        if (item.belt < 0 || item.belt >= static_cast<int>(m_belts.size())) {
            continue;
        }
        const Belt& belt = m_belts[static_cast<size_t>(item.belt)];
        SDL_GPUBuffer* storage = belt.buffer;
        SDL_BindGPUVertexStorageBuffers(pass, 0, &storage, 1);

        // glm is column-major: row r of the rotation is (m[0][r], m[1][r], m[2][r]).
        const glm::mat3& m = item.display_from_ecliptic;
        VertexUniforms vu;
        vu.view_proj = view.view_proj;
        vu.display_x = glm::vec4(m[0][0], m[1][0], m[2][0], 0.0f);
        vu.display_y = glm::vec4(m[0][1], m[1][1], m[2][1], 0.0f);
        vu.display_z = glm::vec4(m[0][2], m[1][2], m[2][2], 0.0f);
        vu.sun = glm::vec4(item.sun_relative, 0.0f);
        vu.params = glm::vec4(item.days_since_epoch, item.point_size_px, sqrt_gm, static_cast<float>(kAuKm));
        vu.viewport = glm::vec4(view.viewport, item.reference_h, 0.0f);
        SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));

        FragmentUniforms fu;
        fu.color = glm::vec4(item.color, item.brightness);
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));

        SDL_DrawGPUPrimitives(pass, belt.count * 6u, 1, 0, 0);
    }
}

} // namespace astraxis
