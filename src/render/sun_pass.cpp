#include "render/sun_pass.hpp"

#include "render/shader.hpp"

#include <shaders/sun.frag.h>
#include <shaders/sun.vert.h>

#include <SDL3/SDL_log.h>
#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr float kMinDiskRadiusPx = 1.2f;
// Total disk "energy" (intensity x pixel area) kept constant when the disk is
// drawn larger than its true size, so a tiny sun does not flicker in and out.
constexpr float kDiskEnergy = 2400.0f;
constexpr float kMaxDiskIntensity = 400.0f;
constexpr float kMinDiskIntensity = 30.0f;
constexpr float kGlowIntensity = 0.8f;
constexpr float kSpriteRadii = 14.0f;
// Point-star scale shared with the background starfield (starfield_pass.cpp).
constexpr float kFaintMag = 6.5f;
constexpr float kFaintIntensity = 0.07f;
constexpr float kMagnitudeGamma = 0.56f;

struct VertexUniforms {
    glm::mat4 view_proj;
    glm::vec4 direction;
    glm::vec4 sprite;
};

struct FragmentUniforms {
    glm::vec4 params;
    glm::vec4 color;
};

} // namespace

bool SunPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;

    SDL_GPUShader* vs = create_shader(device, {.dxil = kSunVertDxil,
                                               .spirv = kSunVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kSunFragDxil,
                                               .spirv = kSunFragSpirv,
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
    // Only where nothing was drawn (depth still at the far value 0).
    info.depth_stencil_state.enable_depth_test = true;
    info.depth_stencil_state.enable_depth_write = false;
    info.depth_stencil_state.compare_op = SDL_GPU_COMPAREOP_GREATER_OR_EQUAL;
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = true;
    info.target_info.depth_stencil_format = format.depth;

    m_pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, fs);
    if (!m_pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Sun pipeline creation failed: %s", SDL_GetError());
        return false;
    }
    return true;
}

void SunPass::shutdown()
{
    if (m_device) {
        SDL_ReleaseGPUGraphicsPipeline(m_device, m_pipeline);
    }
    m_pipeline = nullptr;
    m_device = nullptr;
}

void SunPass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, const Star& star) const
{
    const float true_radius_px = star.angular_radius * view.focal_y * 0.5f * view.viewport.y;
    // Point-star size and brightness from the magnitude (as for catalog stars).
    const float point_size_px = 2.4f + 1.8f * std::clamp((4.5f - star.magnitude) / 6.0f, 0.0f, 1.0f);
    const float point_intensity =
        kFaintIntensity * std::pow(std::pow(10.0f, -0.4f * (star.magnitude - kFaintMag)), kMagnitudeGamma);
    const float radius_px = std::max({true_radius_px, kMinDiskRadiusPx, 0.5f * point_size_px});
    // Small disks keep a constant total brightness; large (nearby) disks keep a
    // constant surface brightness; faint stars stay on the point-star scale.
    const float glare_intensity =
        std::max(kMinDiskIntensity, std::min(kMaxDiskIntensity, kDiskEnergy / (radius_px * radius_px)));
    const float intensity = star.draw_disk ? std::min(point_intensity, glare_intensity) : 0.0f;
    // Glare only for stars much brighter than any night-sky star.
    const float glare = std::clamp((-3.0f - star.magnitude) / 9.0f, 0.0f, 1.0f);
    const float glow = kGlowIntensity * glare * glare;
    if (intensity <= 0.0f && glow <= 0.0f) {
        return;
    }

    VertexUniforms vu;
    vu.view_proj = view.view_proj;
    vu.direction = glm::vec4(star.position, star.finite ? 1.0f : 0.0f);
    vu.sprite = glm::vec4(radius_px * (glow > 0.0f ? kSpriteRadii : 1.5f), 0.0f, view.viewport.x, view.viewport.y);

    FragmentUniforms fu;
    fu.params = glm::vec4(radius_px, intensity, glow, 0.0f);
    fu.color = glm::vec4(star.color, 1.0f);

    SDL_BindGPUGraphicsPipeline(pass, m_pipeline);
    SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));
    SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));
    SDL_DrawGPUPrimitives(pass, 6, 1, 0, 0);
}

} // namespace astraxis
