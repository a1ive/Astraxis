#include "render/pinwheel_pass.hpp"

#include <shaders/pinwheel.vert.h>

#include <glm/mat4x4.hpp>

#include <cmath>

namespace astraxis {

namespace {

struct VertexUniforms {
    glm::mat4 view_proj;
    glm::mat4 view;
    glm::vec4 center;
    glm::vec4 display_x;
    glm::vec4 display_y;
    glm::vec4 display_z;
    glm::vec4 cycle;
    glm::vec4 shape;
    glm::vec4 light;
    glm::vec4 screen;
    std::array<glm::vec4, kPinwheelColorCount> colors;
    std::array<float, kPinwheelRadiusCount> radius;
};

} // namespace

bool PinwheelPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    return m_splats.init(device, format,
                         {.dxil = kPinwheelVertDxil,
                          .spirv = kPinwheelVertSpirv,
                          .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                          .num_storage_buffers = 1,
                          .num_uniform_buffers = 1});
}

bool PinwheelPass::set_grains(const PinwheelGpuGrain* grains, uint32_t count)
{
    return m_splats.set_points(grains, count * static_cast<uint32_t>(sizeof(PinwheelGpuGrain)), count);
}

void PinwheelPass::render(SDL_GPUCommandBuffer* cmd, SplatTargets& targets, const CameraView& view,
                          const PinwheelDrawItem& item, uint32_t width, uint32_t height)
{
    const uint32_t w = SplatPass::target_size(width);
    const uint32_t h = SplatPass::target_size(height);
    // glm is column-major: row r of the rotation is (m[0][r], m[1][r], m[2][r]).
    const glm::mat3& m = item.display_from_icrf;
    VertexUniforms vu;
    vu.view_proj = view.view_proj;
    vu.view = view.view;
    vu.center = glm::vec4(item.center, 0.0f);
    vu.display_x = glm::vec4(m[0][0], m[1][0], m[2][0], 0.0f);
    vu.display_y = glm::vec4(m[0][1], m[1][1], m[2][1], 0.0f);
    vu.display_z = glm::vec4(m[0][2], m[1][2], m[2][2], 0.0f);
    vu.cycle = glm::vec4(item.cycle, item.turns, item.reference_km, item.nucleation_km);
    vu.shape = glm::vec4(item.sigma_along, item.sigma_across, item.wall_sigma_km, item.min_sigma_px);
    const float log_span = std::log(item.radius.back() / item.nucleation_km);
    vu.light = glm::vec4(item.radiance, item.gamma, log_span > 0.0f ? 1.0f / log_span : 0.0f, item.cull_radiance);
    // Splat sizes in the off-screen target's pixels (the projection itself is the scene's).
    vu.screen = glm::vec4(static_cast<float>(w), static_cast<float>(h), 0.5f * static_cast<float>(h) * view.focal_y,
                          item.max_age);
    vu.colors = item.colors;
    vu.radius = item.radius;
    m_splats.render(cmd, targets, &vu, sizeof(vu), width, height);
}

} // namespace astraxis
