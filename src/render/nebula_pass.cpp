#include "render/nebula_pass.hpp"

#include <shaders/nebula.vert.h>

#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

#include <array>
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
    glm::vec4 shape;
    glm::vec4 light;
    glm::vec4 color;
    glm::vec4 screen;
    glm::vec4 layer;
    glm::vec4 thermal;
    std::array<glm::vec4, kNebulaThermalCount> thermal_colors;
};

} // namespace

bool NebulaPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    return m_splats.init(device, format,
                         {.dxil = kNebulaVertDxil,
                          .spirv = kNebulaVertSpirv,
                          .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                          .num_storage_buffers = 1,
                          .num_uniform_buffers = 1});
}

bool NebulaPass::set_grains(const NebulaGpuGrain* grains, uint32_t count)
{
    return m_splats.set_points(grains, count * static_cast<uint32_t>(sizeof(NebulaGpuGrain)), count);
}

void NebulaPass::render(SDL_GPUCommandBuffer* cmd, SplatTargets& targets, const CameraView& view,
                        const std::vector<NebulaDrawItem>& items, uint32_t width, uint32_t height)
{
    // Splat sizes in the off-screen target's pixels (the projection itself is the scene's).
    const float h = static_cast<float>(SplatPass::target_size(height));
    const glm::vec4 screen_size(static_cast<float>(SplatPass::target_size(width)), h, 0.5f * h * view.focal_y, 0.0f);
    std::vector<VertexUniforms> uniforms(items.size());
    std::vector<SplatPass::Draw> draws(items.size());
    for (size_t k = 0; k < items.size(); ++k) {
        const NebulaDrawItem& item = items[k];
        // glm is column-major: row r of the rotation is (m[0][r], m[1][r], m[2][r]).
        const glm::mat3& m = item.display_from_icrf;
        VertexUniforms& vu = uniforms[k];
        vu.view_proj = view.view_proj;
        vu.view = view.view;
        vu.center = glm::vec4(item.center, item.scale);
        vu.display_x = glm::vec4(m[0][0], m[1][0], m[2][0], 0.0f);
        vu.display_y = glm::vec4(m[0][1], m[1][1], m[2][1], 0.0f);
        vu.display_z = glm::vec4(m[0][2], m[1][2], m[2][2], 0.0f);
        vu.shape = glm::vec4(item.patch_sigma_km, item.wall_sigma_km, item.min_sigma_px, item.asymmetry_g);
        vu.light = glm::vec4(item.radiance, item.gamma, item.reference_km, item.cull_radiance);
        vu.color = glm::vec4(item.color, static_cast<float>(item.light));
        vu.screen = screen_size;
        vu.screen.w = static_cast<float>(item.first);
        vu.layer = glm::vec4(item.radius_share, item.column, item.temperature_k, 1.0f / item.temperature_km);
        const double steps = (kNebulaThermalCount - 1) / std::log(m_thermal.hi_k / m_thermal.lo_k);
        vu.thermal = glm::vec4(item.temperature_index, static_cast<float>(std::log(m_thermal.lo_k)),
                               static_cast<float>(steps), 0.0f);
        vu.thermal_colors = m_thermal.colors;
        draws[k] = {&vu, static_cast<uint32_t>(sizeof(VertexUniforms)), item.count};
    }
    m_splats.render(cmd, targets, draws.data(), draws.size(), width, height);
}

} // namespace astraxis
