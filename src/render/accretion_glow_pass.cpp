#include "render/accretion_glow_pass.hpp"

#include <shaders/accretion_glow.vert.h>

#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

namespace astraxis {

namespace {

constexpr float kMinSigmaPx = 0.7f; // in the off-screen target's pixels

struct VertexUniforms {
    glm::mat4 view_proj;
    glm::vec4 screen;
    glm::vec4 to_local[3];
    glm::vec4 center;
    glm::vec4 occluder_centers[kMaxAccretionGlowOccluders];
    glm::vec4 occluder_radii[kMaxAccretionGlowOccluders];
};
static_assert(sizeof(VertexUniforms) == 208, "VertexUniforms must match accretion_glow.vert.hlsl");

} // namespace

bool AccretionGlowPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    return m_splats.init(device, format,
                         {.dxil = kAccretionGlowVertDxil,
                          .spirv = kAccretionGlowVertSpirv,
                          .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                          .num_storage_buffers = 1,
                          .num_uniform_buffers = 1});
}

void AccretionGlowPass::render(SDL_GPUCommandBuffer* cmd, SplatTargets& targets, const CameraView& view,
                               uint32_t width, uint32_t height)
{
    targets.drawn = false;
    if (m_points.empty() || m_items.empty() ||
        !m_splats.stream_points(cmd, m_points.data(), static_cast<uint32_t>(m_points.size() * sizeof(DustSplat)),
                                static_cast<uint32_t>(m_points.size()))) {
        return;
    }
    const uint32_t w = SplatPass::target_size(width);
    const uint32_t h = SplatPass::target_size(height);
    std::vector<VertexUniforms> uniforms(m_items.size());
    std::vector<SplatPass::Draw> draws(m_items.size());
    for (size_t i = 0; i < m_items.size(); ++i) {
        const AccretionGlowItem& item = m_items[i];
        VertexUniforms& vu = uniforms[i];
        vu.view_proj = view.view_proj;
        // Splat sizes in the off-screen target's pixels (the projection itself is the scene's).
        vu.screen = glm::vec4(static_cast<float>(w), static_cast<float>(h),
                              0.5f * static_cast<float>(h) * view.focal_y, kMinSigmaPx);
        // glm is column-major: row r of the rotation is (m[0][r], m[1][r], m[2][r]).
        const glm::mat3& m = item.to_local;
        for (int r = 0; r < 3; ++r) {
            vu.to_local[r] = glm::vec4(m[0][r], m[1][r], m[2][r], 0.0f);
        }
        vu.to_local[0].w = static_cast<float>(item.first);
        vu.center = glm::vec4(item.center, static_cast<float>(item.occluder_count));
        for (int k = 0; k < kMaxAccretionGlowOccluders; ++k) {
            vu.occluder_centers[k] = glm::vec4(item.occluder_centers[k], 0.0f);
            vu.occluder_radii[k] = glm::vec4(item.occluder_radii[k], 1.0f);
        }
        draws[i] = {&vu, static_cast<uint32_t>(sizeof(VertexUniforms)), item.count};
    }
    m_splats.render(cmd, targets, draws.data(), draws.size(), width, height);
}

} // namespace astraxis
