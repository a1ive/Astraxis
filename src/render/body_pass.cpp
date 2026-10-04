#include "render/body_pass.hpp"

#include "core/math.hpp"
#include "render/gpu_buffer.hpp"
#include "render/shader.hpp"
#include "render/texture.hpp"

#include <shaders/body.frag.h>
#include <shaders/body.vert.h>

#include <SDL3/SDL_log.h>

#include <cmath>
#include <cstdint>
#include <vector>

namespace astraxis {

namespace {

constexpr int kSlices = 128; // longitude segments
constexpr int kStacks = 64;  // latitude segments

struct VertexUniforms {
    glm::mat4 model;
    glm::mat4 view_proj;
    glm::mat4 rotation;
    glm::vec4 inv_scale; // w = map u of longitude 0 (turns)
};

struct FragmentUniforms {
    glm::vec4 sun;
    glm::vec4 light2;
    glm::vec4 light2_color;
    glm::vec4 color;
    glm::vec4 params;
    glm::vec4 occluders[kMaxOccluders];
    glm::vec4 ring_center;
    glm::vec4 ring_normal;
    glm::vec4 ring_radii;
    glm::vec4 atmo_rayleigh;
    glm::vec4 atmo_haze;
};

struct BodyVertex {
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 uv;
    float albedo;
};

void build_sphere(std::vector<BodyVertex>& vertices, std::vector<uint32_t>& indices)
{
    // The seam column is duplicated (i = 0 and i = kSlices) so u runs 0..1 without wrapping.
    for (int j = 0; j <= kStacks; ++j) {
        const double lat = -kPi / 2.0 + kPi * j / kStacks;
        for (int i = 0; i <= kSlices; ++i) {
            const double lon = kTwoPi * i / kSlices;
            const glm::vec3 p(static_cast<float>(std::cos(lat) * std::cos(lon)),
                              static_cast<float>(std::cos(lat) * std::sin(lon)), static_cast<float>(std::sin(lat)));
            const glm::vec2 uv(static_cast<float>(i) / kSlices, 1.0f - static_cast<float>(j) / kStacks);
            vertices.push_back({p, p, uv, 1.0f});
        }
    }
    // Counter-clockwise when seen from outside.
    for (int j = 0; j < kStacks; ++j) {
        for (int i = 0; i < kSlices; ++i) {
            const uint32_t a = static_cast<uint32_t>(j * (kSlices + 1) + i);
            const uint32_t b = a + 1;
            const uint32_t c = a + kSlices + 1;
            const uint32_t d = c + 1;
            indices.insert(indices.end(), {a, b, d, a, d, c});
        }
    }
}

} // namespace

bool BodyPass::init(SDL_GPUDevice* device, const SceneTargetFormat& format)
{
    m_device = device;

    SDL_GPUShader* vs = create_shader(device, {.dxil = kBodyVertDxil,
                                               .spirv = kBodyVertSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_VERTEX,
                                               .num_uniform_buffers = 1});
    SDL_GPUShader* fs = create_shader(device, {.dxil = kBodyFragDxil,
                                               .spirv = kBodyFragSpirv,
                                               .stage = SDL_GPU_SHADERSTAGE_FRAGMENT,
                                               .num_samplers = 2,
                                               .num_uniform_buffers = 1});
    if (!vs || !fs) {
        SDL_ReleaseGPUShader(device, vs);
        SDL_ReleaseGPUShader(device, fs);
        return false;
    }

    SDL_GPUVertexBufferDescription vb = {};
    vb.slot = 0;
    vb.pitch = sizeof(BodyVertex);
    vb.input_rate = SDL_GPU_VERTEXINPUTRATE_VERTEX;

    SDL_GPUVertexAttribute attrs[4] = {};
    attrs[0] = {.location = 0, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT3,
                .offset = offsetof(BodyVertex, position)};
    attrs[1] = {.location = 1, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT3,
                .offset = offsetof(BodyVertex, normal)};
    attrs[2] = {.location = 2, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT2,
                .offset = offsetof(BodyVertex, uv)};
    attrs[3] = {.location = 3, .buffer_slot = 0, .format = SDL_GPU_VERTEXELEMENTFORMAT_FLOAT,
                .offset = offsetof(BodyVertex, albedo)};

    SDL_GPUColorTargetDescription color = {};
    color.format = format.color;

    SDL_GPUGraphicsPipelineCreateInfo info = {};
    info.vertex_shader = vs;
    info.fragment_shader = fs;
    info.vertex_input_state.vertex_buffer_descriptions = &vb;
    info.vertex_input_state.num_vertex_buffers = 1;
    info.vertex_input_state.vertex_attributes = attrs;
    info.vertex_input_state.num_vertex_attributes = 4;
    info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
    info.rasterizer_state.fill_mode = SDL_GPU_FILLMODE_FILL;
    info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_BACK;
    info.rasterizer_state.front_face = SDL_GPU_FRONTFACE_COUNTER_CLOCKWISE;
    info.multisample_state.sample_count = format.samples;
    info.depth_stencil_state.enable_depth_test = true;
    info.depth_stencil_state.enable_depth_write = true;
    info.depth_stencil_state.compare_op = SDL_GPU_COMPAREOP_GREATER; // reversed-Z
    info.target_info.color_target_descriptions = &color;
    info.target_info.num_color_targets = 1;
    info.target_info.has_depth_stencil_target = true;
    info.target_info.depth_stencil_format = format.depth;

    m_pipeline = SDL_CreateGPUGraphicsPipeline(device, &info);
    SDL_ReleaseGPUShader(device, vs);
    SDL_ReleaseGPUShader(device, fs);
    if (!m_pipeline) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Body pipeline creation failed: %s", SDL_GetError());
        return false;
    }

    std::vector<BodyVertex> vertices;
    std::vector<uint32_t> indices;
    build_sphere(vertices, indices);
    m_sphere.index_count = static_cast<uint32_t>(indices.size());
    m_sphere.vertices = create_static_buffer(device, SDL_GPU_BUFFERUSAGE_VERTEX, vertices.data(),
                                             static_cast<uint32_t>(vertices.size() * sizeof(BodyVertex)));
    m_sphere.indices = create_static_buffer(device, SDL_GPU_BUFFERUSAGE_INDEX, indices.data(),
                                            static_cast<uint32_t>(indices.size() * sizeof(uint32_t)));

    SDL_GPUSamplerCreateInfo sampler = {};
    sampler.min_filter = SDL_GPU_FILTER_LINEAR;
    sampler.mag_filter = SDL_GPU_FILTER_LINEAR;
    sampler.mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_LINEAR;
    sampler.address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_REPEAT;
    sampler.address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    sampler.address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    sampler.max_anisotropy = 8.0f;
    sampler.enable_anisotropy = true;
    sampler.max_lod = 1000.0f;
    m_sampler = SDL_CreateGPUSampler(device, &sampler);
    m_white = create_solid_texture(device, 255, 255, 255, 255);

    SDL_GPUSamplerCreateInfo profile = {};
    profile.min_filter = SDL_GPU_FILTER_LINEAR;
    profile.mag_filter = SDL_GPU_FILTER_LINEAR;
    profile.mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_LINEAR;
    profile.address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    profile.address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    profile.address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
    profile.max_lod = 1000.0f;
    m_profile_sampler = SDL_CreateGPUSampler(device, &profile);
    m_no_rings = create_profile_texture(device, {glm::vec2(0.0f)});

    return m_sphere.vertices && m_sphere.indices && m_sampler && m_white && m_profile_sampler && m_no_rings;
}

void BodyPass::shutdown()
{
    if (!m_device) {
        return;
    }
    clear_meshes();
    SDL_ReleaseGPUBuffer(m_device, m_sphere.vertices);
    SDL_ReleaseGPUBuffer(m_device, m_sphere.indices);
    m_sphere = {};
    SDL_ReleaseGPUSampler(m_device, m_sampler);
    SDL_ReleaseGPUSampler(m_device, m_profile_sampler);
    SDL_ReleaseGPUTexture(m_device, m_white);
    SDL_ReleaseGPUTexture(m_device, m_no_rings);
    m_sampler = nullptr;
    m_profile_sampler = nullptr;
    m_white = nullptr;
    m_no_rings = nullptr;
    SDL_ReleaseGPUGraphicsPipeline(m_device, m_pipeline);
    m_pipeline = nullptr;
    m_device = nullptr;
}

void BodyPass::clear_meshes()
{
    for (Mesh& mesh : m_meshes) {
        SDL_ReleaseGPUBuffer(m_device, mesh.vertices);
        SDL_ReleaseGPUBuffer(m_device, mesh.indices);
    }
    m_meshes.clear();
}

int BodyPass::add_mesh(std::span<const glm::vec3> positions, std::span<const glm::vec3> normals,
                       std::span<const float> albedo, std::span<const float> map_u,
                       std::span<const uint32_t> indices)
{
    if (positions.empty() || normals.size() != positions.size() || albedo.size() != positions.size() ||
        (!map_u.empty() && map_u.size() != positions.size()) || indices.empty()) {
        return -1;
    }
    // The shader takes the map latitude from the vertex direction (uv.y is unused).
    std::vector<BodyVertex> vertices(positions.size());
    for (size_t k = 0; k < positions.size(); ++k) {
        vertices[k] = {positions[k], normals[k], glm::vec2(map_u.empty() ? 0.0f : map_u[k], 0.0f), albedo[k]};
    }
    Mesh mesh;
    mesh.index_count = static_cast<uint32_t>(indices.size());
    mesh.vertices = create_static_buffer(m_device, SDL_GPU_BUFFERUSAGE_VERTEX, vertices.data(),
                                         static_cast<uint32_t>(vertices.size() * sizeof(BodyVertex)));
    mesh.indices = create_static_buffer(m_device, SDL_GPU_BUFFERUSAGE_INDEX, indices.data(),
                                        static_cast<uint32_t>(indices.size() * sizeof(uint32_t)));
    if (!mesh.vertices || !mesh.indices) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Body mesh creation failed: %s", SDL_GetError());
        SDL_ReleaseGPUBuffer(m_device, mesh.vertices);
        SDL_ReleaseGPUBuffer(m_device, mesh.indices);
        return -1;
    }
    m_meshes.push_back(mesh);
    return static_cast<int>(m_meshes.size()) - 1;
}

void BodyPass::draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, const SunLight& sun,
                    std::span<const BodyDrawItem> items) const
{
    if (items.empty()) {
        return;
    }

    SDL_BindGPUGraphicsPipeline(pass, m_pipeline);
    const Mesh* bound = nullptr;

    for (const BodyDrawItem& item : items) {
        const bool own_mesh = item.mesh >= 0 && item.mesh < static_cast<int>(m_meshes.size());
        const Mesh* mesh = own_mesh ? &m_meshes[static_cast<size_t>(item.mesh)] : &m_sphere;
        if (mesh != bound) {
            SDL_GPUBufferBinding vb = {.buffer = mesh->vertices, .offset = 0};
            SDL_BindGPUVertexBuffers(pass, 0, &vb, 1);
            SDL_GPUBufferBinding ib = {.buffer = mesh->indices, .offset = 0};
            SDL_BindGPUIndexBuffer(pass, &ib, SDL_GPU_INDEXELEMENTSIZE_32BIT);
            bound = mesh;
        }

        VertexUniforms vu;
        vu.model = item.model;
        vu.view_proj = view.view_proj;
        vu.rotation = item.rotation;
        vu.inv_scale = glm::vec4(item.inv_scale, -item.texture_left_lon_deg / 360.0f);
        SDL_PushGPUVertexUniformData(cmd, 0, &vu, sizeof(vu));

        FragmentUniforms fu = {};
        fu.sun = glm::vec4(item.sun_direction, item.sun_angular_radius);
        fu.light2 = glm::vec4(item.light2_direction, item.light2_angular_radius);
        fu.light2_color = glm::vec4(item.light2_color, 0.0f);
        fu.color = glm::vec4(item.color, static_cast<float>(item.style));
        fu.params = glm::vec4(static_cast<float>(item.occluder_count), sun.ambient, item.texture ? 1.0f : 0.0f,
                              item.flip_u ? 1.0f : 0.0f);
        for (int k = 0; k < item.occluder_count && k < kMaxOccluders; ++k) {
            fu.occluders[k] = item.occluders[k];
        }
        fu.ring_center = glm::vec4(item.ring_center, item.ring_profile ? 1.0f : 0.0f);
        fu.ring_normal = glm::vec4(item.ring_normal, 0.0f);
        fu.ring_radii = glm::vec4(item.ring_inner_km, item.ring_outer_km, item.ring_samples, 0.0f);
        if (item.atmosphere) {
            const AtmosphereOptics& a = *item.atmosphere;
            fu.atmo_rayleigh = glm::vec4(a.rayleigh_depth, a.rayleigh_scale_height);
            fu.atmo_haze = glm::vec4(a.haze_attenuation, a.haze_scale_height);
        }
        SDL_PushGPUFragmentUniformData(cmd, 0, &fu, sizeof(fu));

        SDL_GPUTextureSamplerBinding tex[2] = {
            {.texture = item.texture ? item.texture : m_white, .sampler = m_sampler},
            {.texture = item.ring_profile ? item.ring_profile : m_no_rings, .sampler = m_profile_sampler},
        };
        SDL_BindGPUFragmentSamplers(pass, 0, tex, 2);

        SDL_DrawGPUIndexedPrimitives(pass, mesh->index_count, 1, 0, 0, 0);
    }
}

} // namespace astraxis
