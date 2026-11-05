#pragma once

#include "scene/disk.hpp"

#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

#include <array>

namespace astraxis {

inline constexpr int kDiskOccluders = 8;
inline constexpr int kDiskThermalCount = 16;
inline constexpr int kDiskSources = 4;
// The view the shader takes for a thick disk in scattered light, after DiskView.
inline constexpr int kDiskViewThick = 3;

// Layout shared with shaders/disk.hlsli. Lengths are in au, axes are ICRF, and
// positions are relative to the star unless said otherwise.
struct DiskUniforms {
    glm::mat4 ray_to_icrf{1.0f}; // inverse projection/view, rotated to ICRF axes
    glm::vec4 camera{0.0f};      // w: DiskView or kDiskViewThick
    glm::vec4 params{0.0f};      // radiance at the reference (0: nothing to draw), gamma, reference, components
    glm::vec4 range{0.0f};       // ln inner, ln outer, outer, occluders
    glm::vec4 phase_g{0.0f};     // w: a source's least size (radians)
    glm::vec4 phase_w{0.0f};
    glm::vec4 tint{1.0f};        // millimeter view's false color (linear); w: sources
    std::array<glm::vec4, kDiskMaxComponents> axis_x{}; // w: sigma
    std::array<glm::vec4, kDiskMaxComponents> axis_y{}; // w: second sigma
    std::array<glm::vec4, kDiskMaxComponents> axis_z{}; // w: second Gaussian's share
    std::array<glm::vec4, kDiskMaxComponents> color{};  // w: flare radius
    std::array<glm::vec4, kDiskMaxComponents> temperature{}; // K, radius, index, flare index
    std::array<glm::vec4, kDiskOccluders> occluders{};  // relative to the camera; w: radius
    std::array<glm::vec4, kDiskSources> sources{};      // relative to the camera; w: sigma
    std::array<glm::vec4, kDiskSources> source_light{}; // linear radiance at the center
    // Even in ln T over kDiskThermalLoK .. HiK: false color; w: ln of nebula_infrared_emission.
    std::array<glm::vec4, kDiskThermalCount> thermal{};
    // Each component's column per radius (x, y, z, w), or its vertical optical
    // depth in a thick disk.
    std::array<glm::vec4, kDiskRadialCount> table{};
};

// SDL's Vulkan backend binds 4 KiB of each uniform buffer.
static_assert(sizeof(DiskUniforms) <= 4096, "a uniform buffer is bound 4 KiB at a time");

inline constexpr double kDiskThermalLoK = 30.0;
inline constexpr double kDiskThermalHiK = 1000.0;

} // namespace astraxis
