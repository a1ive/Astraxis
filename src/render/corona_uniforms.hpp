#pragma once

#include "scene/corona.hpp"
#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>
#include <array>

namespace astraxis {

inline constexpr int kCoronaOccluders = 8;

// Layout shared by corona.hlsli. Positions and lengths are in solar radii.
struct CoronaUniforms {
    glm::mat4 ray_to_solar{1.0f}; // inverse projection/view, rotated to solar axes
    glm::vec4 camera{0.0f};      // camera relative to Sun, solar axes
    glm::vec4 params{0.0f};     // gain, activity, outer radius, occulter count
    glm::vec4 color{1.0f};      // RGB: linear white-light tint; W: display gamma
    glm::vec4 screen{0.0f};     // full-resolution width, height, unused
    std::array<glm::vec4, kCoronaOccluders> occluders{}; // relative to camera, radius
    std::array<glm::vec4, kCoronaRadialCount> radial{}; // finite-disk scattering A,B
};

} // namespace astraxis
