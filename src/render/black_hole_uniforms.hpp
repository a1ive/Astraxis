#pragma once

#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

namespace astraxis {

// Per-frame inputs of shaders/black_hole.hlsli (see there for the meaning),
// shared by the trace prepass and the composite.
// Everything is in the hole's frame (spin and disk normal along +z), in units of M.
// The camera axes are those of its static observer frame (ObserverFrame), as
// covariant 4-vectors (w = t) and as coordinate directions.
struct BlackHoleUniforms {
    glm::vec4 camera{0.0f};      // xyz = camera position, w = spin a
    glm::vec4 right{0.0f};       // right axis * tan(fov/2) * aspect
    glm::vec4 up{0.0f};          // up axis * tan(fov/2)
    glm::vec4 forward{0.0f};     // forward axis
    glm::vec4 observer{0.0f};    // 4-velocity
    glm::vec4 right_dir{0.0f};   // xyz = direction * tan(fov/2) * aspect, w = tan(fov/2) * aspect
    glm::vec4 up_dir{0.0f};      // xyz = direction * tan(fov/2), w = tan(fov/2)
    glm::vec4 forward_dir{0.0f}; // xyz = direction, w = u^t
    glm::mat4 to_sky{1.0f};      // hole frame -> ICRF
    glm::vec4 radii{0.0f};       // r+, disk inner, disk outer, integration radius
    glm::vec4 params{0.0f};      // b_max, depth, disk brightness, disk time
    glm::vec4 disk{0.0f};        // peak temperature, 1 / profile max, max steps, edge fade start
    glm::vec4 targets{0.0f};     // prepass size, scene size (pixels); set by BlackHolePass::trace
};

} // namespace astraxis
