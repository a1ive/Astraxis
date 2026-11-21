#pragma once

#include "ephem/jet.hpp"

#include <glm/vec4.hpp>

namespace astraxis {

// Per-frame inputs of shaders/jet.frag.hlsl, in the jet frame (ephem/jet.hpp: z along
// the mean axis toward the approaching jet, x toward the nozzle's tilt), units of M.
struct JetUniforms {
    glm::vec4 camera{0.0f};   // xyz = camera position, w = shadow radius (rays passing closer end there)
    glm::vec4 right{0.0f};    // xyz = right axis * tan(fov/2) * aspect
    glm::vec4 up{0.0f};       // xyz = up axis * tan(fov/2)
    glm::vec4 forward{0.0f};  // xyz = forward axis, w = smallest drawn radius per unit distance (radians)
    glm::vec4 shape{0.0f};    // parabola k, cone slope^2, base, length
    glm::vec4 speed{0.0f};    // height of Gamma = 1, Lorentz index, smallest Gamma beta, end of acceleration
    glm::vec4 speed2{0.0f};   // deceleration index, Gamma beta at the end of acceleration, spectral index,
                              // profile area
    glm::vec4 tilt{0.0f};     // half-angle (rad), nozzle end, outer ridge dx/dz and dy/dz
    glm::vec4 tilt_log{0.0f}; // ln of the four tilt heights
    glm::vec4 envelope{0.0f}; // s^2 < A z + B z^2 + C bounds the jet: A, B with the tilt, A, B beyond it
    glm::vec4 envelope2{0.0f}; // x = C
    glm::vec4 envelope3{0.0f}; // approaching jet beyond the first ridge node: A, B, C, that height
    glm::vec4 color{0.0f};    // rgb = color x display gain, w = exposure; drawn as rgb sqrt(w I)
    glm::vec4 inner{0.0f};    // x = first knot start, y = jet radius there, z = knot count, w = ridge nodes
    glm::vec4 knots[kMaxJetKnots] = {}; // start, end, receding and approaching diffuse brightness
    glm::vec4 ridge[kMaxJetKnots] = {}; // approaching ridge: height, (x, y) offset
    glm::vec4 feature_bounds{0.0f}; // z min, z max, cylinder radius squared, count
    glm::vec4 feature_centers[kMaxJetFeatures] = {}; // xyz in the jet frame, w = radius
    glm::vec4 feature_light[(kMaxJetFeatures + 3) / 4] = {}; // packed surface brightnesses
};
static_assert(sizeof(JetUniforms) <= 4096, "SDL's Vulkan backend limits a uniform binding to 4 KiB");

} // namespace astraxis
