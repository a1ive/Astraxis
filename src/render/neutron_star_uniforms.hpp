#pragma once

#include <glm/mat4x4.hpp>
#include <glm/vec4.hpp>

#include <array>

namespace astraxis {

struct NeutronStarUniforms {
    static constexpr int kSamples = 129;
    glm::vec4 right;
    glm::vec4 up;
    glm::vec4 forward;
    glm::mat4 to_sky;
    glm::vec4 observer;   // xyz = outward radial direction, w = distance / radius
    glm::vec4 geometry;   // limb angle, near plane / radius, sin(limb), sky solid angle at view centre
    glm::vec4 surface;    // xyz = linear HDR color, w = blend to the relativistic view
    glm::vec4 hot_spot;   // xyz = magnetic pole, w = strength of the two antipodal spots
    glm::vec4 stars;      // sky grid size, star brightness, Milky Way brightness, sky visibility
    glm::vec4 sky_params; // max star dot radius, magnification cap, angular pixel size
    // x: surface angle vs emission cosine; y: sky deflection vs sqrt of
    // log(angle / limb) / log(pi / limb) (denser sampling near the surface).
    std::array<glm::vec4, kSamples> bending;
};

} // namespace astraxis
