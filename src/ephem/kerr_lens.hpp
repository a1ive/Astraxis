#pragma once

#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

namespace astraxis {

// Images of a point (a body's center) lensed by a Kerr black hole, as seen by a
// static camera: the directions in which light rays from the camera reach it.
// Units G = c = M = 1, Cartesian Kerr-Schild coordinates with the spin along +z
// (as kerr_null.hpp); the rays are traced backward from the camera with the same
// equations and step as the shader, but integrated all the way (no weak-field
// legs), and refined near the target.
//
// A ray is labelled by its coordinate direction at the camera d (the spatial
// part of its contravariant momentum): for any static observer frame there, the
// local direction is n = normalize(inverse(right_dir, up_dir, forward_dir) d), so
// d does not depend on how the camera is turned. The miss of a ray is the
// vector from the target to the ray's closest point; an image is a direction
// whose miss is zero, found by Gauss-Newton on the two angles of d.
struct LensImage {
    glm::dvec3 direction{0.0}; // unit coordinate direction at the camera
    // d miss / d direction as a linear map of tangent perturbations of `direction`
    // (M per radian): kept to start the next solve (warm start).
    glm::dmat3 jacobian{0.0};
    int jacobian_age = 0; // solves since it was last measured by finite differences
    // The image of a small sphere of radius rho (M) is an ellipse with semi-axes
    // rho * angle_per_m[k] (radians of coordinate direction) along axis[k]
    // (unit tangents of `direction`).
    glm::dvec3 axis[2]{glm::dvec3(0.0), glm::dvec3(0.0)};
    double angle_per_m[2] = {0.0, 0.0};
    // Image area over that of the same sphere seen at its coordinate distance
    // without the hole; negative for a mirrored image.
    double magnification = 1.0;
    glm::dvec3 emission{0.0}; // unit coordinate direction of the backward ray at the target
    bool disk_first = false;  // the ray crosses the disk before reaching the target
    double miss = 0.0;        // remaining miss (M)
};

struct LensSettings {
    double spin = 0.0;
    double disk_r_in = 0.0; // disk annulus that can hide a target (r_out = 0: no disk)
    double disk_r_out = 0.0;
    double tolerance_rad = 1e-6;  // accepted angular error of an image
    int jacobian_max_age = 8;     // warm solves before the Jacobian is measured again
    // The secondary image is not looked for while the primary's magnification exceeds 1 by
    // less than this: far from the hole |mu-| = mu+ - 1 (point lens), so it is that faint.
    double min_secondary_magnification = 0.0;
    // Search for a secondary image without a guess (a scan of rays around the hole); off
    // to only refine a known one.
    bool search_secondary = true;
};

// The primary image (out[0], on the target's side of the hole) and the secondary
// (out[1], around the other side), each when found. `guesses` (may be null) are
// the previous solution, out[0..1] of an earlier call with `found`: images are
// refined from them when given, else searched for. Returns which were found
// (bit 0: primary, bit 1: secondary). Higher-order images (rays winding more
// than once around the hole) are not looked for: they are exponentially fainter.
int find_lens_images(const glm::dvec3& camera, const glm::dvec3& target, const LensSettings& settings,
                     const LensImage* guesses, int guesses_found, LensImage out[2]);

} // namespace astraxis
