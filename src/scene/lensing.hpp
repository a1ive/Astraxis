#pragma once

#include "ephem/kerr_lens.hpp"

#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <vector>

namespace astraxis {

class Scene;
struct Body;
struct CameraView;

// Axes of a black hole's own frame (spin and disk normal along +z) in the
// display frame: the frame the ray tracer works in, scaled by M.
glm::dmat3 black_hole_axes(const Scene& scene, const Body& hole);

// One image of a body bent around a black hole, as seen by the camera: where
// and how to draw the body so that it shows there (SceneRenderer), and where its
// label goes (LabelLayout). Only bodies are lensed, not orbit lines.
struct LensedView {
    // Where the body's center is drawn, camera-relative (km, display frame): along
    // the image's direction, at the body's distance, but pulled in front of the
    // hole's depth when the light reaches the body before the disk or the horizon
    // (the ray traced composite then leaves it showing), and pushed behind it when
    // the light crosses the disk on the way (the disk covers it).
    glm::dvec3 offset{0.0};
    // Applied to the body about its center, in this order (display frame): `turn`
    // brings the side the light leaves from to face the camera (the far side, in a
    // secondary image), with its light; `stretch` gives the image's size and shape
    // (an ellipse: arcs near the shadow), identity for an unlensed body.
    glm::dmat3 turn{1.0};
    glm::dmat3 stretch{1.0};
    double angular_radius_per_km = 0.0; // the image's mean angular radius per km of body radius (radians)
    double magnification = 1.0;         // area over that without the hole; negative: mirrored
    bool behind_disk = false;
    LensImage solution; // the solver's state, refined on the next frame
};

struct BodyLens {
    bool lensed = false; // false: the body is drawn where it is
    int found = 0;       // bit 0: the primary image, bit 1: the secondary (find_lens_images)
    int search_wait = 0; // frames until a secondary image not found is searched for again
    LensedView images[2];
};

// For each body orbiting a black hole (through any parent), its primary and
// secondary images from the camera; refines those `lens` already holds (the
// last frame's for this output). A body the camera is right next to, or for which
// no image is found, is not lensed. Drawn bodies are sized as the solution
// gives, so this also applies to bodies too small to see (their labels).
void lens_bodies(const Scene& scene, const CameraView& view, std::vector<BodyLens>& lens);

} // namespace astraxis
