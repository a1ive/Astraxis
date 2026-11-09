#pragma once

namespace astraxis {

// A static observer outside a spherical, opaque star. Distances (including
// mass GM/c²) share a unit. The stellar surface must lie outside r = 3M.
struct SchwarzschildStarLight {
    double mass;
    double radius;
    double observer;

    double limb_angle() const;
    // Angle swept from the observer's radius vector to the surface hit.
    // emission_cosine is measured in the local frame at the surface.
    double surface_angle(double emission_cosine) const;
    // Extra in-plane turn relative to the local ray direction. Only for rays
    // missing the surface. angle is measured from the inward radial direction.
    double sky_deflection(double angle) const;
};

} // namespace astraxis
