#pragma once

#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <vector>

namespace astraxis {

struct ShapeModel;

// One patch of an expanding nebula's shell (see Nebula): where it is at the
// reference time relative to the star at the center (ICRF, km), the shell's
// outward normal there (unit, ICRF) and its share of the dust relative to the
// mean patch.
struct NebulaGrain {
    glm::vec3 position{0.0f};
    float weight = 1.0f;
    glm::vec3 normal{0.0f, 0.0f, 1.0f};
};

// A dusty shell thrown off in one eruption and coasting outward ever since
// (eta Carinae's Homunculus): every part moves at a constant velocity from the
// star, so the shape keeps its proportions and only grows, as
// (t - t_ejection) / (t_reference - t_ejection), "homologous expansion". The
// shell is a triangle mesh (a 3D model made from observations) covered with
// patches of dust that scatter the starlight: each a flat 3D Gaussian on the
// shell, splatted like a pinwheel's (shaders/nebula.vert.hlsl).
//
// The light of a patch: the star's light falls off as 1 / d^2 from the center,
// and the column of dust through the shell as 1 / scale^2 while the same mass
// spreads over a growing surface, so a face-on patch's surface brightness is
// ~ 1 / (d0^2 scale^4) for its distance d0 at the reference time; scattering
// follows the Henyey-Greenstein phase function, forward scattering brightest
// (the side of the shell between the star and the camera). The shell is taken
// to be optically thin: no shadows, the far side seen through the near one.
struct Nebula {
    int center = -1; // the body the shell expands from (and is lit from)
    double t_ejection_tdb = 0.0;
    double t_reference_tdb = 0.0; // when the mesh has its size
    // Gaussian thickness of the shell, across it, at the reference time (km).
    double wall_sigma_km = 0.0;
    // Gaussian size of a patch along the shell at the reference time (km),
    // from the grains' spacing.
    double patch_sigma_km = 0.0;
    // Distance of the farthest patch from the center at the reference time (km).
    double extent_km = 0.0;
    double asymmetry_g = 0.0; // Henyey-Greenstein <cos theta>
    glm::vec3 color{1.0f};    // of the scattered light (sRGB, largest component 1)
    // Display: a face-on patch display_reference_km from the star at the
    // reference time, seen at right angles to the starlight (phase function 1),
    // has the renderer's reference brightness; the range of surface
    // brightness (distance, time and phase function) is raised to display_gamma.
    double display_reference_km = 0.0;
    double display_gamma = 0.5;
    std::vector<NebulaGrain> grains;

    bool enabled() const { return center >= 0 && !grains.empty(); }
};

// The rotation from a shell model's own axes to ICRF, for a nebula seen at sky
// position (ra_deg, dec_deg): the model's +z axis points toward position angle
// axis_pa_deg (from north through east), tilted axis_tilt_deg from the line of
// sight toward the observer (0: straight at the observer, 90: in the plane of
// the sky); the model's direction at azimuth roll_deg around +z (from +x
// toward +y) faces the observer's side of that axis.
glm::dmat3 nebula_model_to_icrf(double ra_deg, double dec_deg, double axis_pa_deg, double axis_tilt_deg,
                                double roll_deg);

// Covers the shell (a mesh in model units, centered on `model_center`) with
// `count` patches: spread over the triangles in proportion to their area and
// placed by a low-discrepancy sequence within each, with the mesh's smoothed
// normals. Positions are scaled by km_per_unit and turned by model_to_icrf.
// Sets nebula.grains and nebula.patch_sigma_km.
void build_nebula_grains(const ShapeModel& shape, const glm::dvec3& model_center, double km_per_unit,
                         const glm::dmat3& model_to_icrf, int count, Nebula& nebula);

// The size at t_tdb relative to the reference time; 0 before the ejection.
double nebula_scale(const Nebula& nebula, double t_tdb);

// Distance of the farthest patch from the center at t_tdb (km).
double nebula_extent_km(const Nebula& nebula, double t_tdb);

// Henyey-Greenstein phase function relative to isotropic scattering, for the
// cosine of the scattering angle (1: light going straight on).
double henyey_greenstein(double g, double cos_theta);

// The exposure the eye would adapt to at t_tdb (a factor on the shell's
// brightness, 0..1): the young, small shell is brighter by scale^(-4 gamma)
// than at the reference time, and is shown dimmed to at most
// kNebulaMaxBrightening times as bright, the older shell as it is.
double nebula_exposure(const Nebula& nebula, double t_tdb);

} // namespace astraxis
