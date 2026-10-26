#pragma once

#include "scene/scene.hpp"

#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <span>
#include <vector>

namespace astraxis {

// A star as a disk on the sky: position (any frame, km), radius and effective
// temperature. A star filling its Roche lobe is its surface instead (its
// facets' temperatures; radius_km then bounds it), with its axes in the frame
// of `position` (columns).
struct StarDisk {
    glm::dvec3 position{0.0};
    double radius_km = 1.0;
    double temperature_k = 5772.0;
    const RocheSurface* roche = nullptr;
    glm::dmat3 axes{1.0};
};

// Planck's spectral radiance at `wavelength_nm` up to a constant factor
// (1 / (exp(hc / lambda k T) - 1)).
double planck_shape(double temperature_k, double wavelength_nm);

// Light of the disks seen from far away in direction `toward_observer` (unit):
// the sum over the disks of B(T) times the limb-darkened area left visible by
// the disks in front, in units of R^2 (km^2) times planck_shape. A Roche
// surface is summed over its facets, hidden where they fall behind a disk in
// front; in front of a disk, it hides its outline on the sky (convex).
double disks_flux(std::span<const StarDisk> disks, const glm::dvec3& toward_observer, double wavelength_nm,
                  double limb_darkening);
// The same with nothing hidden; a Roche surface seen from where it is brightest.
double disks_full_flux(std::span<const StarDisk> disks, double wavelength_nm, double limb_darkening);

// The scene's light curve at t_tdb, in magnitudes fainter than
// LightCurve::reference_flux (>= 0; without a reference, than the stars
// unhidden at t_tdb). Uses icrf_state_at, so the scene need not be at that time.
double light_curve_magnitude(const Scene& scene, double t_tdb);

// A star's disk at t_tdb (ICRF position from the scene origin; a pulsating
// star's radius and temperature then; a Roche-lobe star's surface and axes).
StarDisk star_disk_at(const Scene& scene, int star, double t_tdb);

// Sets LightCurve::reference_flux: each star at its brightest (over its
// pulsation table, if it pulsates), nothing hidden. The loader calls it.
void set_light_curve_reference(Scene& scene);

// Sets scene.light_curve.depth_mag to the deepest dip found over `days` from
// t_tdb (at least 0.1 mag), rounded up to a tenth of a magnitude.
void measure_light_curve_depth(Scene& scene, double t_tdb, double days);

// An observed light curve's magnitude at t_tdb, interpolated between the
// points on either side; NaN outside the points or across a gap (gap_s).
double observed_magnitude(const ObservedLightCurve& curve, double t_tdb);

// As observed_magnitude, but joined across gaps and held at the first and
// last points: the magnitude to light the scene by at any time.
double observed_magnitude_held(const ObservedLightCurve& curve, double t_tdb);

// How much brighter (negative) or fainter than at the reference time the
// curve's stars, and what they light, are shown at t_tdb (magnitudes; 0
// without stars), and the same as a factor on their light.
double observed_magnitude_offset(const ObservedLightCurve& curve, double t_tdb);
double observed_light_factor(const ObservedLightCurve& curve, double t_tdb);

// Plot opacity at a simulated-seconds-per-real-second rate (0 when paused).
// Fades as the whole span passes in 4 to 2 real seconds.
double light_curve_visibility(const LightCurve& curve, double rate);

// The light curve on a fixed grid of times around the clock, for the plot:
// count samples spaced span_days / (count - 1), from `before` (fraction of the
// span) ahead of the current time. Samples already computed are kept as the
// clock moves, so a frame costs only the new ones.
class LightCurveSamples {
public:
    static constexpr int kCount = 241;
    static constexpr double kBefore = 0.75;

    void reset() { m_first = 0; m_magnitudes.clear(); }
    void update(const Scene& scene, double t_tdb);

    // Time of sample i (TDB seconds since J2000) and its magnitude.
    double time(int i) const { return static_cast<double>(m_first + i) * m_step_s; }
    const std::vector<float>& magnitudes() const { return m_magnitudes; }
    bool empty() const { return m_magnitudes.empty(); }

private:
    long long m_first = 0; // index of sample 0 on the grid of m_step_s
    double m_step_s = 0.0;
    std::vector<float> m_magnitudes;
    std::vector<float> m_scratch;
};

} // namespace astraxis
