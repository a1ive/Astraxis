#pragma once

#include <glm/vec3.hpp>

#include <string>
#include <vector>

namespace astraxis {

// A relativistic jet pair from a black hole, drawn in its visible (optical synchrotron)
// light. This is the CPU reference of shaders/jet.frag.hlsl (same formulas).
//
// Lengths are in units of M (GM/c^2). The jet frame has z along the jets' mean axis
// (the precession axis; +z is the jet approaching the Earth), x toward the side the
// nozzle tilts to at the moment. At height z (zeta = |z|):
//   radius      w^2 = k zeta + (c zeta)^2: parabolic near the hole, conical far out;
//   streamlines through the same fraction psi = s / w of the radius;
//   speed       Gamma = (zeta / z1)^index (at least gamma_beta_min in Gamma beta) up to
//               accel_end, then Gamma beta ~ zeta^decel_index;
//   cross-section: a limb-brightened ridge at psi = 0.82 (its outer half maximum at the
//               radius) and a weaker spine;
//   brightness per unit length L(zeta): below the first knot, L ~ w^-(1 + alpha) (a
//               toroidal field B ~ 1/w, electrons at a fixed fraction of its pressure,
//               j ~ B^(3 + alpha)); along the knots, their measured fluxes;
//   Doppler boost D^(2 + alpha), D = 1 / (Gamma (1 - beta v.n)), F_nu ~ nu^-alpha.
// The nozzle (the jet within the tilt range) leans by the precession half-angle: the
// spine at height z lies at z tan(angle tau(zeta)) along x, tau rising from 0 to 1
// between tilt_m[0] and tilt_m[1] and falling back between tilt_m[2] and tilt_m[3]
// (log-smoothstep), so near the hole and far out the jet follows the mean axis.
struct JetKnot {
    std::string name;
    double start_m = 0.0; // deprojected distance range along the jet (M)
    double end_m = 0.0;
    double flux = 0.0;    // measured flux density (any unit: only ratios matter)
    double line_brightness = 0.0; // intrinsic, per unit length (set by prepare_jet)
};

struct JetPrecession {
    double t0_year = 2000.0;            // reference time of Cui et al. 2023 eq. 1 (Julian years)
    double omega_rad_per_year = 0.0;    // precession rate; the nozzle turns by -omega (t - t0) about the axis
    double half_angle_deg = 0.0;        // half-opening angle of the precession cone
    double tilt_m[4] = {10.0, 100.0, 1.0e4, 1.0e5};
};

struct JetParams {
    double parabola_k_m = 15.4;     // w^2 = k zeta + (cone_slope zeta)^2
    double cone_slope = 0.007;
    double lorentz_unity_m = 600.0; // Gamma = 1 at this height
    double lorentz_index = 0.16;
    double gamma_beta_min = 0.4;
    double accel_end_m = 1.0e6;
    double decel_index = 0.0;
    double spectral_index = 0.9;    // alpha in F_nu ~ nu^-alpha
    double base_m = 5.0;            // emission fades in from base_m to 4 base_m
    double length_m = 2.0e7;        // and out over the last 10%
    double display_gain = 1.0;
    std::vector<JetKnot> knots;     // ordered outward, at most kMaxJetKnots
    bool precesses = false;
    JetPrecession precession;
    // Where the jet is on the sky, for the knots' projected positions and the
    // precession angles (the approaching jet is the one toward the Earth).
    double ra_deg = 0.0;
    double dec_deg = 0.0;
    double distance_pc = 0.0;
};

inline constexpr int kMaxJetKnots = 12;

// Cross-section profile (psi = s / w) and its area integral 2 pi int profile psi dpsi.
double jet_profile(double psi);
double jet_profile_area();

double jet_radius(const JetParams& p, double zeta);
// dw / dzeta.
double jet_radius_slope(const JetParams& p, double zeta);
// Gamma beta of the flow at height zeta.
double jet_gamma_beta(const JetParams& p, double zeta);
// Fraction of the precession half-angle the spine leans by at height zeta.
double jet_tilt_fraction(const JetPrecession& p, double zeta);
// Intrinsic brightness per unit length (relative units).
double jet_line_brightness(const JetParams& p, double zeta);
// Doppler factor of the flow on the spine at height zeta seen at angle theta from it.
double jet_doppler(const JetParams& p, double zeta, double theta);

// Fills in the knots' intrinsic line brightnesses (their measured fluxes per unit
// length with the Doppler boost toward the Earth, at angle theta_earth from the axis,
// taken out).
void prepare_jet(JetParams& p, double theta_earth);

// Observer-frame emissivity (relative units: line brightness per M^2) at x (jet frame,
// M) seen from the unit direction to_camera, with the nozzle tilted by tan_tilt (tan of
// the precession half-angle; 0 for none). min_width_m widens the jet (conserving its
// line brightness) where it would be thinner than a pixel.
double jet_emissivity(const JetParams& p, const glm::dvec3& x, const glm::dvec3& to_camera, double tan_tilt,
                      double min_width_m = 0.0);

// Directions at a place on the sky (ICRF unit vectors): north, east and toward the observer.
struct SkyAxes {
    glm::dvec3 north{0.0};
    glm::dvec3 east{0.0};
    glm::dvec3 toward_observer{0.0};
};
SkyAxes sky_axes(double ra_deg, double dec_deg);

// Position angle (radians east of north) of the projection of v on the sky.
double sky_position_angle(const glm::dvec3& v, const SkyAxes& sky);

// The nozzle's direction (ICRF) at `year` (Julian years) for a cone about the unit
// precession axis `axis` (pointing toward the observer's side): Cui et al. 2023 eqs. 1-2,
//   j = R_z(eta_p) R_y(theta) R_z(-omega (t - t0)) R_y(phi_jet) z,
// in their observer frame (x north, y east, z toward the observer), with eta_p and theta
// the axis' position angle and angle from the line of sight.
glm::dvec3 jet_nozzle(const JetPrecession& p, const glm::dvec3& axis, const SkyAxes& sky, double year);

} // namespace astraxis
