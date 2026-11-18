#pragma once

#include "ephem/kerr_null.hpp"

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

namespace astraxis {

// Hot, geometrically thick, radiatively inefficient accretion flow (RIAF) around
// a Kerr black hole, the semi-analytic model of Broderick et al. 2009 (ApJ 697,
// 45, sec. II.1, eqs. 1, 2, 4), seen in thermal synchrotron light. This is the
// CPU reference for the hot flow in shaders/black_hole.hlsli (same formulas).
//
// Positions are Cartesian Kerr-Schild points in units of M (spin along +z);
// r is the Boyer-Lindquist radius, z the height along the spin axis and
// s = sqrt(r^2 - z^2):
//   n_e = n_e0 (r / 2)^-1.1 exp(-z^2 / 2 s^2),   T_e = T_e0 (r / 2)^-0.84,
//   B^2 / 8 pi = n_e m_p c^2 r_S / (12 r) / beta  (toroidal field).
struct RiafParams {
    double n_e0_cm3 = 0.0; // thermal electron density at r = 2 M
    double t_e0_k = 0.0;   // electron temperature at r = 2 M
    double beta = 10.0;    // gas / magnetic pressure
    double outer_m = 50.0; // no emission beyond this radius (M)
};

// The renderer integrates rays only within 0.8 of its 150 M integration radius
// (kRayIntegrationRadius in scene_renderer.cpp): gas farther out is never drawn.
inline constexpr double kHotFlowMaxOuterM = 120.0;

struct RiafPlasma {
    double n_e_cm3 = 0.0;
    double t_e_k = 0.0;
    double b_gauss = 0.0;
};

// All zero beyond outer_m and inside the horizon.
RiafPlasma riaf_plasma(const glm::dvec3& x, double a, const RiafParams& p);

// Modified Bessel function K_2 (x > 0), from the polynomial approximations of
// Abramowitz & Stegun 9.8.1-9.8.8 (as the shader evaluates it).
double bessel_k2(double x);

// Thermal synchrotron emissivity j_nu (erg s^-1 cm^-3 Hz^-1 sr^-1) of electrons
// at temperature t_e_k in a field b_gauss, at pitch angle sin_pitch: the fit of
// Leung, Gammie & Noble 2011 (ApJ 737, 21, eq. 72).
double thermal_synchrotron_emissivity(double nu_hz, double n_e_cm3, double t_e_k, double b_gauss, double sin_pitch);

// Isotropic emission (as Broderick et al. 2009 assume): the angle-averaged sin(pitch) = pi / 4.
inline constexpr double kRiafSinPitch = 0.78539816339744831;

// How much of a hot flow to draw when 10 M spans `ring_pixels` pixels on screen:
// each pixel point-samples one ray, so a ring of a few pixels would flicker as
// the view moves. 0 below 1.5 pixels, 1 above 4.
double riaf_far_view_fade(double ring_pixels);

// Contravariant 4-velocity u^mu of the gas (w = u^t), Kerr-Schild coordinates.
// Outside the ISCO the gas turns at the Keplerian rate of its radius,
// Omega = 1 / (r^3/2 + a), on every shell. Inside it the gas falls with the
// energy -u_t and angular momentum u_phi it had at the ISCO at the same polar
// angle (ballistic plunge at fixed theta), so the flow is continuous there.
glm::dvec4 riaf_four_velocity(const glm::dvec3& x, double a);

// A flare hot spot (GRAVITY 2018, 2023; Wielgus et al. 2022): a Gaussian blob of
// hotter gas on a Keplerian circle in the hole's equatorial plane, moving with the
// flow. It adds optical depth tau0 exp(-d^2 / 2 sigma^2) / (sigma sqrt(2 pi)) per
// unit fluid-frame length (tau0 along a line through its centre) with source
// function T_s, scaled by `strength` (its brightness envelope and fades).
struct HotSpot {
    double radius_m = 9.0;  // orbit radius (Boyer-Lindquist, M)
    double sigma_m = 2.548; // Gaussian sigma (M): 6 r_g across at half maximum
    double phase_rad = 0.0; // azimuth (hole frame, from +x) at "now" at the hole
    double tau0 = 0.1;
    double t_s_k = 0.0;     // source brightness temperature (K)
    double strength = 0.0;  // 0..1 (0: no spot)
};

// Orbital period 2 pi (R^3/2 + a) M / c (distant observer's time).
double hot_spot_period_s(double a, double radius_m, double mass_km);
// Kerr-Schild Cartesian centre at coordinate time t_m (M) after "now":
// x + i y = (R + i a) e^(i phi), z = 0, phi = phase + Omega t, Omega = 1 / (R^3/2 + a).
glm::dvec3 hot_spot_center(const HotSpot& s, double a, double t_m);

struct FlowRay {
    double tb_k = 0.0; // observed brightness temperature (Rayleigh-Jeans), K
    double tau = 0.0;  // total optical depth along the ray
    RayOutcome outcome = RayOutcome::MaxSteps;
    // Diagnostics for tests: the largest source function g T_e met, and g T_e
    // where the optical depth from the camera reaches 1 (0: never).
    double max_source_k = 0.0;
    double tau1_source_k = 0.0;
};

// Radiative transfer along a ray traced backward from a camera at rest (photon
// with observed energy 1, e.g. frame_photon()), as the shader does it: at the end
// of each step (ray_step), the gas emits g T_e (1 - e^-dtau) seen through the
// optical depth in front of it, with g = nu_obs / nu_emit = 1 / (p.u) and
// dtau = alpha_nu(nu_emit) (p.u) h M, alpha_nu = j_nu / B_nu(T_e) (Kirchhoff).
// With a hot spot the ray also accumulates Kerr-Schild coordinate time (trapezoid
// over each step) and sees the spot where it was when the light left it: at
// t = dt_ray + |camera|, "now" at the hole being when light from its centre reaches
// the camera along a straight line.
FlowRay trace_flow_ray(Photon ph, double a, double mass_cm, double nu_obs_hz, const RiafParams& p,
                       double escape_radius, int max_steps = 4000, const HotSpot* spot = nullptr);

} // namespace astraxis
