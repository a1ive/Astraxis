#pragma once

#include "scene/scene.hpp"

namespace astraxis {

// Unit spin axis in ICRF at t_tdb. Geodetic precession follows WT02 Eq. 2
// (Weisberg & Taylor 2002, astro-ph/0205280); otherwise the orbit normal.
// Requires an enabled pulsar with its orbit reference resolved by the loader.
glm::dvec3 pulsar_spin_axis(const Scene& scene, const Body& body, double t_tdb);

// The corotating frame in the current display frame (columns): z along the spin, x the
// spin-perpendicular part of the magnetic axis, which is cos(alpha) z + sin(alpha) x.
// Precession uses simulated time; the turn about the spin uses the slowed display clock.
bool pulsar_corotating_axes(const Scene& scene, const Body& body, double t_tdb, double real_time, glm::dmat3& out);

// Beam directions in the current display frame (the magnetic axis and its opposite).
bool pulsar_beams(const Scene& scene, const Body& body, double t_tdb, double real_time, glm::dvec3 out[2]);

// c / Omega for the true spin period (km).
double light_cylinder_km(double spin_period_s);

// The field of a centered dipole, corotating, inside the light cylinder (Goldreich & Julian
// 1969): lines r = L sin^2(theta) about the magnetic axis, which is inclined by
// `inclination_deg` to the spin (z), in the corotating frame of pulsar_corotating_axes, km.
// In each of `azimuths` magnetic meridian planes, a line is closed when it stays within the
// light cylinder (distance from the spin axis below `light_cylinder`); the last closed one
// just touches it. Closed lines at `closed_fractions` of that L; open lines leave the polar
// caps at `open_factors` times it (and along the magnetic axis), from both poles, to the
// cylinder. w is the fade: 0, rising to 1 as an open line nears the cylinder.
using PulsarFieldLine = Body::Pulsar::FieldLine;
std::vector<PulsarFieldLine> pulsar_field_lines(double inclination_deg, double radius_km, double light_cylinder,
                                                int azimuths = 8);
// The L of the last closed line in the magnetic meridian plane at magnetic azimuth phi
// (0: the plane through the spin axis, on the side the magnetic axis leans to).
double pulsar_last_closed_l(double inclination_deg, double phi, double light_cylinder);

// Illustrative flash from one pole; both arguments are unit display directions.
double pulsar_beam_flash(const Body::Pulsar& pulsar, const glm::dvec3& axis, const glm::dvec3& to_camera);

} // namespace astraxis
