#pragma once

#include "scene/scene.hpp"

namespace astraxis {

// Unit spin axis in ICRF at t_tdb. Geodetic precession follows WT02 Eq. 2
// (Weisberg & Taylor 2002, astro-ph/0205280); otherwise the orbit normal.
// Requires an enabled pulsar with its orbit reference resolved by the loader.
glm::dvec3 pulsar_spin_axis(const Scene& scene, const Body& body, double t_tdb);

// Beam directions in the current display frame. Precession uses simulated
// time; rotation about the spin axis uses the slowed display clock.
bool pulsar_beams(const Scene& scene, const Body& body, double t_tdb, double real_time, glm::dvec3 out[2]);

// Illustrative flash from one pole; both arguments are unit display directions.
double pulsar_beam_flash(const Body::Pulsar& pulsar, const glm::dvec3& axis, const glm::dvec3& to_camera);

} // namespace astraxis
