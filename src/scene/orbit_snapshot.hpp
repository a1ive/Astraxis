#pragma once

#include "ephem/kepler.hpp"

#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <string>
#include <vector>

namespace astraxis {

// An orbit as it was (or will be) at one epoch, drawn whole while an event is
// active ([events.orbit_snapshots]): several of them show an orbit turning over
// millennia (Mercury's perihelion). The elements are osculating ones referred to
// the mean ecliptic and equinox of J2000, as JPL Horizons gives them.
struct OrbitSnapshot {
    double t_tdb = 0.0;
    KeplerElements elements; // radians; the mean anomaly is not used
    std::string label;       // next to the periapsis; may be empty
    // Also the orbit as it would be without general relativity: the periapsis turned
    // back by the relativistic advance since OrbitSnapshots::newtonian_from_tdb.
    bool newtonian = false;
    std::string newtonian_label;
};

struct OrbitSnapshots {
    int body = -1;   // whose orbits (in its orbit color); -1: none
    int center = -1; // the body's parent, which they go around
    double newtonian_from_tdb = 0.0;
    std::vector<OrbitSnapshot> orbits; // oldest first

    bool enabled() const { return body >= 0 && !orbits.empty(); }
};

// From the mean ecliptic and equinox of J2000 to ICRF, with the IAU 1976
// obliquity of 84381.448" that JPL Horizons uses.
glm::dmat3 ecliptic_j2000_to_icrf();

// The periapsis advance per second that general relativity adds to an orbit
// around a mass `gm` (one body, the Schwarzschild term; rad/s): 6 pi GM / (c^2 a
// (1 - e^2)) per orbit, or 3 GM^1.5 / (c^2 a^2.5 (1 - e^2)).
double relativistic_apsidal_rate(double gm_km3_s2, double a_km, double e);

// The snapshot's elements without the relativistic advance since t_from_tdb (in the
// orbit's own plane; the planets' share of the turning stays).
KeplerElements newtonian_elements(const OrbitSnapshot& snapshot, double gm_km3_s2, double t_from_tdb);

// A point of the orbit at true anomaly `nu` (0: the periapsis), in ICRF axes
// relative to the center (km).
glm::dvec3 orbit_point_icrf(const KeplerElements& ecliptic_elements, double nu);

} // namespace astraxis
