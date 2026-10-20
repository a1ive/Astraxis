#pragma once

#include "ephem/ephemeris.hpp"
#include "ephem/kerr_geodesic.hpp"
#include "ephem/kerr_null.hpp"
#include "ephem/motion.hpp"
#include "ephem/visual_orbit.hpp"

#include <glm/mat3x3.hpp>

namespace astraxis {

// A star orbiting a black hole along a Kerr geodesic.
//
// Initial conditions: the Newtonian state of the osculating `elements` at
// `elements_epoch_tdb` (e.g. the GRAVITY S2 elements are osculating at the 2010
// apocentre). The orbit plane is the elements' plane, and the black hole's
// spin axis is taken along the orbit normal (equatorial geodesic) - exact for
// a non-spinning hole, an assumption otherwise.
//
// Trajectories are integrated ahead of time into a knot table (Hermite
// interpolation, as for baked ephemerides), either once over a fixed span or
// as a rolling window that follows the clock (for hypothetical objects with
// short periods, restarting from the elements after large time jumps).
struct KerrOrbitSetup {
    double gm_km3_s2 = 0.0;
    double spin = 0.0; // a / M, in [0, 1)
    VisualOrbit elements;
    double elements_epoch_tdb = 0.0;
    double step_eta = 0.003; // step = eta * sqrt(r^3 / M) (geometric units)
};

class KerrOrbitMotion final : public MotionSource {
public:
    // Fixed span [t_begin, t_end].
    KerrOrbitMotion(const KerrOrbitSetup& setup, double t_begin, double t_end);
    // Rolling window keeping about `window_days` around the requested times.
    static std::unique_ptr<KerrOrbitMotion> rolling(const KerrOrbitSetup& setup, double window_days);

    State eval(double t_tdb) const override;
    bool valid_at(double t_tdb) const override;
    void history_times(double t0, double t1, int max_points, std::vector<double>& out) const override;

    const KerrEquatorial& metric() const { return m_kerr; }
    const EphemerisTable& table() const { return m_table; }
    // Orbit plane in ICRF (columns: P toward the initial periapsis, Q, normal).
    const glm::dmat3& plane() const { return m_plane; }

private:
    struct PlaneState {
        double t = 0.0; // TDB seconds
        glm::dvec2 x{0.0};
        glm::dvec2 u{0.0};
    };

    explicit KerrOrbitMotion(const KerrOrbitSetup& setup);
    PlaneState initial_state(double epoch_tdb) const;
    // Integrates from `from` to time t_end, appending (forward) or prepending
    // (backward) knots; returns the final state.
    PlaneState integrate(PlaneState from, double t_end, std::vector<EphemerisTable::Knot>& out) const;
    EphemerisTable::Knot knot(const PlaneState& s) const;
    void ensure_covered(double t_tdb) const;
    void restart(double t_tdb) const;

    KerrOrbitSetup m_setup;
    KerrEquatorial m_kerr;
    glm::dmat3 m_plane{1.0};
    double m_mass_km = 0.0;

    bool m_rolling = false;
    double m_window_s = 0.0;
    mutable EphemerisTable m_table;
    // Integrator states at the table ends (rolling mode); invalid once trimmed.
    mutable PlaneState m_front;
    mutable PlaneState m_back;
    mutable bool m_front_valid = false;
    mutable bool m_back_valid = false;
};

// A prograde circular orbit in a Kerr hole's equatorial plane, at Boyer-Lindquist
// radius r (units of M; stable from the ISCO outward). Exact, so no integration:
// the coordinate (distant observer's) angular velocity is Omega = 1 / (r^3/2 + a)
// (Bardeen, Press & Teukolsky 1972, eq. 2.16). Placed at the Kerr-Schild
// Cartesian radius sqrt(r^2 + a^2), as the ray tracer draws the disk.
class KerrCircularMotion final : public MotionSource {
public:
    // pole: the spin axis (ICRF); phase_deg: angle from the ascending node of the
    // equator on the ICRF equator (as for the IAU W) at J2000.
    KerrCircularMotion(double gm_km3_s2, double spin, double r_m, const glm::dvec3& pole, double phase_deg);

    State eval(double t_tdb) const override;
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override;

    double period_s() const;

private:
    double m_radius_km = 0.0;
    double m_omega = 0.0; // rad/s, distant observer's time
    double m_phase = 0.0; // rad at J2000
    glm::dvec3 m_e1{1.0, 0.0, 0.0};
    glm::dvec3 m_e2{0.0, 1.0, 0.0};
};

// A body on a Kerr geodesic of any inclination, e.g. a long orbit that dives
// close to the hole out of its equatorial plane. Integrated in Cartesian
// Kerr-Schild coordinates with the ray tracer's Hamiltonian (kerr_null.hpp,
// H = -1/2, RK4 in proper time with steps ~ eta r^3/2) over a fixed span into a
// knot table, like KerrOrbitMotion.
//
// It starts at an apoapsis: at Cartesian distance apoapsis_m (units of M) in
// the direction given by Keplerian angles in the hole's frame (z = spin axis,
// x = ascending node of the hole's equator on the ICRF equator), moving
// perpendicular to the radius in the plane they define. The speed there is
// found by bisection so that the next periapsis is at Boyer-Lindquist radius
// periapsis_m: just above the capture limit, the body "zooms" in and "whirls"
// around the hole before flying out again.
struct KerrApoapsisSetup {
    double gm_km3_s2 = 0.0;
    double spin = 0.0;
    glm::dvec3 pole{0.0, 0.0, 1.0}; // the hole's spin axis (ICRF)
    double apoapsis_m = 0.0;
    double periapsis_m = 0.0;
    double inclination_deg = 0.0;   // to the hole's equator
    double node_deg = 0.0;
    double arg_apoapsis_deg = 0.0;  // from the node, in the orbit plane
    double apoapsis_tdb = 0.0;
    double step_eta = 0.002;
};

class KerrApoapsisMotion final : public MotionSource {
public:
    KerrApoapsisMotion(const KerrApoapsisSetup& setup, double t_begin, double t_end);

    State eval(double t_tdb) const override;
    bool valid_at(double t_tdb) const override;
    void history_times(double t0, double t1, int max_points, std::vector<double>& out) const override;

    // Speed at the apoapsis (units of c) and the particle there.
    double apoapsis_speed() const { return m_speed; }
    const Photon& initial_state() const { return m_start; }
    // Hole frame -> ICRF (columns: x, y, spin axis).
    const glm::dmat3& hole_frame() const { return m_frame; }
    double mass_km() const { return m_mass_km; }

private:
    Photon particle(double speed) const;
    // Lowest Boyer-Lindquist radius of the first inward pass (0: captured).
    double first_periapsis(double speed) const;
    void integrate(Photon s, double t_m, double t_end_m, std::vector<EphemerisTable::Knot>& out) const;

    KerrApoapsisSetup m_setup;
    glm::dmat3 m_frame{1.0};
    glm::dvec3 m_p{1.0, 0.0, 0.0}; // toward the apoapsis (hole frame)
    glm::dvec3 m_q{0.0, 1.0, 0.0}; // direction of motion there
    double m_mass_km = 0.0;
    double m_speed = 0.0;
    Photon m_start;
    EphemerisTable m_table;
};

} // namespace astraxis
