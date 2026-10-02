#pragma once

#include "ephem/ephemeris.hpp"
#include "ephem/kerr_geodesic.hpp"
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

} // namespace astraxis
