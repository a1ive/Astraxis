#pragma once

#include <glm/vec3.hpp>

#include <vector>

namespace astraxis {

// Position (km) and velocity (km/s).
struct State {
    glm::dvec3 position{0.0};
    glm::dvec3 velocity{0.0};
};

// Source of a body's motion relative to its parent, in ICRF axes.
// Display frames (e.g. rotating Sun-Earth) are applied by the scene.
class MotionSource {
public:
    virtual ~MotionSource() = default;

    // State at t_tdb (TDB seconds since J2000).
    virtual State eval(double t_tdb) const = 0;

    // False when the body does not exist at t_tdb (e.g. before launch).
    virtual bool valid_at(double /*t_tdb*/) const { return true; }

    // One revolution of the osculating orbit at t_tdb, starting at the body's
    // current position and going backwards in time, `count` points with the
    // last equal to the first. Returns false if the source has no closed orbit.
    virtual bool sample_orbit(double /*t_tdb*/, int /*count*/, std::vector<glm::dvec3>& /*out*/) const
    {
        return false;
    }

    // Times at which to sample a history trail over [t0, t1] (ascending, both
    // ends included, at most `max_points`). The default is uniform; ephemerides
    // return their knots, which are dense exactly where the path bends.
    virtual void history_times(double t0, double t1, int max_points, std::vector<double>& out) const
    {
        out.clear();
        const int n = max_points < 2 ? 2 : max_points;
        for (int i = 0; i < n; ++i) {
            out.push_back(t0 + (t1 - t0) * static_cast<double>(i) / static_cast<double>(n - 1));
        }
    }
};

// A body at a fixed position relative to its parent (e.g. a distant star).
class FixedMotion final : public MotionSource {
public:
    explicit FixedMotion(const glm::dvec3& position)
        : m_position(position)
    {
    }
    State eval(double /*t_tdb*/) const override { return {m_position, glm::dvec3(0.0)}; }

private:
    glm::dvec3 m_position;
};

} // namespace astraxis
