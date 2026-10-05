#pragma once

#include <glm/vec3.hpp>

#include <memory>
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

// Another motion scaled by a constant. Places both bodies of a pair around
// their barycenter from one relative orbit r (secondary relative to primary):
// the primary at -m2 / (m1 + m2) r, the secondary at m1 / (m1 + m2) r.
class ScaledMotion final : public MotionSource {
public:
    ScaledMotion(std::unique_ptr<MotionSource> inner, double scale)
        : m_inner(std::move(inner))
        , m_scale(scale)
    {
    }
    State eval(double t_tdb) const override
    {
        const State s = m_inner->eval(t_tdb);
        return {s.position * m_scale, s.velocity * m_scale};
    }
    bool valid_at(double t_tdb) const override { return m_inner->valid_at(t_tdb); }
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override
    {
        if (!m_inner->sample_orbit(t_tdb, count, out)) {
            return false;
        }
        for (glm::dvec3& p : out) {
            p *= m_scale;
        }
        return true;
    }
    void history_times(double t0, double t1, int max_points, std::vector<double>& out) const override
    {
        m_inner->history_times(t0, t1, max_points, out);
    }
    const MotionSource& inner() const { return *m_inner; }
    double scale() const { return m_scale; }

private:
    std::unique_ptr<MotionSource> m_inner;
    double m_scale;
};

// Another motion offset by a sibling's (both relative to the same parent):
// for data given relative to that sibling, e.g. a spacecraft orbiting Mercury
// that stays a child of the Sun for its cruise. `anchor` is not owned.
class OffsetMotion final : public MotionSource {
public:
    // `anchors`: the anchor's motion, then those of its ancestors up to (not including) the
    // body's parent; their sum is the anchor relative to the parent.
    OffsetMotion(std::unique_ptr<MotionSource> inner, std::vector<const MotionSource*> anchors)
        : m_inner(std::move(inner))
        , m_anchors(std::move(anchors))
    {
    }
    State eval(double t_tdb) const override
    {
        State s = m_inner->eval(t_tdb);
        for (const MotionSource* anchor : m_anchors) {
            const State a = anchor->eval(t_tdb);
            s.position += a.position;
            s.velocity += a.velocity;
        }
        return s;
    }
    bool valid_at(double t_tdb) const override
    {
        if (!m_inner->valid_at(t_tdb)) {
            return false;
        }
        for (const MotionSource* anchor : m_anchors) {
            if (!anchor->valid_at(t_tdb)) {
                return false;
            }
        }
        return true;
    }
    void history_times(double t0, double t1, int max_points, std::vector<double>& out) const override
    {
        m_inner->history_times(t0, t1, max_points, out);
    }
    const MotionSource& inner() const { return *m_inner; }

private:
    std::unique_ptr<MotionSource> m_inner;
    std::vector<const MotionSource*> m_anchors;
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
