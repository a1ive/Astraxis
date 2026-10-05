#pragma once

#include "ephem/motion.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace astraxis {

// A baked ephemeris (tools/bake/horizons_bake.py): knots of (t, position,
// velocity) in ICRF, relative to a center body, interpolated with cubic
// Hermite splines.
//
// With a reference GM (AXEPH2), the spline interpolates only the deviation
// from the two-body arc that starts at the earlier knot. Orbits around a
// planet (moons, orbiters) then need several times fewer knots.
//
// A table may also have gaps: knots further apart than its maximum spacing
// enclose no data (e.g. a moon kept only around a spacecraft's flybys).
//
// File format (little endian):
//   char[8] "AXEPH1\0\0" or "AXEPH2\0\0"; int32 target NAIF id; int32 center NAIF id;
//   uint32 knot count; uint32 maximum knot spacing in seconds (0: no gaps);
//   AXEPH2 only: float64 reference GM of the center (km^3/s^2)
//   knots: float64 t (TDB s since J2000), x, y, z (km), vx, vy, vz (km/s)
class EphemerisTable {
public:
    struct Knot {
        double t = 0.0;
        glm::dvec3 position{0.0};
        glm::dvec3 velocity{0.0};
    };

    bool load(const std::filesystem::path& path, std::string* error);
    // Tables built in memory (integrated trajectories, tests).
    void set_knots(std::vector<Knot> knots)
    {
        m_knots = std::move(knots);
        m_max_gap = 0.0;
    }
    std::vector<Knot>& mutable_knots() { return m_knots; }

    int target_id() const { return m_target_id; }
    double reference_gm() const { return m_reference_gm; } // 0: plain Hermite
    int center_id() const { return m_center_id; }
    double start() const { return m_knots.front().t; }
    double end() const { return m_knots.back().t; }
    // Whether t has data: inside [start, end] and not in a gap.
    bool covers(double t) const;
    // Start and end of the stretch of data around t (requires covers(t)).
    void segment(double t, double* seg_start, double* seg_end) const;
    double max_gap() const { return m_max_gap; } // 0: no gaps

    // Interpolated state; t is clamped to [start, end].
    State eval(double t) const;
    const std::vector<Knot>& knots() const { return m_knots; }

    // History sampling times over [t0, t1] (see MotionSource::history_times):
    // the table's knots, thinned to at most max_points.
    void history_times(double t0, double t1, int max_points, std::vector<double>& out) const;

private:
    int m_target_id = 0;
    int m_center_id = 0;
    double m_reference_gm = 0.0;
    double m_max_gap = 0.0;
    std::vector<Knot> m_knots;
};

// What an ephemeris-driven body does outside the table's time span.
enum class Extrapolation {
    None,   // not present (e.g. a spacecraft before launch)
    Linear, // continue along the end velocity (e.g. cruising spacecraft)
    Kepler, // two-body orbit from the nearer end of the table (its reference GM), both ways
};

class EphemerisMotion final : public MotionSource {
public:
    // `fallback` (optional) is used outside the table span and takes precedence
    // over `extrapolation` (Kepler extrapolation needs a table with a reference GM); within `blend_s` of either end of the span (or of a
    // stretch between gaps) the position fades between the two, so a less accurate fallback (e.g. mean
    // elements) does not jump. `parent_gm` (km^3/s^2, optional) enables the
    // osculating-orbit trail.
    EphemerisMotion(std::shared_ptr<const EphemerisTable> table, std::unique_ptr<MotionSource> fallback,
                    Extrapolation extrapolation, double parent_gm, double blend_s = 0.0);

    State eval(double t_tdb) const override;
    bool valid_at(double t_tdb) const override;
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override;
    void history_times(double t0, double t1, int max_points, std::vector<double>& out) const override;

    const EphemerisTable& table() const { return *m_table; }
    const MotionSource* fallback() const { return m_fallback.get(); }

private:
    std::shared_ptr<const EphemerisTable> m_table;
    std::unique_ptr<MotionSource> m_fallback;
    Extrapolation m_extrapolation;
    double m_parent_gm;
    double m_blend_s;
};

} // namespace astraxis
