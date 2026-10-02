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
// File format (little endian):
//   char[8] "AXEPH1\0\0"; int32 target NAIF id; int32 center NAIF id;
//   uint32 knot count; uint32 reserved;
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
    void set_knots(std::vector<Knot> knots) { m_knots = std::move(knots); }
    std::vector<Knot>& mutable_knots() { return m_knots; }

    int target_id() const { return m_target_id; }
    int center_id() const { return m_center_id; }
    double start() const { return m_knots.front().t; }
    double end() const { return m_knots.back().t; }
    bool covers(double t) const { return !m_knots.empty() && t >= start() && t <= end(); }

    // Interpolated state; t is clamped to [start, end].
    State eval(double t) const;
    const std::vector<Knot>& knots() const { return m_knots; }

    // History sampling times over [t0, t1] (see MotionSource::history_times):
    // the table's knots, thinned to at most max_points.
    void history_times(double t0, double t1, int max_points, std::vector<double>& out) const;

private:
    int m_target_id = 0;
    int m_center_id = 0;
    std::vector<Knot> m_knots;
};

// What an ephemeris-driven body does outside the table's time span.
enum class Extrapolation {
    None,   // not present (e.g. a spacecraft before launch)
    Linear, // continue along the end velocity (e.g. cruising spacecraft)
};

class EphemerisMotion final : public MotionSource {
public:
    // `fallback` (optional) is used outside the table span and takes precedence
    // over `extrapolation`. `parent_gm` (km^3/s^2, optional) enables the
    // osculating-orbit trail.
    EphemerisMotion(std::shared_ptr<const EphemerisTable> table, std::unique_ptr<MotionSource> fallback,
                    Extrapolation extrapolation, double parent_gm);

    State eval(double t_tdb) const override;
    bool valid_at(double t_tdb) const override;
    bool sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const override;
    void history_times(double t0, double t1, int max_points, std::vector<double>& out) const override;

    const EphemerisTable& table() const { return *m_table; }

private:
    std::shared_ptr<const EphemerisTable> m_table;
    std::unique_ptr<MotionSource> m_fallback;
    Extrapolation m_extrapolation;
    double m_parent_gm;
};

} // namespace astraxis
