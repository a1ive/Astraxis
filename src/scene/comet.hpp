#pragma once

#include "ephem/motion.hpp"

#include <glm/vec3.hpp>

#include <cstdint>
#include <deque>
#include <limits>
#include <vector>

namespace astraxis {

class Scene;

// Haser scale lengths of C2 at 1 au, parent and daughter (km), both growing as
// r_h^2: A'Hearn et al. 1995 (Icarus 118, 223), as tabulated in Gilbert et al.
// 2010 (arXiv:0910.0416, Table 3), with an outflow speed of 1 km/s.
inline constexpr double kC2ParentScaleKm = 2.2e4;
inline constexpr double kC2DaughterScaleKm = 6.6e4;

// Total visual magnitude of a comet seen from 1 au, at `r_au` from the sun:
// m = M1 + 5 log10(delta) + K1 log10(r) with delta = 1 (JPL SBDB's comet
// magnitude parameters).
double comet_total_magnitude(double m1, double k1, double r_au);

// Direction of an ion (type I) tail: along the solar wind as the comet sees it,
// the radial wind at `solar_wind_km_s` minus the comet's own heliocentric
// velocity (the aberration of the tail). Heliocentric position and velocity in
// any one set of axes; the result is a unit vector in the same axes.
glm::dvec3 ion_tail_direction(const glm::dvec3& helio_position, const glm::dvec3& helio_velocity,
                              double solar_wind_km_s);

// A dust grain in the Finson-Probstein model (Finson & Probstein 1968, ApJ 154,
// 327): it leaves the nucleus with the comet's heliocentric state plus an
// ejection velocity and then moves on a two-body orbit around the sun, whose
// attraction radiation pressure reduces by the factor 1 - beta (beta = radiation
// pressure / gravity, larger for smaller grains: Burns, Lamy & Soter 1979,
// Icarus 40, 1). Heliocentric state `dt` seconds after the ejection;
// `universal_anomaly` as for propagate_kepler.
State dust_grain_state(const State& comet_at_ejection, const glm::dvec3& ejection_velocity, double beta,
                       double gm_sun, double dt, double* universal_anomaly = nullptr);

// One grain of a dust tail, as drawn.
struct DustGrain {
    glm::dvec3 position{0.0}; // ICRF, relative to the scene origin
    double weight = 0.0;      // share of the tail's light (see DustTail::update)
    double age_s = 0.0;
    double ejection_km_s = 0.0; // its ejection speed
    double random = 0.0;        // in [0, 1), fixed for the grain (e.g. to thin them out)
};

// The dust tail of one comet as grains emitted at fixed times (slots of
// span / kSlots, each grain at its own time within its slot) over the last
// `dust_tail_days`, for a spread of beta and ejection directions. A grain's
// position at any time follows from its emission state alone, so jumps and
// reversed time need no history; slots are cached while they stay in the span,
// with each grain's last universal anomaly to start the next solution.
class DustTail {
public:
    static constexpr int kSlots = 128;
    static constexpr int kBetas = 16; // strata, log-spaced; each grain draws its beta within one
    static constexpr int kDirections = 4;

    // Grains of comet `body` (a child of the root: the sun, or the planet a
    // scene is centered on) at t. Weights: the production
    // at emission (from the total magnitude law, less the r_h^-2 of sunlight:
    // r_h^-(0.4 K1 - 2)) times the grain size distribution (scattering cross
    // section per log beta ~ beta^0.5 for n(a) ~ a^-3.5, Dohnanyi 1969, JGR 74,
    // 2531), normalized to a sum of 1, then times (r_comet / r_grain)^2 for the
    // sunlight where each grain is. The oldest tenth fades out. Grains are ICRF,
    // so they depend on t alone: a second update at the same t (another output
    // in the same frame, or paused time) keeps them.
    void update(const Scene& scene, int body, double t_tdb);
    const std::vector<DustGrain>& grains() const { return m_grains; }

private:
    struct Grain {
        bool valid = false;
        double emitted_tdb = 0.0;
        double beta = 0.0;
        double production = 0.0;
        double ejection_km_s = 0.0;
        double random = 0.0;
        State start; // heliocentric, at emission, ejection included
        double anomaly = 0.0;
    };
    struct Slot {
        int64_t index = 0;
        std::vector<Grain> grains;
    };
    Slot make_slot(const Scene& scene, int body, int64_t index, double slot_s) const;

    int m_body = -1;
    double m_last_tdb = std::numeric_limits<double>::quiet_NaN(); // of m_grains
    std::deque<Slot> m_slots; // consecutive indices
    std::vector<DustGrain> m_grains;
};

// The tail of a comet orbiting a planet (Shoemaker-Levy 9), after [SEK] Z. Sekanina 1996, in "The Collision of Comet
// Shoemaker-Levy 9 and Jupiter" (Cambridge), pp. 55-80, Sec. 3.4 and 4.4: the fragments were not active. Their tails
// were dust released once, between July and the end of 1992, of 1 to 10 cm and smaller, drifting apart at 0.1 to 0.4 m/s
// and pushed by the sunlight. Relative to the comet a grain follows
//     d'' = G(t) d + beta f(t),
// with G the gravity gradient of the planet and the Sun along the comet's path and f the Sun's gravity (so radiation
// pressure, beta f). That is linear: a grain is d = P(t) v + beta b(t), one P and one b per time of release. They are
// integrated the first time they are needed and kept, so a grain's place is a function of the time alone.
class BoundDustTail {
public:
    static constexpr int kReleases = 8;
    static constexpr int kGrains = 192;
    static constexpr double kBetaMin = 2.0e-5;   // [SEK]: 5e-5 to 2e-4 at 15,000 km from the fragment
    static constexpr double kBetaMax = 2.0e-3;   // Reaches over 100,000 km by January 1994 [SEK Sec. 3.4]
    static constexpr double kSpeedMinKmS = 1.0e-4; // [SEK Sec. 4.4]
    static constexpr double kSpeedMaxKmS = 4.0e-4;

    // Grains of comet `body` at t. Weights: scattering cross section per log beta ~ beta for the n(a) ~ a^-4 of [SEK],
    // normalized to a sum of 1. A grain is left out when its offset is no longer small (the linearization).
    void update(const Scene& scene, int body, double t_tdb);
    const std::vector<DustGrain>& grains() const { return m_grains; }

    // The offset from the comet, ICRF, of a grain released at release_time(release) (after an update or build).
    void build(const Scene& scene, int body);
    glm::dvec3 offset(int release, double beta, const glm::dvec3& velocity, double t_tdb) const;
    double release_time(int release) const { return m_release[static_cast<size_t>(release)]; }
    size_t row_count() const { return m_times.size(); }

private:
    struct Row {
        float p[9]; // column by column
        float b[3];
    };
    struct Seed {
        int release = 0;
        double beta = 0.0;
        glm::dvec3 velocity{0.0};
        double weight = 0.0;
        double random = 0.0;
    };
    const Row& row(int release, size_t index) const { return m_rows[static_cast<size_t>(release) * m_times.size() + index]; }

    int m_body = -1;
    double m_last_tdb = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> m_times;
    std::vector<double> m_release;
    std::vector<size_t> m_first_row; // per release
    std::vector<Row> m_rows;         // kReleases x m_times.size()
    std::vector<Seed> m_seeds;
    std::vector<DustGrain> m_grains;
};

} // namespace astraxis
