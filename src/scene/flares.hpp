#pragma once

#include <cstdint>
#include <optional>
#include <vector>

namespace astraxis {

// Flares of a black hole's hot flow (Sgr A*): a hot spot circling the hole for
// 30-90 minutes (ephem/riaf.hpp HotSpot). Recorded flares come from
// observations; random ones are drawn at an observed average rate from the time
// alone (6-hour slots, each with at most one flare lying wholly inside it), so
// jumps, reversal and pausing always show the same flare at the same moment.
struct FlareRecord {
    double t_peak_tdb = 0.0;
    double duration_s = 3600.0; // full width at half maximum of the brightness envelope
    double radius_m = 9.0;      // orbit radius (M)
    double phase_deg = 0.0;     // spot azimuth at t_peak, in the hole frame of hole_frame_from_pole()
};

struct FlareSettings {
    std::vector<FlareRecord> recorded;
    double random_per_day = 0.0; // 0: no random flares
    uint64_t seed = 0;
};

struct ActiveFlare {
    FlareRecord flare;
    double envelope = 0.0; // 0..1: exp(-4 ln 2 (t - t_peak)^2 / FWHM^2), cut beyond 1.5 FWHM
    bool recorded = false;
};

std::optional<ActiveFlare> active_flare(const FlareSettings& s, double t_tdb);

// Flares at high time warps would flash by: 1 up to 300x (and paused, 0), 0 from 1200x.
double flare_warp_fade(double warp);

} // namespace astraxis
