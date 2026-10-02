#pragma once

namespace astraxis {

// Simulation clock: TDB seconds since J2000, advanced by real time x warp.
struct SimClock {
    double t_tdb = 0.0;
    double warp = 3600.0; // simulated seconds per real second (magnitude)
    bool paused = false;
    bool reverse = false;

    void advance(double real_dt)
    {
        if (!paused) {
            t_tdb += real_dt * warp * (reverse ? -1.0 : 1.0);
        }
    }
};

} // namespace astraxis
