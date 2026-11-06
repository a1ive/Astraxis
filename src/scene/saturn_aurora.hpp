#pragma once

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

namespace astraxis {

// Carbary 2012, doi:10.1029/2012JA017670, equation 2 / table 1.
// Statistical 2007-2009 UV morphology, in solar local time, not System III.
struct SaturnAurora {
    static constexpr double kPeakKm = 1100.0; // Gerard 2009, doi:10.1029/2008GL036554
    static constexpr double kScaleKm = 200.0; // assumed Chapman scale
    static constexpr double kBottomKm = 500.0;
    static constexpr double kTopKm = 3000.0;
    // Visible layers, Dyudina 2015 Fig. 8 and Appendix A.2 (tools/geomag/saturn_visible_colors.py):
    // I/F in RED, GRN, BL1 per unit clear-filter brightness of the "brightest" and "high" boxes,
    // and the high layer's column for its 11 kR against the low layer's 150 kR, edge-on.
    static constexpr glm::dvec3 kVisibleLowRgb{1.000, 0.491, 0.692};
    static constexpr glm::dvec3 kVisibleHighRgb{0.930, 0.420, 1.384};
    static constexpr double kVisibleHighColumn = 0.075;
    // North/south planetocentric colatitude, Gaussian sigma (all radians),
    // dawn enhancement (0..1). Width and brightness interpolation are assumed.
    static glm::dvec4 profile(double local_hours);
    static glm::dvec3 noon(const glm::dvec3& sun_body);
};

} // namespace astraxis
