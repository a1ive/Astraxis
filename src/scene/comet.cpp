#include "scene/comet.hpp"

#include <glm/geometric.hpp>

#include <cmath>

namespace astraxis {

double comet_total_magnitude(double m1, double k1, double r_au)
{
    return m1 + k1 * std::log10(r_au);
}

glm::dvec3 ion_tail_direction(const glm::dvec3& helio_position, const glm::dvec3& helio_velocity,
                              double solar_wind_km_s)
{
    return glm::normalize(solar_wind_km_s * glm::normalize(helio_position) - helio_velocity);
}

} // namespace astraxis
