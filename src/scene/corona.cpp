#include "scene/corona.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <algorithm>
#include <cmath>

namespace astraxis {

double corona_electron_density(double r)
{
    if (r < 1.0 || !std::isfinite(r)) {
        return 0.0;
    }
    // Reginald et al. 2023, Solar Physics 298:73, Section 3.1, p. 8 (Leblanc profile).
    const double q = 1.0 / (r * r);
    return q * (3.3e5 + q * (4.1e6 + q * 8.0e7));
}

glm::dvec2 corona_scattering(double r, double u)
{
    if (r < 1.0 || !std::isfinite(r)) {
        return {};
    }
    // Azimuthally integrated Thomson cross section: 3*sigma_T/8 times
    // [(3/2-mu^2/2) + (3*mu^2/2-1/2)*cos(chi)^2]. Here mu is the
    // incident ray's cosine to the radial direction. Solar surface cosine
    // is sqrt(1-r^2*(1-mu^2)); divide by disk-mean intensity 1-u/3.
    // CODATA 2022 sigma_T, m^2; cm^-3 * km converts to m^-2 by 1e9.
    constexpr double sigma_t = 6.6524587051e-29;
    constexpr int samples = 256;
    const double mu0 = std::sqrt(std::max(0.0, 1.0 - 1.0 / (r * r)));
    const double span = (1.0 / (r * r)) / (1.0 + mu0);
    glm::dvec2 moments(0.0);
    for (int i = 0; i < samples; ++i) {
        const double mu = mu0 + span * (i + 0.5) / samples;
        const double surface_mu = std::sqrt(std::clamp(1.0 - r * r * (1.0 - mu * mu), 0.0, 1.0));
        const double intensity = (1.0 - u + u * surface_mu) / (1.0 - u / 3.0);
        moments += intensity * glm::dvec2(1.5 - 0.5 * mu * mu, 1.5 * mu * mu - 0.5);
    }
    return moments * (span / samples * 3.0 * sigma_t / 8.0 * kSunRadiusKm * 1e9 *
                      corona_electron_density(r));
}

double corona_activity(double t_tdb)
{
    // NASA 2020 "What Will Solar Cycle 25 Look Like?": April 2014 maximum,
    // December 2019 minimum. NASA/NOAA 2024-10-15 announcement: maximum phase.
    // The last anchor represents that phase, not a fitted sunspot-number peak.
    static const double dates[] = {
        tdb_from_utc_jd(jd_from_calendar({2014, 4, 1, 0, 0, 0})),
        tdb_from_utc_jd(jd_from_calendar({2019, 12, 1, 0, 0, 0})),
        tdb_from_utc_jd(jd_from_calendar({2024, 10, 15, 0, 0, 0})),
    };
    if (t_tdb <= dates[0] || t_tdb >= dates[2]) {
        return 1.0;
    }
    const int i = t_tdb < dates[1] ? 0 : 1;
    double f = (t_tdb - dates[i]) / (dates[i + 1] - dates[i]);
    f = f * f * (3.0 - 2.0 * f);
    return i == 0 ? 1.0 - f : f;
}

double corona_visible_fraction(double sun_radius, double occulter_radius, double separation)
{
    if (!(sun_radius > 0.0) || !(occulter_radius > 0.0) || separation >= sun_radius + occulter_radius) {
        return 1.0;
    }
    if (occulter_radius >= sun_radius + separation) {
        return 0.0;
    }
    constexpr int rings = 128;
    double visible = 0.0;
    double total = 0.0;
    for (int i = 0; i < rings; ++i) {
        const double q = (i + 0.5) / rings; // equal-area annuli
        const double r = sun_radius * std::sqrt(q);
        const double intensity = 0.4 + 0.6 * std::sqrt(1.0 - q);
        double covered = 0.0;
        if (separation < 1e-15) {
            covered = r < occulter_radius ? 1.0 : 0.0;
        } else {
            const double c = (r * r + separation * separation - occulter_radius * occulter_radius) /
                             (2.0 * r * separation);
            covered = std::acos(std::clamp(c, -1.0, 1.0)) / kPi;
        }
        visible += intensity * (1.0 - covered);
        total += intensity;
    }
    return visible / total;
}

} // namespace astraxis
