#include "scene/impact.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <glm/geometric.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace astraxis {

namespace {

double clamp01(double x)
{
    return std::clamp(x, 0.0, 1.0);
}

double smooth(double edge0, double edge1, double x)
{
    const double s = clamp01((x - edge0) / (edge1 - edge0));
    return s * s * (3.0 - 2.0 * s);
}

struct Point {
    double x, y;
};

template <size_t N>
double interpolate(const Point (&table)[N], double x)
{
    if (x <= table[0].x) {
        return table[0].y;
    }
    for (size_t i = 1; i < N; ++i) {
        if (x <= table[i].x) {
            return table[i - 1].y + (table[i].y - table[i - 1].y) * (x - table[i - 1].x) / (table[i].x - table[i - 1].x);
        }
    }
    return table[N - 1].y;
}

constexpr double kTwoDaysS = 2.0 * kSecondsPerDay;

} // namespace

const ImpactClassParams& impact_class_params(ImpactClass klass)
{
    // Sizes [HAM95] Table 2 and G's crescent, 13,000 km; what is not from there is a display choice.
    // Flashes where none was measured: about those measured in the class, [CHA95] N for class 3, and
    // [MAR95]'s 1 sigma, 0.5%, for the ones too weak to see. Fireball shares from the light curves' shapes.
    static const ImpactClassParams kC1{13000.0, 1000.0, 1.0, 1.0, 2, 0.0, 0.045, 0.047, 0.85};
    static const ImpactClassParams kC2a{6500.0, 700.0, 0.75, 0.5, 2, 0.0, 0.0, 0.04, 0.7};
    static const ImpactClassParams kC2b{5500.0, 600.0, 0.65, 0.4, 1, 0.0, 0.0, 0.025, 0.6};
    static const ImpactClassParams kC2c{4500.0, 600.0, 0.55, 0.35, 0, 0.0, 0.0, 0.015, 0.2};
    static const ImpactClassParams kC3{0.0, 600.0, 0.45, 0.1, 0, 4.0, 0.0, 0.013, 0.08};
    static const ImpactClassParams kC4{0.0, 0.0, 0.0, 0.02, 0, 0.0, 0.0, 0.005, 0.08};
    switch (klass) {
    case ImpactClass::C1:
        return kC1;
    case ImpactClass::C2a:
        return kC2a;
    case ImpactClass::C2b:
        return kC2b;
    case ImpactClass::C2c:
        return kC2c;
    case ImpactClass::C3:
        return kC3;
    default:
        return kC4;
    }
}

glm::dvec3 impact_plume_axis()
{
    const double tilt = kPlumeTiltDeg * kDegToRad;
    const double azimuth = kPlumeAzimuthDeg * kDegToRad;
    return glm::dvec3(std::sin(tilt) * std::sin(azimuth), std::sin(tilt) * std::cos(azimuth), std::cos(tilt));
}

double impact_plume_vanguard_share()
{
    return std::pow(kPlumeVMaxKmS / kPlumeVMinKmS, -kPlumeAlpha);
}

double impact_plume_top_km(double age_s)
{
    return std::max(0.0, kPlumeVMaxKmS * age_s - 0.5 * kPlumeGravityKmS2 * age_s * age_s);
}

double impact_plume_end_s()
{
    return 2.0 * kPlumeVMaxKmS / kPlumeGravityKmS2;
}

double impact_plume_density(const glm::dvec3& p, double age_s)
{
    if (age_s <= 0.0 || p.z < 0.0) {
        return 0.0;
    }
    // Who is here now left at this velocity: position = v t - g t^2 / 2 z.
    const glm::dvec3 v = (p + glm::dvec3(0.0, 0.0, 0.5 * kPlumeGravityKmS2 * age_s * age_s)) / age_s;
    const double speed = glm::length(v);
    if (speed <= 0.0) {
        return 0.0;
    }
    const double half_cone = 0.5 * kPlumeConeDeg * kDegToRad;
    const double soft = 6.0 * kDegToRad;
    const double cone = smooth(std::cos(half_cone + soft), std::cos(half_cone - soft),
                               glm::dot(v / speed, impact_plume_axis()));
    if (cone <= 0.0) {
        return 0.0;
    }
    const double solid_angle = kTwoPi * (1.0 - std::cos(half_cone));
    // Normalized so that the whole, vanguard included, adds up to 1.
    const double k = std::pow(kPlumeVMinKmS, kPlumeAlpha);
    double f = 0.0;
    if (speed >= kPlumeVMinKmS && speed <= kPlumeVMaxKmS) {
        f = k * kPlumeAlpha * std::pow(speed, -kPlumeAlpha - 3.0) / solid_angle;
    }
    const double width = kPlumeVanguardWidth * kPlumeVMaxKmS;
    const double d = (speed - kPlumeVMaxKmS) / width;
    f += impact_plume_vanguard_share() * std::exp(-0.5 * d * d) /
         (solid_angle * kPlumeVMaxKmS * kPlumeVMaxKmS * width * std::sqrt(kTwoPi));
    return cone * f / (age_s * age_s * age_s);
}

double impact_flash(double age_s, double fireball_share)
{
    if (age_s <= 0.0) {
        return 0.0;
    }
    const double rise = smooth(0.0, 2.0, age_s);
    const double bolide = std::exp(-std::max(age_s - 2.0, 0.0) / 3.0);
    // The stronger the fireball, the longer it lasts: K's 52 s [CHA95] against N's 15.
    const double lasts = 15.0 + 36.0 * fireball_share;
    const double fireball = 1.0 - smooth(0.3 * lasts, lasts, age_s);
    return rise * ((1.0 - fireball_share) * bolide + fireball_share * fireball);
}

double impact_flash_share(const ImpactSite& site)
{
    return site.flash_share > 0.0 ? site.flash_share : impact_class_params(site.klass).flash_share;
}

double impact_flash_temperature_k(double age_s)
{
    // NIMS first saw G about 5 s after the photometer did [CAR95], hence 21 s.
    static const Point kTable[] = {{0.0, std::log(18000.0)}, {5.0, std::log(7600.0)}, {21.0, std::log(2200.0)}, {65.0, std::log(450.0)}};
    return std::exp(interpolate(kTable, age_s));
}

double impact_fireball_height_km(double age_s)
{
    const double t = std::clamp(age_s, 0.0, 40.0);
    return 0.05 * t * t;
}

double impact_fireball_sigma_km(double age_s)
{
    return std::clamp(age_s, 2.0, 60.0);
}

double impact_ring_radius_km(double age_s, bool inner)
{
    return 586.0 + (inner ? 0.189 : 0.454) * age_s;
}

double impact_core_depth(double age_s)
{
    // Days after the impact and the depth in the cores of G and K, [SL98] Table IV; then the
    // value in the band on 1994 Dec 17, by when a site's own pattern has gone.
    static const Point kTable[] = {{0.0, 3.4}, {7.0, 3.4}, {17.0, 2.7}, {41.0, 2.2}, {46.0, 2.25}, {152.0, 0.6}};
    return interpolate(kTable, age_s / kSecondsPerDay);
}

double impact_band_depth(double days)
{
    // 0.6 in December 1994 is [SL98] Table IV, like the rest of the table it is at 890 nm. In visible
    // light the band was detected only until then, with a few darker regions in February 1995 (Sec.
    // 3), so the later values are lower than theirs (0.5, 0.5, 0.3), by choice.
    static const Point kTable[] = {{5.0, 0.0},    {40.0, 0.5},  {154.0, 0.6}, {205.0, 0.25},
                                   {272.0, 0.1}, {710.0, 0.03}, {kImpactBandDays, 0.0}};
    return interpolate(kTable, days);
}

double impact_band_profile(double lat_deg, double days)
{
    // The main part of the band, planetographic -43 to -50 in late August 1994 [SL98] Sec. 3, is -39.4 to
    // -46.1 planetocentric. The debris' limits move at 0.033 deg/day north and 0.004 south (Fig. 4); what
    // spreads out is thin, so that part has a fifth of the depth (a display choice).
    const double main_part = smooth(-48.0, -46.1, lat_deg) * (1.0 - smooth(-39.4, -37.5, lat_deg));
    const double south = -46.1 - 0.004 * days;
    const double north = -39.4 + 0.033 * days;
    const double spread = smooth(south - 2.0, south, lat_deg) * (1.0 - smooth(north, north + 2.0, lat_deg));
    return std::max(main_part, 0.2 * spread);
}

double impact_site_fade(const ImpactClassParams& params, double age_s)
{
    if (params.fade_days > 0.0) {
        return std::exp(-age_s / (params.fade_days * kSecondsPerDay));
    }
    return 1.0 - smooth(46.0 * kSecondsPerDay, 152.0 * kSecondsPerDay, age_s);
}

double impact_wind_km_s(double north_km)
{
    // 4 m/s per 1000 km, with extremes of -10 and 20 m/s [SL98] Table III; east to the north of the
    // sites, west to the south [BEE96] Sec. 3.2.1.
    return std::clamp(4.0e-6 * north_km, -0.010, 0.020);
}

double impact_scar_shape(const ImpactClassParams& params, double east_km, double north_km, double age_s, double seed)
{
    if (age_s <= 0.0) {
        return 0.0;
    }
    const double azimuth = kPlumeAzimuthDeg * kDegToRad; // the same for the plume and the crescent
    const double ax = std::sin(azimuth);
    const double ay = std::cos(azimuth);
    const double r = std::hypot(east_km, north_km);
    const double x = east_km * ax + north_km * ay; // toward the crescent
    const double y = -east_km * ay + north_km * ax;
    const double ejecta = params.ejecta_km + params.expansion_km_s * std::min(age_s, kTwoDaysS);
    // Comes down over the plume's flight: the first falls back after about 6 min [HD01], the most by 17.
    const double landed = smooth(300.0, 1100.0, age_s);
    double sum = 0.0;

    const double wave = 0.025 * params.ejecta_km;
    const double xw = wave > 0.0 ? x + wave * std::sin(y / (3.5 * wave) + seed) : x;
    const double yw = wave > 0.0 ? y + wave * std::sin(x / (2.8 * wave) + 1.7 * seed) : y;

    const double grown = params.ejecta_km > 0.0 ? ejecta / params.ejecta_km : 1.0;
    const double core = params.core_km * (1.0 + (grown - 1.0) * clamp01(age_s / kTwoDaysS));
    double dark = std::exp(-0.5 * (xw * xw + yw * yw) / (core * core));

    if (params.ejecta_km > 0.0) {
        const double e = params.ejecta_km;
        // The streak, from the center toward the crescent [HAM95].
        const double across = yw / (0.035 * e);
        dark = std::max(dark, 0.85 * smooth(-0.02 * e, 0.0, xw) * (1.0 - smooth(0.22 * e, 0.36 * e, xw)) * std::exp(-0.5 * across * across));
        // The crescent: at least 180 deg around, with rays from a point 1000-2000 km toward it [HAM95]. It is optically
        // thin: optical depth about 0.8 at 890 nm against the core's 3.4 ([WEST] Fig. 3, by tools/sl9/debris_optics.py).
        const double x0 = xw - 0.115 * e;
        const double rd = std::hypot(x0, yw);
        const double az = std::atan2(yw, x0);
        const double spread = az / 1.6;
        const double rr = (rd / ejecta - 0.76) / (0.1 * std::max(1.0 - 0.14 * az * az, 0.3));
        const double rays = 1.0 + 0.25 * std::sin(11.0 * az + seed + 2.0 * std::sin(3.0 * az + 1.3 * seed));
        sum += 0.22 * rays * std::exp(-spread * spread - 0.5 * rr * rr);
        // Fainter traces all around.
        const double halo = r / ejecta;
        sum += 0.05 * std::exp(-3.0 * halo * halo);
    }
    sum = (sum + dark) * landed;

    // Rings, sharpest soon after the impact and gone within hours.
    if (params.rings > 0 && age_s < 10.0 * 3600.0) {
        const double fade = std::exp(-age_s / (3.0 * 3600.0));
        const double outer = (r - impact_ring_radius_km(age_s, false)) / 120.0;
        sum += 0.12 * fade * smooth(120.0, 600.0, age_s) * std::exp(-outer * outer);
        if (params.rings > 1) {
            const double inner = (r - impact_ring_radius_km(age_s, true)) / 90.0;
            // Mainly to the northwest of the center [HAM96].
            sum += 0.07 * fade * smooth(120.0, 600.0, age_s) * std::exp(-inner * inner) * smooth(0.0, 0.5, -x / std::max(r, 1.0));
        }
    }
    return sum;
}

namespace {
double debris_air_mass(double mu, double mu0)
{
    return 1.0 / std::max(mu, 0.1) + 1.0 / std::max(mu0, 0.1);
}
} // namespace

glm::dvec3 impact_debris_factor(double depth, double mu, double mu0)
{
    const double through = std::exp(-kDebrisDimming * depth * debris_air_mass(mu, mu0));
    return through + (1.0 - through) * glm::dvec3(kDebrisRatioRed, kDebrisRatioGreen, kDebrisRatioBlue);
}

double impact_seed(const ImpactSite& site)
{
    return std::fmod(site.west_lon_deg * 0.37, kTwoPi);
}

double impact_scar_depth(const std::vector<ImpactSite>& sites, double lat_deg, double west_lon_deg, double t_tdb)
{
    double depth = 0.0;
    double first = std::numeric_limits<double>::infinity();
    const double cos_lat = std::cos(lat_deg * kDegToRad);
    const double circle = kTwoPi * kImpactRadiusKm * cos_lat;
    for (const ImpactSite& site : sites) {
        first = std::min(first, site.t_tdb);
        const double age = t_tdb - site.t_tdb;
        const ImpactClassParams& params = impact_class_params(site.klass);
        if (age <= 0.0 || params.weight <= 0.0) {
            continue;
        }
        const double north = (lat_deg - site.lat_deg) * kDegToRad * kImpactRadiusKm;
        const double east = -std::remainder(west_lon_deg - site.west_lon_deg, 360.0) * kDegToRad * kImpactRadiusKm * cos_lat;
        const double sheared = std::remainder(east - impact_wind_km_s(north) * age, circle);
        depth += params.weight * impact_core_depth(age) * impact_site_fade(params, age) *
                 impact_scar_shape(params, sheared, north, age, impact_seed(site));
    }
    if (std::isfinite(first) && t_tdb > first) {
        const double days = (t_tdb - first) / kSecondsPerDay;
        depth += impact_band_depth(days) * impact_band_profile(lat_deg, days);
    }
    return depth;
}

} // namespace astraxis
