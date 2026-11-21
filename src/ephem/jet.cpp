#include "ephem/jet.hpp"

#include "core/math.hpp"

#include <glm/geometric.hpp>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr double kRidgePsi = 0.82;   // outer half maximum at psi = 0.82 + 1.1774 sigma = 1
constexpr double kRidgeSigma = 0.15;
constexpr double kSpineShare = 0.3;  // the spine's peak relative to the ridge's
constexpr double kSpineSigma = 0.2;
constexpr double kProfileEnd = 1.4;  // nothing beyond 1.4 radii

double smoothstep(double e0, double e1, double x)
{
    const double t = std::clamp((x - e0) / (e1 - e0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

// Soft-edged box over [a, b]: the edges blend over +-30% of its length.
double soft_box(double a, double b, double x)
{
    const double h = 0.3 * (b - a);
    return smoothstep(a - h, a + h, x) * (1.0 - smoothstep(b - h, b + h, x));
}

} // namespace

glm::dvec3 jet_feature_position(const JetAstrometry& a, const JetFeature& f, double year)
{
    const double t = std::clamp(year, f.first_year, f.last_year);
    // The tabulated position and the fitted velocity refer to different epochs.
    const double dt = t - f.first_year;
    const glm::dvec2 mean_motion = f.motion_mas_yr + (0.5 * dt + f.first_year - a.velocity_year) * f.acceleration_mas_yr2;
    const glm::dvec2 xy = f.offset_arcsec + (0.001 * dt) * mean_motion;
    return xy.x * a.along_m + xy.y * a.across_m;
}

double jet_feature_visibility(const JetAstrometry& a, const JetFeature& f, double year)
{
    double v = 1.0;
    if (f.first_year > a.first_year) {
        v *= smoothstep(f.first_year - kJetFeatureRampYears, f.first_year, year);
    }
    if (f.last_year < a.last_year) {
        v *= 1.0 - smoothstep(f.last_year, f.last_year + kJetFeatureRampYears, year);
    }
    return v;
}

double jet_knot_diffuse_share(const JetAstrometry& a, const JetKnot& k, double year)
{
    double total = 0.0, shown = 0.0;
    for (const JetFeature& f : k.features) {
        total += f.weight;
        shown += f.weight * jet_feature_visibility(a, f, year);
    }
    return total > 0.0 ? 1.0 - k.feature_fraction * shown / total : 1.0;
}

double jet_ridge_offset(const JetParams& p, double zeta)
{
    const std::vector<glm::dvec2>& r = p.ridge;
    if (r.empty() || zeta <= r.front().x) {
        return r.empty() ? 0.0 : r.front().y;
    }
    for (size_t i = 1; i < r.size(); ++i) {
        if (zeta < r[i].x) {
            const double u = (zeta - r[i - 1].x) / (r[i].x - r[i - 1].x);
            return r[i - 1].y + u * (r[i].y - r[i - 1].y);
        }
    }
    return r.back().y;
}

double jet_ridge_slope_bound(const JetParams& p)
{
    double k = 0.0;
    for (size_t i = 1; i < p.ridge.size(); ++i) {
        k = std::max(k, std::abs(p.ridge[i].y) / (p.ridge[i].x - p.ridge.front().x));
    }
    return k;
}

double jet_profile(double psi)
{
    if (psi > kProfileEnd) {
        return 0.0;
    }
    const double r = (psi - kRidgePsi) / kRidgeSigma;
    const double c = psi / kSpineSigma;
    return std::exp(-0.5 * r * r) + kSpineShare * std::exp(-0.5 * c * c);
}

double jet_profile_area()
{
    static const double area = [] {
        double sum = 0.0;
        constexpr int kSteps = 4000;
        for (int i = 0; i < kSteps; ++i) {
            const double psi = (i + 0.5) * kProfileEnd / kSteps;
            sum += jet_profile(psi) * psi;
        }
        return kTwoPi * sum * kProfileEnd / kSteps;
    }();
    return area;
}

double jet_radius(const JetParams& p, double zeta)
{
    return std::sqrt(p.parabola_k_m * zeta + p.cone_slope * p.cone_slope * zeta * zeta);
}

double jet_radius_slope(const JetParams& p, double zeta)
{
    return (p.parabola_k_m + 2.0 * p.cone_slope * p.cone_slope * zeta) / (2.0 * std::max(jet_radius(p, zeta), 1e-12));
}

double jet_gamma_beta(const JetParams& p, double zeta)
{
    const auto accelerating = [&](double z) {
        const double gamma = std::pow(std::max(z, 1e-12) / p.lorentz_unity_m, p.lorentz_index);
        return std::max(std::sqrt(std::max(gamma * gamma - 1.0, 0.0)), p.gamma_beta_min);
    };
    if (zeta <= p.accel_end_m) {
        return accelerating(zeta);
    }
    return accelerating(p.accel_end_m) * std::pow(zeta / p.accel_end_m, p.decel_index);
}

double jet_tilt_fraction(const JetPrecession& p, double zeta)
{
    const double l = std::log(std::max(zeta, 1e-12));
    return smoothstep(std::log(p.tilt_m[0]), std::log(p.tilt_m[1]), l) *
           (1.0 - smoothstep(std::log(p.tilt_m[2]), std::log(p.tilt_m[3]), l));
}

double jet_line_brightness(const JetParams& p, double zeta, bool exclude_features)
{
    const double alpha = p.spectral_index;
    if (p.knots.empty()) {
        return std::pow(jet_radius(p, zeta), -(1.0 + alpha));
    }
    const JetKnot& first = p.knots.front();
    const double h = 0.3 * (first.end_m - first.start_m);
    double sum = first.line_brightness * std::pow(jet_radius(p, zeta) / jet_radius(p, first.start_m), -(1.0 + alpha)) *
                 (1.0 - smoothstep(first.start_m - h, first.start_m + h, zeta));
    for (const JetKnot& k : p.knots) {
        const double share = exclude_features ? 1.0 - k.feature_fraction : 1.0;
        sum += share * k.line_brightness * soft_box(k.start_m, k.end_m, zeta);
    }
    return sum;
}

double jet_doppler(const JetParams& p, double zeta, double theta)
{
    const double gb = jet_gamma_beta(p, zeta);
    const double gamma = std::sqrt(1.0 + gb * gb);
    return 1.0 / (gamma - gb * std::cos(theta));
}

void prepare_jet(JetParams& p, double theta_earth)
{
    const double boost = 2.0 + p.spectral_index;
    double first = 0.0;
    for (JetKnot& k : p.knots) {
        const double mid = 0.5 * (k.start_m + k.end_m);
        k.line_brightness = k.flux / (k.end_m - k.start_m) / std::pow(jet_doppler(p, mid, theta_earth), boost);
        if (first == 0.0) {
            first = k.line_brightness;
        }
    }
    for (JetKnot& k : p.knots) {
        k.line_brightness /= first; // the first knot at 1
        double total_weight = 0.0;
        for (const JetFeature& f : k.features) {
            total_weight += f.weight;
        }
        for (JetFeature& f : k.features) {
            f.flux = k.flux * k.feature_fraction / first * (f.weight / total_weight);
        }
    }
}

double jet_emissivity(const JetParams& p, const glm::dvec3& x, const glm::dvec3& to_camera, double tan_tilt,
                      double min_width_m, glm::dvec2 outer_slope, glm::dvec2 across_m)
{
    const double zeta = std::abs(x.z);
    if (zeta < p.base_m || zeta > p.length_m) {
        return 0.0;
    }
    const double side = x.z >= 0.0 ? 1.0 : -1.0;
    const double lean = p.precesses ? std::tan(std::atan(tan_tilt) * jet_tilt_fraction(p.precession, zeta)) : 0.0;
    const double split = p.precesses ? p.precession.tilt_m[3] : 0.0;
    const glm::dvec2 ridge = side > 0.0 ? jet_ridge_offset(p, zeta) * across_m : glm::dvec2(0.0);
    const glm::dvec2 lateral =
        glm::dvec2(x.x - x.z * lean, x.y) - side * std::max(zeta - split, 0.0) * outer_slope - ridge;
    const double s = glm::length(lateral);
    const double w = jet_radius(p, zeta);
    const double width = std::sqrt(w * w + min_width_m * min_width_m);
    const double psi = s / width;
    if (psi > kProfileEnd) {
        return 0.0;
    }
    const double fade = smoothstep(p.base_m, 4.0 * p.base_m, zeta) * (1.0 - smoothstep(0.9 * p.length_m, p.length_m, zeta));
    const double cross = jet_profile(psi) / (jet_profile_area() * width * width);

    // Along the local streamline: the leaning axis, spreading as the radius grows.
    glm::dvec3 dir = side * glm::dvec3(lean, 0.0, 1.0);
    if (zeta > split) {
        dir += side * glm::dvec3(outer_slope, 0.0);
    }
    if (s > 0.0) {
        dir += glm::dvec3(lateral / s, 0.0) * (s / w) * jet_radius_slope(p, zeta);
    }
    dir = glm::normalize(dir);
    const double gb = jet_gamma_beta(p, zeta);
    const double gamma = std::sqrt(1.0 + gb * gb);
    const double doppler = 1.0 / (gamma - gb * glm::dot(dir, to_camera));
    return fade * jet_line_brightness(p, zeta, side > 0.0) * cross * std::pow(doppler, 2.0 + p.spectral_index);
}

SkyAxes sky_axes(double ra_deg, double dec_deg)
{
    const double a = ra_deg * kDegToRad;
    const double d = dec_deg * kDegToRad;
    SkyAxes out;
    out.toward_observer = -glm::dvec3(std::cos(d) * std::cos(a), std::cos(d) * std::sin(a), std::sin(d));
    out.east = glm::dvec3(-std::sin(a), std::cos(a), 0.0);
    out.north = glm::dvec3(-std::sin(d) * std::cos(a), -std::sin(d) * std::sin(a), std::cos(d));
    return out;
}

double sky_position_angle(const glm::dvec3& v, const SkyAxes& sky)
{
    return std::atan2(glm::dot(v, sky.east), glm::dot(v, sky.north));
}

glm::dvec3 jet_nozzle(const JetPrecession& p, const glm::dvec3& axis, const SkyAxes& sky, double year)
{
    const double theta = std::acos(std::clamp(glm::dot(axis, sky.toward_observer), -1.0, 1.0));
    const double eta = sky_position_angle(axis, sky);
    const double phi = p.half_angle_deg * kDegToRad;
    const double turn = -p.omega_rad_per_year * (year - p.t0_year);
    // Precession frame, then R_y(theta), then R_z(eta_p).
    const glm::dvec3 jp(std::sin(phi) * std::cos(turn), std::sin(phi) * std::sin(turn), std::cos(phi));
    const glm::dvec3 jy(jp.x * std::cos(theta) + jp.z * std::sin(theta), jp.y,
                        -jp.x * std::sin(theta) + jp.z * std::cos(theta));
    const glm::dvec3 jo(jy.x * std::cos(eta) - jy.y * std::sin(eta), jy.x * std::sin(eta) + jy.y * std::cos(eta), jy.z);
    return glm::normalize(jo.x * sky.north + jo.y * sky.east + jo.z * sky.toward_observer);
}

} // namespace astraxis
