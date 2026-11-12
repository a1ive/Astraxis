#include "scene/magnetosphere.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <locale>
#include <sstream>

namespace astraxis {

namespace {

constexpr double kInfinity = std::numeric_limits<double>::infinity();
// B^2 / (2 mu0) in nPa for B in nT.
constexpr double kMagneticPressureNpaPerNt2 = 1e-9 / (2.0 * 1.25663706212e-6);

// Joy et al. 2002, Figure 3: each coefficient of their equation 2 is c0 + c1 P^(-1/4) for
// A, B, C and c0 + c1 P for D, E, F, with P in nPa and lengths in units of 120 Jupiter radii.
constexpr double kJoyScale = 120.0;
constexpr double kJoyMagnetopause[6][2] = {{-0.134, 0.488},  {-0.581, -0.225}, {-0.186, -0.016},
                                           {-0.014, 0.096},  {-0.814, -0.811}, {-0.050, 0.168}};
constexpr double kJoyBowShock[6][2] = {{-1.107, 1.591}, {-0.566, -0.812}, {0.048, -0.059},
                                       {0.077, -0.038}, {-0.874, -0.299}, {-0.055, 0.124}};

BoundarySurface joy_surface(const double table[6][2], double pressure_npa)
{
    BoundarySurface s;
    s.shape = BoundarySurface::Shape::Quadric;
    s.scale = kJoyScale;
    const double p = std::max(pressure_npa, 1e-4);
    const double root = std::pow(p, -0.25);
    for (int i = 0; i < 6; ++i) {
        s.coefficients[i] = table[i][0] + table[i][1] * (i < 3 ? root : p);
    }
    return s;
}

// Half-width of a flared surface at x = -downtail, without the cylinder.
double flared_half_width(double standoff, double flaring, double downtail)
{
    // x(theta) falls monotonically from the nose; look for it between 90 deg and 180.
    double low = 0.5 * kPi;
    double high = kPi - 1e-9;
    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (low + high);
        const double r = standoff * std::pow(2.0 / (1.0 + std::cos(mid)), flaring);
        (r * std::cos(mid) > -downtail ? low : high) = mid;
    }
    const double r = standoff * std::pow(2.0 / (1.0 + std::cos(low)), flaring);
    return r * std::sin(low);
}

} // namespace

bool SolarWindRecord::load(const std::filesystem::path& file, std::string* error)
{
    const auto failed = [&](const std::string& what) {
        if (error) {
            *error = file.string() + ": " + what;
        }
        return false;
    };
    std::ifstream in(file, std::ios::binary);
    if (!in) {
        return failed("cannot open");
    }
    const std::vector<char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    uint16_t version = 0;
    uint16_t first_year = 0;
    uint32_t count = 0;
    if (bytes.size() < 12 || std::memcmp(bytes.data(), "AXSW", 4) != 0) {
        return failed("not a solar wind record (AXSW)");
    }
    std::memcpy(&version, bytes.data() + 4, 2);
    std::memcpy(&first_year, bytes.data() + 6, 2);
    std::memcpy(&count, bytes.data() + 8, 4);
    static_assert(sizeof(Hour) == 4);
    if (version != 1 || count == 0 || bytes.size() != 12 + sizeof(Hour) * static_cast<size_t>(count)) {
        return failed("unsupported version or wrong size");
    }
    m_hours.resize(count);
    std::memcpy(m_hours.data(), bytes.data() + 12, sizeof(Hour) * m_hours.size());
    // January 1 of the first year, 00:00 UTC (Gregorian: days before it since
    // JD 1721425.5, 0001-01-01).
    const double y = first_year - 1.0;
    m_start_jd_utc = 1721425.5 + 365.0 * y + std::floor(y / 4.0) - std::floor(y / 100.0) + std::floor(y / 400.0);
    return true;
}

double SolarWindRecord::hour(double t_tdb) const
{
    return (utc_jd_from_tdb(t_tdb) - m_start_jd_utc) * 24.0;
}

double SolarWindRecord::start_tdb() const
{
    return tdb_from_utc_jd(m_start_jd_utc);
}

double SolarWindRecord::end_tdb() const
{
    return tdb_from_utc_jd(m_start_jd_utc + static_cast<double>(m_hours.size()) / 24.0);
}

SolarWind SolarWindRecord::mean(double t0_tdb, double t1_tdb, const SolarWind& fallback) const
{
    if (t1_tdb < t0_tdb) {
        std::swap(t0_tdb, t1_tdb);
    }
    const double x0 = hour(t0_tdb);
    const double x1 = std::max(hour(t1_tdb), x0 + 1e-6);
    const double n = static_cast<double>(m_hours.size());
    // Sums over the measured hours within the span, and how much of it they are.
    double sum[3]{}, measured[3]{};
    const size_t last = static_cast<size_t>(std::clamp(std::ceil(x1), 0.0, n));
    for (size_t k = static_cast<size_t>(std::clamp(std::floor(x0), 0.0, n)); k < last; ++k) {
        const double weight = std::min(x1, k + 1.0) - std::max(x0, static_cast<double>(k));
        const Hour& h = m_hours[k];
        if (h.pressure != 0xffff) {
            sum[0] += weight * h.pressure * 0.01;
            measured[0] += weight;
        }
        if (h.bz != -128) {
            sum[1] += weight * h.bz * 0.5;
            measured[1] += weight;
        }
        if (h.speed != 0xff) {
            sum[2] += weight * (150.0 + 8.0 * h.speed);
            measured[2] += weight;
        }
    }
    const double span = x1 - x0;
    SolarWind wind = fallback;
    // Never zero, which would put the boundary at infinity: the record has hours of 0.00 nPa
    // (in 1979).
    wind.pressure_npa = std::max((sum[0] + (span - measured[0]) * fallback.pressure_npa) / span, 0.02);
    wind.bz_nt = (sum[1] + (span - measured[1]) * fallback.bz_nt) / span;
    wind.speed_km_s = (sum[2] + (span - measured[2]) * fallback.speed_km_s) / span;
    return wind;
}

double SolarWindRecord::coverage(double t0_tdb, double t1_tdb) const
{
    if (t1_tdb < t0_tdb) {
        std::swap(t0_tdb, t1_tdb);
    }
    const double x0 = hour(t0_tdb);
    const double x1 = std::max(hour(t1_tdb), x0 + 1e-6);
    const double n = static_cast<double>(m_hours.size());
    double measured = 0.0;
    const size_t last = static_cast<size_t>(std::clamp(std::ceil(x1), 0.0, n));
    for (size_t k = static_cast<size_t>(std::clamp(std::floor(x0), 0.0, n)); k < last; ++k) {
        if (m_hours[k].pressure != 0xffff) {
            measured += std::min(x1, k + 1.0) - std::max(x0, static_cast<double>(k));
        }
    }
    return measured / (x1 - x0);
}

bool BoundaryTable::load(const std::filesystem::path& file, std::string* error)
{
    const auto failed = [&](const std::string& what) {
        if (error) {
            *error = file.string() + ": " + what;
        }
        return false;
    };
    std::ifstream in(file);
    if (!in) {
        return failed("cannot open");
    }
    std::string line;
    std::getline(in, line);
    std::istringstream header(line);
    header.imbue(std::locale::classic());
    std::string magic;
    double step_deg = 0.0;
    if (!(header >> magic >> step_deg) || magic != "AXBOUNDARY1" || !(step_deg > 0.0 && step_deg <= 15.0)) {
        return failed("not a boundary table (AXBOUNDARY1)");
    }
    m_step = step_deg * kDegToRad;
    int found = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '\r') {
            continue;
        }
        std::istringstream row(line);
        row.imbue(std::locale::classic());
        std::string surface, state;
        int azimuth = -1;
        Profile profile;
        if (!(row >> surface >> state >> azimuth >> profile.flaring) || (surface != "shock" && surface != "magnetopause") ||
            (state != "earth_like" && state != "pole_on") || azimuth < 0 || azimuth > 270 || azimuth % 90 != 0 ||
            !(profile.flaring > 0.0 && profile.flaring < 2.0)) {
            return failed("bad line: " + line.substr(0, 40));
        }
        for (double r = 0.0; row >> r;) {
            if (!(r > 0.0)) {
                return failed("distances must be positive");
            }
            profile.radii.push_back(r);
        }
        Profile& slot = m_profiles[surface == "shock" ? Shock : Magnetopause][state == "pole_on" ? 1 : 0][azimuth / 90];
        if (profile.radii.size() < 4 || !slot.radii.empty()) {
            return failed("a profile is too short or given twice");
        }
        slot = std::move(profile);
        ++found;
    }
    if (found != 16) {
        return failed("needs the shock and the magnetopause in two states at four azimuths");
    }
    return true;
}

double BoundaryTable::at(const Profile& profile, double theta) const
{
    const double x = theta / m_step;
    const size_t last = profile.radii.size() - 1;
    if (x < static_cast<double>(last)) {
        const size_t i = static_cast<size_t>(x);
        return profile.radii[i] + (x - static_cast<double>(i)) * (profile.radii[i + 1] - profile.radii[i]);
    }
    const double facing = 1.0 + std::cos(theta);
    if (!(facing > 1e-9)) {
        return kInfinity;
    }
    return profile.radii[last] * std::pow((1.0 + std::cos(static_cast<double>(last) * m_step)) / facing, profile.flaring);
}

double BoundaryTable::radius(Surface surface, double blend, const glm::dvec3& d,
                             const BoundaryAzimuth& azimuth) const
{
    const double theta = std::acos(std::clamp(d.x, -1.0, 1.0));
    double r[4];
    for (int k = 0; k < 4; ++k) {
        r[k] = glm::mix(at(m_profiles[surface][0][k], theta), at(m_profiles[surface][1][k], theta), blend);
    }
    if (!std::isfinite(r[0] + r[1] + r[2] + r[3])) {
        return kInfinity;
    }
    // Evaluate the harmonics at the azimuth before rotation, then reflect their north-south term.
    const double across = std::hypot(d.y, d.z);
    const double c = across > 1e-12 ? (azimuth.rotation.x * d.y + azimuth.rotation.y * d.z) / across : 1.0;
    const double s = across > 1e-12 ? (azimuth.rotation.x * d.z - azimuth.rotation.y * d.y) / across : 0.0;
    const double mean = 0.25 * (r[0] + r[1] + r[2] + r[3]);
    const double variation = 0.5 * (r[0] - r[2]) * c + 0.5 * (r[1] - r[3]) * azimuth.north_south * s +
                             0.25 * (r[0] - r[1] + r[2] - r[3]) * (c * c - s * s);
    return mean + azimuth.strength * variation;
}

double BodyMagnetosphere::pressure_for_standoff(double radii) const
{
    return surface_field_nt * surface_field_nt * kMagneticPressureNpaPerNt2 / std::pow(radii, compressibility);
}

double BoundarySurface::radius(const glm::dvec3& d) const
{
    switch (shape) {
    case Shape::Flared: {
        if (!(1.0 + d.x > 1e-12)) {
            return kInfinity;
        }
        const double r = standoff * std::pow(2.0 / (1.0 + d.x), flaring);
        if (flare_end > 0.0 && r * d.x < -flare_end) {
            const double s = std::sqrt(std::max(1.0 - d.x * d.x, 0.0));
            return s > 1e-12 ? tail_radius / s : kInfinity;
        }
        return r;
    }
    case Shape::Conic: {
        if (focus == 0.0) {
            const double denominator = 1.0 + eccentricity * d.x;
            return denominator > 1e-12 ? standoff * (1.0 + eccentricity) / denominator : kInfinity;
        }
        // |q| + e q.x = l with q = r d - (focus, 0, 0): a quadratic in r.
        const double l = (standoff - focus) * (1.0 + eccentricity);
        const double m = l + eccentricity * focus;
        const double a = 1.0 - eccentricity * eccentricity * d.x * d.x;
        const double b = 2.0 * d.x * (m * eccentricity - focus);
        const double c = focus * focus - m * m; // negative: the center is inside
        if (std::abs(a) < 1e-12) {
            return b > 1e-12 ? -c / b : kInfinity;
        }
        const double discriminant = b * b - 4.0 * a * c;
        if (discriminant < 0.0) {
            return kInfinity;
        }
        // The root on the branch around the focus: l - e q.x >= 0.
        for (const double r : {(-b + std::sqrt(discriminant)) / (2.0 * a), (-b - std::sqrt(discriminant)) / (2.0 * a)}) {
            if (r > 0.0 && m - eccentricity * r * d.x >= 0.0) {
                return r;
            }
        }
        return kInfinity;
    }
    case Shape::Tabulated:
        return table ? scale * table->radius(which, blend, d, azimuth) : kInfinity;
    case Shape::Quadric: {
        // a r^2 + b r + A = 0 along the direction, with A > 0: the positive root nearest the planet.
        const double* c = coefficients;
        const double a = c[2] * d.x * d.x + c[4] * d.y * d.y + c[5] * d.x * d.y - d.z * d.z;
        const double b = c[1] * d.x + c[3] * d.y;
        const double discriminant = b * b - 4.0 * a * c[0];
        if (!(c[0] > 0.0) || discriminant < 0.0) {
            return kInfinity;
        }
        const double denominator = -b + std::sqrt(discriminant);
        return denominator > 1e-12 ? scale * 2.0 * c[0] / denominator : kInfinity;
    }
    }
    return kInfinity;
}

bool BoundarySurface::contains(const glm::dvec3& point) const
{
    const double r = glm::length(point);
    return r == 0.0 || r < radius(point / r);
}

double BoundarySurface::tail_half_width(double downtail) const
{
    double low = 0.5 * kPi;
    double high = kPi - 1e-9;
    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (low + high);
        const double r = radius({std::cos(mid), std::sin(mid), 0.0});
        (r * std::cos(mid) > -downtail ? low : high) = mid;
    }
    return radius({std::cos(low), std::sin(low), 0.0}) * std::sin(low);
}

BoundarySurface magnetopause(const BodyMagnetosphere& body, const SolarWind& wind)
{
    const double pressure = std::max(wind.pressure_npa, 1e-6);
    BoundarySurface s;
    switch (body.model) {
    case MagnetosphereModel::Joy2002:
        return joy_surface(kJoyMagnetopause, pressure);
    case MagnetosphereModel::Tabulated:
        s.shape = BoundarySurface::Shape::Tabulated;
        s.table = body.boundaries;
        s.which = BoundaryTable::Magnetopause;
        s.blend = wind.phase;
        s.azimuth = wind.azimuth;
        s.scale = body.standoff_scale;
        s.standoff = s.nose();
        return s;
    case MagnetosphereModel::Shue1998:
        // Shue et al. 1998, equations 10 and 11.
        s.standoff = (10.22 + 1.29 * std::tanh(0.184 * (wind.bz_nt + 8.14))) * std::pow(pressure, -1.0 / 6.6);
        s.flaring = (0.58 - 0.007 * wind.bz_nt) * (1.0 + 0.024 * std::log(pressure));
        break;
    case MagnetosphereModel::PowerLaw:
        s.standoff = wind.standoff_radii > 0.0    ? wind.standoff_radii
                     : body.standoff_at_1npa > 0.0 ? body.standoff_at_1npa * std::pow(pressure, -1.0 / body.compressibility)
                                                   : body.standoff_radii;
        s.flaring = body.flaring + body.flaring_slope * pressure;
        break;
    case MagnetosphereModel::Heliosphere:
        s.standoff = wind.standoff_radii > 0.0 ? wind.standoff_radii : body.standoff_radii;
        s.flaring = 0.5; // 1 / cos(theta / 2) = (2 / (1 + cos theta))^(1/2)
        break;
    case MagnetosphereModel::PressureBalance:
    case MagnetosphereModel::None:
        s.standoff = std::pow(body.surface_field_nt * body.surface_field_nt * kMagneticPressureNpaPerNt2 / pressure,
                              1.0 / body.compressibility);
        s.flaring = body.flaring;
        break;
    }
    if (body.flare_end_radii > 0.0) {
        s.flare_end = body.flare_end_radii;
        s.tail_radius = flared_half_width(s.standoff, s.flaring, s.flare_end);
    }
    return s;
}

BoundarySurface bow_shock(const BodyMagnetosphere& body, const SolarWind& wind)
{
    if (body.model == MagnetosphereModel::Joy2002) {
        return joy_surface(kJoyBowShock, wind.pressure_npa);
    }
    BoundarySurface s;
    if (body.model == MagnetosphereModel::Tabulated) {
        s.shape = BoundarySurface::Shape::Tabulated;
        s.table = body.boundaries;
        s.which = BoundaryTable::Shock;
        s.blend = wind.phase;
        // The simulated shock's flaring also reflects the IMF. Keep its azimuth in NSO.
        s.azimuth.north_south = wind.azimuth.north_south;
        s.scale = body.shock_scale;
        s.standoff = s.nose();
        return s;
    }
    s.shape = BoundarySurface::Shape::Conic;
    const bool given = body.model == MagnetosphereModel::Heliosphere || body.model == MagnetosphereModel::PowerLaw;
    if (given && wind.shock_radii > 0.0) {
        s.standoff = wind.shock_radii;
    } else if (body.model == MagnetosphereModel::PowerLaw && body.shock_at_1npa > 0.0) {
        s.standoff = body.shock_at_1npa * std::pow(std::max(wind.pressure_npa, 1e-6), -1.0 / body.shock_compressibility);
    } else {
        s.standoff = body.shock_standoff_ratio * magnetopause(body, wind).standoff;
    }
    s.eccentricity = body.shock_eccentricity;
    s.focus = body.shock_focus_radii;
    return s;
}

double dipole_phase(const glm::dvec3& dipole, const glm::dvec3& spin, const glm::dvec3& upstream)
{
    // The cosine between dipole and wind is `along` + `across` cos(turn): pole-on is the end
    // of that range farther from zero.
    const double along = glm::dot(dipole, spin) * glm::dot(upstream, spin);
    const double across = std::sqrt(std::max((1.0 - std::pow(glm::dot(dipole, spin), 2.0)) * (1.0 - std::pow(glm::dot(upstream, spin), 2.0)), 0.0));
    const double cosine = (along < 0.0 ? -1.0 : 1.0) * glm::dot(dipole, upstream);
    if (across < 1e-9) {
        return std::abs(cosine) > std::sqrt(0.5) ? 1.0 : 0.0;
    }
    return std::clamp((cosine - (std::abs(along) - across)) / (2.0 * across), 0.0, 1.0);
}

glm::dmat3 solar_wind_axes(const glm::dvec3& from_sun, const glm::dvec3& velocity, double wind_speed_km_s)
{
    const glm::dvec3 outward = glm::normalize(from_sun);
    glm::dvec3 north = glm::cross(from_sun, velocity);
    if (!(glm::length(north) > 0.0)) {
        north = glm::cross(outward, std::abs(outward.z) < 0.9 ? glm::dvec3(0.0, 0.0, 1.0) : glm::dvec3(1.0, 0.0, 0.0));
    }
    north = glm::normalize(north);
    // The wind as the planet sees it.
    const glm::dvec3 flow = wind_speed_km_s > 0.0 ? wind_speed_km_s * outward - velocity : outward;
    const glm::dvec3 x = -glm::normalize(flow);
    const glm::dvec3 y = glm::normalize(glm::cross(north, x));
    return glm::dmat3(x, y, glm::cross(x, y));
}

void boundary_wireframe(const BoundarySurface& surface, double tail, int meridians,
                        std::vector<std::vector<glm::dvec4>>& lines)
{
    lines.clear();
    const double nose = surface.nose();
    if (!(nose > 0.0) || !std::isfinite(nose) || meridians < 3 || !(tail > 0.0)) {
        return;
    }
    const auto direction = [](double theta, double phi) {
        return glm::dvec3(std::cos(theta), std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi));
    };
    // Angles from the axis, in even steps of length along the dusk meridian, and a few
    // more toward the tail for the meridians that reach the cut later.
    constexpr int kSteps = 56;
    constexpr int kRingEvery = 8;
    const double step = (2.0 * nose + tail) / kSteps;
    std::vector<double> thetas{0.0};
    for (int i = 0; i < 4 * kSteps; ++i) {
        const double theta = thetas.back();
        const double r = surface.radius(direction(theta, 0.0));
        if (!std::isfinite(r) || r * std::cos(theta) < -tail) {
            break;
        }
        const double ahead = surface.radius(direction(std::min(theta + 1e-4, kPi - 1e-9), 0.0));
        const double slope = std::isfinite(ahead) ? (ahead - r) / 1e-4 : 1e12;
        const double next = theta + std::clamp(step / std::sqrt(r * r + slope * slope), 1e-5, 6.0 * kDegToRad);
        if (next >= kPi - 1e-6) {
            break;
        }
        thetas.push_back(next);
    }
    for (int i = 0; i < 6; ++i) {
        thetas.push_back(thetas.back() + 0.5 * (kPi - thetas.back()));
    }

    const int ring_points = std::max(48, 4 * meridians);
    std::vector<glm::dvec4> ring;
    const auto fade = [&](double x) { return std::clamp(-x / tail, 0.0, 1.0); };
    // The point of a ring's azimuth at one angle from the axis; false beyond the cut.
    const auto point_at = [&](double theta, double phi, glm::dvec4& out) {
        const glm::dvec3 d = direction(theta, phi);
        const double r = surface.radius(d);
        if (!std::isfinite(r) || r * d.x < -tail) {
            return false;
        }
        out = glm::dvec4(r * d, fade(r * d.x));
        return true;
    };
    // The angle at which an azimuth meets the cut, between an angle inside and one beyond.
    const auto cut_at = [&](double inside, double beyond, double phi) {
        for (int i = 0; i < 50; ++i) {
            const double mid = 0.5 * (inside + beyond);
            glm::dvec4 unused;
            (point_at(mid, phi, unused) ? inside : beyond) = mid;
        }
        glm::dvec4 p;
        point_at(inside, phi, p);
        return p;
    };
    const auto trace = [&](double phi, size_t count, std::vector<glm::dvec4>& line) {
        line.clear();
        for (size_t i = 0; i < count; ++i) {
            glm::dvec4 p;
            if (!point_at(thetas[i], phi, p)) {
                if (i > 0) {
                    line.push_back(cut_at(thetas[i - 1], thetas[i], phi));
                }
                return false;
            }
            line.push_back(p);
        }
        return true;
    };

    for (int m = 0; m < meridians; ++m) {
        std::vector<glm::dvec4> line;
        trace(kTwoPi * m / meridians, thetas.size(), line);
        if (line.size() >= 2) {
            lines.push_back(std::move(line));
        }
    }
    // Rings at even steps, as long as the whole ring is short of the cut, then the cut itself.
    std::vector<glm::dvec4> scratch;
    for (size_t i = kRingEvery / 2; i < thetas.size(); i += kRingEvery) {
        ring.clear();
        bool whole = true;
        for (int k = 0; k <= ring_points && whole; ++k) {
            glm::dvec4 p;
            whole = point_at(thetas[i], kTwoPi * k / ring_points, p);
            ring.push_back(p);
        }
        if (!whole) {
            break;
        }
        lines.push_back(ring);
    }
    ring.clear();
    for (int k = 0; k <= ring_points; ++k) {
        if (!trace(kTwoPi * (k % ring_points) / ring_points, thetas.size(), scratch) && !scratch.empty()) {
            ring.push_back(scratch.back());
        }
    }
    if (ring.size() == static_cast<size_t>(ring_points) + 1) {
        lines.push_back(ring);
    }
}

} // namespace astraxis
