#include "scene/aurora.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>

namespace astraxis {

namespace {

constexpr double kIntervalsPerDay = 8.0;

double clamp_kp(double kp)
{
    return std::clamp(kp, 0.0, 9.0);
}

// Sigernes et al. (2011), Table 1: for each boundary, the coefficients b0..b3
// (rows) of A0, A1, A2, A3, tau1, tau2, tau3 (columns) in their eq. 3,
// x = b0 + b1 log|AL| + b2 log^2|AL| + b3 log^3|AL|. (The poleward b0 of A2
// is printed "-4-44".)
using StarkovTable = double[4][7];
constexpr StarkovTable kPoleward = {
    {-0.07, -10.06, -4.44, -3.77, -6.61, 6.37, -4.48},
    {24.54, 19.83, 7.47, 7.90, 10.17, -1.10, 10.16},
    {-12.53, -9.33, -3.01, -4.73, -5.80, 0.34, -5.87},
    {2.15, 1.24, 0.25, 0.91, 1.19, -0.38, 0.98},
};
constexpr StarkovTable kEquatorward = {
    {1.61, -9.59, -12.07, -6.56, -2.22, -23.98, -20.07},
    {23.21, 17.78, 17.49, 11.44, 1.50, 42.79, 36.67},
    {-10.97, -7.20, -7.96, -6.73, -0.58, -26.96, -24.20},
    {2.03, 0.96, 1.15, 1.31, 0.08, 5.56, 5.11},
};
constexpr StarkovTable kDiffuse = {
    {3.44, -2.41, -0.74, -2.12, -1.68, 8.69, 8.61},
    {29.77, 7.89, 3.94, 3.24, -2.48, -20.73, -5.34},
    {-16.38, -4.32, -3.09, -1.67, 1.58, 13.03, -1.36},
    {3.35, 0.87, 0.72, 0.31, -0.28, -2.14, 0.76},
};

// Eq. 2: theta = A0 + A1 cos[15 (t + tau1)] + A2 cos[15 (2t + tau2)] + A3 cos[15 (3t + tau3)],
// angles in degrees (15 per hour), with the coefficients of eq. 3.
double starkov_colatitude(const StarkovTable& table, double log_al, double t)
{
    double c[7];
    for (int j = 0; j < 7; ++j) {
        c[j] = table[0][j] + log_al * (table[1][j] + log_al * (table[2][j] + log_al * table[3][j]));
    }
    const double k = 15.0 * kDegToRad;
    return c[0] + c[1] * std::cos(k * (t + c[4])) + c[2] * std::cos(k * (2.0 * t + c[5])) +
           c[3] * std::cos(k * (3.0 * t + c[6]));
}

// IGRF-14 dipole coefficients (nT): epoch, g10, g11, h11 (igrf14coeffs.txt,
// https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt; DGRF from 1945
// to 2020, IGRF otherwise), and the 2025-30 secular variation (nT/yr).
struct DipoleEpoch {
    double year, g10, g11, h11;
};
constexpr DipoleEpoch kIgrfDipole[] = {
    {1900.0, -31543.0, -2298.0, 5922.0},       {1905.0, -31464.0, -2298.0, 5909.0},
    {1910.0, -31354.0, -2297.0, 5898.0},       {1915.0, -31212.0, -2306.0, 5875.0},
    {1920.0, -31060.0, -2317.0, 5845.0},       {1925.0, -30926.0, -2318.0, 5817.0},
    {1930.0, -30805.0, -2316.0, 5808.0},       {1935.0, -30715.0, -2306.0, 5812.0},
    {1940.0, -30654.0, -2292.0, 5821.0},       {1945.0, -30594.0, -2285.0, 5810.0},
    {1950.0, -30554.0, -2250.0, 5815.0},       {1955.0, -30500.0, -2215.0, 5820.0},
    {1960.0, -30421.0, -2169.0, 5791.0},       {1965.0, -30334.0, -2119.0, 5776.0},
    {1970.0, -30220.0, -2068.0, 5737.0},       {1975.0, -30100.0, -2013.0, 5675.0},
    {1980.0, -29992.0, -1956.0, 5604.0},       {1985.0, -29873.0, -1905.0, 5500.0},
    {1990.0, -29775.0, -1848.0, 5406.0},       {1995.0, -29692.0, -1784.0, 5306.0},
    {2000.0, -29619.4, -1728.2, 5186.1},       {2005.0, -29554.63, -1669.05, 5077.99},
    {2010.0, -29496.57, -1586.42, 4944.26},    {2015.0, -29441.46, -1501.77, 4795.99},
    {2020.0, -29403.41, -1451.37, 4653.35},    {2025.0, -29350.0, -1410.3, 4545.5},
};
constexpr DipoleEpoch kIgrfSecular = {2030.0, 12.6, 10.0, -21.5};

// Linear RGB (I/F) of 1 kR of each line, white-balanced to sunlight
// (tools/geomag/aurora_colors.py).
constexpr glm::dvec3 kVioletPerKr{2.1513e-12, 0.0, 9.0337e-12};
constexpr glm::dvec3 kBluePerKr{1.3230e-10, 0.0, 6.0553e-10};
constexpr glm::dvec3 kGreenPerKr{1.7361e-09, 6.3813e-09, 0.0};

} // namespace

const glm::dvec3 kRedPerKr{4.4807e-09, 0.0, 2.9771e-10};

bool KpRecord::load(const std::filesystem::path& file, std::string* error)
{
    std::ifstream in(file, std::ios::binary);
    if (!in) {
        *error = "cannot open " + file.string();
        return false;
    }
    const std::vector<char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    uint16_t version = 0;
    uint16_t first_year = 0;
    uint32_t count = 0;
    if (bytes.size() < 12 || std::memcmp(bytes.data(), "AXKP", 4) != 0) {
        *error = file.string() + ": not a Kp record (AXKP)";
        return false;
    }
    std::memcpy(&version, bytes.data() + 4, 2);
    std::memcpy(&first_year, bytes.data() + 6, 2);
    std::memcpy(&count, bytes.data() + 8, 4);
    if (version != 1 || count == 0 || bytes.size() != 12 + static_cast<size_t>(count)) {
        *error = file.string() + ": unsupported version or wrong size";
        return false;
    }
    m_thirds.assign(bytes.begin() + 12, bytes.end());
    std::vector<uint8_t> valid;
    for (uint8_t& v : m_thirds) {
        if (v > 27) {
            v = 255; // missing: replaced by the median below
        } else {
            valid.push_back(v);
        }
    }
    if (valid.empty()) {
        *error = file.string() + ": no valid Kp";
        return false;
    }
    std::nth_element(valid.begin(), valid.begin() + static_cast<std::ptrdiff_t>(valid.size() / 2), valid.end());
    const uint8_t median = valid[valid.size() / 2];
    m_median = median / 3.0;
    m_prefix.assign(m_thirds.size() + 1, 0.0);
    for (size_t k = 0; k < m_thirds.size(); ++k) {
        if (m_thirds[k] == 255) {
            m_thirds[k] = median;
        }
        m_prefix[k + 1] = m_prefix[k] + m_thirds[k];
    }
    // January 1 of the first year, 00:00 UTC (Gregorian: days before it since
    // JD 1721425.5, 0001-01-01).
    const double y = first_year - 1.0;
    m_start_jd_utc = 1721425.5 + 365.0 * y + std::floor(y / 4.0) - std::floor(y / 100.0) + std::floor(y / 400.0);
    return true;
}

double KpRecord::interval(double t_tdb) const
{
    return (utc_jd_from_tdb(t_tdb) - m_start_jd_utc) * kIntervalsPerDay;
}

double KpRecord::start_tdb() const
{
    return tdb_from_utc_jd(m_start_jd_utc);
}

double KpRecord::end_tdb() const
{
    return tdb_from_utc_jd(m_start_jd_utc + static_cast<double>(m_thirds.size()) / kIntervalsPerDay);
}

double KpRecord::at(double t_tdb) const
{
    const double x = std::floor(interval(t_tdb));
    if (m_thirds.empty() || x < 0.0 || x >= static_cast<double>(m_thirds.size())) {
        return m_median;
    }
    return m_thirds[static_cast<size_t>(x)] / 3.0;
}

double KpRecord::cumulative(double x) const
{
    const double n = static_cast<double>(m_thirds.size());
    const double median = 3.0 * m_median;
    if (x <= 0.0) {
        return x * median;
    }
    if (x >= n) {
        return m_prefix.back() + (x - n) * median;
    }
    const double k = std::floor(x);
    const size_t i = static_cast<size_t>(k);
    return m_prefix[i] + (x - k) * m_thirds[i];
}

double KpRecord::mean(double t0_tdb, double t1_tdb) const
{
    if (m_thirds.empty()) {
        return m_median;
    }
    if (t1_tdb < t0_tdb) {
        std::swap(t0_tdb, t1_tdb);
    }
    const double x0 = interval(t0_tdb);
    const double x1 = interval(t1_tdb);
    if (x1 - x0 < 1e-9) {
        return at(t0_tdb);
    }
    return (cumulative(x1) - cumulative(x0)) / (x1 - x0) / 3.0;
}

double starkov_al_nt(double kp)
{
    kp = clamp_kp(kp);
    return 18.0 - 12.3 * kp + 27.2 * kp * kp - 2.0 * kp * kp * kp;
}

AuroralBoundaries starkov_boundaries(double kp, double mlt_hours)
{
    const double log_al = std::log10(starkov_al_nt(kp));
    return {
        .diffuse_deg = starkov_colatitude(kDiffuse, log_al, mlt_hours),
        .equatorward_deg = starkov_colatitude(kEquatorward, log_al, mlt_hours),
        .poleward_deg = starkov_colatitude(kPoleward, log_al, mlt_hours),
    };
}

double electron_hemispheric_power_gw(double kp)
{
    kp = clamp_kp(kp);
    return 2.5103 + kp * (9.9207 + kp * (-2.1825 + kp * 0.3033));
}

glm::dvec3 geomagnetic_north_pole(double t_tdb)
{
    const double year = 2000.0 + t_tdb / (kDaysPerJulianYear * kSecondsPerDay);
    constexpr size_t n = std::size(kIgrfDipole);
    DipoleEpoch d;
    if (year <= kIgrfDipole[0].year) {
        d = kIgrfDipole[0];
    } else if (year >= kIgrfDipole[n - 1].year) {
        const DipoleEpoch& last = kIgrfDipole[n - 1];
        const double dt = std::min(year, kIgrfSecular.year) - last.year;
        d = {year, last.g10 + kIgrfSecular.g10 * dt, last.g11 + kIgrfSecular.g11 * dt,
             last.h11 + kIgrfSecular.h11 * dt};
    } else {
        size_t k = 0;
        while (kIgrfDipole[k + 1].year < year) {
            ++k;
        }
        const DipoleEpoch& a = kIgrfDipole[k];
        const DipoleEpoch& b = kIgrfDipole[k + 1];
        const double f = (year - a.year) / (b.year - a.year);
        d = {year, a.g10 + f * (b.g10 - a.g10), a.g11 + f * (b.g11 - a.g11), a.h11 + f * (b.h11 - a.h11)};
    }
    // The dipole moment points along -(g11, h11, g10); with g10 < 0 that is
    // north, toward the north geomagnetic pole.
    return glm::normalize(glm::dvec3(-d.g11, -d.h11, -d.g10));
}

double magnetic_local_time(const glm::dvec3& n, const glm::dvec3& pole, const glm::dvec3& sun)
{
    glm::dvec3 noon = sun - glm::dot(sun, pole) * pole;
    noon = glm::normalize(noon);
    const glm::dvec3 dusk = glm::cross(pole, noon);
    return 12.0 + std::atan2(glm::dot(n, dusk), glm::dot(n, noon)) * (12.0 / kPi);
}

AuroraLineYields aurora_line_yields(double mean_energy_kev)
{
    // Table 2: kR per erg cm^-2 s^-1 at 1 keV and at 10 keV (fixed energy flux).
    const double f = std::clamp(std::log10(mean_energy_kev), 0.0, 1.0);
    const auto interpolate = [f](double at_1kev, double at_10kev) {
        return std::exp((1.0 - f) * std::log(at_1kev) + f * std::log(at_10kev));
    };
    return {
        .violet = interpolate(0.189, 0.066),
        .blue = interpolate(0.072, 0.025),
        .green = interpolate(0.325, 0.445),
        .red = interpolate(0.525, 0.082),
    };
}

double red_over_blue(double blue_kr)
{
    // 2 at 0.1 kR, 0.4 at 3 kR: exponent ln(0.2) / ln(30).
    const double exponent = std::log(0.2) / std::log(30.0);
    return 2.0 * std::pow(std::clamp(blue_kr, 0.1, 3.0) / 0.1, exponent);
}

AuroraEmission aurora_emission(double kp)
{
    AuroraEmission e;
    e.kp = clamp_kp(kp);
    e.hemispheric_power_gw = electron_hemispheric_power_gw(e.kp);

    // Areas (m^2) of the oval and the diffuse aurora at the flux altitude.
    const double radius_m = (6371.0 + AuroraEmission::kFluxAltitudeKm) * 1000.0;
    const double dphi = kTwoPi / AuroraEmission::kOvalSamples;
    double oval_area = 0.0;
    double diffuse_area = 0.0;
    for (int k = 0; k < AuroraEmission::kOvalSamples; ++k) {
        const AuroralBoundaries b = starkov_boundaries(e.kp, 24.0 * k / AuroraEmission::kOvalSamples);
        const double diffuse = std::max(b.diffuse_deg, b.equatorward_deg) * kDegToRad;
        const double equatorward = std::max(b.equatorward_deg, b.poleward_deg) * kDegToRad;
        const double poleward = b.poleward_deg * kDegToRad;
        e.oval[static_cast<size_t>(k)] = glm::vec3(diffuse, equatorward, poleward);
        oval_area += radius_m * radius_m * dphi * (std::cos(poleward) - std::cos(equatorward));
        diffuse_area += radius_m * radius_m * dphi * (std::cos(equatorward) - std::cos(diffuse));
    }
    // 1 W m^-2 = 1e3 erg cm^-2 s^-1.
    const double watts = e.hemispheric_power_gw * 1e9;
    e.oval_flux = 1e3 * watts / (oval_area + AuroraEmission::kDiffuseFlux * diffuse_area);

    const AuroraLineYields y = aurora_line_yields(AuroraEmission::kMeanEnergyKev);
    e.blue_kr = e.oval_flux * y.blue;
    e.low_rgb = e.oval_flux * y.violet * kVioletPerKr + e.blue_kr * (kBluePerKr + kGreenOverBlue * kGreenPerKr);
    return e;
}

} // namespace astraxis
