#include "scene/hr_diagram.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace astraxis {

namespace {

// [TO96] Tables 1 and 2 at Z = 0.02 (each coefficient's constant term).
constexpr double kLzAlpha = 0.39704170;
constexpr double kLzBeta = 8.52762600;
constexpr double kLzGamma = 0.00025546;
constexpr double kLzDelta = 5.43288900;
constexpr double kLzEpsilon = 5.56357900;
constexpr double kLzZeta = 0.78866060;
constexpr double kLzEta = 0.00586685;
constexpr double kRzTheta = 1.71535900;
constexpr double kRzIota = 6.59778800;
constexpr double kRzKappa = 10.08855000;
constexpr double kRzLambda = 1.01249500;
constexpr double kRzMu = 0.07490166;
constexpr double kRzNu = 0.01077422;
constexpr double kRzXi = 3.08223400;
constexpr double kRzOmicron = 17.84778000;
constexpr double kRzPi = 0.00022582;

// [HPT00] Appendix at Z = 0.02: zeta = log(Z / 0.02) = 0, so a_n is its alpha
// column, with the Appendix's limits applied (those that hold at Z = 0.02).
constexpr double kSigma = -1.6989700043360187; // log10(0.02)
constexpr double kZ = 0.02;
constexpr double kHydrogen = 0.76 - 3.0 * kZ; // X, as in SSE
constexpr double kMassHook = 1.0185;          // M_hook at zeta = 0

// t_BGB and t_hook.
constexpr double a1 = 1.593890e3, a2 = 2.706708e3, a3 = 1.466143e2, a4 = 4.141960e-2, a5 = 3.426349e-1;
constexpr double a6 = 1.949814e1, a7 = 4.903830, a8 = 5.212154e-2, a9 = 1.312179, a10 = 8.073972e-1;
// L_TMS (a11 = a11' a14, a12 = a12' a14).
constexpr double a14 = 3.858911e3;
constexpr double a11 = 1.031538 * a14, a12 = 1.043715 * a14, a13 = 7.859573e2, a15 = 2.888720e2, a16 = 7.196580;
// R_TMS (a18 = a18' a20, a19 = a19' a20); c1 is Z independent.
constexpr double a20 = 2.652091e1;
constexpr double a18 = 2.187715e-1 * a20, a19 = 1.466440 * a20, a21 = 1.472103, a22 = 3.071048;
constexpr double a23 = 2.617890, a24 = 1.075567e-2, a25 = 1.476246, a26 = 5.502535;
constexpr double c1 = -8.672073e-2;
// alpha_L.
constexpr double a45 = 2.321400e-1, a46 = 1.163659e-2, a47 = 1.048020e-2, a48 = 1.555590;
// beta_L.
constexpr double a54 = 3.855707e-1, a55 = 3.579064e-1, a56 = 9.587587e-1;
// alpha_R.
constexpr double a58 = 4.907546e-1, a59 = 4.537070, a60 = 1.796220, a61 = 2.256216;
constexpr double a65 = 1.564231e-3, a67 = 5.210157;
// beta_R.
constexpr double a69 = 1.071489, a70 = 7.108492e-1, a71 = 3.478514, a73 = 3.969331e-3;
// gamma.
constexpr double a76c = 1.192334e-2, a77c = -1.668868e-1, a78c = 7.615495e-1, a79c = 9.409838;

double pow10(double x)
{
    return std::pow(10.0, x);
}

struct Coefficients {
    double a17, a49, a50, a51, a52, a53, a57;
    double a62, a63, a64, a66, a68, alpha_r_at_a66;
    double a72, a74, a75, a76, a77, a78, a79, a80, a81;
};

Coefficients make_coefficients()
{
    Coefficients c{};
    c.a17 = pow10(std::max(0.097 - 0.1072 * (kSigma + 3.0),
                           std::max(0.097, std::min(0.1461, 0.1461 + 0.1237 * (kSigma + 2.0)))));
    c.a49 = std::max(9.77e-2, 0.145);
    c.a50 = std::min(0.24, 0.306);
    c.a51 = std::min(0.33, 0.3625);
    c.a52 = std::min(std::max(1.1064, 0.9), 1.0); // Z > 0.01
    c.a53 = std::min(std::max(1.19, 1.0), 1.1);   // Z > 0.01
    c.a57 = std::max(0.6355, std::max(1.25, std::min(1.4, 1.5135)));
    c.a62 = std::max(0.065, 8.43e-2);
    c.a63 = 7.36e-2;
    c.a64 = std::max(0.091, std::min(0.121, 0.136));
    c.a66 = std::max(1.477, std::min(1.6, -0.308));
    c.a66 = std::max(0.8, std::min(0.8, c.a66));
    c.a68 = std::max(0.9, std::min(1.116, 1.0));
    c.alpha_r_at_a66 = a58 * std::pow(c.a66, a60) / (a59 + std::pow(c.a66, a61));
    if (c.a68 > c.a66) {
        c.a64 = c.alpha_r_at_a66;
        c.a68 = c.a66;
    }
    c.a72 = std::max(9.132108e-1, 0.95); // Z > 0.01
    c.a74 = std::max(1.4, std::min(1.6, 1.6));
    c.a75 = std::max(std::max(1.0, std::min(8.109e-1, 1.27)), 0.6355);
    c.a76 = std::max(a76c, -0.1015564);
    c.a77 = std::max(-0.3868776, std::min(0.0, a77c));
    c.a78 = std::max(0.0, std::min(a78c, 7.454));
    c.a79 = std::min(a79c, std::max(2.0, -13.3));
    c.a80 = std::max(0.0585542, -0.2711);
    c.a81 = std::min(1.5, std::max(0.4, 2.493));
    return c;
}

const Coefficients& coefficients()
{
    static const Coefficients c = make_coefficients();
    return c;
}

// Times in Myr, as in [HPT00].
double t_bgb_myr(double m)
{
    return (a1 + a2 * std::pow(m, 4.0) + a3 * std::pow(m, 5.5) + std::pow(m, 7.0)) /
           (a4 * m * m + a5 * std::pow(m, 7.0));
}

double t_ms_myr(double m)
{
    const double mu = std::max(0.5, 1.0 - 0.01 * std::max(a6 / std::pow(m, a7), a8 + a9 / std::pow(m, a10)));
    const double x = 0.95; // max(0.95, min(0.95 - 0.03 (zeta + 0.30103), 0.99)) at zeta = 0
    return std::max(mu, x) * t_bgb_myr(m);
}

double l_tms(double m)
{
    return (a11 * std::pow(m, 3.0) + a12 * std::pow(m, 4.0) + a13 * std::pow(m, a16 + 1.8)) /
           (a14 + a15 * std::pow(m, 5.0) + std::pow(m, a16));
}

double r_tms(double m)
{
    const Coefficients& c = coefficients();
    const auto low = [](double mm) { return (a18 + a19 * std::pow(mm, a21)) / (a20 + std::pow(mm, a22)); };
    const auto high = [](double mm) {
        return (c1 * mm * mm * mm + a23 * std::pow(mm, a26) + a24 * std::pow(mm, a26 + 1.5)) /
               (a25 + std::pow(mm, 5.0));
    };
    const double m_star = c.a17 + 0.1;
    if (m <= c.a17) {
        return std::max(1.5 * zams_star(m).radius, low(m)); // the guard for extrapolation to low masses
    }
    if (m >= m_star) {
        return high(m);
    }
    return low(c.a17) + (high(m_star) - low(c.a17)) / 0.1 * (m - c.a17);
}

double alpha_l(double m)
{
    const Coefficients& c = coefficients();
    const auto formula = [](double mm) {
        return (a45 + a46 * std::pow(mm, a48)) / (std::pow(mm, 0.4) + a47 * std::pow(mm, 1.9));
    };
    if (m >= 2.0) {
        return formula(m);
    }
    if (m <= 0.5) {
        return c.a49;
    }
    if (m <= 0.7) {
        return c.a49 + (0.3 - c.a49) / 0.2 * (m - 0.5);
    }
    if (m <= c.a52) {
        return 0.3 + (c.a50 - 0.3) / (c.a52 - 0.7) * (m - 0.7);
    }
    if (m <= c.a53) {
        return c.a50 + (c.a51 - c.a50) / (c.a53 - c.a52) * (m - c.a52);
    }
    return c.a51 + (formula(2.0) - c.a51) / (2.0 - c.a53) * (m - c.a53);
}

double beta_l(double m)
{
    const Coefficients& c = coefficients();
    double beta = std::max(0.0, a54 - a55 * std::pow(m, a56));
    if (m > c.a57 && beta > 0.0) {
        const double b = a54 - a55 * std::pow(c.a57, a56);
        beta = std::max(0.0, b - 10.0 * b * (m - c.a57));
    }
    return beta;
}

double alpha_r(double m)
{
    const Coefficients& c = coefficients();
    if (m <= 0.5) {
        return c.a62;
    }
    if (m <= 0.65) {
        return c.a62 + (c.a63 - c.a62) / 0.15 * (m - 0.5);
    }
    if (m <= c.a68) {
        return c.a63 + (c.a64 - c.a63) / (c.a68 - 0.65) * (m - 0.65);
    }
    if (m <= c.a66) {
        return c.a64 + (c.alpha_r_at_a66 - c.a64) / (c.a66 - c.a68) * (m - c.a68);
    }
    if (m <= a67) {
        return a58 * std::pow(m, a60) / (a59 + std::pow(m, a61));
    }
    return a58 * std::pow(a67, a60) / (a59 + std::pow(a67, a61)) + a65 * (m - a67);
}

double beta_r(double m)
{
    const Coefficients& c = coefficients();
    const auto formula = [](double mm) { return a69 * std::pow(mm, 3.5) / (a70 + std::pow(mm, a71)); };
    double beta;
    if (m <= 1.0) {
        beta = 1.06;
    } else if (m <= c.a74) {
        beta = 1.06 + (c.a72 - 1.06) / (c.a74 - 1.0) * (m - 1.0);
    } else if (m <= 2.0) {
        beta = c.a72 + (formula(2.0) - c.a72) / (2.0 - c.a74) * (m - c.a74);
    } else if (m <= 16.0) {
        beta = formula(m);
    } else {
        beta = formula(16.0) + a73 * (m - 16.0);
    }
    return beta - 1.0;
}

double gamma_r(double m)
{
    const Coefficients& c = coefficients();
    if (m > c.a75 + 0.1) {
        return 0.0;
    }
    double b = std::max(0.0, c.a76 + c.a77 * std::pow(1.0 - c.a78, c.a79));
    double gamma;
    if (m <= 1.0) {
        gamma = c.a76 + c.a77 * std::pow(std::fabs(m - c.a78), c.a79);
    } else if (m <= c.a75) {
        gamma = b + (c.a80 - b) * std::pow((m - 1.0) / (c.a75 - 1.0), c.a81);
    } else {
        if (c.a75 > 1.0) {
            b = c.a80;
        }
        gamma = b - 10.0 * b * (m - c.a75);
    }
    return std::max(gamma, 0.0);
}

double segment_distance(HrPoint p, HrPoint a, HrPoint b)
{
    const double dx = b.log_t - a.log_t;
    const double dy = b.log_l - a.log_l;
    const double len2 = dx * dx + dy * dy;
    double s = len2 > 0.0 ? ((p.log_t - a.log_t) * dx + (p.log_l - a.log_l) * dy) / len2 : 0.0;
    s = std::clamp(s, 0.0, 1.0);
    return std::hypot(p.log_t - (a.log_t + s * dx), p.log_l - (a.log_l + s * dy));
}

HrPoint hr_point(const StellarModelPoint& s)
{
    return {std::log10(effective_temperature_k(s.luminosity, s.radius)), std::log10(s.luminosity)};
}

} // namespace

StellarModelPoint zams_star(double m)
{
    const double sq = std::sqrt(m);
    const double l = (kLzAlpha * std::pow(m, 5.0) * sq + kLzBeta * std::pow(m, 11.0)) /
                     (kLzGamma + std::pow(m, 3.0) + kLzDelta * std::pow(m, 5.0) + kLzEpsilon * std::pow(m, 7.0) +
                      kLzZeta * std::pow(m, 8.0) + kLzEta * std::pow(m, 9.0) * sq);
    const double r = ((kRzTheta * m * m + kRzIota * std::pow(m, 6.0)) * sq + kRzKappa * std::pow(m, 11.0) +
                      (kRzLambda + kRzMu * sq) * std::pow(m, 19.0)) /
                     (kRzNu + kRzXi * m * m + (kRzOmicron * std::pow(m, 8.0) + std::pow(m, 18.0) + kRzPi * std::pow(m, 19.0)) * sq);
    return {l, r};
}

double main_sequence_lifetime_yr(double m)
{
    return t_ms_myr(m) * 1e6;
}

StellarModelPoint main_sequence_star(double m, double age_yr)
{
    const StellarModelPoint zams = zams_star(m);
    const double tau = std::clamp(age_yr / main_sequence_lifetime_yr(m), 0.0, 1.0);

    constexpr double kEta = 10.0; // the exponent for Z > 0.0009
    const double al = alpha_l(m);
    const double bl = beta_l(m);
    const double lx = std::log10(l_tms(m) / zams.luminosity);
    const double log_l = al * tau + bl * std::pow(tau, kEta) + (lx - al - bl) * tau * tau;

    const double ar = alpha_r(m);
    const double br = beta_r(m);
    const double g = gamma_r(m);
    const double rx = std::log10(r_tms(m) / zams.radius);
    const double log_r =
        ar * tau + br * std::pow(tau, 10.0) + g * std::pow(tau, 40.0) + (rx - ar - br - g) * tau * tau * tau;

    double r = zams.radius * pow10(log_r);
    if (m < kMassHook) {
        // Low-mass stars are partly degenerate.
        r = std::max(r, 0.0258 * std::pow(1.0 + kHydrogen, 5.0 / 3.0) / std::cbrt(m));
    }
    return {zams.luminosity * pow10(log_l), r};
}

double white_dwarf_radius(double m)
{
    constexpr double kChandrasekhar = 1.44;
    const double x = std::pow(kChandrasekhar / m, 2.0 / 3.0) - std::pow(m / kChandrasekhar, 2.0 / 3.0);
    return 0.0115 * std::sqrt(std::max(x, 0.0));
}

double effective_temperature_k(double l, double r)
{
    return 5772.0 * std::pow(l / (r * r), 0.25);
}

bool HrBand::contains(HrPoint p) const
{
    bool inside = false;
    for (size_t i = 0, j = outline.size() - 1; i < outline.size(); j = i++) {
        const HrPoint& a = outline[i];
        const HrPoint& b = outline[j];
        if ((a.log_l > p.log_l) != (b.log_l > p.log_l)) {
            const double t = a.log_t + (p.log_l - a.log_l) * (b.log_t - a.log_t) / (b.log_l - a.log_l);
            if (p.log_t < t) {
                inside = !inside;
            }
        }
    }
    return inside;
}

double HrBand::distance(HrPoint p) const
{
    if (outline.size() < 3 || contains(p)) {
        return 0.0;
    }
    double d = std::numeric_limits<double>::infinity();
    for (size_t i = 0, j = outline.size() - 1; i < outline.size(); j = i++) {
        d = std::min(d, segment_distance(p, outline[j], outline[i]));
    }
    return d;
}

HrBand main_sequence_band(double age_yr)
{
    constexpr int kSamples = 160;
    constexpr double kLowMass = 0.1;
    constexpr double kHighMass = 100.0;
    HrBand band;
    std::vector<HrPoint> upper;
    band.outline.reserve(2 * (kSamples + 1));
    upper.reserve(kSamples + 1);
    for (int i = 0; i <= kSamples; ++i) {
        const double m = kLowMass * std::pow(kHighMass / kLowMass, static_cast<double>(i) / kSamples);
        band.outline.push_back(hr_point(zams_star(m)));
        upper.push_back(hr_point(main_sequence_star(m, age_yr)));
    }
    band.outline.insert(band.outline.end(), upper.rbegin(), upper.rend());
    return band;
}

HrBand white_dwarf_band(double hot_k, double cool_k)
{
    const auto at = [](double t, double r) {
        const double l = r * r * std::pow(t / 5772.0, 4.0);
        return HrPoint{std::log10(t), std::log10(l)};
    };
    const double light = white_dwarf_radius(0.2);
    const double heavy = white_dwarf_radius(1.2);
    // Lines of constant radius are straight on the diagram; the lower one first.
    return HrBand{{at(hot_k, heavy), at(cool_k, heavy), at(cool_k, light), at(hot_k, light)}};
}

} // namespace astraxis
