#include "scene/flares.hpp"

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr double kSlotS = 6.0 * 3600.0;
constexpr double kCutFwhm = 1.5; // the envelope is ~0.002 there
constexpr double kMinDurationS = 1800.0;
constexpr double kMaxDurationS = 5400.0;
constexpr double kMinRadiusM = 7.0;
constexpr double kMaxRadiusM = 11.0;

uint64_t splitmix64(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Uniform in [0, 1) from the slot, the seed and a channel.
double unit(int64_t slot, uint64_t seed, uint64_t channel)
{
    const uint64_t h = splitmix64(static_cast<uint64_t>(slot) ^ splitmix64(seed * 8 + channel));
    return static_cast<double>(h >> 11) * (1.0 / 9007199254740992.0);
}

double envelope(const FlareRecord& f, double t)
{
    const double x = (t - f.t_peak_tdb) / f.duration_s;
    return std::abs(x) > kCutFwhm ? 0.0 : std::exp(-4.0 * std::log(2.0) * x * x);
}

double half_window(const FlareRecord& f)
{
    return kCutFwhm * f.duration_s;
}

std::optional<FlareRecord> random_flare(const FlareSettings& s, int64_t slot)
{
    if (unit(slot, s.seed, 0) >= s.random_per_day * kSlotS / 86400.0) {
        return std::nullopt;
    }
    FlareRecord f;
    f.duration_s = kMinDurationS + (kMaxDurationS - kMinDurationS) * unit(slot, s.seed, 1);
    f.radius_m = kMinRadiusM + (kMaxRadiusM - kMinRadiusM) * unit(slot, s.seed, 2);
    f.phase_deg = 360.0 * unit(slot, s.seed, 3);
    // Wholly inside its slot, so random flares never overlap.
    const double margin = half_window(f);
    f.t_peak_tdb = static_cast<double>(slot) * kSlotS + margin + (kSlotS - 2.0 * margin) * unit(slot, s.seed, 4);
    // Recorded flares take precedence: a random one overlapping any is dropped whole.
    for (const FlareRecord& r : s.recorded) {
        if (std::abs(r.t_peak_tdb - f.t_peak_tdb) < half_window(r) + margin) {
            return std::nullopt;
        }
    }
    return f;
}

} // namespace

std::optional<ActiveFlare> active_flare(const FlareSettings& s, double t_tdb)
{
    std::optional<ActiveFlare> best;
    for (const FlareRecord& r : s.recorded) {
        const double e = envelope(r, t_tdb);
        if (e > 0.0 && (!best || e > best->envelope)) {
            best = ActiveFlare{r, e, true};
        }
    }
    if (best || s.random_per_day <= 0.0) {
        return best;
    }
    if (const std::optional<FlareRecord> f = random_flare(s, static_cast<int64_t>(std::floor(t_tdb / kSlotS)))) {
        const double e = envelope(*f, t_tdb);
        if (e > 0.0) {
            return ActiveFlare{*f, e, false};
        }
    }
    return std::nullopt;
}

double flare_warp_fade(double warp)
{
    const double t = std::clamp((std::abs(warp) - 300.0) / 900.0, 0.0, 1.0);
    return 1.0 - t * t * (3.0 - 2.0 * t);
}

} // namespace astraxis
