#include "core/time.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>

namespace astraxis {

namespace {

// IERS Leap_Second.dat (https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat),
// "updated through IERS Bulletin 72 issued in July 2026, expires 28 June 2027".
// Each entry: MJD (UTC) from which TAI - UTC takes the given value.
struct LeapEntry {
    double mjd;
    double tai_minus_utc;
};

constexpr LeapEntry kLeapSeconds[] = {
    {41317.0, 10}, {41499.0, 11}, {41683.0, 12}, {42048.0, 13}, {42413.0, 14}, {42778.0, 15},
    {43144.0, 16}, {43509.0, 17}, {43874.0, 18}, {44239.0, 19}, {44786.0, 20}, {45151.0, 21},
    {45516.0, 22}, {46247.0, 23}, {47161.0, 24}, {47892.0, 25}, {48257.0, 26}, {48804.0, 27},
    {49169.0, 28}, {49534.0, 29}, {50083.0, 30}, {50630.0, 31}, {51179.0, 32}, {53736.0, 33},
    {54832.0, 34}, {56109.0, 35}, {57204.0, 36}, {57754.0, 37},
};

constexpr double kMjdOffset = 2400000.5;

} // namespace

double tai_minus_utc(double jd_utc)
{
    const double mjd = jd_utc - kMjdOffset;
    double value = kLeapSeconds[0].tai_minus_utc;
    for (const LeapEntry& e : kLeapSeconds) {
        if (mjd >= e.mjd) {
            value = e.tai_minus_utc;
        } else {
            break;
        }
    }
    return value;
}

double tdb_from_utc_jd(double jd_utc)
{
    return (jd_utc - kJ2000Jd) * kSecondsPerDay + tai_minus_utc(jd_utc) + kTtMinusTai;
}

double utc_jd_from_tdb(double t_tdb)
{
    // TAI - UTC is piecewise constant: two fixed-point steps settle it.
    double jd_utc = kJ2000Jd + (t_tdb - 69.184) / kSecondsPerDay;
    for (int i = 0; i < 2; ++i) {
        jd_utc = kJ2000Jd + (t_tdb - tai_minus_utc(jd_utc) - kTtMinusTai) / kSecondsPerDay;
    }
    return jd_utc;
}

CalendarDateTime calendar_from_jd(double jd)
{
    const double shifted = jd + 0.5;
    const int64_t z = static_cast<int64_t>(std::floor(shifted));
    const double f = shifted - static_cast<double>(z);

    int64_t a = z;
    if (z >= 2299161) {
        const int64_t alpha = static_cast<int64_t>(std::floor((static_cast<double>(z) - 1867216.25) / 36524.25));
        a = z + 1 + alpha - static_cast<int64_t>(std::floor(static_cast<double>(alpha) / 4.0));
    }
    const int64_t b = a + 1524;
    const int64_t c = static_cast<int64_t>(std::floor((static_cast<double>(b) - 122.1) / 365.25));
    const int64_t d = static_cast<int64_t>(std::floor(365.25 * static_cast<double>(c)));
    const int64_t e = static_cast<int64_t>(std::floor(static_cast<double>(b - d) / 30.6001));

    CalendarDateTime out;
    out.day = static_cast<int>(b - d - static_cast<int64_t>(std::floor(30.6001 * static_cast<double>(e))));
    out.month = static_cast<int>(e < 14 ? e - 1 : e - 13);
    out.year = static_cast<int>(out.month > 2 ? c - 4716 : c - 4715);

    int64_t seconds_of_day = static_cast<int64_t>(std::floor(f * kSecondsPerDay + 1e-6));
    if (seconds_of_day > 86399) {
        seconds_of_day = 86399;
    }
    out.hour = static_cast<int>(seconds_of_day / 3600);
    out.minute = static_cast<int>((seconds_of_day / 60) % 60);
    out.second = static_cast<int>(seconds_of_day % 60);
    return out;
}

double jd_from_calendar(const CalendarDateTime& c)
{
    int y = c.year;
    int m = c.month;
    if (m <= 2) {
        y -= 1;
        m += 12;
    }
    // Gregorian calendar from 1582-10-15 on.
    const bool gregorian = c.year > 1582 || (c.year == 1582 && (c.month > 10 || (c.month == 10 && c.day >= 15)));
    double b = 0.0;
    if (gregorian) {
        const double a = std::floor(y / 100.0);
        b = 2.0 - a + std::floor(a / 4.0);
    }
    const double day_fraction = (c.hour * 3600.0 + c.minute * 60.0 + c.second) / kSecondsPerDay;
    return std::floor(365.25 * (y + 4716)) + std::floor(30.6001 * (m + 1)) + c.day + day_fraction + b - 1524.5;
}

bool parse_utc(const char* text, double* t_tdb)
{
    CalendarDateTime c;
    const int fields = std::sscanf(text, "%d-%d-%d %d:%d:%d", &c.year, &c.month, &c.day, &c.hour, &c.minute,
                                   &c.second);
    if (fields != 3 && fields != 5 && fields != 6) {
        return false;
    }
    if (c.month < 1 || c.month > 12 || c.day < 1 || c.day > 31 || c.hour < 0 || c.hour > 23 || c.minute < 0 ||
        c.minute > 59 || c.second < 0 || c.second > 60) {
        return false;
    }
    *t_tdb = tdb_from_utc_jd(jd_from_calendar(c));
    return true;
}

std::string format_utc(double t_tdb)
{
    const CalendarDateTime c = calendar_from_jd(utc_jd_from_tdb(t_tdb));
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d UTC", c.year, c.month, c.day, c.hour, c.minute,
                  c.second);
    return buf;
}

} // namespace astraxis
