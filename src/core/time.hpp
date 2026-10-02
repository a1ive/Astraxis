#pragma once

#include <string>

namespace astraxis {

// Simulation time is TDB seconds since J2000.0 (JD 2451545.0 TDB).

inline constexpr double kSecondsPerDay = 86400.0;
inline constexpr double kDaysPerJulianYear = 365.25;
inline constexpr double kDaysPerJulianCentury = 36525.0;
inline constexpr double kJ2000Jd = 2451545.0;
inline constexpr double kUnixEpochJd = 2440587.5;

// TT - TAI (exact, by definition). TDB - TT (periodic, < 2 ms) is ignored.
inline constexpr double kTtMinusTai = 32.184;

// TAI - UTC in seconds at a UTC Julian date, from the IERS leap-second table
// (Leap_Second.dat, IERS Bulletin C). Before 1972-01-01 (when UTC used rubber
// seconds) the 1972 value of 10 s is used as an approximation.
double tai_minus_utc(double jd_utc);

// Simulation time from a UTC Julian date, and back. Leap seconds are handled;
// times inside a leap second (23:59:60) map to 23:59:59.
double tdb_from_utc_jd(double jd_utc);
double utc_jd_from_tdb(double t_tdb);

inline double jd_from_tdb(double t_tdb)
{
    return kJ2000Jd + t_tdb / kSecondsPerDay;
}

// Julian epoch (decimal year, e.g. 2018.379 = J2018.379) <-> simulation time.
inline double tdb_from_julian_year(double year)
{
    return (year - 2000.0) * kDaysPerJulianYear * kSecondsPerDay;
}

inline double julian_centuries_from_tdb(double t_tdb)
{
    return t_tdb / (kSecondsPerDay * kDaysPerJulianCentury);
}

// Unix time counts UTC days of exactly 86400 s (leap seconds are not counted).
inline double tdb_from_unix_utc(double unix_seconds)
{
    return tdb_from_utc_jd(unix_seconds / kSecondsPerDay + kUnixEpochJd);
}

struct CalendarDateTime {
    int year = 0;
    int month = 0;
    int day = 0;
    int hour = 0;
    int minute = 0;
    int second = 0;
};

// Gregorian calendar date (Julian calendar before 1582-10-15), truncated to
// whole seconds. Algorithm: Meeus, "Astronomical Algorithms", 2nd ed., ch. 7.
CalendarDateTime calendar_from_jd(double jd);

// Inverse of calendar_from_jd (Meeus ch. 7).
double jd_from_calendar(const CalendarDateTime& c);

// "YYYY-MM-DD HH:MM:SS UTC" for a simulation time.
std::string format_utc(double t_tdb);

// Parses "YYYY-MM-DD", "YYYY-MM-DD HH:MM" or "YYYY-MM-DD HH:MM:SS" (UTC).
// Returns false on malformed input.
bool parse_utc(const char* text, double* t_tdb);

} // namespace astraxis
