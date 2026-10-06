// Tests for core/: calendar and time scales, star colors, the Galactic frame.

#include "test_util.hpp"

#include "core/color.hpp"
#include "core/math.hpp"
#include "core/time.hpp"
#include "scene/star_catalog.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace astraxis;

namespace {

void test_calendar()
{
    CalendarDateTime c = calendar_from_jd(kJ2000Jd);
    check(c.year == 2000 && c.month == 1 && c.day == 1 && c.hour == 12 && c.minute == 0 && c.second == 0,
          "J2000 = 2000-01-01 12:00:00");

    c = calendar_from_jd(2299160.5);
    check(c.year == 1582 && c.month == 10 && c.day == 15, "Gregorian reform start 1582-10-15");

    c = calendar_from_jd(kUnixEpochJd);
    check(c.year == 1970 && c.month == 1 && c.day == 1 && c.hour == 0, "Unix epoch 1970-01-01");

    // 2026-10-02 00:00:00 UTC = Unix 1790899200.
    const double t = tdb_from_unix_utc(1790899200.0);
    check(format_utc(t) == "2026-10-02 00:00:00 UTC", "Unix -> UTC round trip");

    double parsed = 0.0;
    check(parse_utc("2026-10-02 00:00:00", &parsed) && std::abs(parsed - t) < 1e-3, "parse_utc", parsed - t);
    check(parse_utc("1979-03-05", &parsed) && format_utc(parsed) == "1979-03-05 00:00:00 UTC", "parse_utc date only");
    check(parse_utc("1066-10-14 09:30", &parsed) && format_utc(parsed) == "1066-10-14 09:30:00 UTC",
          "parse_utc Julian calendar");
    check(!parse_utc("2026-13-01", &parsed), "parse_utc rejects bad month");

    // Leap seconds: TDB - UTC = TAI-UTC + 32.184 s.
    check(std::abs(tdb_from_utc_jd(2451545.0) - (32.0 + 32.184)) < 1e-6, "TDB-UTC at J2000 (TAI-UTC 32 s)");
    parse_utc("1979-03-05 12:00:00", &parsed);
    const double jd_1979 = jd_from_calendar({1979, 3, 5, 12, 0, 0});
    check(std::abs(parsed - ((jd_1979 - kJ2000Jd) * kSecondsPerDay + 18.0 + 32.184)) < 1e-6,
          "TDB-UTC in 1979 (TAI-UTC 18 s)");
    check(format_utc(parsed) == "1979-03-05 12:00:00 UTC", "1979 UTC round trip");
    // Across the 2016-12-31 leap second the UTC clock advances 1 s less than TDB.
    double a = 0.0;
    double b = 0.0;
    parse_utc("2016-12-31 23:59:00", &a);
    parse_utc("2017-01-01 00:01:00", &b);
    check(std::abs((b - a) - 121.0) < 1e-4, "leap second at 2017-01-01", b - a);
    check(!parse_utc("hello", &parsed), "parse_utc rejects garbage");
}

void test_star_catalog_and_color()
{
    std::vector<CatalogStar> stars;
    std::string error;
    check(load_star_catalog(ASTRAXIS_ASSET_DIR "/stars/bsc5.csv", stars, &error), "BSC5 loads");
    check(stars.size() > 9000, "BSC5 star count", static_cast<double>(stars.size()));
    double brightest = 99.0;
    for (const CatalogStar& s : stars) {
        brightest = std::min(brightest, s.vmag);
    }
    check(std::abs(brightest + 1.46) < 1e-9, "brightest star is Sirius (-1.46)", brightest);

    // The Sun (B-V = 0.65) is ~5800 K; blackbody colors go from red to blue.
    const double t_sun = temperature_from_bv(0.65);
    check(t_sun > 5600.0 && t_sun < 6000.0, "Sun temperature from B-V", t_sun);
    const glm::dvec3 cool = blackbody_linear_srgb(3000.0);
    const glm::dvec3 hot = blackbody_linear_srgb(15000.0);
    check(cool.r > cool.b && hot.b > hot.r, "blackbody color ordering");
    const glm::dvec3 d65ish = blackbody_linear_srgb(6500.0);
    check(d65ish.r > 0.9 && d65ish.g > 0.9 && d65ish.b > 0.85, "6500 K is near white", d65ish.b);
}

void test_galactic_frame()
{
    // The Galactic centre direction (l = b = 0) in ICRS: RA 266.405, Dec -28.936 deg
    // (Hipparcos definition, ESA SP-1200 vol. 1, sec. 1.5.3).
    const glm::dvec3 gc = galactic_to_icrf() * glm::dvec3(1.0, 0.0, 0.0);
    double ra = std::atan2(gc.y, gc.x) * kRadToDeg;
    ra = ra < 0.0 ? ra + 360.0 : ra;
    const double dec = std::asin(gc.z) * kRadToDeg;
    check(std::abs(ra - 266.405) < 0.01 && std::abs(dec + 28.936) < 0.01, "Galactic centre direction", ra);
}

} // namespace

void run_core_tests()
{
    test_calendar();
    test_star_catalog_and_color();
    test_galactic_frame();
}
