// Shared helpers for the self-checking tests.

#include "test_util.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "scene/scene_loader.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace astraxis;

namespace {

int g_failures = 0;

} // namespace

void check(bool ok, const char* what, double value)
{
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s (value = %.9g)\n", what, value);
    }
}

int failure_count()
{
    return g_failures;
}

Scene load_scene_or_die(const char* file)
{
    Scene scene;
    std::string error;
    if (!load_scene_file(std::string(ASTRAXIS_ASSET_DIR "/scenes/") + file, scene, &error)) {
        std::printf("FAIL: cannot load %s: %s\n", file, error.c_str());
        std::exit(1);
    }
    return scene;
}

Scene load_jupiter()
{
    Scene scene;
    std::string error;
    if (!load_scene_file(ASTRAXIS_ASSET_DIR "/scenes/jupiter.toml", scene, &error)) {
        std::printf("FAIL: cannot load jupiter.toml: %s\n", error.c_str());
        std::exit(1);
    }
    return scene;
}

double tdb_from_jd_tdb(double jd)
{
    return (jd - kJ2000Jd) * kSecondsPerDay;
}

glm::dvec3 unit_toward(double ra_deg, double dec_deg)
{
    const double ra = ra_deg * kDegToRad;
    const double dec = dec_deg * kDegToRad;
    return {std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec)};
}
