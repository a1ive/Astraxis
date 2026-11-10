// Shared helpers for the self-checking tests.

#include "test_util.hpp"

#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#define STBI_ONLY_PNG
#include <stb_image.h>

#include "core/math.hpp"
#include "core/time.hpp"
#include "scene/scene_loader.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

using namespace astraxis;

namespace {

int g_failures = 0;

} // namespace

double TestImage::at(double lat_deg, double east_lon_deg, int c) const
{
    const double u = (east_lon_deg + 180.0) / 360.0;
    const int x = static_cast<int>((u - std::floor(u)) * width) % width;
    const int y = std::clamp(static_cast<int>((90.0 - lat_deg) / 180.0 * height), 0, height - 1);
    return rgb[(static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)) * 3 +
               static_cast<size_t>(c)] / 255.0;
}

TestImage load_test_image(const char* asset_path)
{
    TestImage image;
    int channels = 0;
    const std::string path = (std::filesystem::path(ASTRAXIS_ASSET_DIR) / asset_path).string();
    stbi_uc* pixels = stbi_load(path.c_str(), &image.width, &image.height, &channels, 3);
    if (pixels) {
        image.rgb.assign(pixels, pixels + static_cast<size_t>(image.width) * static_cast<size_t>(image.height) * 3);
        stbi_image_free(pixels);
    }
    return image;
}

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

void check_model_sky(const Scene& scene, const char* what)
{
    const std::filesystem::path map = std::filesystem::path(ASTRAXIS_ASSET_DIR) / scene.sky.milky_way;
    check(scene.sky.milky_way.starts_with("sky/") && std::filesystem::exists(map), what);
    check(scene.sky.stars.size() > 500, what, static_cast<double>(scene.sky.stars.size()));
    check(scene.sky.milky_way_brightness > 0.0 && scene.sky.milky_way_brightness < 10.0, what,
          scene.sky.milky_way_brightness);
    const int sun = scene.find("Sun");
    if (sun >= 0) {
        const glm::dvec3 sun_pc = scene.icrf_state_at(sun, 0.0).position / kParsecKm;
        const double miss = glm::length(scene.sky.viewer_pc + sun_pc) / glm::length(scene.sky.viewer_pc);
        check(miss < 1e-6, what, miss);
    }
    double brightest = 99.0;
    for (const CatalogStar& s : scene.sky.stars) {
        brightest = std::min(brightest, s.vmag);
    }
    std::printf("info: %s: %zu stars, brightest V = %.2f, Milky Way brightness %.3f\n", what, scene.sky.stars.size(),
                brightest, scene.sky.milky_way_brightness);
}
