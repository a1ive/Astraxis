// Shared helpers for the self-checking tests (no framework). Each *_tests.cpp
// runs its domain's tests from run_<domain>_tests(); test_main.cpp calls them all.

#pragma once

#include "scene/scene.hpp"

#include <glm/vec3.hpp>

#include <vector>

// Records and prints a failure unless ok.
void check(bool ok, const char* what, double value = 0.0);
int failure_count();

// Loads assets/scenes/<file>; exits on failure.
astraxis::Scene load_scene_or_die(const char* file);
astraxis::Scene load_jupiter();

double tdb_from_jd_tdb(double jd);

// Unit vector in ICRF toward the given RA/Dec.
glm::dvec3 unit_toward(double ra_deg, double dec_deg);

// A scene far from the Sun shows a model sky (tools/sky/make_galaxy_sky.py):
// its map and stars exist, the map's brightness comes from the stars file, and
// the viewer sits opposite the scene's "Sun" body, if it has one.
void check_model_sky(const astraxis::Scene& scene, const char* what);

// An image of the assets (JPEG or PNG) as 8-bit RGB, row after row from the top; empty if
// it cannot be read.
struct TestImage {
    int width = 0;
    int height = 0;
    std::vector<unsigned char> rgb;
    // Channel c (0..2) at east longitude and latitude on a 2:1 map whose left edge is 180 W,
    // as 0..1.
    double at(double lat_deg, double east_lon_deg, int c) const;
};
TestImage load_test_image(const char* asset_path);

void run_settings_tests();
void run_core_tests();
void run_ephem_tests();
void run_planet_scene_tests();
void run_shape_tests();
void run_appearance_tests();
void run_mission_tests();
void run_stellar_system_tests();
void run_relativity_tests();
void run_scene_tests();
