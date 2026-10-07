// Shared helpers for the self-checking tests (no framework). Each *_tests.cpp
// runs its domain's tests from run_<domain>_tests(); test_main.cpp calls them all.

#pragma once

#include "scene/scene.hpp"

#include <glm/vec3.hpp>

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
