// Minimal self-checking tests for core / ephem / scene and the settings file
// (no framework).

#include "test_util.hpp"

#include <cstdio>

int main()
{
    run_settings_tests();
    run_core_tests();
    run_ephem_tests();
    run_planet_scene_tests();
    run_shape_tests();
    run_appearance_tests();
    run_mission_tests();
    run_stellar_system_tests();
    run_relativity_tests();
    run_scene_tests();

    const int failures = failure_count();
    if (failures == 0) {
        std::printf("All tests passed.\n");
        return 0;
    }
    std::printf("%d failure(s).\n", failures);
    return 1;
}
