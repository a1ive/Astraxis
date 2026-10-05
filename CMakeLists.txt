cmake_minimum_required(VERSION 3.28)

project(Astraxis VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(VCLTL)
include(Dependencies)
include(Shaders)

function(astraxis_set_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive- /utf-8)
        # Allow portable C library calls (sscanf, ...) without the *_s variants.
        target_compile_definitions(${target} PRIVATE _CRT_SECURE_NO_WARNINGS)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    target_compile_definitions(${target} PUBLIC GLM_FORCE_DEPTH_ZERO_TO_ONE GLM_FORCE_RIGHT_HANDED)
endfunction()

# --- Simulation library (core / ephem / scene): no SDL or GPU dependencies ---

add_library(astraxis_sim STATIC
    src/core/color.cpp
    src/core/color.hpp
    src/core/math.hpp
    src/core/sim_clock.hpp
    src/core/time.cpp
    src/core/time.hpp
    src/ephem/ephemeris.cpp
    src/ephem/ephemeris.hpp
    src/ephem/kerr_geodesic.cpp
    src/ephem/kerr_geodesic.hpp
    src/ephem/kerr_null.cpp
    src/ephem/kerr_null.hpp
    src/ephem/kerr_orbit.cpp
    src/ephem/kerr_orbit.hpp
    src/ephem/nbody.cpp
    src/ephem/nbody.hpp
    src/ephem/post_keplerian.cpp
    src/ephem/post_keplerian.hpp
    src/ephem/visual_orbit.cpp
    src/ephem/visual_orbit.hpp
    src/ephem/jpl_approx_planet.cpp
    src/ephem/jpl_approx_planet.hpp
    src/ephem/kepler.cpp
    src/ephem/kepler.hpp
    src/ephem/mean_element_orbit.cpp
    src/ephem/mean_element_orbit.hpp
    src/ephem/motion.hpp
    src/scene/camera.cpp
    src/scene/camera.hpp
    src/scene/camera_director.cpp
    src/scene/camera_director.hpp
    src/scene/comet.cpp
    src/scene/comet.hpp
    src/scene/label_layout.cpp
    src/scene/label_layout.hpp
    src/scene/scene.cpp
    src/scene/scene.hpp
    src/scene/scene_loader.cpp
    src/scene/scene_loader.hpp
    src/scene/shape_model.cpp
    src/scene/shape_model.hpp
    src/scene/simulation.cpp
    src/scene/simulation.hpp
    src/scene/star_catalog.cpp
    src/scene/star_catalog.hpp
)
target_include_directories(astraxis_sim PUBLIC src)
target_link_libraries(astraxis_sim PUBLIC glm::glm-header-only PRIVATE tomlplusplus::tomlplusplus)
astraxis_set_warnings(astraxis_sim)

# --- Shaders ---------------------------------------------------------------

astraxis_add_shaders(astraxis_shaders
    shaders/atmosphere.vert.hlsl
    shaders/atmosphere.frag.hlsl
    shaders/beam.vert.hlsl
    shaders/beam.frag.hlsl
    shaders/belt.vert.hlsl
    shaders/belt.frag.hlsl
    shaders/bloom_downsample.frag.hlsl
    shaders/black_hole.frag.hlsl
    shaders/black_hole_composite.frag.hlsl
    shaders/milky_way.frag.hlsl
    shaders/bloom_upsample.frag.hlsl
    shaders/body.vert.hlsl
    shaders/body.frag.hlsl
    shaders/comet.frag.hlsl
    shaders/dust.vert.hlsl
    shaders/dust.frag.hlsl
    shaders/composite.frag.hlsl
    shaders/fullscreen.vert.hlsl
    shaders/plume.vert.hlsl
    shaders/plume.frag.hlsl
    shaders/orbit_line.vert.hlsl
    shaders/orbit_line.frag.hlsl
    shaders/ring.vert.hlsl
    shaders/ring.frag.hlsl
    shaders/starfield.vert.hlsl
    shaders/starfield.frag.hlsl
    shaders/sun.vert.hlsl
    shaders/sun.frag.hlsl
)

# --- Executable ------------------------------------------------------------

add_executable(astraxis
    src/app/main.cpp
    src/app/app.cpp
    src/app/app.hpp
    src/app/app_ui.cpp
    src/platform/paths.cpp
    src/platform/paths.hpp
    src/platform/window.cpp
    src/platform/window.hpp
    src/render/atmosphere_pass.cpp
    src/render/atmosphere_pass.hpp
    src/render/beam_pass.cpp
    src/render/beam_pass.hpp
    src/render/belt_pass.cpp
    src/render/belt_pass.hpp
    src/render/black_hole_pass.cpp
    src/render/black_hole_pass.hpp
    src/render/body_pass.cpp
    src/render/body_pass.hpp
    src/render/comet_pass.cpp
    src/render/comet_pass.hpp
    src/render/dust_pass.cpp
    src/render/dust_pass.hpp
    src/render/gpu_buffer.cpp
    src/render/gpu_buffer.hpp
    src/render/gpu_device.cpp
    src/render/gpu_device.hpp
    src/render/orbit_line_pass.cpp
    src/render/orbit_line_pass.hpp
    src/render/plume_pass.cpp
    src/render/plume_pass.hpp
    src/render/post_process.cpp
    src/render/post_process.hpp
    src/render/render_output.cpp
    src/render/render_output.hpp
    src/render/scene_targets.cpp
    src/render/scene_targets.hpp
    src/render/ring_pass.cpp
    src/render/ring_pass.hpp
    src/render/shader.cpp
    src/render/shader.hpp
    src/render/texture.cpp
    src/render/texture.hpp
    src/render/starfield_pass.cpp
    src/render/starfield_pass.hpp
    src/render/sun_pass.cpp
    src/render/sun_pass.hpp
    src/view/scene_renderer.cpp
    src/view/scene_renderer.hpp
)

target_link_libraries(astraxis PRIVATE
    astraxis_sim
    SDL3::SDL3
    imgui
    stb
    astraxis_shaders
)
astraxis_set_warnings(astraxis)

if(WIN32)
    target_sources(astraxis PRIVATE res/astraxis.rc)
endif()

# Runtime assets (scenes, star catalog, textures) live next to the executable.
add_custom_command(TARGET astraxis POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/assets" "$<TARGET_FILE_DIR:astraxis>/assets"
    COMMENT "Copying assets"
)

if(MSVC)
    # No console window in Release; keep it in Debug for log output.
    set_target_properties(astraxis PROPERTIES
        WIN32_EXECUTABLE $<CONFIG:Release>
        VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT astraxis)
endif()

# --- Tests -----------------------------------------------------------------

enable_testing()

add_executable(astraxis_tests tests/sim_tests.cpp)
target_link_libraries(astraxis_tests PRIVATE astraxis_sim)
target_compile_definitions(astraxis_tests PRIVATE ASTRAXIS_ASSET_DIR="${CMAKE_CURRENT_SOURCE_DIR}/assets")
astraxis_set_warnings(astraxis_tests)
add_test(NAME astraxis_tests COMMAND astraxis_tests)
