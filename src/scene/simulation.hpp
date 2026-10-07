#pragma once

#include "core/sim_clock.hpp"
#include "scene/camera.hpp"
#include "scene/camera_director.hpp"
#include "scene/scene.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// The simulation as seen on one output: the camera view at the output's size
// and what depends on it. Each output (window, monitor) has its own.
struct OutputView {
    CameraView camera;
    std::vector<float> body_fades; // per body: label/marker/orbit visibility (satellite_fades)
    double content_scale = 1.0;    // output pixels per point
};

// What is shown, independent of where and how: the scene, its clock, the
// camera and the auto tour, and the animation clocks of the effects. A host
// (window, screensaver, wallpaper) drives it with real time and input, and
// SceneRenderer draws it.
class Simulation {
public:
    // Loads a scene file and starts at the scene's own time and view; a running
    // tour goes on in the new scene. On failure nothing changes.
    bool load_scene(const std::filesystem::path& file, std::string* error);

    // Advances the clock, the scene, the tour and the camera by `real_dt` seconds.
    void update(double real_dt);
    // The view on an output of `width` x `height` pixels with `content_scale`
    // pixels per point, as of the last update().
    void compute_view(int width, int height, double content_scale, OutputView& out) const;

    void reset_to_now();
    void set_focus(int body); // stops the tour
    void set_frame(int frame);
    void jump_to_event(size_t index);
    // Views the focus from (slightly off) the sun's direction.
    void look_from_sun_side();

    void start_tour();
    void stop_tour();
    bool touring() const { return m_director.active(); }

    const Scene& scene() const { return m_scene; }
    SimClock& clock() { return m_clock; }
    const SimClock& clock() const { return m_clock; }
    OrbitCamera& camera() { return m_camera; }
    const OrbitCamera& camera() const { return m_camera; }
    const CameraDirector& director() const { return m_director; }

    // Accretion-disk animation clock (units of M), and the real seconds the
    // clock has run (pulsar beams, plume and comet animation); both wrapped to
    // keep float precision.
    double disk_time() const { return m_disk_time; }
    double animation_time() const { return m_animation_time; }
    // Real seconds since the scene was loaded, paused or not (the title card).
    double scene_age() const { return m_scene_age; }
    // The event last jumped to (-1: none since the scene was loaded) and the
    // real seconds since (its caption).
    int last_event() const { return m_last_event; }
    double event_age() const { return m_event_age; }
    // Moves the caption's clock (the host shows or dismisses the caption).
    void set_event_age(double age) { m_event_age = age; }

private:
    // Applies (or, when the shot ends, undoes) a time warp asked for by the tour.
    void apply_tour_warp();

    Scene m_scene;
    SimClock m_clock;
    OrbitCamera m_camera;
    CameraDirector m_director;

    double m_tour_warp = 0.0;  // warp currently imposed by the tour (0 = none)
    double m_saved_warp = 0.0; // the user's warp, restored afterwards
    double m_disk_time = 0.0;
    double m_animation_time = 0.0;
    double m_scene_age = 0.0;
    int m_last_event = -1;
    double m_event_age = 0.0;
};

} // namespace astraxis
