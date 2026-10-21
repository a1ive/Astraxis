#pragma once

#include "core/sim_clock.hpp"
#include "scene/camera.hpp"
#include "scene/camera_director.hpp"
#include "scene/light_curve.hpp"
#include "scene/scene.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// When an event's caption is shown, in real seconds after the jump: it waits
// for the camera's flight, fades in, stays for a reading time that grows with
// its length, and fades out.
struct CaptionTimes {
    double fade_in = 0.0;  // starts fading in
    double shown = 0.0;    // fully shown
    double fade_out = 0.0; // starts fading out
    double gone = 0.0;
};
CaptionTimes caption_times(const std::string& caption);

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
//
// The tour's event shots go through the scene's events in time order: the
// picture fades to black, the clock jumps to shortly before the event (or,
// for an event coming up anyway, just keeps going), the warp eases down to the
// event's own as its time comes (for a flyby: shortly before closest
// approach), the camera frames the flyby, and the event's caption is shown
// while the shot lasts.
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
    // Ends the active event's trails, also when `frame` is already active.
    void set_frame(int frame);
    // Jumps to an event and frames it as the tour does (across the flyby).
    void jump_to_event(size_t index);
    // The closest approach of the event's `flyby` body to its focus, within
    // 300 s at `warp` (simulated seconds per second) after the event's time;
    // none if the event names no flyby.
    Flyby event_flyby(size_t index, double warp) const;
    // Views the focus from (slightly off) the sun's direction.
    void look_from_sun_side();

    // `seed` 0: a random tour.
    void start_tour(uint32_t seed = 0);
    void stop_tour();
    bool touring() const { return m_director.active(); }
    // The event the tour is showing or heading for (-1: none).
    int tour_event() const { return m_tour_event; }
    // The clock has reached it (its caption is up), rather than being on the way.
    bool tour_event_arrived() const { return m_tour_event >= 0 && m_event_phase == EventPhase::Hold; }
    // Brightness of the picture, 0..1: the tour fades to black to move the clock.
    double fade() const;
    // Brightness of the starry sky, 0..1: it fades out while a rotating display frame
    // would spin it past faster than the eye can follow (turning with Janus at days
    // per second).
    double sky_visibility() const;

    const Scene& scene() const { return m_scene; }
    SimClock& clock() { return m_clock; }
    const SimClock& clock() const { return m_clock; }
    OrbitCamera& camera() { return m_camera; }
    const OrbitCamera& camera() const { return m_camera; }
    const CameraDirector& director() const { return m_director; }
    // The scene's light curve around the clock (empty without one).
    const LightCurveSamples& light_curve() const { return m_light_curve; }

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
    enum class EventPhase {
        FadeOut,  // the previous shot fades to black, then the clock jumps
        Approach, // the clock runs up to the event, slowing down
        Hold,     // at the event's warp, with its caption
    };

    // Applies (or, when the shot ends, undoes) a time warp asked for by the tour.
    void apply_tour_warp();
    // The warp the user chose (the tour's own warp put aside).
    double user_warp() const { return m_tour_warp > 0.0 ? m_saved_warp : m_clock.warp; }

    // Event shots of the tour.
    bool can_show_event(const SceneEvent& e) const;
    // The next event to show, in time order: the first not yet shown in this
    // round after the last one the tour showed (at the start of the tour: after
    // the clock), else the earliest not shown; -1 if none. An event shown at
    // weeks per second runs the clock past the next ones, which are still shown.
    int next_tour_event() const;
    // Simulated time before the event at which the approach starts.
    double approach_lead(int index) const;
    // Close enough ahead to be reached without moving the clock.
    bool event_imminent(int index) const;
    void begin_tour_event();
    void jump_to_tour_event();
    void aim_at_tour_event(double duration);
    void run_tour_event(double real_dt);
    void end_tour_event();
    // The camera shot of event `index`: from where the event says, else
    // across the path of `flyby` (if any) turned toward the sun, else from
    // the default angles. The scene must be at about the event's time.
    CameraShot event_shot(size_t index, const Flyby& flyby) const;

    Scene m_scene;
    SimClock m_clock;
    OrbitCamera m_camera;
    CameraDirector m_director;
    LightCurveSamples m_light_curve;

    double m_tour_warp = 0.0;  // warp currently imposed by the tour (0 = none)
    double m_saved_warp = 0.0; // the user's warp, restored afterwards

    int m_tour_event = -1;
    EventPhase m_event_phase = EventPhase::Hold;
    double m_event_warp = 0.0;     // the event's own warp
    double m_event_max_warp = 0.0; // the approach starts at this warp
    double m_event_shot_warp = 0.0; // the warp now
    Flyby m_event_flyby;           // the pass the event shows, if any
    double m_event_arrival = 0.0;  // when the approach ends: the event's time, or shortly before the flyby
    std::vector<bool> m_events_shown; // in this round of the tour
    int m_last_tour_event = -1;       // the last event the tour showed (-1: none since it started)
    double m_fade = 1.0;              // linear; fade() eases it
    double m_fade_target = 1.0;
    double m_sky = 1.0; // linear; sky_visibility() eases it
    double m_disk_time = 0.0;
    double m_animation_time = 0.0;
    double m_scene_age = 0.0;
    int m_last_event = -1;
    double m_event_age = 0.0;
};

} // namespace astraxis
