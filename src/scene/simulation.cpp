#include "scene/simulation.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "scene/jovian_plasma_sheet.hpp"
#include "scene/satellite_aurora.hpp"
#include "scene/scene_loader.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <utility>

namespace astraxis {

namespace {

// A moon's label, marker and orbit fade out as its orbit shrinks on screen
// from kSatelliteShowPx to kSatelliteHidePx (points, not pixels).
constexpr double kSatelliteHidePx = 10.0;
constexpr double kSatelliteShowPx = 28.0;
constexpr double kDiskMaxWarp = 60.0; // cap on the disk animation speed (ISCO period ~92 M)

// Event caption timing, real seconds after the jump: the camera flies for 2 s
// (Simulation::jump_to_event); the caption stays for a leisurely reading time.
constexpr double kCaptionDelay = 1.5;
constexpr double kCaptionFadeIn = 0.8;
constexpr double kCaptionFadeOut = 1.5;
constexpr double kCaptionWordsPerSecond = 3.0;
constexpr double kCaptionExtraSeconds = 4.0;
constexpr double kCaptionMinSeconds = 8.0;
constexpr double kCaptionMaxSeconds = 45.0;

// Event shots of the tour.
constexpr double kEventFadeSeconds = 1.2; // to black, and back after the jump
// On the approach the warp is the simulated time still to go over this many
// real seconds: it falls exponentially, by e every kApproachSeconds...
constexpr double kApproachSeconds = 3.0;
// ...from at most kLeadWarpRatio times the event's warp (a flyby at 1 min/s
// starts 3.2 hours before, at about an hour per second, and takes 12 s to
// slow down), but no faster than kMaxLeadWarp: an event that is itself shown
// at days per second is simply jumped to.
constexpr double kLeadWarpRatio = 64.0;
constexpr double kMaxLeadWarp = 86400.0;
// An event coming up within this many real seconds at the user's warp is
// waited for (the current shot is cut short if it would run past it).
constexpr double kMaxWaitSeconds = 20.0;
constexpr double kEventTransitionSeconds = 5.0; // camera flight to an event reached without a jump
// The shot holds the event for its caption, and if a spacecraft passes
// close by, until a little after its closest approach.
constexpr double kMinEventSeconds = 25.0;
constexpr double kMaxEventSeconds = 100.0;
constexpr double kAfterFlybySeconds = 15.0;
// The approach of a flyby ends this long (real seconds at the event's warp)
// before closest approach, so that it comes while the caption is up. Events
// are often set an hour or so before it.
constexpr double kBeforeFlybySeconds = 20.0;
// An event's flyby is looked for this long (real seconds at the event's
// warp) after the event's time: some are set hours before closest approach.
constexpr double kFlybySearchSeconds = 300.0;
// A flyby is seen across its path, turned this far (tan) toward the sun so
// that the body is mostly lit.
constexpr double kFlybySunTilt = 0.6;
constexpr double kEventDrift = 0.004; // rad/s

// The sky turns with a rotating display frame: fully shown up to the first rate
// (deg per real second; the Earth-Moon frame at 2 days per second turns 26),
// gone from the second, easing over kSkyFadeSeconds.
constexpr double kSkyTurnShownDegPerS = 45.0;
constexpr double kSkyTurnHiddenDegPerS = 180.0;
constexpr double kSkyFadeSeconds = 0.8;
// A pinwheel's spiral turning faster than the eye can follow (at 60 frames per
// second, 4 turns per second jump 24 deg a frame) fades out the same way.
constexpr double kPinwheelShownTurnsPerS = 1.0;
constexpr double kPinwheelHiddenTurnsPerS = 4.0;
// The dust's exposure follows pinwheel_exposure with this time constant (real
// seconds), in log brightness.
constexpr double kPinwheelExposureSeconds = 1.5;

double smoothstep(double x)
{
    x = std::clamp(x, 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

double event_warp(const SceneEvent& e, double user_warp)
{
    return e.warp > 0.0 ? e.warp : user_warp;
}

double lead_warp(double warp)
{
    return std::max(warp, std::min(warp * kLeadWarpRatio, kMaxLeadWarp));
}

// Above this warp (2 h/s, Io's eclipse takes about a second) its aurora keeps the level it
// has between eclipses: the search below would run nearly every frame, for a flicker.
constexpr double kIoEclipseMaxWarp = 7200.0;
constexpr double kIoEclipseSearchStep = 300.0; // s: much less than the 2 hours of an eclipse

// Whether a moon's center is in the shadow of its parent, a sphere of the equatorial radius
// (as the shader's lighting), without the penumbra.
bool in_parent_umbra(const Scene& scene, int moon, double t)
{
    const Body& body = scene.bodies[static_cast<size_t>(moon)];
    const glm::dvec3 center = scene.icrf_state_at(moon, t).position;
    const glm::dvec3 to_sun = glm::normalize(scene.sun_icrf_state_at(t).position - center);
    const glm::dvec3 to_parent = scene.icrf_state_at(body.parent, t).position - center;
    return glm::dot(to_parent, to_sun) > 0.0 &&
           glm::length(glm::cross(to_parent, to_sun)) < scene.bodies[static_cast<size_t>(body.parent)].equatorial_radius_km;
}

// The time between low and high (in step's seconds apart) where the umbra test turns to `to`.
double refine_umbra_crossing(const Scene& scene, int moon, double low, double high, bool to)
{
    for (int k = 0; k < 24; ++k) {
        const double middle = 0.5 * (low + high);
        (in_parent_umbra(scene, moon, middle) == to ? high : low) = middle;
    }
    return high;
}

} // namespace

CaptionTimes caption_times(const std::string& caption)
{
    size_t words = 0;
    bool in_word = false;
    for (const char c : caption) {
        const bool space = c == ' ' || c == '\n';
        if (!space && !in_word) {
            ++words;
        }
        in_word = !space;
    }
    const double reading = std::clamp(static_cast<double>(words) / kCaptionWordsPerSecond + kCaptionExtraSeconds,
                                      kCaptionMinSeconds, kCaptionMaxSeconds);
    CaptionTimes t;
    t.fade_in = kCaptionDelay;
    t.shown = t.fade_in + kCaptionFadeIn;
    t.fade_out = t.shown + reading;
    t.gone = t.fade_out + kCaptionFadeOut;
    return t;
}

bool Simulation::load_scene(const std::filesystem::path& file, std::string* error)
{
    Scene scene;
    if (!load_scene_file(file, scene, error)) {
        return false;
    }
    m_scene = std::move(scene);
    m_rotor_turns.assign(m_scene.bodies.size(), 0.0);
    m_io = -1;
    m_io_eclipse_valid = false;
    for (size_t k = 0; k < m_scene.bodies.size(); ++k) {
        if (m_scene.bodies[k].aurora.model == AuroraModel::Io && m_scene.bodies[k].parent >= 0) {
            m_io = static_cast<int>(k);
        }
    }
    m_scene_age = 0.0;
    m_last_event = -1;
    m_events_shown.assign(m_scene.events.size(), false);

    const bool touring = m_director.active();
    stop_tour();
    if (m_scene.view.start_now) {
        reset_to_now();
    } else {
        m_clock.t_tdb = m_scene.view.start_tdb;
    }
    m_clock.warp = m_scene.view.warp;
    m_scene.update(m_clock.t_tdb);
    m_exposure_t_tdb = m_clock.t_tdb;
    // The plot's range: the deepest eclipse within a few spans of the start.
    measure_light_curve_depth(m_scene, m_clock.t_tdb, 4.0 * m_scene.light_curve.span_days);
    m_light_curve.reset();
    m_light_curve.update(m_scene, m_clock.t_tdb);

    m_camera = OrbitCamera{};
    m_camera.set_up_axis(m_scene.up_axis());
    const int focus = m_scene.view.focus;
    m_camera.focus(focus, m_scene);
    if (m_scene.view.distance_km > 0.0) {
        m_camera.fly_to({focus, m_scene.view.yaw_rad, m_scene.view.pitch_rad, m_scene.view.distance_km, false}, 0.0, m_scene);
    }
    m_camera.update(0.0, m_scene);
    if (m_scene.pinwheel.enabled()) {
        m_pinwheel_exposure = astraxis::pinwheel_exposure(m_scene.pinwheel, m_clock.t_tdb, pinwheel_camera_km());
    }
    if (m_scene.view.initial_event >= 0) {
        jump_to_event(static_cast<size_t>(m_scene.view.initial_event));
        m_camera.update(0.0, m_scene);
    }
    if (touring) {
        start_tour();
    }
    return true;
}

void Simulation::update(double real_dt)
{
    // Set from outside since the last frame (an event, the time slider): a jump.
    bool jumped = m_clock.t_tdb != m_exposure_t_tdb;
    m_clock.advance(real_dt);
    const double t_run = m_clock.t_tdb;
    m_scene.update(m_clock.t_tdb);
    m_scene_age += real_dt;
    m_event_age += real_dt;

    const double turn = m_clock.paused ? 0.0 : m_scene.frame_turn_rate(m_clock.t_tdb) * m_clock.warp * kRadToDeg;
    const double sky_target =
        std::clamp((kSkyTurnHiddenDegPerS - turn) / (kSkyTurnHiddenDegPerS - kSkyTurnShownDegPerS), 0.0, 1.0);
    const double sky_step = real_dt / kSkyFadeSeconds;
    m_sky = sky_target < m_sky ? std::max(sky_target, m_sky - sky_step) : std::min(sky_target, m_sky + sky_step);
    const Pinwheel& pinwheel = m_scene.pinwheel;
    const double spiral_turns =
        m_clock.paused || !pinwheel.enabled() ? 0.0 : std::abs(m_clock.warp) / pinwheel.period_s;
    const double pinwheel_target = std::clamp(
        (kPinwheelHiddenTurnsPerS - spiral_turns) / (kPinwheelHiddenTurnsPerS - kPinwheelShownTurnsPerS), 0.0, 1.0);
    m_pinwheel = pinwheel_target < m_pinwheel ? std::max(pinwheel_target, m_pinwheel - sky_step)
                                              : std::min(pinwheel_target, m_pinwheel + sky_step);

    // Accretion-disk animation: follows the clock, but capped so that the inner
    // disk does not blur at high time warp.
    for (const Body& body : m_scene.bodies) {
        if (body.kind == BodyKind::BlackHole && body.gm_km3_s2 > 0.0) {
            const double t_m = body.gm_km3_s2 / std::pow(kSpeedOfLightKmS, 3.0); // seconds per M
            const double rate = m_clock.paused ? 0.0 : std::min(m_clock.warp, kDiskMaxWarp);
            m_disk_time = std::fmod(m_disk_time + real_dt * rate * (m_clock.reverse ? -1.0 : 1.0) / t_m, 1.0e6);
            break;
        }
    }

    if (!m_clock.paused) {
        m_animation_time = std::fmod(m_animation_time + real_dt, 1.0e6);
    }
    // Rotors: at their own rate, but no faster on screen than kRotorMaxTurnsPerSecond.
    for (size_t i = 0; i < m_scene.bodies.size() && i < m_rotor_turns.size(); ++i) {
        const Body::Rotor& rotor = m_scene.bodies[i].rotor;
        if (rotor.shape && !m_clock.paused) {
            const double rate = std::min(m_clock.warp / rotor.period_s, kRotorMaxTurnsPerSecond);
            const double turns = m_rotor_turns[i] + real_dt * rate * (m_clock.reverse ? -1.0 : 1.0);
            m_rotor_turns[i] = turns - std::floor(turns);
        }
    }

    // The aurora: Kp over a window that grows with the warp, and its animation clock.
    bool has_aurora = false;
    for (const Body& body : m_scene.bodies) {
        has_aurora = has_aurora || bool(body.aurora.kp) || bool(body.aurora.jovian);
        if (body.aurora.kp) {
            const double window = std::max(3.0 * 3600.0, std::abs(m_clock.warp) * kAuroraKpWindowSeconds);
            m_aurora_kp = body.aurora.kp->mean(m_clock.t_tdb - 0.5 * window, m_clock.t_tdb + 0.5 * window);
            break;
        }
    }
    if (has_aurora && !m_clock.paused) {
        const double rate = std::min(std::abs(m_clock.warp), kAuroraMaxWarp);
        m_aurora_time = std::fmod(m_aurora_time + real_dt * rate * (m_clock.reverse ? -1.0 : 1.0), 1.0e5);
    }
    update_io_eclipse();

    if (m_director.active()) {
        const int next = next_tour_event();
        const bool imminent = next >= 0 && event_imminent(next);
        m_director.set_events(next >= 0, imminent);
        // An event coming up during the shot: cut the shot short rather than run past it.
        if (imminent && m_tour_event < 0 &&
            m_scene.events[static_cast<size_t>(next)].t_tdb - m_clock.t_tdb < m_clock.warp * m_director.time_left()) {
            m_director.hold_event(0.0);
        }
    }
    m_director.set_time_warp(m_clock.paused ? 0.0 : user_warp());
    m_director.update(real_dt, m_scene, m_camera);
    if (m_director.take_event_request()) {
        begin_tour_event();
    } else if (m_tour_event >= 0 && m_director.current_shot() != CameraDirector::ShotKind::Event) {
        end_tour_event();
    }
    run_tour_event(real_dt);
    apply_tour_warp();
    m_camera.update(real_dt, m_scene);

    // Use the final clock after tour jumps and warp changes, as the plot does.
    // Keep cached samples while hidden; catch up before the plot reappears.
    if (light_curve_visibility(m_scene.light_curve, m_clock.paused ? 0.0 : m_clock.warp) > 0.0) {
        m_light_curve.update(m_scene, m_clock.t_tdb);
    }

    // The eye adapts to the dust over a few seconds; across a jump of the
    // clock (the tour's, behind a fade to black) it is already adapted.
    if (pinwheel.enabled()) {
        jumped = jumped || m_clock.t_tdb != t_run;
        const double target = astraxis::pinwheel_exposure(pinwheel, m_clock.t_tdb, pinwheel_camera_km());
        if (jumped) {
            m_pinwheel_exposure = target;
        } else {
            const double k = 1.0 - std::exp(-real_dt / kPinwheelExposureSeconds);
            m_pinwheel_exposure *= std::pow(target / m_pinwheel_exposure, k);
        }
    }
    m_exposure_t_tdb = m_clock.t_tdb;
}

double Simulation::pinwheel_camera_km() const
{
    const Body& center = m_scene.bodies[static_cast<size_t>(m_scene.pinwheel.center)];
    return glm::length(m_camera.position() - center.world_position);
}

double Simulation::rotor_turns(int body) const
{
    return body >= 0 && static_cast<size_t>(body) < m_rotor_turns.size() ? m_rotor_turns[static_cast<size_t>(body)]
                                                                           : 0.0;
}

double Simulation::fade() const
{
    return smoothstep(m_fade);
}

double Simulation::sky_visibility() const
{
    return smoothstep(m_sky);
}

void Simulation::update_io_eclipse()
{
    if (m_io < 0 || std::abs(m_clock.warp) > kIoEclipseMaxWarp) {
        m_io_eclipse_valid = false;
        return;
    }
    const double t = m_clock.t_tdb;
    if (m_io_eclipse_valid && t >= m_io_eclipse.ingress && t < m_io_next_ingress) {
        return;
    }
    m_io_eclipse_valid = false;
    // Back to the last ingress, then forward through the egress to the next ingress. A
    // period is under two days; the search gives up on a moon that is not eclipsed.
    constexpr int kSteps = static_cast<int>(2.0 * 86400.0 / kIoEclipseSearchStep);
    double ingress = t;
    bool inside = in_parent_umbra(m_scene, m_io, t);
    bool found = false;
    for (int k = 0; k < kSteps && !found; ++k) {
        const bool earlier = in_parent_umbra(m_scene, m_io, ingress - kIoEclipseSearchStep);
        if (inside && !earlier) {
            ingress = refine_umbra_crossing(m_scene, m_io, ingress - kIoEclipseSearchStep, ingress, true);
            found = true;
        } else {
            ingress -= kIoEclipseSearchStep;
            inside = earlier;
        }
    }
    if (!found) {
        return;
    }
    // The first step that is out of the shadow, then the first that is back in it.
    double time = ingress + kIoEclipseSearchStep;
    while (in_parent_umbra(m_scene, m_io, time)) {
        time += kIoEclipseSearchStep;
    }
    const double egress = refine_umbra_crossing(m_scene, m_io, time - kIoEclipseSearchStep, time, false);
    for (int k = 0; k < kSteps && !in_parent_umbra(m_scene, m_io, time); ++k) {
        time += kIoEclipseSearchStep;
    }
    m_io_next_ingress = refine_umbra_crossing(m_scene, m_io, time - kIoEclipseSearchStep, time, true);
    m_io_eclipse = {ingress, egress};
    m_io_eclipse_valid = t < m_io_next_ingress;
}

bool Simulation::io_eclipse(MoonEclipse& eclipse) const
{
    eclipse = m_io_eclipse;
    return m_io_eclipse_valid;
}

MoonAuroraState Simulation::moon_aurora(int body) const
{
    MoonAuroraState state;
    if (body < 0 || static_cast<size_t>(body) >= m_scene.bodies.size()) {
        return state;
    }
    const Body& moon = m_scene.bodies[static_cast<size_t>(body)];
    const bool placed = glm::dot(moon.in_parent_axes, moon.in_parent_axes) > 0.0;
    if (moon.aurora.model == AuroraModel::Ganymede && placed) {
        const double steady = plasma_sheet_steadiness(m_clock.warp);
        const auto shares = ganymede_sheet_shares(jovian_magnetic_latitude_deg(moon.in_parent_axes));
        state.map_shares = glm::dvec4(1.0 - steady, steady * shares[0], steady * shares[1], steady * shares[2]);
    }
    if ((moon.aurora.model == AuroraModel::Europa || moon.aurora.model == AuroraModel::Callisto) && placed) {
        const bool europa = moon.aurora.model == AuroraModel::Europa;
        const double scale = europa ? kEuropaPlasmaScaleRj : kCallistoPlasmaScaleRj;
        const double now = plasma_sheet_density(centrifugal_equator_height_rj(moon.in_parent_axes), scale);
        const double mean = plasma_sheet_mean_density(glm::length(moon.in_parent_axes), scale);
        const double steady = plasma_sheet_steadiness(m_clock.warp);
        state.plasma = std::lerp(mean, now, steady) / (europa ? europa_reference_density() : 1.0);
        if (europa) {
            // A long exposure of a sheet going by evens out above and below.
            state.field = glm::normalize(glm::mix(glm::dvec3(0.0, 0.0, 1.0), jovian_field_line_at_moon(moon.in_parent_axes), steady));
            state.sheet = -steady * jovian_magnetic_latitude_deg(moon.in_parent_axes) / kJovianDipoleTiltDeg;
        }
    }
    if (moon.aurora.model == AuroraModel::Io) {
        // Between eclipses, and when they go by too fast to follow, as before an eclipse.
        const double t = m_clock.t_tdb;
        state.sodium = (m_io_eclipse_valid ? io_sodium_level(m_io_eclipse, t) : 1.0) / kIoSodiumReferenceLevel;
        state.so2 = (m_io_eclipse_valid ? io_so2_level(m_io_eclipse, t) : 1.0) / io_so2_reference_level();
    }
    if (moon.aurora.model == AuroraModel::Mars) {
        state.mars_limb_rayleigh = MarsAurora::limb_rayleigh_at(m_clock.t_tdb);
    }
    return state;
}

double Simulation::pinwheel_visibility() const
{
    return smoothstep(m_pinwheel);
}

void Simulation::compute_view(int width, int height, double content_scale, OutputView& out) const
{
    CameraView& view = out.camera;
    view = m_camera.view(width, height);
    // Stars are ICRF directions: rotate them into the display frame.
    view.sky_view_proj = view.view_proj * glm::mat4(glm::mat3(glm::transpose(m_scene.current_transform().axes)));
    out.content_scale = content_scale;

    const double px_per_radian = 0.5 * view.viewport.y * view.focal_y;
    satellite_fades(m_scene, view.position, px_per_radian, kSatelliteHidePx * content_scale,
                    kSatelliteShowPx * content_scale, out.body_fades);
    const int target = m_camera.target();
    if (target >= 0 && target < static_cast<int>(out.body_fades.size())) {
        out.body_fades[static_cast<size_t>(target)] = 1.0f; // the focus always shows
    }
    lens_bodies(m_scene, view, out.lens);
}

void Simulation::reset_to_now()
{
    const auto since_epoch = std::chrono::system_clock::now().time_since_epoch();
    m_clock.t_tdb = tdb_from_unix_utc(std::chrono::duration<double>(since_epoch).count());
}

void Simulation::set_focus(int body)
{
    if (body >= 0 && body < static_cast<int>(m_scene.bodies.size())) {
        stop_tour();
        if (body != m_camera.target()) {
            m_camera.focus(body, m_scene);
        }
    }
}

void Simulation::set_frame(int frame)
{
    if (frame == m_scene.active_frame()) {
        m_scene.set_active_event(-1); // picked again: its own trails, not the event's
        return;
    }
    m_scene.set_active_frame(frame);
    m_scene.update(m_clock.t_tdb);
    m_camera.set_up_axis(m_scene.up_axis());
}

void Simulation::jump_to_event(size_t index)
{
    const SceneEvent& e = m_scene.events[index];
    stop_tour();
    m_last_event = static_cast<int>(index);
    m_event_age = 0.0;
    m_clock.t_tdb = e.t_tdb;
    if (e.warp > 0.0) {
        m_clock.warp = e.warp;
    }
    set_frame(e.frame >= 0 ? e.frame : 0); // events without a frame use the root frame
    m_scene.set_active_event(static_cast<int>(index)); // its trails, until the frame changes
    m_scene.update(m_clock.t_tdb);
    if (e.focus >= 0) {
        m_camera.set_drift(0.0);
        m_camera.fly_to(event_shot(index, event_flyby(index, m_clock.warp)), 2.0, m_scene);
    }
}

Flyby Simulation::event_flyby(size_t index, double warp) const
{
    const SceneEvent& e = m_scene.events[index];
    if (e.focus < 0 || e.flyby < 0) {
        return {};
    }
    return find_flyby(m_scene, e.focus, e.flyby, e.t_tdb, e.t_tdb + warp * kFlybySearchSeconds);
}

CameraShot Simulation::event_shot(size_t index, const Flyby& flyby) const
{
    const SceneEvent& e = m_scene.events[index];
    CameraShot shot;
    shot.target = e.focus;
    shot.observer = e.observer;
    shot.observer_offset_km = e.observer_offset_km;
    shot.aim_offset_km = e.aim_offset_km;
    shot.observer_body_fixed = e.observer_body_fixed;
    shot.fov_deg = e.fov_deg;
    shot.distance = e.distance_km > 0.0 ? e.distance_km : m_camera.default_distance(e.focus, m_scene);
    shot.yaw = -1.2;
    shot.pitch = 0.3;
    if (e.from_orbit_normal >= 0) {
        m_camera.angles_for_direction(m_scene.orbit_normal(e.from_orbit_normal), &shot.yaw, &shot.pitch);
    } else if (e.from_earth) {
        m_camera.angles_for_direction(m_scene.current_transform().direction_to_display(-m_scene.sky.viewer_pc),
                                     &shot.yaw, &shot.pitch);
    } else if (e.from_body >= 0 || e.from_sun) {
        const glm::dvec3 from = e.from_sun ? m_scene.sun_position()
                                           : m_scene.bodies[static_cast<size_t>(e.from_body)].world_position;
        const glm::dvec3 to_body = from - m_scene.bodies[static_cast<size_t>(e.focus)].world_position;
        m_camera.angles_for_direction(to_body, &shot.yaw, &shot.pitch);
        shot.yaw += e.phase_deg * kDegToRad;
    } else if (flyby.body >= 0 && flyby.normal_icrf != glm::dvec3(0.0)) {
        // A flyby close-up: across the path, so that it shows as the curve it
        // is, and turned toward the sun.
        glm::dvec3 dir = m_scene.current_transform().direction_to_display(flyby.normal_icrf);
        const glm::dvec3 to_light =
            m_scene.light_position(e.focus) - m_scene.bodies[static_cast<size_t>(e.focus)].world_position;
        if (glm::dot(to_light, to_light) > 0.0) {
            const glm::dvec3 sun = glm::normalize(to_light);
            dir = glm::normalize((glm::dot(dir, sun) < 0.0 ? -dir : dir) + kFlybySunTilt * sun);
        }
        m_camera.angles_for_direction(dir, &shot.yaw, &shot.pitch);
    }
    if (!std::isnan(e.elevation_deg)) {
        shot.pitch = e.elevation_deg * kDegToRad;
    }
    return shot;
}

void Simulation::look_from_sun_side()
{
    stop_tour();
    // Slightly off the sun line, so moon shadows on the planet are not hidden behind the moons.
    const glm::dvec3 target = m_scene.bodies[static_cast<size_t>(m_camera.target())].world_position;
    m_camera.look_from(m_scene.sun_position() - target);
    m_camera.rotate(0.35, 0.1);
}

void Simulation::start_tour(uint32_t seed)
{
    if (seed == 0) {
        seed = static_cast<uint32_t>(std::chrono::steady_clock::now().time_since_epoch().count()) | 1u;
    }
    m_last_tour_event = -1; // a new round; the story picks up from the clock
    m_events_shown.assign(m_scene.events.size(), false);
    const int next = next_tour_event();
    m_director.set_events(next >= 0, next >= 0 && event_imminent(next));
    m_director.set_time_warp(m_clock.paused ? 0.0 : user_warp());
    m_director.start(m_scene, m_camera, seed);
    if (m_director.take_event_request()) {
        begin_tour_event();
        run_tour_event(0.0);
    }
    apply_tour_warp();
}

void Simulation::stop_tour()
{
    if (m_director.active()) {
        m_director.stop(m_camera);
        end_tour_event();
        m_fade = 1.0;
        apply_tour_warp();
    }
}

bool Simulation::can_show_event(const SceneEvent& e) const
{
    return e.focus >= 0 && e.focus < static_cast<int>(m_scene.bodies.size());
}

int Simulation::next_tour_event() const
{
    if (m_clock.paused || m_clock.reverse) {
        return -1; // events are shown going forward
    }
    // Events in order of time, then of the file (several may share a time).
    const auto before = [&](int a, int b) {
        const double ta = m_scene.events[static_cast<size_t>(a)].t_tdb;
        const double tb = m_scene.events[static_cast<size_t>(b)].t_tdb;
        return ta < tb || (ta == tb && a < b);
    };
    const auto after_cursor = [&](int k) {
        return m_last_tour_event >= 0 ? before(m_last_tour_event, k)
                                      : m_scene.events[static_cast<size_t>(k)].t_tdb > m_clock.t_tdb;
    };
    // Among the events not shown in this round; failing that (a new round), among all.
    for (const bool new_round : {false, true}) {
        int ahead = -1;
        int earliest = -1;
        for (size_t i = 0; i < m_scene.events.size(); ++i) {
            const int k = static_cast<int>(i);
            if (!can_show_event(m_scene.events[i]) ||
                (!new_round && i < m_events_shown.size() && m_events_shown[i])) {
                continue;
            }
            if (after_cursor(k) && (ahead < 0 || before(k, ahead))) {
                ahead = k;
            }
            if (earliest < 0 || before(k, earliest)) {
                earliest = k;
            }
        }
        if (earliest >= 0) {
            return ahead >= 0 ? ahead : earliest;
        }
    }
    return -1;
}

double Simulation::approach_lead(int index) const
{
    return lead_warp(event_warp(m_scene.events[static_cast<size_t>(index)], user_warp())) * kApproachSeconds;
}

bool Simulation::event_imminent(int index) const
{
    const SceneEvent& e = m_scene.events[static_cast<size_t>(index)];
    if (m_clock.paused || m_clock.reverse || (e.frame >= 0 ? e.frame : 0) != m_scene.active_frame()) {
        return false;
    }
    const double ahead = e.t_tdb - m_clock.t_tdb;
    return ahead > 0.0 && ahead <= std::max(approach_lead(index), user_warp() * kMaxWaitSeconds);
}

void Simulation::begin_tour_event()
{
    const int index = next_tour_event();
    if (index < 0) {
        end_tour_event();
        m_director.skip_event(m_scene, m_camera);
        return;
    }
    m_events_shown.resize(m_scene.events.size(), false);
    if (m_events_shown[static_cast<size_t>(index)]) {
        m_events_shown.assign(m_scene.events.size(), false); // a new round
    }
    m_events_shown[static_cast<size_t>(index)] = true;
    m_last_tour_event = index;

    const SceneEvent& e = m_scene.events[static_cast<size_t>(index)];
    const bool imminent = event_imminent(index);
    m_tour_event = index;
    m_event_warp = event_warp(e, user_warp());
    // Coming up anyway: only ever slowed down from the user's warp.
    m_event_max_warp = imminent ? std::max(m_event_warp, user_warp()) : lead_warp(m_event_warp);
    m_event_shot_warp = m_clock.warp;

    // The pass the shot shows, if one comes close while it lasts.
    m_event_flyby = event_flyby(static_cast<size_t>(index), m_event_warp);
    m_event_arrival = e.t_tdb;
    if (m_event_flyby.body >= 0) {
        m_event_arrival = std::max(e.t_tdb, m_event_flyby.t_tdb - m_event_warp * kBeforeFlybySeconds);
    }

    if (imminent) {
        m_event_phase = EventPhase::Approach;
        m_scene.set_active_event(index);
        aim_at_tour_event(kEventTransitionSeconds);
    } else {
        m_event_phase = EventPhase::FadeOut;
        m_fade_target = 0.0;
    }
}

void Simulation::jump_to_tour_event()
{
    const SceneEvent& e = m_scene.events[static_cast<size_t>(m_tour_event)];
    m_clock.t_tdb = m_event_arrival - approach_lead(m_tour_event);
    set_frame(e.frame >= 0 ? e.frame : 0);
    m_scene.set_active_event(m_tour_event);
    m_scene.update(m_clock.t_tdb);
    aim_at_tour_event(0.0);
    m_event_phase = EventPhase::Approach;
    m_fade_target = 1.0;
}

void Simulation::aim_at_tour_event(double duration)
{
    const SceneEvent& e = m_scene.events[static_cast<size_t>(m_tour_event)];
    const CameraShot shot = event_shot(static_cast<size_t>(m_tour_event), m_event_flyby);
    double drift = (m_tour_event % 2 == 0 ? 1.0 : -1.0) * kEventDrift;
    if (e.observer >= 0 || e.from_earth || e.from_orbit_normal >= 0 || (e.from_body < 0 && !e.from_sun && m_event_flyby.body >= 0)) {
        drift = 0.0; // preserve the requested sightline
    }
    m_camera.set_drift(drift);
    m_camera.fly_to(shot, duration, m_scene);
}

void Simulation::run_tour_event(double real_dt)
{
    const double step = real_dt / kEventFadeSeconds;
    m_fade = m_fade_target < m_fade ? std::max(m_fade_target, m_fade - step) : std::min(m_fade_target, m_fade + step);
    if (m_tour_event < 0) {
        return;
    }
    const SceneEvent& e = m_scene.events[static_cast<size_t>(m_tour_event)];
    switch (m_event_phase) {
    case EventPhase::FadeOut:
        if (m_fade <= 0.0) {
            jump_to_tour_event();
        }
        break;
    case EventPhase::Approach: {
        if (m_clock.paused || m_clock.reverse) {
            m_director.hold_event(0.0); // stopped or turned back by the user: go on with the tour
            break;
        }
        const double remaining = m_event_arrival - m_clock.t_tdb;
        if (remaining > 0.0) {
            m_event_shot_warp = std::clamp(remaining / kApproachSeconds, m_event_warp, m_event_max_warp);
            break;
        }
        // Arrived: show the caption, and hold until it is read and any flyby is over.
        m_event_phase = EventPhase::Hold;
        m_event_shot_warp = m_event_warp;
        m_last_event = m_tour_event;
        m_event_age = 0.0;
        double hold = std::clamp(caption_times(e.caption).gone + 1.0, kMinEventSeconds, kMaxEventSeconds);
        if (m_event_flyby.body >= 0) {
            hold = std::max(hold, (m_event_flyby.t_tdb - m_event_arrival) / m_event_warp + kAfterFlybySeconds);
        }
        m_director.hold_event(hold);
        break;
    }
    case EventPhase::Hold:
        break;
    }
}

void Simulation::end_tour_event()
{
    m_tour_event = -1;
    m_event_flyby = Flyby{};
    m_fade_target = 1.0;
}

void Simulation::apply_tour_warp()
{
    double want = m_director.shot_warp();
    if (m_tour_event >= 0) {
        if (m_event_phase == EventPhase::FadeOut) {
            return; // the previous shot's warp until the cut
        }
        want = m_event_shot_warp;
    }
    if (want == m_tour_warp) {
        return;
    }
    if (want > 0.0) {
        if (m_tour_warp <= 0.0) {
            m_saved_warp = m_clock.warp;
        }
        m_clock.warp = want;
    } else {
        m_clock.warp = m_saved_warp;
    }
    m_tour_warp = want;
}

} // namespace astraxis
