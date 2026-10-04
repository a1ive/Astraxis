#pragma once

#include "core/sim_clock.hpp"
#include "platform/window.hpp"
#include "render/beam_pass.hpp"
#include "render/belt_pass.hpp"
#include "render/black_hole_pass.hpp"
#include "render/body_pass.hpp"
#include "render/orbit_line_pass.hpp"
#include "render/post_process.hpp"
#include "render/renderer.hpp"
#include "render/ring_pass.hpp"
#include "render/starfield_pass.hpp"
#include "render/sun_pass.hpp"
#include "scene/camera.hpp"
#include "scene/camera_director.hpp"
#include "scene/scene.hpp"

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

union SDL_Event;

namespace astraxis {

// Command line: --scene <file stem> (default solar_system), --event <n> (1-based,
// in the order of the Events list) to start at that event.
struct LaunchOptions {
    std::string scene = "solar_system";
    int event = 0; // 0 = none
};

class App {
public:
    bool init(const LaunchOptions& options);
    void run();
    void shutdown();

private:
    void handle_events();
    void handle_key(const SDL_Event& event);
    void update(double real_dt);
    void build_ui();
    void build_control_panel();
    void build_info_panel(); // details of the focused body
    void build_labels();
    void render();

    void build_body_items();
    void build_orbit_lines();
    void build_belt_items();
    void reset_to_now();
    void set_focus(int body);

    // Scenes and their GPU resources.
    bool load_scene(size_t index);
    void release_body_textures();
    // Renders the background (stars, Milky Way) into the lensing cube map,
    // creating it on first use.
    void update_sky_cube();
    // Applies (or, when the shot ends, undoes) a time warp asked for by the tour.
    void apply_tour_warp();
    void create_sky_cube();
    // Traces the (first visible) black hole off-screen; composited in the scene pass.
    void trace_black_hole(SDL_GPUCommandBuffer* cmd, uint32_t width, uint32_t height);
    void set_frame(int frame);
    void jump_to_event(size_t index);

    // Auto tour (idle mode).
    void start_tour();
    void stop_tour();
    // Any user input; `takes_camera` stops the tour (camera control returns to the user).
    void note_activity(bool takes_camera);
    bool panel_visible() const;

    Window m_window;
    Renderer m_renderer;
    StarfieldPass m_starfield;
    BodyPass m_bodies;
    RingPass m_rings;
    OrbitLinePass m_orbits;
    BeltPass m_belts;
    SunPass m_sun;
    BeamPass m_beams;
    BlackHolePass m_black_hole;
    SDL_GPUTexture* m_sky_cube = nullptr; // starfield cube map for lensing (created on demand)
    SDL_GPUTextureFormat m_sky_cube_format = SDL_GPU_TEXTUREFORMAT_INVALID;
    uint32_t m_sky_cube_size = 0;
    SDL_GPUTexture* m_milky_way = nullptr; // the scene's Milky Way map, if any
    double m_tour_warp = 0.0;  // warp currently imposed by the tour (0 = none)
    double m_saved_warp = 0.0; // the user's warp, restored afterwards
    double m_disk_time = 0.0;             // accretion-disk animation clock (units of M)
    double m_pulsar_time = 0.0;           // real seconds the clock has run (pulsar beam sweep)
    std::vector<CatalogStar> m_catalog;   // the catalog sky (seen from Earth)
    bool m_scene_sky_stars = false;       // the starfield shows the scene's own stars
    std::vector<BeamDrawItem> m_beam_items;
    PostProcess m_post;
    PostSettings m_post_settings;

    std::filesystem::path m_asset_dir;
    std::vector<std::filesystem::path> m_scene_files;
    size_t m_scene_index = 0;
    std::string m_scene_error;
    Scene m_scene;
    std::vector<SDL_GPUTexture*> m_body_textures; // per body, may be null
    std::vector<SDL_GPUTexture*> m_ring_textures; // per body: ring profile, may be null
    std::vector<int> m_body_meshes;               // per body: BodyPass mesh index (-1: ellipsoid)
    OrbitCamera m_camera;
    CameraDirector m_director;
    SimClock m_clock;
    CameraView m_view;

    bool m_running = false;
    bool m_imgui_ready = false;
    bool m_dragging = false;
    uint64_t m_last_counter = 0;

    // View options.
    bool m_show_ui = true;
    bool m_show_orbits = true;
    bool m_show_belts = true;
    bool m_show_labels = true;
    bool m_show_info = true;
    bool m_show_demo = false;
    bool m_auto_tour = true;     // start the tour after kIdleSeconds without input
    double m_idle_time = 0.0;    // seconds since the last input
    double m_pointer_idle = 0.0; // seconds since the mouse last moved
    bool m_cursor_hidden = false;
    float m_star_brightness = 1.0f;
    float m_line_width = 1.6f;

    // "Go to date" input.
    char m_goto_text[32] = "";
    bool m_goto_error = false;

    // Per-frame scratch.
    std::vector<BodyDrawItem> m_body_items;
    std::vector<RingDrawItem> m_ring_items;
    std::vector<BeltDrawItem> m_belt_items;
    std::vector<int> m_belt_ids;          // per scene belt: BeltPass index (-1: not uploaded)
    std::vector<float> m_belt_reference_h; // per scene belt: median absolute magnitude
    std::vector<double> m_belt_radius_km;  // per scene belt: median semi-major axis
    std::vector<glm::dvec3> m_trail_points;
    std::vector<float> m_trail_fades;
    std::vector<float> m_body_fades;  // per body: label/marker/orbit visibility (satellite_fades)
    std::vector<float> m_label_alpha; // per body: displayed label opacity, eased toward its target
    std::vector<glm::vec4> m_line_points;
};

} // namespace astraxis
