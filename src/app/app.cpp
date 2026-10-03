#include "app/app.hpp"

#include "core/color.hpp"
#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/kerr_null.hpp"
#include "platform/paths.hpp"
#include "render/texture.hpp"
#include "scene/scene_loader.hpp"
#include "scene/star_catalog.hpp"

#include <SDL3/SDL.h>
#include <glm/gtc/matrix_transform.hpp>
#include <imgui.h>
#include <imgui_impl_sdl3.h>
#include <imgui_impl_sdlgpu3.h>

#include <algorithm>
#include <system_error>

namespace astraxis {

namespace {

// Upper bound on one frame's real time step, so that dragging the window or
// hitting a breakpoint does not cause a large jump in the simulation.
constexpr double kMaxRealDt = 0.1;

constexpr double kMouseRotateSpeed = 0.005; // radians per pixel
constexpr double kIdleSeconds = 60.0;       // auto tour starts after this long without input
constexpr double kPointerHideSeconds = 3.0; // during the tour, hide panel and cursor after this
constexpr int kTrailPoints = 2048;
constexpr double kMaxWarp = 1e8; // ~3.2 years per second
constexpr float kOrbitOpacity = 0.55f;
constexpr float kAmbient = 0.012f;
// Surface radiance of star spheres (linear HDR, times the blackbody color).
constexpr float kStarSurface = 6.0f;
// Apparent magnitude of the Sun at 1 au (V ~ -26.74); used for bolometric-ish
// apparent brightness m = kSunMagAt1Au - 2.5 log10(L / d_au^2).
constexpr double kSunMagAt1Au = -26.74;
// Body pass styles beyond SurfaceStyle (see body.frag.hlsl).
constexpr int kStyleStar = 2;
constexpr int kStyleBlackHole = 3;

// Black hole rendering.
constexpr float kRayIntegrationRadius = 150.0f; // M; weak-field deflection beyond
constexpr float kRayMaxSteps = 400.0f;
constexpr float kLensEdgeFade = 0.6f;     // fade the lensed region over its outer 40%
constexpr double kDiskMaxWarp = 60.0;     // cap on the disk animation speed (ISCO period ~92 M)

double apparent_magnitude(double luminosity_solar, double distance_km)
{
    const double d_au = distance_km / kAuKm;
    return kSunMagAt1Au - 2.5 * std::log10(luminosity_solar / (d_au * d_au));
}

#ifdef NDEBUG
constexpr bool kGpuDebug = false;
#else
constexpr bool kGpuDebug = true;
#endif

glm::vec3 to_render(const glm::dvec3& world, const glm::dvec3& camera)
{
    return glm::vec3(world - camera);
}

} // namespace

bool App::init(const LaunchOptions& options)
{
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_Init failed: %s", SDL_GetError());
        return false;
    }

    if (!m_window.create("Astraxis", 1280, 720)) {
        return false;
    }
    if (!m_renderer.init(m_window.handle(), kGpuDebug)) {
        return false;
    }
    m_asset_dir = asset_directory();
    SDL_LogInfo(SDL_LOG_CATEGORY_APPLICATION, "Assets: %s", m_asset_dir.string().c_str());

    std::vector<CatalogStar> catalog;
    std::string catalog_error;
    if (!load_star_catalog(m_asset_dir / "stars" / "bsc5.csv", catalog, &catalog_error)) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Star catalog unavailable (%s); using procedural stars",
                    catalog_error.c_str());
    }

    const SceneTargetFormat& format = m_renderer.scene_format();
    if (!m_starfield.init(m_renderer.device(), format, catalog) || !m_bodies.init(m_renderer.device(), format) ||
        !m_rings.init(m_renderer.device(), format) ||
        !m_orbits.init(m_renderer.device(), format) || !m_sun.init(m_renderer.device(), format) ||
        !m_black_hole.init(m_renderer.device(), format) ||
        !m_post.init(m_renderer.device(), format.color, m_renderer.swapchain_format())) {
        return false;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.IniFilename = nullptr; // no imgui.ini in the working directory

    const float scale = m_window.content_scale();
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 6.0f;
    style.Colors[ImGuiCol_WindowBg].w = 0.78f;
    style.ScaleAllSizes(scale);
    style.FontScaleDpi = scale;

    ImGui_ImplSDL3_InitForSDLGPU(m_window.handle());
    ImGui_ImplSDLGPU3_InitInfo init_info = {};
    init_info.Device = m_renderer.device();
    init_info.ColorTargetFormat = m_renderer.swapchain_format();
    init_info.MSAASamples = SDL_GPU_SAMPLECOUNT_1;
    ImGui_ImplSDLGPU3_Init(&init_info);
    m_imgui_ready = true;

    std::error_code ec;
    for (const auto& entry : std::filesystem::directory_iterator(m_asset_dir / "scenes", ec)) {
        if (entry.path().extension() == ".toml") {
            m_scene_files.push_back(entry.path());
        }
    }
    std::sort(m_scene_files.begin(), m_scene_files.end());
    bool found = false;
    for (size_t i = 0; i < m_scene_files.size(); ++i) {
        if (m_scene_files[i].stem() == options.scene) {
            m_scene_index = i;
            found = true;
        }
    }
    if (!found) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Scene '%s' not found", options.scene.c_str());
    }

    reset_to_now();
    if (m_scene_files.empty() || !load_scene(m_scene_index)) {
        const std::string message = m_scene_files.empty()
                                        ? "No scenes found in " + (m_asset_dir / "scenes").string()
                                        : m_scene_error;
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "%s", message.c_str());
        SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Astraxis", message.c_str(), m_window.handle());
        return false;
    }
    if (options.event > 0) {
        if (static_cast<size_t>(options.event) <= m_scene.events.size()) {
            jump_to_event(static_cast<size_t>(options.event - 1));
        } else {
            SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Scene has no event %d", options.event);
        }
    }

    m_last_counter = SDL_GetPerformanceCounter();
    m_running = true;
    return true;
}

void App::run()
{
    while (m_running) {
        handle_events();
        if (!m_running) {
            break;
        }

        const uint64_t now = SDL_GetPerformanceCounter();
        const double real_dt = static_cast<double>(now - m_last_counter) / SDL_GetPerformanceFrequency();
        m_last_counter = now;

        if (m_window.is_minimized()) {
            SDL_Delay(50);
            continue;
        }

        update(std::min(real_dt, kMaxRealDt));
        build_ui();
        render();
    }
}

void App::shutdown()
{
    if (m_renderer.device()) {
        SDL_WaitForGPUIdle(m_renderer.device());
    }
    if (m_imgui_ready) {
        ImGui_ImplSDL3_Shutdown();
        ImGui_ImplSDLGPU3_Shutdown();
        ImGui::DestroyContext();
        m_imgui_ready = false;
    }
    release_body_textures();
    if (m_sky_cube) {
        SDL_ReleaseGPUTexture(m_renderer.device(), m_sky_cube);
        m_sky_cube = nullptr;
    }
    if (m_milky_way) {
        SDL_ReleaseGPUTexture(m_renderer.device(), m_milky_way);
        m_milky_way = nullptr;
    }
    m_black_hole.shutdown();
    m_post.shutdown();
    m_sun.shutdown();
    m_orbits.shutdown();
    m_rings.shutdown();
    m_bodies.shutdown();
    m_starfield.shutdown();
    m_renderer.shutdown();
    m_window.destroy();
    SDL_Quit();
}

void App::handle_events()
{
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        ImGui_ImplSDL3_ProcessEvent(&event);
        const ImGuiIO& io = ImGui::GetIO();

        switch (event.type) {
        case SDL_EVENT_QUIT:
            m_running = false;
            break;
        case SDL_EVENT_WINDOW_CLOSE_REQUESTED:
            if (event.window.windowID == SDL_GetWindowID(m_window.handle())) {
                m_running = false;
            }
            break;
        case SDL_EVENT_KEY_DOWN:
            handle_key(event);
            break;
        case SDL_EVENT_MOUSE_BUTTON_DOWN:
            if (!io.WantCaptureMouse &&
                (event.button.button == SDL_BUTTON_LEFT || event.button.button == SDL_BUTTON_RIGHT)) {
                m_dragging = true;
                note_activity(true);
            } else {
                note_activity(false);
            }
            break;
        case SDL_EVENT_MOUSE_BUTTON_UP:
            if (event.button.button == SDL_BUTTON_LEFT || event.button.button == SDL_BUTTON_RIGHT) {
                m_dragging = false;
            }
            break;
        case SDL_EVENT_MOUSE_MOTION:
            m_pointer_idle = 0.0;
            if (m_dragging) {
                m_camera.rotate(-event.motion.xrel * kMouseRotateSpeed, event.motion.yrel * kMouseRotateSpeed);
                note_activity(true);
            } else {
                note_activity(false);
            }
            break;
        case SDL_EVENT_MOUSE_WHEEL:
            if (!io.WantCaptureMouse) {
                m_camera.zoom(event.wheel.y);
                note_activity(true);
            } else {
                note_activity(false);
            }
            break;
        default:
            break;
        }
    }
}

void App::handle_key(const SDL_Event& event)
{
    // WantCaptureKeyboard is set whenever an ImGui window has nav focus,
    // so only let ImGui keep keys while a text field is being edited.
    if (ImGui::GetIO().WantTextInput) {
        return;
    }
    const SDL_Keycode key = event.key.key;
    if (key == SDLK_ESCAPE) {
        m_running = false;
        return;
    }
    if (event.key.repeat) {
        return;
    }
    note_activity(false);

    switch (key) {
    case SDLK_A:
        if (m_director.active()) {
            stop_tour();
        } else {
            start_tour();
        }
        break;
    case SDLK_SPACE:
        m_clock.paused = !m_clock.paused;
        break;
    case SDLK_R:
        m_clock.reverse = !m_clock.reverse;
        break;
    case SDLK_N:
        reset_to_now();
        break;
    case SDLK_H:
    case SDLK_F1:
        m_show_ui = !m_show_ui;
        break;
    case SDLK_O:
        m_show_orbits = !m_show_orbits;
        break;
    case SDLK_L:
        m_show_labels = !m_show_labels;
        break;
    case SDLK_LEFTBRACKET:
        m_clock.warp = std::max(1.0, m_clock.warp / 2.0);
        break;
    case SDLK_RIGHTBRACKET:
        m_clock.warp = std::min(kMaxWarp, m_clock.warp * 2.0);
        break;
    default:
        if (key >= SDLK_1 && key <= SDLK_9) {
            set_focus(static_cast<int>(key - SDLK_1));
        }
        break;
    }
}

void App::reset_to_now()
{
    SDL_Time now_ns = 0;
    SDL_GetCurrentTime(&now_ns);
    m_clock.t_tdb = tdb_from_unix_utc(static_cast<double>(now_ns) * 1e-9);
}

void App::set_focus(int body)
{
    if (body >= 0 && body < static_cast<int>(m_scene.bodies.size())) {
        stop_tour();
        if (body != m_camera.target()) {
            m_camera.focus(body, m_scene);
        }
    }
}

bool App::load_scene(size_t index)
{
    Scene scene;
    if (!load_scene_file(m_scene_files[index], scene, &m_scene_error)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "%s", m_scene_error.c_str());
        return false;
    }
    m_scene_error.clear();

    if (m_renderer.device()) {
        SDL_WaitForGPUIdle(m_renderer.device());
    }
    release_body_textures();
    m_scene = std::move(scene);
    m_scene_index = index;

    m_body_textures.assign(m_scene.bodies.size(), nullptr);
    for (size_t i = 0; i < m_scene.bodies.size(); ++i) {
        const Body& body = m_scene.bodies[i];
        if (body.texture.empty()) {
            continue;
        }
        std::string error;
        m_body_textures[i] = load_texture_srgb(m_renderer.device(), m_asset_dir / body.texture, &error);
        if (!m_body_textures[i]) {
            SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Texture for %s unavailable: %s", body.name.c_str(),
                        error.c_str());
        }
    }

    m_starfield.set_milky_way(nullptr, 0, 0.0f);
    if (m_milky_way) {
        SDL_ReleaseGPUTexture(m_renderer.device(), m_milky_way);
        m_milky_way = nullptr;
    }
    if (!m_scene.sky.milky_way.empty()) {
        std::string error;
        uint32_t map_width = 0;
        m_milky_way = load_texture_srgb(m_renderer.device(), m_asset_dir / m_scene.sky.milky_way, &error, &map_width);
        if (m_milky_way) {
            m_starfield.set_milky_way(m_milky_way, map_width, static_cast<float>(m_scene.sky.milky_way_brightness));
        } else {
            SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Milky Way map unavailable: %s", error.c_str());
        }
    }

    for (const Body& body : m_scene.bodies) {
        if (body.kind == BodyKind::BlackHole) {
            update_sky_cube(); // the sky differs per scene (Milky Way or not)
            break;
        }
    }

    const bool touring = m_director.active();
    stop_tour();
    if (m_scene.view.start_now) {
        reset_to_now();
    } else {
        m_clock.t_tdb = m_scene.view.start_tdb;
    }
    m_clock.warp = m_scene.view.warp;
    m_scene.update(m_clock.t_tdb);

    m_camera = OrbitCamera{};
    m_camera.set_up_axis(m_scene.up_axis());
    const int focus = m_scene.view.focus;
    m_camera.focus(focus, m_scene);
    if (m_scene.view.distance_km > 0.0) {
        m_camera.fly_to({focus, m_scene.view.yaw_rad, m_scene.view.pitch_rad, m_scene.view.distance_km, false}, 0.0, m_scene);
    }
    if (touring) {
        start_tour();
    }
    return true;
}

void App::set_frame(int frame)
{
    if (frame == m_scene.active_frame()) {
        return;
    }
    m_scene.set_active_frame(frame);
    m_scene.update(m_clock.t_tdb);
    m_camera.set_up_axis(m_scene.up_axis());
}

void App::jump_to_event(size_t index)
{
    const SceneEvent& e = m_scene.events[index];
    stop_tour();
    m_clock.t_tdb = e.t_tdb;
    if (e.warp > 0.0) {
        m_clock.warp = e.warp;
    }
    set_frame(e.frame >= 0 ? e.frame : 0); // events without a frame use the root frame
    m_scene.update(m_clock.t_tdb);
    if (e.focus >= 0) {
        const double distance =
            e.distance_km > 0.0 ? e.distance_km : m_camera.default_distance(e.focus, m_scene);
        double yaw = -1.2;
        double pitch = 0.3;
        if (e.from_orbit_normal >= 0) {
            m_camera.angles_for_direction(m_scene.orbit_normal(e.from_orbit_normal), &yaw, &pitch);
        } else if (e.from_body >= 0) {
            const glm::dvec3 to_body = m_scene.bodies[static_cast<size_t>(e.from_body)].world_position -
                                       m_scene.bodies[static_cast<size_t>(e.focus)].world_position;
            m_camera.angles_for_direction(to_body, &yaw, &pitch);
        }
        m_camera.fly_to({e.focus, yaw, pitch, distance, false}, 2.0, m_scene);
    }
}

void App::update_sky_cube()
{
    if (!m_sky_cube) {
        create_sky_cube();
        if (!m_sky_cube) {
            return;
        }
    }
    const uint32_t size = m_sky_cube_size;

    // Match the on-screen stars: same angular size where possible, same total
    // flux (intensity x angular area). Screen pixel ~ 2 tan(fov/2) / height;
    // cube texel at a face centre = 2 / size. The Milky Way is a radiance and
    // needs no such scaling.
    int width = 0;
    int height = 0;
    SDL_GetWindowSizeInPixels(m_window.handle(), &width, &height);
    const double screen_px = 2.0 * std::tan(22.5 * kDegToRad) / std::max(height, 1);
    const double texel = 2.0 / size;
    const double size_scale = std::max(0.8, screen_px / texel);
    const double flux_scale = std::pow(screen_px / (size_scale * texel), 2.0);
    m_starfield.render_cubemap(m_sky_cube, size, m_sky_cube_format, static_cast<float>(m_star_brightness * flux_scale),
                               static_cast<float>(size_scale));
    m_black_hole.set_sky(m_sky_cube);
}

void App::create_sky_cube()
{
    SDL_GPUDevice* device = m_renderer.device();
    const SDL_GPUTextureUsageFlags usage = SDL_GPU_TEXTUREUSAGE_SAMPLER | SDL_GPU_TEXTUREUSAGE_COLOR_TARGET;
    SDL_GPUTextureFormat format = SDL_GPU_TEXTUREFORMAT_R11G11B10_UFLOAT;
    uint32_t size = 2048;
    if (!SDL_GPUTextureSupportsFormat(device, format, SDL_GPU_TEXTURETYPE_CUBE, usage)) {
        format = SDL_GPU_TEXTUREFORMAT_R16G16B16A16_FLOAT;
        size = 1024;
    }
    SDL_GPUTextureCreateInfo info = {};
    info.type = SDL_GPU_TEXTURETYPE_CUBE;
    info.format = format;
    info.usage = usage;
    info.width = size;
    info.height = size;
    info.layer_count_or_depth = 6;
    info.num_levels = 1;
    m_sky_cube = SDL_CreateGPUTexture(device, &info);
    if (!m_sky_cube) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Sky cube map creation failed: %s", SDL_GetError());
        return;
    }
    m_sky_cube_format = format;
    m_sky_cube_size = size;
    SDL_LogInfo(SDL_LOG_CATEGORY_GPU, "Sky cube map %u px", size);
}

void App::trace_black_hole(SDL_GPUCommandBuffer* cmd, uint32_t width, uint32_t height)
{
    if (!m_black_hole.ready()) {
        return;
    }
    const FrameTransform& frame = m_scene.current_transform();
    // Camera basis (display frame) from the rotation-only view matrix rows.
    const glm::dvec3 right(m_view.view[0][0], m_view.view[1][0], m_view.view[2][0]);
    const glm::dvec3 up(m_view.view[0][1], m_view.view[1][1], m_view.view[2][1]);
    const glm::dvec3 forward(-m_view.view[0][2], -m_view.view[1][2], -m_view.view[2][2]);
    const double tan_half = 1.0 / m_view.focal_y;
    const double aspect = m_view.viewport.x / m_view.viewport.y;

    for (const Body& body : m_scene.bodies) {
        if (body.kind != BodyKind::BlackHole || !body.visible) {
            continue;
        }
        const double m_km = body.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
        const double a = body.spin;

        // Hole frame: z along the spin / disk axis.
        const glm::dvec3 z = glm::normalize(frame.direction_to_display(body.pole));
        const glm::dvec3 helper = std::abs(z.z) < 0.9 ? glm::dvec3(0.0, 0.0, 1.0) : glm::dvec3(1.0, 0.0, 0.0);
        const glm::dvec3 x = glm::normalize(glm::cross(helper, z));
        const glm::dmat3 hole(x, glm::cross(z, x), z); // hole -> display
        const glm::dmat3 to_hole = glm::transpose(hole);

        const glm::dvec3 rel = m_view.position - body.world_position;
        const double distance = glm::length(rel);

        BlackHoleUniforms u;
        u.camera = glm::vec4(glm::vec3(to_hole * rel / m_km), static_cast<float>(a));
        u.right = glm::vec4(glm::vec3(to_hole * right * (tan_half * aspect)), 0.0f);
        u.up = glm::vec4(glm::vec3(to_hole * up * tan_half), 0.0f);
        u.forward = glm::vec4(glm::vec3(to_hole * forward), 0.0f);
        u.to_sky = glm::mat4(glm::mat3(frame.axes * hole));

        const double r_in = isco_radius(a);
        const bool disk = body.disk_outer_m > r_in;
        u.radii = glm::vec4(static_cast<float>(horizon_radius(a)), static_cast<float>(r_in),
                            static_cast<float>(disk ? body.disk_outer_m : 0.0), kRayIntegrationRadius);

        // Lensing region: where the weak-field deflection 4/b exceeds half a pixel.
        const double pixel = 2.0 * tan_half / m_view.viewport.y;
        const double b_max = std::max(4.0 / (0.5 * pixel), 1.2 * kRayIntegrationRadius);
        const float depth = m_view.proj[3][2] / static_cast<float>(distance); // reversed-Z: near / distance

        // Peak of the profile r^-3/4 (1 - sqrt(r_in/r))^1/4 is at r = 49/36 r_in.
        const double r_peak = 49.0 / 36.0 * r_in;
        const double profile_max = std::pow(r_peak, -0.75) * std::pow(1.0 - std::sqrt(r_in / r_peak), 0.25);
        u.params = glm::vec4(static_cast<float>(b_max), depth, disk ? static_cast<float>(body.disk_brightness) : 0.0f,
                             static_cast<float>(m_disk_time));
        u.disk = glm::vec4(static_cast<float>(body.disk_temperature_k), static_cast<float>(1.0 / profile_max),
                           kRayMaxSteps, kLensEdgeFade);
        m_black_hole.trace(cmd, width, height, u);
        return; // one hole per scene: the trace targets hold a single result
    }
}

void App::release_body_textures()
{
    for (SDL_GPUTexture* texture : m_body_textures) {
        if (texture) {
            SDL_ReleaseGPUTexture(m_renderer.device(), texture);
        }
    }
    m_body_textures.clear();
}

void App::start_tour()
{
    m_director.set_time_warp(m_clock.paused ? 0.0 : m_clock.warp);
    m_director.start(m_scene, m_camera, static_cast<uint32_t>(SDL_GetTicksNS()));
    apply_tour_warp();
    m_pointer_idle = kPointerHideSeconds; // hide the panel right away
}

void App::apply_tour_warp()
{
    const double want = m_director.shot_warp();
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

void App::stop_tour()
{
    if (m_director.active()) {
        m_director.stop(m_camera);
        apply_tour_warp();
    }
}

void App::note_activity(bool takes_camera)
{
    m_idle_time = 0.0;
    if (takes_camera) {
        stop_tour();
    }
}

bool App::panel_visible() const
{
    return m_show_ui && (!m_director.active() || m_pointer_idle < kPointerHideSeconds);
}

void App::update(double real_dt)
{
    m_clock.advance(real_dt);
    m_scene.update(m_clock.t_tdb);

    // Accretion-disk animation: follows the clock, but capped so that the inner
    // disk does not blur at high time warp; wrapped to keep float precision.
    for (const Body& body : m_scene.bodies) {
        if (body.kind == BodyKind::BlackHole && body.gm_km3_s2 > 0.0) {
            const double t_m = body.gm_km3_s2 / std::pow(kSpeedOfLightKmS, 3.0); // seconds per M
            const double rate = m_clock.paused ? 0.0 : std::min(m_clock.warp, kDiskMaxWarp);
            m_disk_time = std::fmod(m_disk_time + real_dt * rate * (m_clock.reverse ? -1.0 : 1.0) / t_m, 1.0e6);
            break;
        }
    }

    m_idle_time += real_dt;
    m_pointer_idle += real_dt;
    if (m_auto_tour && !m_director.active() && m_idle_time > kIdleSeconds) {
        start_tour();
    }
    m_director.set_time_warp(m_clock.paused ? 0.0 : m_clock.warp);
    m_director.update(real_dt, m_scene, m_camera);
    apply_tour_warp();
    m_camera.update(real_dt, m_scene);

    // During the tour, hide the cursor while the mouse is still.
    const bool hide_cursor = m_director.active() && m_pointer_idle >= kPointerHideSeconds;
    if (hide_cursor != m_cursor_hidden) {
        m_cursor_hidden = hide_cursor;
        if (hide_cursor) {
            SDL_HideCursor();
        } else {
            SDL_ShowCursor();
        }
    }

    int width = 0;
    int height = 0;
    SDL_GetWindowSizeInPixels(m_window.handle(), &width, &height);
    m_view = m_camera.view(width, height);
    // Stars are ICRF directions: rotate them into the display frame.
    m_view.sky_view_proj =
        m_view.view_proj * glm::mat4(glm::mat3(glm::transpose(m_scene.current_transform().axes)));
}

void App::build_body_items()
{
    m_body_items.clear();
    m_ring_items.clear();
    const glm::dvec3& cam = m_view.position;
    const int star = m_scene.star_index();
    const double sun_radius = star >= 0 ? m_scene.bodies[static_cast<size_t>(star)].equatorial_radius_km : kSunRadiusKm;

    // Bodies that can cast shadows.
    std::vector<int> casters;
    for (size_t k = 0; k < m_scene.bodies.size(); ++k) {
        const Body& b = m_scene.bodies[k];
        if (b.visible && b.kind == BodyKind::Planet) {
            casters.push_back(static_cast<int>(k));
        }
    }

    for (size_t i = 0; i < m_scene.bodies.size(); ++i) {
        const Body& body = m_scene.bodies[i];
        const bool drawable = body.kind == BodyKind::Planet || body.kind == BodyKind::Star ||
                              body.kind == BodyKind::BlackHole;
        if (!body.visible || !drawable) {
            continue;
        }
        const glm::vec3 radii(static_cast<float>(body.equatorial_radius_km),
                              static_cast<float>(body.equatorial_radius_km),
                              static_cast<float>(body.polar_radius_km));

        BodyDrawItem item;
        glm::dmat3 orientation = body.orientation;
        if (i < m_body_textures.size() && m_body_textures[i]) {
            // Turn the mesh so that u = 0 sits at the map's left-edge longitude
            // (east longitude; a west-positive map is mirrored via flip_u).
            item.texture = m_body_textures[i];
            item.flip_u = body.texture_west_positive;
            const double left_east = body.texture_west_positive ? -body.texture_left_lon_deg
                                                                : body.texture_left_lon_deg;
            orientation = orientation * rotation_z(left_east * kDegToRad);
        }
        item.rotation = glm::mat4(glm::mat3(orientation));
        item.model = glm::translate(glm::mat4(1.0f), to_render(body.world_position, cam)) * item.rotation *
                     glm::scale(glm::mat4(1.0f), radii);
        item.inv_scale = 1.0f / radii;
        item.color = body.color;
        item.style = static_cast<int>(body.style);
        if (body.kind == BodyKind::Star) {
            item.style = kStyleStar;
            item.color = glm::vec3(blackbody_linear_srgb(body.temperature_k)) * kStarSurface;
            item.texture = nullptr;
            m_body_items.push_back(item);
            continue;
        }
        if (body.kind == BodyKind::BlackHole) {
            if (m_black_hole.ready()) {
                continue; // ray traced (shadow, lensing, disk) by trace_black_hole()
            }
            item.style = kStyleBlackHole;
            item.texture = nullptr;
            m_body_items.push_back(item);
            continue;
        }

        // Bodies orbiting a black hole with an accretion disk are lit by the disk,
        // which outshines the scene's stars there; bodies around a star (e.g. a
        // planet of Proxima in a multiple system) by that star; the rest by the
        // scene's star / the sun.
        glm::dvec3 light = m_scene.sun_position();
        double light_radius = sun_radius;
        const Body* parent = body.parent >= 0 ? &m_scene.bodies[static_cast<size_t>(body.parent)] : nullptr;
        if (parent && parent->kind == BodyKind::BlackHole && parent->disk_outer_m > 0.0) {
            light = parent->world_position;
            light_radius = parent->disk_outer_m * parent->gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
        } else if (const int star_k = m_scene.lighting_star(static_cast<int>(i)); star_k >= 0) {
            light = m_scene.bodies[static_cast<size_t>(star_k)].world_position;
            light_radius = m_scene.bodies[static_cast<size_t>(star_k)].equatorial_radius_km;
        }
        const glm::dvec3 to_sun = light - body.world_position;
        const double sun_distance = glm::length(to_sun);
        item.sun_direction = glm::vec3(to_sun / sun_distance);
        item.sun_angular_radius = static_cast<float>(light_radius / sun_distance);

        if (!body.rings.bands.empty()) {
            // Everything in the body-fixed frame; the ring plane is its equator.
            const glm::dmat3 to_body = glm::transpose(body.orientation);
            RingDrawItem ring;
            ring.model = glm::translate(glm::mat4(1.0f), to_render(body.world_position, cam)) *
                         glm::mat4(glm::mat3(body.orientation));
            ring.sun_direction = glm::vec3(to_body * (to_sun / sun_distance));
            ring.sun_angular_radius = item.sun_angular_radius;
            ring.camera = glm::vec3(to_body * (cam - body.world_position));
            ring.color = body.rings.color;
            ring.gain = static_cast<float>(body.rings.gain);
            ring.phase_g = static_cast<float>(body.rings.phase_g);
            ring.equatorial_radius = static_cast<float>(body.equatorial_radius_km);
            ring.polar_radius = static_cast<float>(body.polar_radius_km);
            for (const RingBand& band : body.rings.bands) {
                if (ring.band_count < kMaxRingBands) {
                    ring.bands[ring.band_count++] =
                        glm::vec4(static_cast<float>(band.inner_km), static_cast<float>(band.outer_km),
                                  static_cast<float>(band.optical_depth), static_cast<float>(band.thickness_km));
                }
            }
            m_ring_items.push_back(ring);
        }

        // The nearest other bodies are the only plausible shadow casters.
        std::sort(casters.begin(), casters.end(), [&](int a, int b) {
            return glm::length(m_scene.bodies[static_cast<size_t>(a)].world_position - body.world_position) <
                   glm::length(m_scene.bodies[static_cast<size_t>(b)].world_position - body.world_position);
        });
        for (int k : casters) {
            if (k == static_cast<int>(i) || item.occluder_count >= kMaxOccluders) {
                continue;
            }
            const Body& other = m_scene.bodies[static_cast<size_t>(k)];
            item.occluders[item.occluder_count++] = glm::vec4(to_render(other.world_position, cam),
                                                              static_cast<float>(other.equatorial_radius_km));
        }
        m_body_items.push_back(item);
    }
}

void App::build_orbit_lines()
{
    m_orbits.begin();
    if (!m_show_orbits) {
        return;
    }
    const glm::dvec3& cam = m_view.position;

    for (size_t i = 0; i < m_scene.bodies.size(); ++i) {
        m_scene.trail(static_cast<int>(i), m_clock.t_tdb, kTrailPoints, m_trail_points, m_trail_fades);
        if (m_trail_points.size() < 2) {
            continue;
        }
        m_line_points.clear();
        for (size_t k = 0; k < m_trail_points.size(); ++k) {
            m_line_points.emplace_back(to_render(m_trail_points[k], cam), m_trail_fades[k]);
        }
        m_orbits.add_line(m_line_points, glm::vec4(m_scene.bodies[i].orbit_color, kOrbitOpacity));
    }

    // Lines of apsides, one per periapsis within the trail: from the periapsis
    // through the parent out to the apoapsis distance on the far side, fading
    // with age like the trail. They turn by the apsidal advance each orbit. The
    // part over the parent (or a black hole's disk) is left out.
    for (size_t i = 0; i < m_scene.bodies.size(); ++i) {
        const Body& body = m_scene.bodies[i];
        if (!body.mark_periapsides || !body.visible) {
            continue;
        }
        const Body& parent = m_scene.bodies[static_cast<size_t>(body.parent)];
        double gap = parent.equatorial_radius_km * 1.5;
        if (parent.kind == BodyKind::BlackHole) {
            const double m_km = parent.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
            gap = std::max(parent.disk_outer_m, 6.0) * m_km; // the disk, or the shadow (~5.2 M)
        }
        const double span = body.trail_history_days * kSecondsPerDay;
        const int index = static_cast<int>(i);
        const FrameTransform& transform = m_scene.current_transform();
        const glm::vec3 color = glm::mix(body.orbit_color, glm::vec3(1.0f), 0.5f);
        const std::vector<double> times = m_scene.periapsis_times(index, m_clock.t_tdb - span, m_clock.t_tdb);
        for (size_t k = 0; k < times.size(); ++k) {
            const double tp = times[k];
            // Apoapsis distance: half an anomalistic period later (or earlier, for the latest).
            const double half = 0.5 * (k + 1 < times.size() ? times[k + 1] - tp : (k > 0 ? tp - times[k - 1] : 0.0));
            if (half <= 0.0) {
                continue;
            }
            const glm::dvec3 center = m_scene.icrf_state_at(body.parent, tp).position;
            const glm::dvec3 peri = m_scene.icrf_state_at(index, tp).position - center;
            const double r_peri = glm::length(peri);
            const double r_apo = glm::length(m_scene.icrf_state_at(index, tp + half).position -
                                             m_scene.icrf_state_at(body.parent, tp + half).position);
            const glm::dvec3 dir = peri / r_peri;
            const float age = static_cast<float>((m_clock.t_tdb - tp) / span);
            // Signed distances along dir: periapsis side, then apoapsis side.
            const double segments[2][2] = {{r_peri, gap}, {-gap, -r_apo}};
            for (const auto& seg : segments) {
                if (std::max(std::abs(seg[0]), std::abs(seg[1])) <= gap) {
                    continue; // the whole side lies within the gap
                }
                m_line_points.clear();
                m_line_points.emplace_back(to_render(transform.to_display(center + dir * seg[0]), cam), age);
                m_line_points.emplace_back(to_render(transform.to_display(center + dir * seg[1]), cam), age);
                m_orbits.add_line(m_line_points, glm::vec4(color, kOrbitOpacity));
            }
        }
    }
}

void App::render()
{
    ImDrawData* draw_data = ImGui::GetDrawData();

    Frame frame;
    if (!m_renderer.begin_frame(frame)) {
        return;
    }

    if (frame.swapchain && frame.width > 0 && frame.height > 0) {
        build_body_items();
        build_orbit_lines();

        // Uploads (copy passes) must precede the render passes.
        m_orbits.upload(frame.cmd);
        ImGui_ImplSDLGPU3_PrepareDrawData(draw_data, frame.cmd);
        trace_black_hole(frame.cmd, frame.width, frame.height);

        SDL_GPURenderPass* pass = m_renderer.begin_scene_pass(frame, {0.0f, 0.0f, 0.0f, 1.0f});
        m_starfield.draw(frame.cmd, pass, m_view, m_star_brightness);

        SunLight sun;
        sun.ambient = kAmbient;
        m_bodies.draw(frame.cmd, pass, m_view, sun, m_body_items);
        m_rings.draw(frame.cmd, pass, m_view, m_ring_items);
        m_black_hole.composite(pass);

        // Stars: sprites at their real positions (depth-tested against bodies).
        // Without a star body, the sun is drawn as a direction at infinity.
        const float viewport_half = 0.5f * m_view.viewport.y * m_view.focal_y;
        bool any_star = false;
        for (const Body& body : m_scene.bodies) {
            if (body.kind != BodyKind::Star || !body.visible) {
                continue;
            }
            any_star = true;
            const glm::dvec3 to_star = body.world_position - m_view.position;
            const double distance = glm::length(to_star);
            SunPass::Star s;
            s.position = glm::vec3(to_star);
            s.finite = true;
            s.angular_radius = static_cast<float>(body.equatorial_radius_km / distance);
            s.color = glm::vec3(blackbody_linear_srgb(body.temperature_k));
            s.magnitude = static_cast<float>(apparent_magnitude(body.luminosity_solar, distance));
            s.draw_disk = s.angular_radius * viewport_half < 3.0f; // otherwise the sphere shows it
            m_sun.draw(frame.cmd, pass, m_view, s);
        }
        if (!any_star) {
            const glm::dvec3 to_sun = m_scene.sun_position() - m_view.position;
            const double distance = glm::length(to_sun);
            SunPass::Star s;
            s.position = glm::vec3(to_sun / distance);
            s.finite = false;
            s.angular_radius = static_cast<float>(kSunRadiusKm / distance);
            s.color = glm::vec3(blackbody_linear_srgb(5772.0)); // IAU 2015 B3 nominal solar Teff
            s.magnitude = static_cast<float>(apparent_magnitude(1.0, distance));
            m_sun.draw(frame.cmd, pass, m_view, s);
        }

        m_orbits.draw(frame.cmd, pass, m_view, m_line_width);
        SDL_EndGPURenderPass(pass);

        m_post.run(frame.cmd, m_renderer.hdr_texture(), frame.width, frame.height, frame.swapchain,
                   m_post_settings);

        SDL_GPURenderPass* overlay = m_renderer.begin_overlay_pass(frame);
        ImGui_ImplSDLGPU3_RenderDrawData(draw_data, frame.cmd, overlay);
        SDL_EndGPURenderPass(overlay);
    }

    m_renderer.end_frame(frame);
}

} // namespace astraxis
