#include "app/app.hpp"

#include "core/color.hpp"
#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/kerr_null.hpp"
#include "platform/paths.hpp"
#include "render/texture.hpp"
#include "scene/comet.hpp"
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
// A moon's label, marker and orbit fade out as its orbit shrinks on screen
// from kSatelliteShowPx to kSatelliteHidePx (points, not pixels).
constexpr double kSatelliteHidePx = 10.0;
constexpr double kSatelliteShowPx = 28.0;
// Belts keep their configured point brightness while their median radius spans at
// least this many points on screen.
constexpr double kBeltFullRadiusPx = 250.0;
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
// Pulsar beams: drawn as long as the camera is far (shorter, and the far end
// could pass behind the camera), and the star flashes when a beam sweeps over
// the camera.
constexpr double kBeamLengthPerDistance = 0.9;
constexpr float kBeamIntensity = 0.6f;
constexpr double kBeamFlashGain = 300.0;

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

// The two beam directions (display frame) of a pulsar after `real_time`
// seconds of running clock: the magnetic axis, inclined to the spin axis and
// turning about it once per display period. Returns false if not a pulsar.
bool pulsar_beams(const Scene& scene, const Body& body, double real_time, glm::dvec3 out[2])
{
    if (!body.pulsar.enabled || body.pulsar.spin_axis_orbit_of < 0) {
        return false;
    }
    const glm::dvec3 spin = scene.orbit_normal(body.pulsar.spin_axis_orbit_of);
    const glm::dvec3 ref = std::abs(spin.z) < 0.9 ? glm::dvec3(0.0, 0.0, 1.0) : glm::dvec3(1.0, 0.0, 0.0);
    const glm::dvec3 e1 = glm::normalize(glm::cross(spin, ref));
    const glm::dvec3 e2 = glm::cross(spin, e1);
    const double phase = kTwoPi * real_time / body.pulsar.display_period_s;
    const double tilt = body.pulsar.magnetic_inclination_deg * kDegToRad;
    out[0] = std::cos(tilt) * spin + std::sin(tilt) * (std::cos(phase) * e1 + std::sin(phase) * e2);
    out[1] = -out[0];
    return true;
}

glm::vec3 to_render(const glm::dvec3& world, const glm::dvec3& camera)
{
    return glm::vec3(world - camera);
}

// Plumes (illustrative): base width of a jet, width of an umbrella's canopy
// shell (in the ballistic parameter q, see plume.frag.hlsl), a geyser tail's
// width at its e-folding length, and how far the boxes reach (in scale
// heights / tail lengths).
constexpr double kJetBaseKm = 2.0;
constexpr double kUmbrellaShell = 0.25;
constexpr double kTailWidthKm = 6.0;
constexpr double kJetReach = 4.0;
constexpr double kTailReach = 3.0;

// Comets (illustrative brightness). The total magnitude law gives each coma's
// apparent flux, which is compressed like the catalog stars' (flux^gamma, see
// starfield_pass.cpp) by scaling the luminosity with distance: brightness
// relative to a comet of kCometReferenceMag. Seen from closer than
// kComaNearScales daughter scale lengths, the brightness stays as it is there.
// The ion tail's luminosity is a share of the coma's, raised for tails much
// larger than the coma (and lowered for smaller ones) so that their surface
// brightness relative to the coma is compressed with the same gamma: Halley's
// tail is 10^4 times the area of Giacobini-Zinner's for its light. Comae and
// tails smaller than a few pixels are widened (at the same total flux) so
// that distant comets show as soft spots rather than flickering sub-pixel
// cores.
constexpr double kCometMagnitudeGamma = 0.56;
constexpr double kCometReferenceMag = 5.0;
constexpr double kComaGain = 3.0e-5;  // radiance x solid angle of a reference comet
constexpr double kIonTailShare = 3.0; // tail luminosity / coma luminosity, for a tail of the coma's size
constexpr double kComaNearScales = 10.0;
constexpr double kComaReach = 8.0;    // box half-size, daughter scale lengths
constexpr double kIonTailReach = 3.0; // box length, e-folding lengths
constexpr double kIonTailOpening = 0.03; // tan of the tail's half-opening angle
constexpr double kComaMinPx = 2.0;    // daughter scale length on screen, at least
constexpr double kIonTailMinPx = 1.0; // tail radius on screen, at least
// Dust (illustrative): shares of the coma's luminosity (each tail grain's raised
// with its distance d from the nucleus, by (d / Ld)^(2 (1 - gamma)) beyond the
// gas coma's daughter scale length Ld: the tail's surface brightness falls off
// with distance as compressed as the star magnitudes are), the
// Henyey-Greenstein asymmetry of the scattering (forward-throwing), and the
// splats' growth: the spread of ejection speeds and sizes each grain stands
// for (its own ejection speed x age).
constexpr double kDustComaShare = 0.7;
constexpr double kDustTailShare = 1.0;
constexpr double kDustAsymmetry = 0.5;
constexpr double kDustSpread = 1.0; // splat sigma growth / (ejection speed x age)
constexpr double kDustMinPx = 1.0;  // splat sigma on screen, at least
// Splats near the camera grow to cover the screen, and thousands of them cost
// far too much fill rate. Beyond kDustThinPx (sigma on screen) a grain is kept
// with probability (kDustThinPx / sigma)^2, by its fixed random number, and
// brightened by the inverse, which keeps the light on average; between
// kDustFadePx and kDustMaxPx they fade out (a faint haze; the analytic dust
// coma still shows the cloud).
constexpr double kDustThinPx = 8.0;
constexpr double kDustFadePx = 150.0;
constexpr double kDustMaxPx = 300.0;

// Henyey-Greenstein phase function (mean 1 over the sphere); `cos_angle` between
// the incoming sunlight and the direction toward the camera (1 = forward).
double dust_phase(double cos_angle)
{
    const double g = kDustAsymmetry;
    return (1.0 - g * g) / std::pow(1.0 + g * g - 2.0 * g * cos_angle, 1.5);
}

glm::dvec3 lat_lon_direction(const glm::dvec2& lat_lon_deg)
{
    const double lat = lat_lon_deg.x * kDegToRad;
    const double lon = lat_lon_deg.y * kDegToRad;
    return glm::dvec3(std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat));
}

glm::dvec3 surface_point_toward(const Body& body, const glm::dvec3& dir)
{
    const double lat = std::asin(std::clamp(dir.z, -1.0, 1.0)) / kDegToRad;
    const double lon = std::atan2(dir.y, dir.x) / kDegToRad;
    return body_surface_point(body, lat, lon);
}

// A plume's local frame in the body-fixed frame: origin at the surface point
// in direction `dir`, z along the surface normal, x horizontal toward
// azimuth `x_azimuth_deg` (from north through east; 90 = east), y = z cross x.
struct PlumeFrame {
    glm::dvec3 origin{0.0};
    glm::dmat3 axes{1.0}; // columns x, y, z
};

PlumeFrame plume_frame(const Body& body, const glm::dvec3& dir, double x_azimuth_deg)
{
    PlumeFrame f;
    f.origin = surface_point_toward(body, dir);
    const glm::dvec3 up = body_surface_normal(body, f.origin);
    glm::dvec3 east = glm::cross(glm::dvec3(0.0, 0.0, 1.0), up);
    east = glm::length(east) > 1e-9 ? glm::normalize(east) : glm::dvec3(0.0, 1.0, 0.0);
    const glm::dvec3 north = glm::cross(up, east);
    const double a = x_azimuth_deg * kDegToRad;
    const glm::dvec3 x = std::sin(a) * east + std::cos(a) * north;
    f.axes = glm::dmat3(x, glm::cross(up, x), up);
    return f;
}

// Grows an item's box (local frame) to hold a point widened by `radius`.
void enclose(PlumeDrawItem& item, const glm::dvec3& p, double radius, bool first)
{
    const glm::vec3 lo(p - radius);
    const glm::vec3 hi(p + radius);
    item.box_min = first ? lo : glm::min(item.box_min, lo);
    item.box_max = first ? hi : glm::max(item.box_max, hi);
}

// The draw items of one plume (a geyser has two: column and tail).
// `occluders`: camera-relative centers (xyz) and radii (w).
void add_plume_items(const Body& body, const Plume& plume, const glm::dvec3& camera, const glm::dvec3& sun_body,
                     float sun_angular_radius, const glm::vec4* occluders, int occluder_count,
                     std::vector<PlumeDrawItem>& out)
{
    auto make = [&](PlumeDrawItem::Type type, const PlumeFrame& f) {
        const glm::dmat3 rot = body.orientation * f.axes; // local -> display
        const glm::dmat3 to_local = glm::transpose(rot);
        const glm::dvec3 vent = body.world_position + body.orientation * f.origin;
        PlumeDrawItem item;
        item.type = type;
        item.model = glm::translate(glm::mat4(1.0f), glm::vec3(vent - camera)) * glm::mat4(glm::mat3(rot));
        item.to_local = glm::mat3(to_local);
        item.camera_local = glm::vec3(to_local * (camera - vent));
        item.planet_center = glm::vec3(glm::transpose(f.axes) * -f.origin);
        item.planet_radius = static_cast<float>(glm::length(f.origin));
        item.sun_direction = glm::vec3(glm::transpose(f.axes) * sun_body);
        item.sun_angular_radius = sun_angular_radius;
        for (int k = 0; k < occluder_count && item.occluder_count < kMaxPlumeOccluders; ++k) {
            const glm::dvec3 center = glm::dvec3(glm::vec3(occluders[k])) + camera;
            item.occluders[item.occluder_count++] = glm::vec4(glm::vec3(to_local * (center - vent)), occluders[k].w);
        }
        item.albedo = glm::vec3(plume.albedo);
        item.g = static_cast<float>(plume.g);
        return item;
    };
    const glm::dvec3 dir = lat_lon_direction(plume.lat_lon_deg);
    const double r_planet = glm::length(surface_point_toward(body, dir));
    const double height = plume.height_km;
    const double tau = plume.optical_depth;
    switch (plume.type) {
    case Plume::Type::Umbrella: {
        PlumeDrawItem item = make(PlumeDrawItem::Umbrella, plume_frame(body, dir, 90.0));
        const double rc = plume.radius_km;
        item.shape = glm::vec4(static_cast<float>(height), static_cast<float>(rc), static_cast<float>(kUmbrellaShell), 0.0f);
        // Vertical optical depth through the canopy above the vent: alt = H (1 - q^2)
        // there, so the integral of exp(-q / w) d(alt) is 2 H w^2 (w << 1).
        item.density = static_cast<float>(tau / (2.0 * height * kUmbrellaShell * kUmbrellaShell));
        item.box_min = glm::vec3(-1.05 * rc, -1.05 * rc, -rc * rc / (2.0 * r_planet) - 0.05 * height);
        item.box_max = glm::vec3(1.05 * rc, 1.05 * rc, 1.1 * height);
        out.push_back(item);
        break;
    }
    case Plume::Type::Jets: {
        const glm::dvec3 end = lat_lon_direction(plume.end_lat_lon_deg);
        const PlumeFrame f = plume_frame(body, glm::normalize(dir + end), 90.0);
        PlumeDrawItem item = make(PlumeDrawItem::Jets, f);
        const double spread = std::tan(plume.spread_deg * kDegToRad);
        const double top = kJetReach * height;
        const double top_width = 2.0 * (kJetBaseKm + top * spread);
        item.shape = glm::vec4(static_cast<float>(height), static_cast<float>(spread), static_cast<float>(kJetBaseKm), 0.0f);
        item.density = static_cast<float>(tau / (std::sqrt(kPi) * kJetBaseKm)); // across a jet at its base
        item.jet_count = std::min(plume.count, kMaxPlumeDrawJets);
        for (int k = 0; k < item.jet_count; ++k) {
            const double t = item.jet_count > 1 ? static_cast<double>(k) / (item.jet_count - 1) : 0.5;
            const glm::dvec3 vent = surface_point_toward(body, glm::normalize(glm::mix(dir, end, t)));
            const glm::dvec3 pos = glm::transpose(f.axes) * (vent - f.origin);
            const glm::dvec3 axis = glm::transpose(f.axes) * body_surface_normal(body, vent);
            item.jet_positions[k] = glm::vec4(glm::vec3(pos), 0.0f);
            item.jet_axes[k] = glm::vec4(glm::vec3(axis), 0.0f);
            enclose(item, pos, 2.0 * kJetBaseKm, k == 0);
            enclose(item, pos + top * axis, top_width, false);
        }
        out.push_back(item);
        break;
    }
    case Plume::Type::Geyser: {
        const double rc = plume.radius_km;
        const double density = tau / (std::sqrt(kPi) * rc); // across the column
        PlumeDrawItem column = make(PlumeDrawItem::Column, plume_frame(body, dir, 90.0));
        column.shape = glm::vec4(static_cast<float>(height), static_cast<float>(rc), 0.0f, 0.0f);
        column.density = static_cast<float>(density);
        column.box_min = glm::vec3(-3.0 * rc, -3.0 * rc, -rc);
        column.box_max = glm::vec3(3.0 * rc, 3.0 * rc, height + rc);
        out.push_back(column);
        if (plume.tail_km > 0.0) {
            // Along +x, toward the tail's azimuth.
            PlumeDrawItem tail = make(PlumeDrawItem::Tail, plume_frame(body, dir, plume.tail_azimuth_deg));
            const double length = plume.tail_km;
            const double widening = std::max(kTailWidthKm - rc, 0.0) / length;
            const double reach = kTailReach * length;
            const double end_width = 3.0 * (rc + reach * widening);
            tail.shape = glm::vec4(static_cast<float>(height), static_cast<float>(rc), static_cast<float>(length),
                                   static_cast<float>(widening));
            tail.density = static_cast<float>(density);
            tail.box_min = glm::vec3(-rc, -end_width, height - end_width - reach * reach / (2.0 * r_planet));
            tail.box_max = glm::vec3(reach, end_width, height + end_width);
            out.push_back(tail);
        }
        break;
    }
    }
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

    std::vector<CatalogStar>& catalog = m_catalog;
    std::string catalog_error;
    if (!load_star_catalog(m_asset_dir / "stars" / "bsc5.csv", catalog, &catalog_error)) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Star catalog unavailable (%s); using procedural stars",
                    catalog_error.c_str());
    }

    const SceneTargetFormat& format = m_renderer.scene_format();
    if (!m_starfield.init(m_renderer.device(), format, catalog) || !m_bodies.init(m_renderer.device(), format) ||
        !m_rings.init(m_renderer.device(), format) || !m_atmospheres.init(m_renderer.device(), format) ||
        !m_plumes.init(m_renderer.device(), format) || !m_comets.init(m_renderer.device(), format) ||
        !m_dust.init(m_renderer.device(), format) ||
        !m_orbits.init(m_renderer.device(), format) || !m_belts.init(m_renderer.device(), format) ||
        !m_sun.init(m_renderer.device(), format) || !m_beams.init(m_renderer.device(), format) ||
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
    m_beams.shutdown();
    m_sun.shutdown();
    m_orbits.shutdown();
    m_belts.shutdown();
    m_rings.shutdown();
    m_atmospheres.shutdown();
    m_plumes.shutdown();
    m_comets.shutdown();
    m_dust.shutdown();
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
    case SDLK_I:
        m_show_info = !m_show_info;
        break;
    case SDLK_M:
        m_show_atmospheres = !m_show_atmospheres;
        break;
    case SDLK_P:
        m_show_plumes = !m_show_plumes;
        break;
    case SDLK_C:
        m_show_comets = !m_show_comets;
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

    m_body_fades.assign(m_scene.bodies.size(), 1.0f); // until the next update()
    m_dust_tails.assign(m_scene.bodies.size(), DustTail{});
    m_belts.clear();
    m_belt_ids.clear();
    m_belt_reference_h.clear();
    m_belt_radius_km.clear();
    for (const SceneBelt& belt : m_scene.belts) {
        m_belt_ids.push_back(m_belts.add_belt(belt.elements.data(), static_cast<uint32_t>(belt.size())));
        // Weight 1 at the median H: the typical point keeps the configured brightness.
        std::vector<float> h;
        h.reserve(belt.size());
        for (size_t k = 0; k < belt.size(); ++k) {
            h.push_back(belt.elements[k * SceneBelt::kStride + 6]);
        }
        std::nth_element(h.begin(), h.begin() + static_cast<std::ptrdiff_t>(h.size() / 2), h.end());
        m_belt_reference_h.push_back(h.empty() ? 0.0f : h[h.size() / 2]);
        std::vector<float> a;
        a.reserve(belt.size());
        for (size_t k = 0; k < belt.size(); ++k) {
            a.push_back(belt.elements[k * SceneBelt::kStride]);
        }
        std::nth_element(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(a.size() / 2), a.end());
        m_belt_radius_km.push_back(a.empty() ? 0.0 : a[a.size() / 2] * kAuKm);
    }
    m_label_alpha.clear();
    m_body_textures.assign(m_scene.bodies.size(), nullptr);
    m_ring_textures.assign(m_scene.bodies.size(), nullptr);
    m_bodies.clear_meshes();
    m_body_meshes.assign(m_scene.bodies.size(), -1);
    for (size_t i = 0; i < m_scene.bodies.size(); ++i) {
        const Body& body = m_scene.bodies[i];
        if (body.shape) {
            m_body_meshes[i] = m_bodies.add_mesh(body.shape->positions, body.shape->normals, body.shape->albedo,
                                                 body.shape->map_u, body.shape->indices);
        }
        if (!body.rings.profile.empty()) {
            m_ring_textures[i] = create_profile_texture(m_renderer.device(), body.rings.profile);
        }
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

    // A scene may bring its own stars (e.g. a globular cluster seen from inside).
    if (!m_scene.sky.stars.empty()) {
        m_starfield.set_stars(m_scene.sky.stars);
        m_scene_sky_stars = true;
    } else if (m_scene_sky_stars) {
        m_starfield.set_stars(m_catalog);
        m_scene_sky_stars = false;
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
            yaw += e.phase_deg * kDegToRad;
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
    for (auto* textures : {&m_body_textures, &m_ring_textures}) {
        for (SDL_GPUTexture* texture : *textures) {
            if (texture) {
                SDL_ReleaseGPUTexture(m_renderer.device(), texture);
            }
        }
        textures->clear();
    }
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

    if (!m_clock.paused) {
        m_pulsar_time = std::fmod(m_pulsar_time + real_dt, 1.0e6);
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

    const double px_per_radian = 0.5 * m_view.viewport.y * m_view.focal_y;
    const double scale = m_window.content_scale();
    satellite_fades(m_scene, m_view.position, px_per_radian, kSatelliteHidePx * scale, kSatelliteShowPx * scale,
                    m_body_fades);
    const int target = m_camera.target();
    if (target >= 0 && target < static_cast<int>(m_body_fades.size())) {
        m_body_fades[static_cast<size_t>(target)] = 1.0f; // the focus always shows
    }
}

void App::build_body_items()
{
    m_body_items.clear();
    m_ring_items.clear();
    m_atmosphere_items.clear();
    m_atmosphere_optics.resize(m_scene.bodies.size());
    m_plume_items.clear();
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
        // An atmosphere with an opaque deck shows the deck (in the body's
        // color) instead of the surface.
        const Atmosphere& atmosphere = body.atmosphere;
        const bool show_atmosphere = m_show_atmospheres && atmosphere.enabled;
        const double deck_km = show_atmosphere ? atmosphere.deck_altitude_km : 0.0;
        const glm::dvec3 radii_km = glm::dvec3(body.equatorial_radius_km, body.equatorial_radius_b_km,
                                               body.polar_radius_km) + deck_km;
        const glm::vec3 radii(radii_km);

        BodyDrawItem item;
        if (i < m_body_textures.size() && m_body_textures[i] && deck_km <= 0.0) {
            // East longitude of the map's left edge (a west-positive map is mirrored via flip_u).
            item.texture = m_body_textures[i];
            item.flip_u = body.texture_west_positive;
            item.texture_left_lon_deg = static_cast<float>(body.texture_west_positive ? -body.texture_left_lon_deg
                                                                                       : body.texture_left_lon_deg);
        }
        // A shape model is in km already; the ellipsoid is a scaled unit sphere.
        item.mesh = i < m_body_meshes.size() ? m_body_meshes[i] : -1;
        const glm::vec3 scale = item.mesh >= 0 ? glm::vec3(1.0f) : radii;
        item.rotation = glm::mat4(glm::mat3(body.orientation));
        item.model = glm::translate(glm::mat4(1.0f), to_render(body.world_position, cam)) * item.rotation *
                     glm::scale(glm::mat4(1.0f), scale);
        item.inv_scale = 1.0f / scale;
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
        // which outshines the scene's stars there; the rest by the brightest
        // star where they are (e.g. Proxima for its planets in the alpha Cen
        // scene) and a second one if it adds light (the other sun of a
        // circumbinary planet), or else by the sun.
        glm::dvec3 light = m_scene.sun_position();
        double light_radius = sun_radius;
        const Body* parent = body.parent >= 0 ? &m_scene.bodies[static_cast<size_t>(body.parent)] : nullptr;
        StarLight stars[kMaxStarLights];
        int star_count = 0;
        if (parent && parent->kind == BodyKind::BlackHole && parent->disk_outer_m > 0.0) {
            light = parent->world_position;
            light_radius = parent->disk_outer_m * parent->gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
        } else if (star_count = m_scene.lighting_stars(static_cast<int>(i), stars, kMaxStarLights); star_count > 0) {
            const Body& brightest = m_scene.bodies[static_cast<size_t>(stars[0].star)];
            light = brightest.world_position;
            light_radius = brightest.equatorial_radius_km;
        }
        const glm::dvec3 to_sun = light - body.world_position;
        const double sun_distance = glm::length(to_sun);
        item.sun_direction = glm::vec3(to_sun / sun_distance);
        item.sun_angular_radius = static_cast<float>(light_radius / sun_distance);
        if (star_count > 1) {
            // Colors relative to the brightest star, which looks white (the eye
            // adapts to it); the luminance follows the bolometric flux.
            const Body& first = m_scene.bodies[static_cast<size_t>(stars[0].star)];
            const Body& second = m_scene.bodies[static_cast<size_t>(stars[1].star)];
            const glm::dvec3 to_second = second.world_position - body.world_position;
            const double second_distance = glm::length(to_second);
            const glm::dvec3 tint = blackbody_linear_srgb(second.temperature_k) /
                                     glm::max(blackbody_linear_srgb(first.temperature_k), glm::dvec3(1e-3));
            const double luminance = glm::dot(tint, glm::dvec3(0.2126, 0.7152, 0.0722));
            item.light2_direction = glm::vec3(to_second / second_distance);
            item.light2_angular_radius = static_cast<float>(second.equatorial_radius_km / second_distance);
            item.light2_color = glm::vec3(tint * (stars[1].relative_flux / luminance));
        }

        if (show_atmosphere) {
            // In units of the base radius (the shell's unit sphere).
            const double base_km = radii_km.x;
            AtmosphereOptics& optics = m_atmosphere_optics[i];
            optics.rayleigh_depth = glm::vec3(atmosphere.rayleigh_depth);
            optics.rayleigh_scale_height = static_cast<float>(atmosphere.rayleigh_scale_height_km / base_km);
            optics.haze_attenuation =
                glm::vec3(atmosphere.haze_depth * (1.0 - atmosphere.haze_albedo * atmosphere.haze_g));
            optics.haze_scale_height = static_cast<float>(atmosphere.haze_scale_height_km / base_km);
            optics.haze_scattering = glm::vec3(atmosphere.haze_depth * atmosphere.haze_albedo);
            optics.haze_g = static_cast<float>(atmosphere.haze_g);
            item.atmosphere = &optics;

            const glm::dmat3 linear = body.orientation * glm::dmat3(glm::dvec3(radii_km.x, 0.0, 0.0),
                                                                    glm::dvec3(0.0, radii_km.y, 0.0),
                                                                    glm::dvec3(0.0, 0.0, radii_km.z));
            const glm::dmat3 to_unit = glm::inverse(linear);
            AtmosphereDrawItem shell;
            shell.model = glm::translate(glm::mat4(1.0f), to_render(body.world_position, cam)) *
                          glm::mat4(glm::mat3(linear));
            shell.to_unit = glm::mat3(to_unit);
            shell.camera_unit = glm::vec3(to_unit * (cam - body.world_position));
            shell.top = static_cast<float>(1.0 + atmosphere.height_km / base_km);
            shell.sun_direction = item.sun_direction;
            shell.sun_angular_radius = item.sun_angular_radius;
            shell.optics = optics;
            shell.gain = static_cast<float>(atmosphere.gain);
            m_atmosphere_items.push_back(shell);
        }

        SDL_GPUTexture* ring_texture = i < m_ring_textures.size() ? m_ring_textures[i] : nullptr;
        if (ring_texture) {
            // The rings shade the planet...
            item.ring_profile = ring_texture;
            item.ring_center = to_render(body.world_position, cam);
            item.ring_normal = glm::vec3(body.orientation[2]);
            item.ring_inner_km = static_cast<float>(body.rings.inner_km);
            item.ring_outer_km = static_cast<float>(body.rings.outer_km);
            item.ring_samples = static_cast<float>(body.rings.profile.size());

            // ...and are drawn in the body-fixed frame, where their plane is the equator.
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
            ring.inner_km = static_cast<float>(body.rings.inner_km);
            ring.outer_km = static_cast<float>(body.rings.outer_km);
            ring.profile = ring_texture;
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
        if (m_show_plumes) {
            const glm::dvec3 sun_body = glm::transpose(body.orientation) * glm::dvec3(item.sun_direction);
            for (const Plume& plume : body.plumes) {
                add_plume_items(body, plume, cam, sun_body, item.sun_angular_radius, item.occluders,
                                item.occluder_count, m_plume_items);
            }
        }
        m_body_items.push_back(item);
    }
}

void App::build_belt_items()
{
    m_belt_items.clear();
    if (!m_show_belts) {
        return;
    }
    // J2000 ecliptic -> ICRF (SBDB elements use the IAU 1976 obliquity, 84381.448") -> display.
    const double eps = 84381.448 / 3600.0 * kDegToRad;
    const glm::dmat3 icrf_from_ecliptic(glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(0.0, std::cos(eps), std::sin(eps)),
                                        glm::dvec3(0.0, -std::sin(eps), std::cos(eps)));
    const glm::dmat3 display_from_ecliptic = glm::transpose(m_scene.current_transform().axes) * icrf_from_ecliptic;
    const glm::dvec3 sun = m_scene.sun_position() - m_view.position;
    const double px_per_radian = 0.5 * m_view.viewport.y * m_view.focal_y;
    for (size_t k = 0; k < m_scene.belts.size(); ++k) {
        const SceneBelt& belt = m_scene.belts[k];
        BeltDrawItem item;
        item.belt = k < m_belt_ids.size() ? m_belt_ids[k] : -1;
        item.sun_relative = glm::vec3(sun);
        item.display_from_ecliptic = glm::mat3(display_from_ecliptic);
        item.days_since_epoch = static_cast<float>((m_clock.t_tdb - belt.epoch_tdb) / kSecondsPerDay);
        item.point_size_px = static_cast<float>(belt.point_size_px * m_window.content_scale());
        item.color = belt.color;
        // Once the belt shrinks below kBeltFullRadiusPx on screen, dim each point with
        // the belt's apparent area: its surface brightness stays the same instead of
        // tens of thousands of points piling up into a glare around the Sun.
        const double radius_px = k < m_belt_radius_km.size() ? m_belt_radius_km[k] / glm::length(sun) * px_per_radian : 0.0;
        const double area = std::min(1.0, radius_px / (kBeltFullRadiusPx * m_window.content_scale()));
        item.brightness = static_cast<float>(belt.brightness * area * area);
        item.reference_h = k < m_belt_reference_h.size() ? m_belt_reference_h[k] : 0.0f;
        m_belt_items.push_back(item);
    }
}

void App::build_comet_items()
{
    m_comet_items.clear();
    std::vector<DustSplat>& splats = m_dust.splats();
    splats.clear();
    if (!m_show_comets) {
        return;
    }
    const glm::dvec3& cam = m_view.position;
    const double px_per_radian = 0.5 * m_view.viewport.y * m_view.focal_y;
    const FrameTransform& transform = m_scene.current_transform();
    const State sun = m_scene.sun_icrf_state_at(m_clock.t_tdb);
    for (size_t i = 0; i < m_scene.bodies.size(); ++i) {
        const Body& body = m_scene.bodies[i];
        const Body::Comet& comet = body.comet;
        if (!comet.enabled || !body.visible) {
            continue;
        }
        const State state = m_scene.icrf_state_at(static_cast<int>(i), m_clock.t_tdb);
        const glm::dvec3 helio = state.position - sun.position;
        const double r_au = glm::length(helio) / kAuKm;
        const glm::dvec3 tail_axis = glm::normalize(
            transform.direction_to_display(ion_tail_direction(helio, state.velocity - sun.velocity, comet.solar_wind_km_s)));

        // Local frame: x along the tail.
        const glm::dvec3 helper = std::abs(tail_axis.z) < 0.9 ? glm::dvec3(0.0, 0.0, 1.0) : glm::dvec3(1.0, 0.0, 0.0);
        const glm::dvec3 y = glm::normalize(glm::cross(helper, tail_axis));
        const glm::dmat3 rot(tail_axis, y, glm::cross(tail_axis, y)); // local -> display
        const glm::dmat3 to_local = glm::transpose(rot);
        const glm::dvec3 rel = cam - body.world_position;
        const double distance = glm::length(rel);

        // Scale lengths grow as r_h^2 (comet.hpp).
        double lp = kC2ParentScaleKm * r_au * r_au;
        double ld = kC2DaughterScaleKm * r_au * r_au;
        const double widen = std::max(1.0, kComaMinPx / (ld / distance * px_per_radian));
        lp *= widen;
        ld *= widen;

        const double seen_from_km = std::max(distance, kComaNearScales * ld);
        const double magnitude =
            comet_total_magnitude(comet.m1, comet.k1, r_au) + 5.0 * std::log10(seen_from_km / kAuKm);
        const double flux = std::pow(10.0, -0.4 * kCometMagnitudeGamma * (magnitude - kCometReferenceMag));
        const double luminosity = kComaGain * flux * seen_from_km * seen_from_km;

        CometDrawItem base;
        base.model = glm::translate(glm::mat4(1.0f), to_render(body.world_position, cam)) * glm::mat4(glm::mat3(rot));
        base.to_local = glm::mat3(to_local);
        base.camera_local = glm::vec3(to_local * rel);
        // Rays stop at the nucleus and the bodies nearest the camera.
        const double nucleus = std::cbrt(body.equatorial_radius_km * body.equatorial_radius_b_km * body.polar_radius_km);
        base.occluders[base.occluder_count++] = glm::vec4(0.0f, 0.0f, 0.0f, static_cast<float>(nucleus));
        std::vector<std::pair<double, size_t>> near;
        for (size_t k = 0; k < m_scene.bodies.size(); ++k) {
            const Body& other = m_scene.bodies[k];
            if (k != i && other.visible && other.kind == BodyKind::Planet) {
                near.emplace_back(glm::length(other.world_position - cam), k);
            }
        }
        std::sort(near.begin(), near.end());
        for (const auto& entry : near) {
            if (base.occluder_count >= kMaxCometOccluders) {
                break;
            }
            const Body& other = m_scene.bodies[entry.second];
            base.occluders[base.occluder_count++] =
                glm::vec4(glm::vec3(to_local * (other.world_position - body.world_position)),
                          static_cast<float>(other.equatorial_radius_km));
        }

        CometDrawItem coma = base;
        coma.type = CometDrawItem::Coma;
        coma.shape = glm::vec3(static_cast<float>(lp), static_cast<float>(ld), static_cast<float>(nucleus));
        coma.color = comet.coma_color;
        coma.brightness = static_cast<float>(luminosity);
        coma.box_min = glm::vec3(static_cast<float>(-kComaReach * ld));
        coma.box_max = glm::vec3(static_cast<float>(kComaReach * ld));
        m_comet_items.push_back(coma);

        const double length = comet.ion_tail_length_km;
        const double tail_area = 0.5 * comet.ion_tail_diameter_km * length / (ld * ld / (widen * widen));
        const double tail_share = kIonTailShare * std::pow(tail_area, 1.0 - kCometMagnitudeGamma);
        const double w0 = std::max(0.5 * comet.ion_tail_diameter_km, kIonTailMinPx * distance / px_per_radian);
        const double reach = kIonTailReach * length;
        const double end_width = 3.0 * (w0 + reach * kIonTailOpening);
        CometDrawItem tail = base;
        tail.type = CometDrawItem::IonTail;
        tail.shape = glm::vec3(static_cast<float>(w0), static_cast<float>(kIonTailOpening), static_cast<float>(length));
        tail.color = comet.ion_tail_color;
        tail.brightness = static_cast<float>(tail_share * luminosity);
        tail.box_min = glm::vec3(static_cast<float>(-2.0 * w0), static_cast<float>(-end_width), static_cast<float>(-end_width));
        tail.box_max = glm::vec3(static_cast<float>(reach), static_cast<float>(end_width), static_cast<float>(end_width));
        m_comet_items.push_back(tail);

        // Dust coma: grains ejected sunward at v0 sqrt(beta / r_h) turn back at
        // v^2 / (2 beta g) = v0^2 r_h au^2 / (2 GM), whatever beta (fountain model).
        const double gm_sun = m_scene.bodies[static_cast<size_t>(body.parent)].gm_km3_s2;
        const double v0 = comet.dust_speed_km_s;
        const double dust_scale = std::max(v0 * v0 * r_au * kAuKm * kAuKm / (2.0 * gm_sun),
                                           kComaMinPx * distance / px_per_radian);
        const glm::dvec3 sun_display = m_scene.sun_position();
        const double cos_nucleus =
            glm::dot(glm::normalize(body.world_position - sun_display), glm::normalize(cam - body.world_position));
        CometDrawItem dust_coma = base;
        dust_coma.type = CometDrawItem::Coma;
        dust_coma.shape = glm::vec3(static_cast<float>(1e-3 * dust_scale), static_cast<float>(dust_scale),
                                    static_cast<float>(nucleus));
        dust_coma.color = comet.dust_color;
        dust_coma.brightness = static_cast<float>(kDustComaShare * luminosity * dust_phase(cos_nucleus));
        dust_coma.box_min = glm::vec3(static_cast<float>(-kComaReach * dust_scale));
        dust_coma.box_max = glm::vec3(static_cast<float>(kComaReach * dust_scale));
        m_comet_items.push_back(dust_coma);

        // Dust tail: Finson-Probstein grains (comet.hpp).
        DustTail& dust = m_dust_tails[i];
        dust.update(m_scene, static_cast<int>(i), m_clock.t_tdb);
        const double ld_physical = ld / widen;
        const double tail_luminosity = kDustTailShare * luminosity;
        for (const DustGrain& g : dust.grains()) {
            const glm::dvec3 d = g.position - state.position;
            const double boost =
                std::pow(std::max(glm::dot(d, d) / (ld_physical * ld_physical), 1.0), 1.0 - kCometMagnitudeGamma);
            const glm::dvec3 world = transform.to_display(g.position);
            const double cos_angle = glm::dot(glm::normalize(world - sun_display), glm::normalize(cam - world));
            const double sigma = dust_scale + kDustSpread * g.ejection_km_s * g.age_s;
            const double sigma_px = sigma / glm::length(world - cam) * px_per_radian;
            if (sigma_px >= kDustMaxPx) {
                continue;
            }
            const double keep = std::min(1.0, kDustThinPx * kDustThinPx / (sigma_px * sigma_px));
            if (g.random >= keep) {
                continue;
            }
            const double near_fade =
                (1.0 - std::clamp((sigma_px - kDustFadePx) / (kDustMaxPx - kDustFadePx), 0.0, 1.0)) / keep;
            DustSplat splat;
            splat.position = to_render(world, cam);
            splat.sigma_km = static_cast<float>(sigma);
            splat.color = comet.dust_color;
            splat.luminosity =
                static_cast<float>(tail_luminosity * boost * near_fade * g.weight * dust_phase(cos_angle));
            splats.push_back(splat);
        }
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
        // An orbit fades with its body; a travelled path (e.g. a spacecraft's
        // tour) stays while its host system is shown.
        const int host = m_scene.satellite_host(static_cast<int>(i));
        const float fade = m_scene.bodies[i].trail == TrailMode::History
                               ? (host >= 0 ? m_body_fades[static_cast<size_t>(host)] : 1.0f)
                               : m_body_fades[i];
        if (fade <= 0.0f) {
            continue;
        }
        m_scene.trail(static_cast<int>(i), m_clock.t_tdb, kTrailPoints, m_trail_points, m_trail_fades);
        if (m_trail_points.size() < 2) {
            continue;
        }
        m_line_points.clear();
        for (size_t k = 0; k < m_trail_points.size(); ++k) {
            m_line_points.emplace_back(to_render(m_trail_points[k], cam), m_trail_fades[k]);
        }
        m_orbits.add_line(m_line_points, glm::vec4(m_scene.bodies[i].orbit_color, kOrbitOpacity * fade));
    }

    // Lines of apsides, one per periapsis within the trail: from the periapsis
    // through the parent out to the apoapsis distance on the far side, fading
    // with age like the trail. They turn by the apsidal advance each orbit. The
    // part over the parent (or a black hole's disk) is left out.
    for (size_t i = 0; i < m_scene.bodies.size(); ++i) {
        const Body& body = m_scene.bodies[i];
        if (!body.mark_periapsides || !body.visible || m_body_fades[i] <= 0.0f) {
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
                m_orbits.add_line(m_line_points, glm::vec4(color, kOrbitOpacity * m_body_fades[i]));
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
        build_comet_items();

        // Uploads (copy passes) must precede the render passes.
        m_orbits.upload(frame.cmd);
        m_dust.upload(frame.cmd);
        ImGui_ImplSDLGPU3_PrepareDrawData(draw_data, frame.cmd);
        trace_black_hole(frame.cmd, frame.width, frame.height);

        SDL_GPURenderPass* pass = m_renderer.begin_scene_pass(frame, {0.0f, 0.0f, 0.0f, 1.0f});
        m_starfield.draw(frame.cmd, pass, m_view, m_star_brightness);

        SunLight sun;
        sun.ambient = kAmbient;
        m_bodies.draw(frame.cmd, pass, m_view, sun, m_body_items);
        m_atmospheres.draw(frame.cmd, pass, m_view, m_atmosphere_items);
        m_plumes.draw(frame.cmd, pass, m_view, static_cast<float>(m_pulsar_time), m_plume_items);
        m_comets.draw(frame.cmd, pass, m_view, static_cast<float>(m_pulsar_time), m_comet_items);
        m_dust.draw(frame.cmd, pass, m_view, static_cast<float>(kDustMinPx));
        m_rings.draw(frame.cmd, pass, m_view, m_ring_items);
        build_belt_items();
        m_belts.draw(frame.cmd, pass, m_view, m_belt_items);
        m_black_hole.composite(pass);

        // Stars: sprites at their real positions (depth-tested against bodies).
        // Without a star body, the sun is drawn as a direction at infinity.
        const float viewport_half = 0.5f * m_view.viewport.y * m_view.focal_y;
        bool any_star = false;
        m_beam_items.clear();
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
            glm::dvec3 beams[2];
            if (pulsar_beams(m_scene, body, m_pulsar_time, beams)) {
                const glm::dvec3 to_camera = -to_star / distance;
                const double half_angle = body.pulsar.beam_half_angle_deg * kDegToRad;
                double flash = 0.0;
                for (const glm::dvec3& axis : beams) {
                    const double off = std::acos(std::clamp(glm::dot(axis, to_camera), -1.0, 1.0)) / half_angle;
                    flash += std::exp(-off * off);
                    BeamDrawItem beam;
                    beam.apex = glm::vec3(to_star);
                    beam.axis = glm::vec3(axis);
                    beam.length_km = static_cast<float>(kBeamLengthPerDistance * distance);
                    beam.half_angle_rad = static_cast<float>(half_angle);
                    beam.color = body.pulsar.beam_color;
                    beam.intensity = kBeamIntensity;
                    m_beam_items.push_back(beam);
                }
                s.magnitude -= static_cast<float>(2.5 * std::log10(1.0 + kBeamFlashGain * flash));
            }
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

        m_beams.draw(frame.cmd, pass, m_view, m_beam_items);
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
