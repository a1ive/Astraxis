#include "app/app.hpp"

#include "app/label_overlay.hpp"
#include "app/scene_list.hpp"
#include "platform/paths.hpp"

#include <SDL3/SDL.h>
#include <imgui.h>
#include <imgui_impl_sdl3.h>
#include <imgui_impl_sdlgpu3.h>

#include <algorithm>

namespace astraxis {

namespace {

// Upper bound on one frame's real time step, so that dragging the window or
// hitting a breakpoint does not cause a large jump in the simulation.
constexpr double kMaxRealDt = 0.1;

constexpr double kMouseRotateSpeed = 0.005; // radians per pixel
constexpr double kIdleSeconds = 60.0;       // auto tour starts after this long without input
constexpr double kPointerHideSeconds = 3.0; // during the tour, hide panel and cursor after this
constexpr double kMaxWarp = 1e8; // ~3.2 years per second

#ifdef NDEBUG
constexpr bool kGpuDebug = false;
#else
constexpr bool kGpuDebug = true;
#endif

} // namespace

bool App::init(const LaunchOptions& options)
{
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_Init failed: %s", SDL_GetError());
        return false;
    }

    m_config_path = options.config.empty() ? config_path() : options.config;
    std::string config_error;
    m_config_ok = load_settings(m_config_path, SettingsSection::Window, m_settings, &config_error);
    if (!m_config_ok) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Settings not read (and will not be saved): %s",
                    config_error.c_str());
    }
    m_view_options = m_settings.view;
    m_show_labels = m_settings.labels;
    m_show_info = m_settings.info;
    m_auto_tour = m_settings.auto_tour;
    m_fps_limit = options.fps.value_or(m_settings.fps);

    m_start_display = options.display ? find_display(*options.display, "")
                                      : find_display(m_settings.display, m_settings.display_name);
    m_start_fullscreen = options.fullscreen.value_or(m_settings.fullscreen);
    if (!m_window.create("Astraxis", 1280, 720, m_start_display, m_start_fullscreen)) {
        return false;
    }
    if (!m_gpu.init(kGpuDebug) || !m_output.init(m_gpu, m_window.handle())) {
        return false;
    }
    m_asset_dir = asset_directory();
    SDL_LogInfo(SDL_LOG_CATEGORY_APPLICATION, "Assets: %s", m_asset_dir.string().c_str());

    if (!m_scene_renderer.init(m_gpu, m_output.swapchain_format(), m_asset_dir)) {
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
    m_title_font = load_fonts(m_asset_dir);

    ImGui_ImplSDL3_InitForSDLGPU(m_window.handle());
    ImGui_ImplSDLGPU3_InitInfo init_info = {};
    init_info.Device = m_gpu.device();
    init_info.ColorTargetFormat = m_output.swapchain_format();
    init_info.MSAASamples = SDL_GPU_SAMPLECOUNT_1;
    ImGui_ImplSDLGPU3_Init(&init_info);
    m_imgui_ready = true;

    m_scene_files = list_scene_files(m_asset_dir);
    const std::vector<SceneEntry> entries = list_scenes(m_asset_dir);
    for (const std::filesystem::path& file : m_scene_files) {
        const std::string stem = file.stem().string();
        const auto entry = std::find_if(entries.begin(), entries.end(),
                                        [&](const SceneEntry& e) { return e.stem == stem; });
        m_scene_names.push_back(entry != entries.end() ? entry->name : stem);
    }
    const std::string scene = options.scene.value_or(m_settings.scene);
    const int found = find_scene(m_scene_files, scene);
    if (found >= 0) {
        m_scene_index = static_cast<size_t>(found);
    } else {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Scene '%s' not found", scene.c_str());
    }

    if (m_scene_files.empty() || !load_scene(m_scene_index)) {
        const std::string message = m_scene_files.empty()
                                        ? "No scenes found in " + (m_asset_dir / "scenes").string()
                                        : m_scene_error;
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "%s", message.c_str());
        SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Astraxis", message.c_str(), m_window.handle());
        return false;
    }
    if (options.event > 0) {
        if (static_cast<size_t>(options.event) <= m_sim.scene().events.size()) {
            m_sim.jump_to_event(static_cast<size_t>(options.event - 1));
        } else {
            SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Scene has no event %d", options.event);
        }
    }

    m_start_scene = m_scene_files[m_scene_index].stem().string();
    m_last_counter = SDL_GetPerformanceCounter();
    m_frame_start_ns = SDL_GetTicksNS();
    m_running = true;
    m_initialized = true;
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

        if (m_fps_limit > 0) {
            const uint64_t period = SDL_NS_PER_SECOND / static_cast<uint64_t>(m_fps_limit);
            const uint64_t elapsed = SDL_GetTicksNS() - m_frame_start_ns;
            if (elapsed < period) {
                SDL_DelayPrecise(period - elapsed);
            }
        }
        m_frame_start_ns = SDL_GetTicksNS();
    }
}

void App::shutdown()
{
    if (m_initialized) {
        save_settings();
        m_initialized = false;
    }
    if (m_gpu.device()) {
        SDL_WaitForGPUIdle(m_gpu.device());
    }
    if (m_imgui_ready) {
        ImGui_ImplSDL3_Shutdown();
        ImGui_ImplSDLGPU3_Shutdown();
        ImGui::DestroyContext();
        m_imgui_ready = false;
    }
    m_scene_renderer.shutdown();
    m_output.shutdown();
    m_gpu.shutdown();
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
                m_sim.camera().rotate(-event.motion.xrel * kMouseRotateSpeed, event.motion.yrel * kMouseRotateSpeed);
                note_activity(true);
            } else {
                note_activity(false);
            }
            break;
        case SDL_EVENT_MOUSE_WHEEL:
            if (!io.WantCaptureMouse) {
                m_sim.camera().zoom(event.wheel.y);
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

    SimClock& clock = m_sim.clock();
    switch (key) {
    case SDLK_A:
        if (m_sim.touring()) {
            m_sim.stop_tour();
        } else {
            start_tour();
        }
        break;
    case SDLK_SPACE:
        clock.paused = !clock.paused;
        break;
    case SDLK_R:
        clock.reverse = !clock.reverse;
        break;
    case SDLK_N:
        m_sim.reset_to_now();
        break;
    case SDLK_F11:
        m_window.set_fullscreen(!m_window.is_fullscreen());
        break;
    case SDLK_H:
    case SDLK_F1:
        m_show_ui = !m_show_ui;
        break;
    case SDLK_O:
        m_view_options.orbits = !m_view_options.orbits;
        break;
    case SDLK_L:
        m_show_labels = !m_show_labels;
        break;
    case SDLK_I:
        m_show_info = !m_show_info;
        break;
    case SDLK_M:
        m_view_options.atmospheres = !m_view_options.atmospheres;
        break;
    case SDLK_P:
        m_view_options.plumes = !m_view_options.plumes;
        break;
    case SDLK_C:
        m_view_options.comets = !m_view_options.comets;
        break;
    case SDLK_LEFTBRACKET:
        clock.warp = std::max(1.0, clock.warp / 2.0);
        break;
    case SDLK_RIGHTBRACKET:
        clock.warp = std::min(kMaxWarp, clock.warp * 2.0);
        break;
    default:
        if (key >= SDLK_1 && key <= SDLK_9) {
            m_sim.set_focus(static_cast<int>(key - SDLK_1));
        }
        break;
    }
}

bool App::load_scene(size_t index)
{
    const bool touring = m_sim.touring();
    if (!m_sim.load_scene(m_scene_files[index], &m_scene_error)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "%s", m_scene_error.c_str());
        return false;
    }
    m_scene_error.clear();
    m_scene_index = index;
    m_labels.reset();
    if (touring) {
        m_pointer_idle = kPointerHideSeconds; // the tour goes on: hide the panel right away
    }

    int width = 0;
    int height = 0;
    SDL_GetWindowSizeInPixels(m_window.handle(), &width, &height);
    // The frame drawn next must not see the previous scene's bodies.
    m_sim.compute_view(width, height, m_window.content_scale(), m_view);
    m_scene_renderer.load_scene(m_sim.scene(), m_view_options, static_cast<uint32_t>(std::max(height, 0)));
    return true;
}

void App::save_settings()
{
    if (!m_config_ok) {
        return;
    }
    Settings settings = m_settings;
    settings.view = m_view_options;
    settings.labels = m_show_labels;
    settings.info = m_show_info;
    settings.auto_tour = m_auto_tour;
    // The command line's overrides stay out of the file unless changed in the app.
    const std::string scene = m_scene_files[m_scene_index].stem().string();
    if (scene != m_start_scene) {
        settings.scene = scene;
    }
    const uint32_t display = m_window.display();
    if (display != m_start_display) {
        settings.display = display_index(display);
        settings.display_name = display_name(display);
    }
    if (m_window.is_fullscreen() != m_start_fullscreen) {
        settings.fullscreen = m_window.is_fullscreen();
    }

    std::string error;
    if (!astraxis::save_settings(m_config_path, SettingsSection::Window, settings, &error)) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Settings not saved: %s", error.c_str());
    }
}

void App::start_tour()
{
    m_sim.start_tour();
    m_pointer_idle = kPointerHideSeconds; // hide the panel right away
}

void App::note_activity(bool takes_camera)
{
    m_idle_time = 0.0;
    if (takes_camera) {
        m_sim.stop_tour();
    }
}

bool App::panel_visible() const
{
    return m_show_ui && (!m_sim.touring() || m_pointer_idle < kPointerHideSeconds);
}

void App::update(double real_dt)
{
    m_idle_time += real_dt;
    m_pointer_idle += real_dt;
    if (m_auto_tour && !m_sim.touring() && m_idle_time > kIdleSeconds) {
        start_tour();
    }
    m_sim.update(real_dt);

    // During the tour, hide the cursor while the mouse is still.
    const bool hide_cursor = m_sim.touring() && m_pointer_idle >= kPointerHideSeconds;
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
    m_sim.compute_view(width, height, m_window.content_scale(), m_view);
}

void App::render()
{
    ImDrawData* draw_data = ImGui::GetDrawData();

    Frame frame;
    if (!m_output.begin_frame(frame)) {
        return;
    }

    if (frame.swapchain && frame.width > 0 && frame.height > 0) {
        ImGui_ImplSDLGPU3_PrepareDrawData(draw_data, frame.cmd); // a copy pass: before the render passes
        m_scene_renderer.render(frame, m_sim, m_view, m_view_options);

        SDL_GPURenderPass* overlay = RenderOutput::begin_overlay_pass(frame);
        ImGui_ImplSDLGPU3_RenderDrawData(draw_data, frame.cmd, overlay);
        SDL_EndGPURenderPass(overlay);
    }

    m_output.end_frame(frame);
}

} // namespace astraxis
