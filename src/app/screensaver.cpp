#include "app/screensaver.hpp"

#include "app/label_overlay.hpp"
#include "app/scene_list.hpp"
#include "platform/paths.hpp"

#include <SDL3/SDL.h>
#include <imgui.h>
#include <imgui_impl_sdlgpu3.h>

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

constexpr double kMaxRealDt = 0.1;         // as in App
constexpr uint64_t kGraceNs = 300'000'000; // input right after the start is ignored
constexpr float kMouseSlack = 8.0f;        // pixels the mouse may drift before it counts
constexpr int kPreviewFps = 30;            // the preview is tiny; spare the GPU
constexpr uint64_t kLogPeriodNs = 10'000'000'000;

#ifdef NDEBUG
constexpr bool kGpuDebug = false;
#else
constexpr bool kGpuDebug = true;
#endif

} // namespace

bool Screensaver::init(const ScreensaverOptions& options)
{
    m_options = options;

    // Let the display power down while the screensaver runs (SDL otherwise
    // keeps it on), and keep every display's window up when another one
    // takes the focus.
    SDL_SetHint(SDL_HINT_VIDEO_ALLOW_SCREENSAVER, "1");
    SDL_SetHint(SDL_HINT_VIDEO_MINIMIZE_ON_FOCUS_LOSS, "0");
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_Init failed: %s", SDL_GetError());
        return false;
    }

    const std::filesystem::path config = options.config.empty() ? config_path() : options.config;
    Settings settings;
    std::string config_error;
    if (!load_settings(config, SettingsSection::Screensaver, settings, &config_error)) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Settings not read: %s", config_error.c_str());
    }
    m_view_options = settings.view;
    m_show_labels = settings.labels && !options.preview;
    m_fps_limit = options.preview ? kPreviewFps : options.fps.value_or(settings.fps);

    if (!m_gpu.init(kGpuDebug)) {
        return false;
    }
    int display = options.display.value_or(0);
    if (!options.display && settings.display > 0) {
        display = display_index(find_display(settings.display, settings.display_name));
    }
    if (!create_outputs(display)) {
        return false;
    }

    const std::filesystem::path asset_dir = asset_directory();
    const Output& first = *m_outputs.front();
    const SDL_GPUTextureFormat format = first.preview ? first.offscreen.format() : first.output.swapchain_format();
    if (!m_scene_renderer.init(m_gpu, format, asset_dir)) {
        return false;
    }

    const std::vector<std::filesystem::path> scenes = list_scene_files(asset_dir);
    const std::string name = options.scene.value_or(settings.scene);
    int index = find_scene(scenes, name);
    if (index < 0) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Scene '%s' not found", name.c_str());
        index = find_scene(scenes, Settings{}.scene);
    }
    std::string scene_error = "No scenes found in " + (asset_dir / "scenes").string();
    if (index < 0 || !m_sim.load_scene(scenes[static_cast<size_t>(index)], &scene_error)) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "%s", scene_error.c_str());
        if (!options.preview) {
            SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Astraxis", scene_error.c_str(), nullptr);
        }
        return false;
    }
    int max_height = 0;
    for (const auto& out : m_outputs) {
        int width = 0;
        int height = 0;
        if (out->preview) {
            out->window.child_size(width, height);
        } else {
            SDL_GetWindowSizeInPixels(out->window.handle(), &width, &height);
        }
        max_height = std::max(max_height, height);
    }
    m_scene_renderer.load_scene(m_sim.scene(), m_view_options, static_cast<uint32_t>(max_height));
    m_sim.start_tour();

    // ImGui only draws the labels: no platform backend (there is no input to
    // pass on), one frame per output.
    if (m_show_labels) {
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGui::GetIO().IniFilename = nullptr;
        ImGui_ImplSDLGPU3_InitInfo init_info = {};
        init_info.Device = m_gpu.device();
        init_info.ColorTargetFormat = format;
        init_info.MSAASamples = SDL_GPU_SAMPLECOUNT_1;
        ImGui_ImplSDLGPU3_Init(&init_info);
        m_imgui_ready = true;
    }

    if (!options.preview) {
        SDL_HideCursor();
    }
    m_start_ns = SDL_GetTicksNS();
    SDL_LogVerbose(SDL_LOG_CATEGORY_APPLICATION, "Ready after %.2f s",
                   static_cast<double>(m_start_ns) / SDL_NS_PER_SECOND);
    m_frame_start_ns = m_start_ns;
    m_log_start_ns = m_start_ns;
    m_last_counter = SDL_GetPerformanceCounter();
    m_running = true;
    return true;
}

bool Screensaver::create_outputs(int display)
{
    if (m_options.preview) {
        // No swapchain: the frames are read back and drawn with GDI.
        auto out = std::make_unique<Output>();
        out->preview = true;
        if (!out->window.create_child(m_options.parent) ||
            !out->offscreen.init(m_gpu, SDL_GPU_TEXTUREFORMAT_B8G8R8A8_UNORM)) {
            out->window.destroy();
            return false;
        }
        m_outputs.push_back(std::move(out));
        return true;
    }

    // One window per display; `display` (1-based) alone shows the scene, or all if 0.
    int count = 0;
    SDL_DisplayID* displays = SDL_GetDisplays(&count);
    for (int i = 0; displays && i < count; ++i) {
        auto out = std::make_unique<Output>();
        out->scene = display <= 0 || display > count || display == i + 1;
        if (!out->window.create_cover("Astraxis", displays[i])) {
            continue;
        }
        if (!out->output.init(m_gpu, out->window.handle())) {
            out->window.destroy();
            continue;
        }
        m_outputs.push_back(std::move(out));
    }
    SDL_free(displays);
    if (m_outputs.empty()) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "No screensaver window could be created");
        return false;
    }
    // The first window takes the keyboard.
    SDL_RaiseWindow(m_outputs.front()->window.handle());
    return true;
}

void Screensaver::run()
{
    while (m_running) {
        handle_events();
        if (m_options.preview && !m_outputs.front()->window.parent_alive()) {
            quit("preview window gone");
        }
        if (!m_running) {
            break;
        }

        const uint64_t now = SDL_GetPerformanceCounter();
        m_real_dt = std::min(static_cast<double>(now - m_last_counter) / SDL_GetPerformanceFrequency(), kMaxRealDt);
        m_last_counter = now;

        m_sim.update(m_real_dt);
        render();
        log_frame_rates();

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

void Screensaver::shutdown()
{
    if (m_gpu.device()) {
        SDL_WaitForGPUIdle(m_gpu.device());
    }
    if (m_imgui_ready) {
        ImGui_ImplSDLGPU3_Shutdown();
        ImGui::DestroyContext();
        m_imgui_ready = false;
    }
    m_scene_renderer.shutdown();
    for (auto& out : m_outputs) {
        out->output.shutdown();
        out->offscreen.shutdown();
        out->window.destroy();
    }
    m_outputs.clear();
    m_gpu.shutdown();
    SDL_ShowCursor();
    SDL_Quit();
}

void Screensaver::handle_events()
{
    // Input that arrives with the windows (focus changes, the click that
    // started a preview) does not count.
    const bool armed = !m_options.preview && SDL_GetTicksNS() - m_start_ns > kGraceNs;

    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
        case SDL_EVENT_QUIT:
        case SDL_EVENT_WINDOW_CLOSE_REQUESTED:
        case SDL_EVENT_WINDOW_DESTROYED:
            quit("window closed");
            break;
        case SDL_EVENT_KEY_DOWN:
            if (armed) {
                quit("key pressed");
            }
            break;
        case SDL_EVENT_MOUSE_BUTTON_DOWN:
        case SDL_EVENT_MOUSE_WHEEL:
            if (armed) {
                quit("mouse button or wheel");
            }
            break;
        case SDL_EVENT_MOUSE_MOTION: {
            // In desktop coordinates: motion events are relative to whichever
            // window the mouse is over.
            float x = 0.0f;
            float y = 0.0f;
            SDL_GetGlobalMouseState(&x, &y);
            if (!m_have_mouse) {
                m_have_mouse = true;
                m_mouse_x = x;
                m_mouse_y = y;
            } else if (armed && std::fabs(x - m_mouse_x) + std::fabs(y - m_mouse_y) > kMouseSlack) {
                quit("mouse moved");
            }
            break;
        }
        default:
            break;
        }
    }
}

void Screensaver::quit(const char* reason)
{
    if (m_running) {
        SDL_LogInfo(SDL_LOG_CATEGORY_APPLICATION, "Screensaver ends: %s", reason);
    }
    m_running = false;
}

void Screensaver::render()
{
    // Only the first output waits for its swapchain: waiting on each would
    // wait for each display's vertical blank in turn.
    for (size_t i = 0; i < m_outputs.size(); ++i) {
        render_output(*m_outputs[i], i == 0);
    }
}

void Screensaver::render_output(Output& out, bool wait)
{
    Frame frame;
    if (out.preview) {
        int width = 0;
        int height = 0;
        out.window.child_size(width, height);
        if (!out.offscreen.begin_frame(frame, static_cast<uint32_t>(width), static_cast<uint32_t>(height))) {
            return;
        }
        if (frame.swapchain) {
            draw_scene(out, frame);
        }
        if (const uint8_t* pixels = out.offscreen.end_frame(frame)) {
            out.window.draw_child(pixels, width, height);
            ++out.frames;
        }
        return;
    }

    if (!out.output.begin_frame(frame, wait)) {
        return;
    }
    if (frame.swapchain && frame.width > 0 && frame.height > 0) {
        ++out.frames;
        if (out.scene) {
            draw_scene(out, frame);
        } else {
            SDL_GPUColorTargetInfo color = {};
            color.texture = frame.swapchain;
            color.clear_color = {0.0f, 0.0f, 0.0f, 1.0f};
            color.load_op = SDL_GPU_LOADOP_CLEAR;
            color.store_op = SDL_GPU_STOREOP_STORE;
            SDL_EndGPURenderPass(SDL_BeginGPURenderPass(frame.cmd, &color, 1, nullptr));
        }
    }
    out.output.end_frame(frame);
}

void Screensaver::draw_scene(Output& out, const Frame& frame)
{
    const float scale = out.window.content_scale();
    m_sim.compute_view(static_cast<int>(frame.width), static_cast<int>(frame.height), scale, out.view);

    ImDrawData* draw_data = nullptr;
    if (m_show_labels) {
        // As imgui_impl_sdl3 would set it up for this window.
        int w = static_cast<int>(frame.width);
        int h = static_cast<int>(frame.height);
        if (out.window.handle()) {
            SDL_GetWindowSize(out.window.handle(), &w, &h);
        }
        ImGuiIO& io = ImGui::GetIO();
        io.DisplaySize = ImVec2(static_cast<float>(w), static_cast<float>(h));
        io.DisplayFramebufferScale = ImVec2(w > 0 ? static_cast<float>(frame.width) / w : 1.0f,
                                            h > 0 ? static_cast<float>(frame.height) / h : 1.0f);
        io.DeltaTime = static_cast<float>(std::max(m_real_dt, 1e-4));
        ImGui::GetStyle().FontScaleDpi = scale;
        ImGui_ImplSDLGPU3_NewFrame(); // creates the samplers on first use
        ImGui::NewFrame();
        draw_labels(out.labels, m_sim.scene(), out.view, m_sim.camera().target());
        ImGui::Render();
        draw_data = ImGui::GetDrawData();
        ImGui_ImplSDLGPU3_PrepareDrawData(draw_data, frame.cmd); // a copy pass: before the render passes
    }

    m_scene_renderer.render(frame, m_sim, out.view, m_view_options);

    if (draw_data) {
        SDL_GPURenderPass* overlay = RenderOutput::begin_overlay_pass(frame);
        ImGui_ImplSDLGPU3_RenderDrawData(draw_data, frame.cmd, overlay);
        SDL_EndGPURenderPass(overlay);
    }
}

void Screensaver::log_frame_rates()
{
    // SDL_LOGGING=app=verbose shows it.
    const uint64_t now = SDL_GetTicksNS();
    if (now - m_log_start_ns < kLogPeriodNs) {
        return;
    }
    const double seconds = static_cast<double>(now - m_log_start_ns) / SDL_NS_PER_SECOND;
    for (size_t i = 0; i < m_outputs.size(); ++i) {
        SDL_LogVerbose(SDL_LOG_CATEGORY_APPLICATION, "Output %zu: %.1f FPS", i,
                       static_cast<double>(m_outputs[i]->frames) / seconds);
        m_outputs[i]->frames = 0;
    }
    m_log_start_ns = now;
}

} // namespace astraxis
