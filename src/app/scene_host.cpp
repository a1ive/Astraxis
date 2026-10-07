#include "app/scene_host.hpp"

#include "app/label_overlay.hpp"
#include "app/scene_list.hpp"
#include "platform/paths.hpp"

#include <SDL3/SDL.h>
#include <imgui.h>
#include <imgui_impl_sdlgpu3.h>

#include <algorithm>

namespace astraxis {

namespace {

constexpr double kMaxRealDt = 0.1; // as in App
constexpr uint64_t kLogPeriodNs = 10'000'000'000;

#ifdef NDEBUG
constexpr bool kGpuDebug = false;
#else
constexpr bool kGpuDebug = true;
#endif

} // namespace

bool SceneHost::init_gpu()
{
    return m_gpu.init(kGpuDebug);
}

bool SceneHost::add_output(std::unique_ptr<Output> out)
{
    const bool ok = out->gdi ? out->offscreen.init(m_gpu, SDL_GPU_TEXTUREFORMAT_B8G8R8A8_UNORM)
                             : out->output.init(m_gpu, out->window.handle());
    if (!ok) {
        out->window.destroy();
        return false;
    }
    out->output.set_render_scale(m_render_scale);
    m_outputs.push_back(std::move(out));
    return true;
}

void SceneHost::set_render_scale(float scale)
{
    m_render_scale = scale;
    for (auto& out : m_outputs) {
        out->output.set_render_scale(scale);
    }
}

bool SceneHost::any_active() const
{
    for (const auto& out : m_outputs) {
        if (out->active) {
            return true;
        }
    }
    return false;
}

void SceneHost::remove_outputs()
{
    if (m_gpu.device()) {
        SDL_WaitForGPUIdle(m_gpu.device());
    }
    for (auto& out : m_outputs) {
        out->output.shutdown();
        out->offscreen.shutdown();
        out->window.destroy();
    }
    m_outputs.clear();
}

bool SceneHost::start(const Settings& settings, const std::optional<std::string>& scene, bool labels, bool title,
                      std::string* error)
{
    m_view_options = settings.view;
    m_show_labels = labels;
    m_show_title = title;

    // All outputs share one swapchain format (SDR), and the offscreen one is
    // the same: B8G8R8A8_UNORM.
    const std::filesystem::path asset_dir = asset_directory();
    const Output& first = *m_outputs.front();
    const SDL_GPUTextureFormat format = first.gdi ? first.offscreen.format() : first.output.swapchain_format();
    if (!m_scene_renderer.init(m_gpu, format, asset_dir)) {
        *error = "The renderer could not be set up";
        return false;
    }
    m_renderer_ready = true;

    const std::vector<std::filesystem::path> scenes = list_scene_files(asset_dir);
    const std::string name = scene.value_or(settings.scene);
    int index = find_scene(scenes, name);
    if (index < 0) {
        SDL_LogWarn(SDL_LOG_CATEGORY_APPLICATION, "Scene '%s' not found", name.c_str());
        index = find_scene(scenes, Settings{}.scene);
    }
    *error = "No scenes found in " + (asset_dir / "scenes").string();
    if (index < 0 || !m_sim.load_scene(scenes[static_cast<size_t>(index)], error)) {
        return false;
    }
    error->clear();

    int max_height = 0;
    for (const auto& out : m_outputs) {
        int width = 0;
        int height = 0;
        if (out->gdi) {
            out->window.child_size(width, height);
        } else {
            SDL_GetWindowSizeInPixels(out->window.handle(), &width, &height);
        }
        max_height = std::max(max_height, height);
    }
    // The lensing sky's stars are sized for the scene's height.
    m_scene_renderer.load_scene(m_sim.scene(), m_view_options,
                                static_cast<uint32_t>(static_cast<float>(max_height) * m_render_scale));
    m_sim.start_tour();

    // ImGui only draws the labels and the title: no platform backend (there
    // is no input to pass on), one frame per output.
    if (m_show_labels || m_show_title) {
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGui::GetIO().IniFilename = nullptr;
        // The default font is chosen (pixel or vector) for the size it is drawn at.
        ImGui::GetStyle().FontScaleDpi = m_outputs.front()->window.content_scale();
        m_fonts = load_fonts(asset_dir);
        ImGui_ImplSDLGPU3_InitInfo init_info = {};
        init_info.Device = m_gpu.device();
        init_info.ColorTargetFormat = format;
        init_info.MSAASamples = SDL_GPU_SAMPLECOUNT_1;
        ImGui_ImplSDLGPU3_Init(&init_info);
        m_imgui_ready = true;
    }

    const uint64_t now = SDL_GetTicksNS();
    SDL_LogVerbose(SDL_LOG_CATEGORY_APPLICATION, "Ready after %.2f s", static_cast<double>(now) / SDL_NS_PER_SECOND);
    m_log_start_ns = now;
    reset_clock();
    return true;
}

void SceneHost::reset_clock()
{
    m_last_counter = SDL_GetPerformanceCounter();
    m_frame_start_ns = SDL_GetTicksNS();
}

void SceneHost::frame()
{
    const uint64_t now = SDL_GetPerformanceCounter();
    m_real_dt = std::min(static_cast<double>(now - m_last_counter) / SDL_GetPerformanceFrequency(), kMaxRealDt);
    m_last_counter = now;

    m_sim.update(m_real_dt);
    // Every output waits for its swapchain. Only waiting for the first one
    // skipped a third of the frames of a slower one (4K next to 1080p: 60 and
    // 40 FPS, against 56 and 56); the outputs then run at the slowest's pace.
    for (auto& out : m_outputs) {
        if (out->active) {
            render_output(*out);
        }
    }
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

void SceneHost::shutdown()
{
    remove_outputs();
    if (m_imgui_ready) {
        ImGui_ImplSDLGPU3_Shutdown();
        ImGui::DestroyContext();
        m_imgui_ready = false;
    }
    if (m_renderer_ready) {
        m_scene_renderer.shutdown();
        m_renderer_ready = false;
    }
    m_gpu.shutdown();
}

void SceneHost::render_output(Output& out)
{
    Frame frame;
    if (out.gdi) {
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

    if (!out.output.begin_frame(frame)) {
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

void SceneHost::draw_scene(Output& out, const Frame& frame)
{
    // The scene's pixels per point: the display's, times the render scale.
    const float scale = out.window.content_scale();
    const double scene_scale = static_cast<double>(scale) * frame.width / std::max(frame.output_width, 1u);
    m_sim.compute_view(static_cast<int>(frame.width), static_cast<int>(frame.height), scene_scale, out.view);

    ImDrawData* draw_data = nullptr;
    const bool title = m_show_title && scene_title_visible(m_sim.scene_age());
    if (m_show_labels || title) {
        // As imgui_impl_sdl3 would set it up for this window; the labels are
        // drawn at the output's full resolution.
        int w = static_cast<int>(frame.output_width);
        int h = static_cast<int>(frame.output_height);
        if (out.window.handle()) {
            SDL_GetWindowSize(out.window.handle(), &w, &h);
        }
        ImGuiIO& io = ImGui::GetIO();
        io.DisplaySize = ImVec2(static_cast<float>(w), static_cast<float>(h));
        io.DisplayFramebufferScale = ImVec2(w > 0 ? static_cast<float>(frame.output_width) / w : 1.0f,
                                            h > 0 ? static_cast<float>(frame.output_height) / h : 1.0f);
        io.DeltaTime = static_cast<float>(std::max(m_real_dt, 1e-4));
        ImGui::GetStyle().FontScaleDpi = scale;
        ImGui_ImplSDLGPU3_NewFrame(); // creates the samplers on first use
        ImGui::NewFrame();
        if (m_show_labels) {
            draw_labels(out.labels, m_sim.scene(), out.view, m_sim.camera().target());
        }
        if (title) {
            draw_scene_title(m_fonts, m_sim.scene().name, m_sim.scene_age());
        }
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

void SceneHost::log_frame_rates()
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
