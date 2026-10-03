// ImGui panels and on-screen body labels.

#include "app/app.hpp"

#include "core/time.hpp"

#include <glm/glm.hpp>
#include <imgui.h>
#include <imgui_impl_sdl3.h>
#include <imgui_impl_sdlgpu3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iterator>

namespace astraxis {

namespace {

constexpr float kMarkerRadius = 2.5f;     // points; drawn when the body is smaller than this on screen
constexpr float kLabelPickRadius = 10.0f; // points
constexpr float kLabelPadding = 2.0f;     // points around each label when checking overlaps
constexpr float kLabelFadeSeconds = 0.3f; // labels ease in/out over this long

// "1 s = 2.5 h" style description of a time warp.
void format_warp(double warp, bool reverse, char* buf, size_t size)
{
    const char* unit = "s";
    double value = warp;
    if (warp >= kDaysPerJulianYear * kSecondsPerDay) {
        value = warp / (kDaysPerJulianYear * kSecondsPerDay);
        unit = "yr";
    } else if (warp >= kSecondsPerDay) {
        value = warp / kSecondsPerDay;
        unit = "d";
    } else if (warp >= 3600.0) {
        value = warp / 3600.0;
        unit = "h";
    } else if (warp >= 60.0) {
        value = warp / 60.0;
        unit = "min";
    }
    std::snprintf(buf, size, "1 s = %s%.3g %s", reverse ? "-" : "", value, unit);
}

ImU32 to_imgui_color(const glm::vec3& c, float alpha)
{
    return ImGui::ColorConvertFloat4ToU32(ImVec4(c.r, c.g, c.b, alpha));
}

// True if the segment from the camera (origin) to `target` passes through a sphere.
bool occluded(const glm::dvec3& target, const glm::dvec3& center, double radius)
{
    const double dist = glm::length(target);
    const glm::dvec3 dir = target / dist;
    const double along = glm::dot(center, dir);
    if (along <= 0.0 || along >= dist) {
        return false;
    }
    const double perp2 = glm::dot(center, center) - along * along;
    return perp2 < radius * radius;
}

} // namespace

void App::build_ui()
{
    ImGui_ImplSDLGPU3_NewFrame();
    ImGui_ImplSDL3_NewFrame();
    ImGui::NewFrame();

    if (m_show_labels) {
        build_labels();
    }
    if (panel_visible()) {
        build_control_panel();
    }
    if (m_show_demo) {
        ImGui::ShowDemoWindow(&m_show_demo);
    }

    ImGui::Render();
}

void App::build_control_panel()
{
    const ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(12.0f, 12.0f), ImGuiCond_FirstUseEver);
    ImGui::Begin("Astraxis", nullptr, ImGuiWindowFlags_AlwaysAutoResize);

    // --- Scene ---
    if (m_scene_files.size() > 1) {
        const std::string current = m_scene_files[m_scene_index].stem().string();
        ImGui::SetNextItemWidth(180.0f * ImGui::GetStyle().FontScaleDpi);
        if (ImGui::BeginCombo("Scene", current.c_str())) {
            for (size_t i = 0; i < m_scene_files.size(); ++i) {
                const std::string name = m_scene_files[i].stem().string();
                if (ImGui::Selectable(name.c_str(), i == m_scene_index) && i != m_scene_index) {
                    load_scene(i);
                }
            }
            ImGui::EndCombo();
        }
    }
    if (!m_scene_error.empty()) {
        ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.4f, 1.0f), "%s", m_scene_error.c_str());
    }

    // --- Time ---
    ImGui::TextUnformatted(format_utc(m_clock.t_tdb).c_str());

    float warp = static_cast<float>(m_clock.warp);
    char warp_label[64];
    format_warp(m_clock.warp, m_clock.reverse, warp_label, sizeof(warp_label));
    ImGui::SetNextItemWidth(260.0f * ImGui::GetStyle().FontScaleDpi);
    if (ImGui::SliderFloat("##warp", &warp, 1.0f, 1e8f, warp_label, ImGuiSliderFlags_Logarithmic)) {
        m_clock.warp = std::clamp(static_cast<double>(warp), 1.0, 1e8);
    }

    if (ImGui::Button(m_clock.paused ? "Resume" : "Pause")) {
        m_clock.paused = !m_clock.paused;
    }
    ImGui::SameLine();
    ImGui::Checkbox("Reverse", &m_clock.reverse);
    ImGui::SameLine();
    if (ImGui::Button("Now")) {
        reset_to_now();
    }

    struct Preset {
        const char* label;
        double warp;
    };
    static constexpr Preset kPresets[] = {
        {"1x", 1.0}, {"1 h/s", 3600.0}, {"1 d/s", 86400.0}, {"1 mo/s", 30.436875 * 86400.0},
        {"1 yr/s", 365.25 * 86400.0},
    };
    for (size_t i = 0; i < std::size(kPresets); ++i) {
        if (i > 0) {
            ImGui::SameLine();
        }
        if (ImGui::Button(kPresets[i].label)) {
            m_clock.warp = kPresets[i].warp;
        }
    }

    ImGui::SetNextItemWidth(180.0f * ImGui::GetStyle().FontScaleDpi);
    const bool entered = ImGui::InputTextWithHint("##goto", "YYYY-MM-DD HH:MM", m_goto_text, sizeof(m_goto_text),
                                                  ImGuiInputTextFlags_EnterReturnsTrue);
    ImGui::SameLine();
    if (ImGui::Button("Go to") || entered) {
        double t = 0.0;
        m_goto_error = !parse_utc(m_goto_text, &t);
        if (!m_goto_error) {
            m_clock.t_tdb = t;
        }
    }
    if (m_goto_error) {
        ImGui::SameLine();
        ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.4f, 1.0f), "invalid date");
    }

    // --- View ---
    ImGui::SeparatorText("View");
    const float combo_width = 200.0f * ImGui::GetStyle().FontScaleDpi;
    ImGui::SetNextItemWidth(combo_width);
    if (ImGui::BeginCombo("Focus", m_scene.bodies[static_cast<size_t>(m_camera.target())].name.c_str())) {
        for (size_t i = 0; i < m_scene.bodies.size(); ++i) {
            const Body& b = m_scene.bodies[i];
            ImGui::BeginDisabled(!b.visible);
            char label[96];
            std::snprintf(label, sizeof(label), "%zu  %s", i + 1, b.name.c_str());
            if (ImGui::Selectable(label, m_camera.target() == static_cast<int>(i))) {
                set_focus(static_cast<int>(i));
            }
            ImGui::EndDisabled();
        }
        ImGui::EndCombo();
    }
    if (m_scene.frames.size() > 1) {
        ImGui::SetNextItemWidth(combo_width);
        if (ImGui::BeginCombo("Frame", m_scene.frame().name.c_str())) {
            for (size_t i = 0; i < m_scene.frames.size(); ++i) {
                if (ImGui::Selectable(m_scene.frames[i].name.c_str(), m_scene.active_frame() == static_cast<int>(i))) {
                    set_frame(static_cast<int>(i));
                }
            }
            ImGui::EndCombo();
        }
    }
    if (!m_scene.events.empty()) {
        ImGui::SetNextItemWidth(combo_width);
        if (ImGui::BeginCombo("Events", "Jump to...")) {
            for (size_t i = 0; i < m_scene.events.size(); ++i) {
                if (ImGui::Selectable(m_scene.events[i].name.c_str())) {
                    jump_to_event(i);
                }
            }
            ImGui::EndCombo();
        }
    }
    if (ImGui::Button("View from sun side")) {
        stop_tour();
        // Slightly off the sun line, so moon shadows on the planet are not hidden behind the moons.
        const glm::dvec3 target = m_scene.bodies[static_cast<size_t>(m_camera.target())].world_position;
        m_camera.look_from(m_scene.sun_position() - target);
        m_camera.rotate(0.35, 0.1);
    }

    ImGui::SeparatorText("Auto tour");
    if (m_director.active()) {
        ImGui::Text("%s  (next in %.0f s)", CameraDirector::shot_name(m_director.current_shot()),
                    std::max(0.0, m_director.time_left()));
        ImGui::SameLine();
        if (ImGui::Button("Stop")) {
            stop_tour();
        }
    } else if (ImGui::Button("Start tour")) {
        start_tour();
    }
    ImGui::SameLine();
    ImGui::Checkbox("Start when idle (60 s)", &m_auto_tour);

    ImGui::SeparatorText("Display");
    ImGui::Checkbox("Orbits", &m_show_orbits);
    ImGui::SameLine();
    ImGui::Checkbox("Labels", &m_show_labels);
    if (!m_scene.belts.empty()) {
        ImGui::SameLine();
        ImGui::Checkbox("Belts", &m_show_belts);
    }
    ImGui::SameLine();
    ImGui::Checkbox("ImGui demo", &m_show_demo);
    ImGui::SetNextItemWidth(160.0f * ImGui::GetStyle().FontScaleDpi);
    ImGui::SliderFloat("Stars", &m_star_brightness, 0.0f, 2.0f, "%.2f");
    ImGui::SetNextItemWidth(160.0f * ImGui::GetStyle().FontScaleDpi);
    ImGui::SliderFloat("Line width", &m_line_width, 0.5f, 4.0f, "%.1f px");
    ImGui::SetNextItemWidth(160.0f * ImGui::GetStyle().FontScaleDpi);
    ImGui::SliderFloat("Exposure", &m_post_settings.exposure, 0.1f, 8.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
    ImGui::SetNextItemWidth(160.0f * ImGui::GetStyle().FontScaleDpi);
    ImGui::SliderFloat("Bloom", &m_post_settings.bloom_strength, 0.0f, 0.2f, "%.3f");

    ImGui::Separator();
    ImGui::TextDisabled("%.0f FPS  |  %s", io.Framerate, m_renderer.driver_name());
    ImGui::TextDisabled("Drag: rotate   Wheel: zoom   Double-click label: focus");
    ImGui::TextDisabled("Space pause  R reverse  [ ] warp  N now  1-9 focus  A tour  O/L  H hide");

    ImGui::End();
}

void App::build_labels()
{
    const ImGuiIO& io = ImGui::GetIO();
    ImDrawList* draw = ImGui::GetBackgroundDrawList();
    const glm::dvec3& cam = m_view.position;
    const float font_size = ImGui::GetFontSize();
    const bool want_pick = ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left) && !io.WantCaptureMouse;

    auto project = [&](const glm::dvec3& rel, ImVec2* screen, float* w) {
        const glm::vec4 clip = m_view.view_proj * glm::vec4(glm::vec3(rel), 1.0f);
        if (clip.w <= 0.0f) {
            return false;
        }
        *screen = ImVec2((clip.x / clip.w * 0.5f + 0.5f) * io.DisplaySize.x,
                         (0.5f - clip.y / clip.w * 0.5f) * io.DisplaySize.y);
        *w = clip.w;
        return true;
    };

    // Lagrange points and other markers: a small diamond and a name.
    for (const Marker& marker : m_scene.markers) {
        ImVec2 s;
        float w = 0.0f;
        if (!project(marker.world_position - cam, &s, &w)) {
            continue;
        }
        const float r = 4.0f;
        draw->AddQuad(ImVec2(s.x, s.y - r), ImVec2(s.x + r, s.y), ImVec2(s.x, s.y + r), ImVec2(s.x - r, s.y),
                      IM_COL32(180, 200, 255, 170), 1.2f);
        draw->AddText(ImVec2(s.x + r + 3.0f, s.y - font_size * 0.5f), IM_COL32(180, 200, 255, 170),
                      marker.name.c_str());
    }

    // Body labels are placed by priority: the focus, stars, bodies orbiting a
    // star, spacecraft, then moons, larger bodies first. A label that would
    // overlap one already placed is left out; labels ease in and out so that
    // they do not flicker as bodies pass each other.
    struct Candidate {
        int body;
        int rank;
        ImVec2 screen;
        float radius_px;
    };
    std::vector<Candidate> candidates;
    const int focus = m_camera.target();
    m_label_alpha.resize(m_scene.bodies.size(), 0.0f);
    for (size_t i = 0; i < m_scene.bodies.size(); ++i) {
        const Body& body = m_scene.bodies[i];
        if (!body.visible || body.kind == BodyKind::Barycenter || m_body_fades[i] <= 0.0f) {
            m_label_alpha[i] = 0.0f;
            continue;
        }
        const glm::dvec3 rel = body.world_position - cam;

        bool hidden = false;
        for (size_t k = 0; k < m_scene.bodies.size() && !hidden; ++k) {
            const Body& other = m_scene.bodies[k];
            if (k != i && other.visible && other.kind != BodyKind::Spacecraft &&
                other.kind != BodyKind::Barycenter) {
                hidden = occluded(rel, other.world_position - cam, other.equatorial_radius_km);
            }
        }
        ImVec2 screen;
        float w = 0.0f;
        if (hidden || !project(rel, &screen, &w)) {
            m_label_alpha[i] = 0.0f;
            continue;
        }
        const float radius_px = static_cast<float>(body.equatorial_radius_km) / w * m_view.focal_y * 0.5f *
                                io.DisplaySize.y;
        const int host = m_scene.satellite_host(static_cast<int>(i));
        int rank = 4; // moons
        if (static_cast<int>(i) == focus) {
            rank = 0;
        } else if (body.kind == BodyKind::Star) {
            rank = 1;
        } else if (body.kind == BodyKind::Spacecraft) {
            rank = 3;
        } else if (host < 0 || m_scene.bodies[static_cast<size_t>(host)].kind == BodyKind::Star) {
            rank = 2;
        }
        candidates.push_back({static_cast<int>(i), rank, screen, radius_px});
    }
    std::sort(candidates.begin(), candidates.end(), [&](const Candidate& a, const Candidate& b) {
        if (a.rank != b.rank) {
            return a.rank < b.rank;
        }
        return m_scene.bodies[static_cast<size_t>(a.body)].equatorial_radius_km >
               m_scene.bodies[static_cast<size_t>(b.body)].equatorial_radius_km;
    });

    std::vector<ImVec4> placed; // label rectangles: min x, min y, max x, max y
    const float step = io.DeltaTime / kLabelFadeSeconds;
    for (const Candidate& c : candidates) {
        const size_t i = static_cast<size_t>(c.body);
        const Body& body = m_scene.bodies[i];
        const float fade = m_body_fades[i];

        const float offset = std::max(c.radius_px, kMarkerRadius) * 0.7071f + 4.0f;
        const ImVec2 text_pos(c.screen.x + offset, c.screen.y - offset - font_size * 0.5f);
        const ImVec2 text_size = ImGui::CalcTextSize(body.name.c_str());
        const ImVec4 rect(text_pos.x - kLabelPadding, text_pos.y - kLabelPadding,
                          text_pos.x + text_size.x + kLabelPadding, text_pos.y + text_size.y + kLabelPadding);
        bool overlaps = false;
        for (const ImVec4& r : placed) {
            if (rect.x < r.z && r.x < rect.z && rect.y < r.w && r.y < rect.w) {
                overlaps = true;
                break;
            }
        }
        if (!overlaps && fade > 0.05f) {
            placed.push_back(rect);
        }
        float& eased = m_label_alpha[i];
        eased = std::clamp(eased + (overlaps ? -step : step), 0.0f, 1.0f);
        const float alpha = eased * fade;

        // Spacecraft are never drawn as bodies, so they always get a marker.
        if (c.radius_px < kMarkerRadius || body.kind == BodyKind::Spacecraft) {
            draw->AddCircleFilled(c.screen, kMarkerRadius,
                                  to_imgui_color(glm::mix(body.color, glm::vec3(1.0f), 0.3f), 0.95f * fade));
        }
        if (alpha > 0.0f) {
            draw->AddText(text_pos, IM_COL32(220, 225, 235, static_cast<int>(190.0f * alpha)), body.name.c_str());
        }

        if (want_pick && alpha > 0.3f) {
            const float dx = io.MousePos.x - c.screen.x;
            const float dy = io.MousePos.y - c.screen.y;
            if (std::sqrt(dx * dx + dy * dy) < std::max(c.radius_px, kLabelPickRadius)) {
                set_focus(c.body);
            }
        }
    }
}

} // namespace astraxis
