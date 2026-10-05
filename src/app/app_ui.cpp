// ImGui panels and on-screen body labels.

#include "app/app.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <glm/glm.hpp>
#include <imgui.h>
#include <imgui_impl_sdl3.h>
#include <imgui_impl_sdlgpu3.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <iterator>

namespace astraxis {

namespace {

constexpr float kLabelPickRadius = 10.0f; // points

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

// Integer with thousands separators ("384,400").
void format_grouped(double value, char* buf, size_t size)
{
    char digits[32];
    std::snprintf(digits, sizeof(digits), "%.0f", std::fabs(value));
    const size_t n = std::strlen(digits);
    size_t out = 0;
    if (value < 0.0 && out + 1 < size) {
        buf[out++] = '-';
    }
    for (size_t i = 0; i < n && out + 1 < size; ++i) {
        if (i > 0 && (n - i) % 3 == 0 && out + 2 < size) {
            buf[out++] = ',';
        }
        buf[out++] = digits[i];
    }
    buf[out] = '\0';
}

// km below 0.1 au, then au, then light years.
void format_distance(double km, char* buf, size_t size)
{
    const double light_year_km = kSpeedOfLightKmS * kDaysPerJulianYear * kSecondsPerDay;
    if (km < 10.0) {
        std::snprintf(buf, size, "%.3g km", km);
    } else if (km < 0.1 * kAuKm) {
        char grouped[32];
        format_grouped(km, grouped, sizeof(grouped));
        std::snprintf(buf, size, "%s km", grouped);
    } else if (km < 0.1 * light_year_km) {
        std::snprintf(buf, size, "%.4g au", km / kAuKm);
    } else {
        std::snprintf(buf, size, "%.4g ly", km / light_year_km);
    }
}

void format_duration(double seconds, char* buf, size_t size)
{
    const double year = kDaysPerJulianYear * kSecondsPerDay;
    if (seconds < 60.0) {
        std::snprintf(buf, size, "%.3g s", seconds);
    } else if (seconds < 3600.0) {
        std::snprintf(buf, size, "%.3g min", seconds / 60.0);
    } else if (seconds < 2.0 * kSecondsPerDay) {
        std::snprintf(buf, size, "%.3g h", seconds / 3600.0);
    } else if (seconds < 2.0 * year) {
        std::snprintf(buf, size, "%.4g d", seconds / kSecondsPerDay);
    } else {
        std::snprintf(buf, size, "%.4g yr", seconds / year);
    }
}

// One "label  value" row of an info panel table.
void info_row(const char* label, const char* fmt, ...) IM_FMTARGS(2);
void info_row(const char* label, const char* fmt, ...)
{
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::TextDisabled("%s", label);
    ImGui::TableNextColumn();
    va_list args;
    va_start(args, fmt);
    ImGui::TextV(fmt, args);
    va_end(args);
}

bool begin_info_table()
{
    return ImGui::BeginTable("##rows", 2, ImGuiTableFlags_SizingFixedFit);
}

enum class TransportIcon { Rewind, Play, PlayBackward, Pause, FastForward };

// A square button with a drawn icon (the default font has no symbols);
// `lit` shows it as selected.
bool transport_button(const char* id, TransportIcon icon, float size, bool lit)
{
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const bool pressed = ImGui::InvisibleButton(id, ImVec2(size, size));
    const ImGuiCol bg = ImGui::IsItemActive()    ? ImGuiCol_ButtonActive
                        : ImGui::IsItemHovered() ? ImGuiCol_ButtonHovered
                        : lit                    ? ImGuiCol_ButtonActive
                                                 : ImGuiCol_Button;
    ImDrawList* draw = ImGui::GetWindowDrawList();
    draw->AddRectFilled(p, ImVec2(p.x + size, p.y + size), ImGui::GetColorU32(bg), ImGui::GetStyle().FrameRounding);

    const ImU32 fg = ImGui::GetColorU32(ImGuiCol_Text);
    const float cx = p.x + size * 0.5f;
    const float cy = p.y + size * 0.5f;
    const float r = size * 0.22f;
    auto triangle = [&](float tip_x, float base_x) {
        draw->AddTriangleFilled(ImVec2(base_x, cy - r), ImVec2(tip_x, cy), ImVec2(base_x, cy + r), fg);
    };
    switch (icon) {
    case TransportIcon::Rewind:
        triangle(cx - r, cx);
        triangle(cx, cx + r);
        break;
    case TransportIcon::FastForward:
        triangle(cx + r, cx);
        triangle(cx, cx - r);
        break;
    case TransportIcon::Play:
        triangle(cx + r, cx - r * 0.7f);
        break;
    case TransportIcon::PlayBackward:
        triangle(cx - r, cx + r * 0.7f);
        break;
    case TransportIcon::Pause:
        draw->AddRectFilled(ImVec2(cx - r * 0.7f, cy - r), ImVec2(cx - r * 0.2f, cy + r), fg);
        draw->AddRectFilled(ImVec2(cx + r * 0.2f, cy - r), ImVec2(cx + r * 0.7f, cy + r), fg);
        break;
    }
    return pressed;
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
        build_time_bar();
        if (m_show_info) {
            build_info_panel();
        }
    }
    if (m_show_demo) {
        ImGui::ShowDemoWindow(&m_show_demo);
    }

    ImGui::Render();
}

void App::build_control_panel()
{
    const Scene& scene = m_sim.scene();
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

    // --- Time (pause, direction and warp are on the time bar) ---
    ImGui::SeparatorText("Time");
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
            m_sim.clock().warp = kPresets[i].warp;
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
            m_sim.clock().t_tdb = t;
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
    if (ImGui::BeginCombo("Focus", scene.bodies[static_cast<size_t>(m_sim.camera().target())].name.c_str())) {
        for (size_t i = 0; i < scene.bodies.size(); ++i) {
            const Body& b = scene.bodies[i];
            ImGui::BeginDisabled(!b.visible);
            char label[96];
            std::snprintf(label, sizeof(label), "%zu  %s", i + 1, b.name.c_str());
            if (ImGui::Selectable(label, m_sim.camera().target() == static_cast<int>(i))) {
                m_sim.set_focus(static_cast<int>(i));
            }
            ImGui::EndDisabled();
        }
        ImGui::EndCombo();
    }
    if (scene.frames.size() > 1) {
        ImGui::SetNextItemWidth(combo_width);
        if (ImGui::BeginCombo("Frame", scene.frame().name.c_str())) {
            for (size_t i = 0; i < scene.frames.size(); ++i) {
                if (ImGui::Selectable(scene.frames[i].name.c_str(), scene.active_frame() == static_cast<int>(i))) {
                    m_sim.set_frame(static_cast<int>(i));
                }
            }
            ImGui::EndCombo();
        }
    }
    if (!scene.events.empty()) {
        ImGui::SetNextItemWidth(combo_width);
        if (ImGui::BeginCombo("Events", "Jump to...")) {
            for (size_t i = 0; i < scene.events.size(); ++i) {
                if (ImGui::Selectable(scene.events[i].name.c_str())) {
                    m_sim.jump_to_event(i);
                }
            }
            ImGui::EndCombo();
        }
    }
    if (ImGui::Button("View from sun side")) {
        m_sim.look_from_sun_side();
    }

    ImGui::SeparatorText("Auto tour");
    if (m_sim.director().active()) {
        ImGui::Text("%s  (next in %.0f s)", CameraDirector::shot_name(m_sim.director().current_shot()),
                    std::max(0.0, m_sim.director().time_left()));
        ImGui::SameLine();
        if (ImGui::Button("Stop")) {
            m_sim.stop_tour();
        }
    } else if (ImGui::Button("Start tour")) {
        start_tour();
    }
    ImGui::SameLine();
    ImGui::Checkbox("Start when idle (60 s)", &m_auto_tour);

    ImGui::SeparatorText("Display");
    ImGui::Checkbox("Orbits", &m_view_options.orbits);
    ImGui::SameLine();
    ImGui::Checkbox("Labels", &m_show_labels);
    ImGui::SameLine();
    ImGui::Checkbox("Info", &m_show_info);
    if (!scene.belts.empty()) {
        ImGui::SameLine();
        ImGui::Checkbox("Belts", &m_view_options.belts);
    }
    // Second row, so the panel stays narrow.
    bool row_started = false;
    if (std::any_of(scene.bodies.begin(), scene.bodies.end(),
                    [](const Body& b) { return b.atmosphere.enabled; })) {
        ImGui::Checkbox("Atmospheres", &m_view_options.atmospheres);
        ImGui::SetItemTooltip("Off: the surface under the clouds (Venus: Magellan radar; Titan: Cassini ISS) (M)");
        row_started = true;
    }
    if (std::any_of(scene.bodies.begin(), scene.bodies.end(), [](const Body& b) { return !b.plumes.empty(); })) {
        if (row_started) {
            ImGui::SameLine();
        }
        ImGui::Checkbox("Plumes", &m_view_options.plumes);
        ImGui::SetItemTooltip("Volcanic and cryovolcanic plumes, always erupting (P)");
        row_started = true;
    }
    if (std::any_of(scene.bodies.begin(), scene.bodies.end(), [](const Body& b) { return b.comet.enabled; })) {
        if (row_started) {
            ImGui::SameLine();
        }
        ImGui::Checkbox("Comets", &m_view_options.comets);
        ImGui::SetItemTooltip("Comae (C2) and ion tails (C)");
        row_started = true;
    }
    if (row_started) {
        ImGui::SameLine();
    }
    ImGui::Checkbox("ImGui demo", &m_show_demo);
    ImGui::SetNextItemWidth(160.0f * ImGui::GetStyle().FontScaleDpi);
    ImGui::SliderFloat("Stars", &m_view_options.star_brightness, 0.0f, 2.0f, "%.2f");
    ImGui::SetNextItemWidth(160.0f * ImGui::GetStyle().FontScaleDpi);
    ImGui::SliderFloat("Line width", &m_view_options.line_width, 0.5f, 4.0f, "%.1f px");
    ImGui::SetNextItemWidth(160.0f * ImGui::GetStyle().FontScaleDpi);
    ImGui::SliderFloat("Exposure", &m_view_options.post.exposure, 0.1f, 8.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
    ImGui::SetNextItemWidth(160.0f * ImGui::GetStyle().FontScaleDpi);
    ImGui::SliderFloat("Bloom", &m_view_options.post.bloom_strength, 0.0f, 0.2f, "%.3f");

    ImGui::Separator();
    ImGui::TextDisabled("%.0f FPS  |  %s", io.Framerate, m_gpu.driver_name());
    ImGui::TextDisabled("Drag: rotate   Wheel: zoom");
    ImGui::TextDisabled("Double-click label: focus");
    ImGui::TextDisabled("Space pause  R reverse  [ ] warp  N now");
    ImGui::TextDisabled("1-9 focus  A tour  O/L/I/M/P/C  H hide");
    ImGui::TextDisabled("F11 fullscreen  Esc quit");

    ImGui::End();
}

void App::build_time_bar()
{
    SimClock& clock = m_sim.clock();
    const ImGuiIO& io = ImGui::GetIO();
    const ImGuiStyle& style = ImGui::GetStyle();
    const float margin = 12.0f * style.FontScaleDpi;
    const float width = std::min(640.0f * style.FontScaleDpi, io.DisplaySize.x - 2.0f * margin);

    // Bottom center, like a video player's control bar.
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y - margin), ImGuiCond_Always,
                            ImVec2(0.5f, 1.0f));
    ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0.0f), ImVec2(width, FLT_MAX));
    ImGui::SetNextWindowBgAlpha(0.8f);
    ImGui::Begin("##time_bar", nullptr,
                 ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings |
                     ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav | ImGuiWindowFlags_NoMove);

    // Backward / play-pause / forward, the date, and "Now" at the right.
    const float size = ImGui::GetFrameHeight() * 1.8f;
    const bool running = !clock.paused;
    if (transport_button("##backward", TransportIcon::Rewind, size, running && clock.reverse)) {
        clock.reverse = true;
        clock.paused = false;
    }
    ImGui::SetItemTooltip("Run backward (R)");
    ImGui::SameLine();
    const TransportIcon play = running         ? TransportIcon::Pause
                               : clock.reverse ? TransportIcon::PlayBackward
                                               : TransportIcon::Play;
    if (transport_button("##play", play, size, false)) {
        clock.paused = !clock.paused;
    }
    ImGui::SetItemTooltip(running ? "Pause (Space)" : "Resume (Space)");
    ImGui::SameLine();
    if (transport_button("##forward", TransportIcon::FastForward, size, running && !clock.reverse)) {
        clock.reverse = false;
        clock.paused = false;
    }
    ImGui::SetItemTooltip("Run forward (R)");

    const float row_y = ImGui::GetCursorPosY() - size - style.ItemSpacing.y;
    ImGui::SameLine(0.0f, style.ItemSpacing.x * 3.0f);
    ImGui::SetCursorPosY(row_y + (size - ImGui::GetTextLineHeight()) * 0.5f);
    ImGui::TextUnformatted(format_utc(clock.t_tdb).c_str());

    const float now_width = ImGui::CalcTextSize("Now").x + style.FramePadding.x * 4.0f;
    ImGui::SameLine(ImGui::GetContentRegionMax().x - now_width);
    ImGui::SetCursorPosY(row_y);
    if (ImGui::Button("Now", ImVec2(now_width, size))) {
        m_sim.reset_to_now();
    }
    ImGui::SetItemTooltip("Jump to the present (N)");

    // Time warp across the whole width, taller than usual.
    char warp_label[64];
    format_warp(clock.warp, clock.reverse, warp_label, sizeof(warp_label));
    float warp = static_cast<float>(clock.warp);
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(style.FramePadding.x, style.FramePadding.y * 2.5f));
    ImGui::PushStyleVar(ImGuiStyleVar_GrabMinSize, style.GrabMinSize * 2.0f);
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::SliderFloat("##bar_warp", &warp, 1.0f, 1e8f, warp_label, ImGuiSliderFlags_Logarithmic)) {
        clock.warp = std::clamp(static_cast<double>(warp), 1.0, 1e8);
    }
    ImGui::PopStyleVar(2);

    ImGui::End();
}

void App::build_info_panel()
{
    const Scene& scene = m_sim.scene();
    const int focus = m_sim.camera().target();
    if (focus < 0 || focus >= static_cast<int>(scene.bodies.size())) {
        return;
    }
    const Body& body = scene.bodies[static_cast<size_t>(focus)];
    const ImGuiIO& io = ImGui::GetIO();
    const float scale = ImGui::GetStyle().FontScaleDpi;

    // Anchored at the middle of the right edge; purely informational, so it
    // never takes the mouse from the camera.
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x - 12.0f * scale, io.DisplaySize.y * 0.5f), ImGuiCond_Always,
                            ImVec2(1.0f, 0.5f));
    ImGui::SetNextWindowSizeConstraints(ImVec2(240.0f * scale, 0.0f), ImVec2(FLT_MAX, FLT_MAX));
    ImGui::SetNextWindowBgAlpha(0.8f); // labels are drawn behind it
    ImGui::Begin("##info", nullptr,
                 ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings |
                     ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav | ImGuiWindowFlags_NoMove |
                     ImGuiWindowFlags_NoMouseInputs);

    ImGui::TextUnformatted(body.name.c_str());

    // What it orbits; bodies at rest (e.g. the Sun far from an exoplanet system) orbit nothing.
    const Scene::OrbitCenter center = body.visible ? scene.orbit_center(focus) : Scene::OrbitCenter{};
    State relative;
    if (center.body >= 0) {
        const State s = scene.icrf_state_at(focus, scene.time());
        const State c = scene.icrf_state_at(center.body, scene.time());
        relative = {s.position - c.position, s.velocity - c.velocity};
    }
    const bool orbiting = center.body >= 0 && glm::length(relative.velocity) > 1e-6;
    const char* kind = nullptr;
    switch (body.kind) {
    case BodyKind::Star:
        kind = body.pulsar.enabled ? "Pulsar" : "Star";
        break;
    case BodyKind::Spacecraft:
        kind = "Spacecraft";
        break;
    case BodyKind::BlackHole:
        kind = "Black hole";
        break;
    case BodyKind::Barycenter:
        kind = "Barycenter";
        break;
    case BodyKind::Ghost:
        kind = "Comparison orbit";
        break;
    case BodyKind::Planet:
        break;
    }
    const char* host_name = orbiting ? scene.bodies[static_cast<size_t>(center.body)].name.c_str() : "";
    if (kind && orbiting) {
        ImGui::TextDisabled("%s, orbiting %s", kind, host_name);
    } else if (kind) {
        ImGui::TextDisabled("%s", kind);
    } else if (orbiting) {
        ImGui::TextDisabled("Orbiting %s", host_name);
    }
    if (!body.visible) {
        ImGui::TextDisabled("Not present at this time");
    }

    char buf[64];
    char buf2[64];

    // --- Physical ---
    const double gm = body.body_gm_km3_s2;
    const double a = body.equatorial_radius_km;
    const double b = body.equatorial_radius_b_km;
    const double c = body.polar_radius_km;
    const double mean_radius = std::cbrt(a * b * c);
    const bool has_surface = body.kind == BodyKind::Planet || body.kind == BodyKind::Star;
    if (has_surface || gm > 0.0) {
        ImGui::SeparatorText("Physical");
        if (begin_info_table()) {
            if (has_surface) {
                if (body.kind == BodyKind::Star) {
                    std::snprintf(buf2, sizeof(buf2), "%#.3g R_sun", mean_radius / kSunRadiusKm);
                } else {
                    std::snprintf(buf2, sizeof(buf2), "%#.3g R_earth", mean_radius / kEarthRadiusKm);
                }
                if (a != b) {
                    info_row("Radii", "%.4g x %.4g x %.4g km", a, b, c);
                    info_row("", "%s (mean)", buf2);
                } else if (a != c) {
                    format_grouped(a, buf, sizeof(buf));
                    info_row("Radius", "%s km equatorial", buf);
                    format_grouped(c, buf, sizeof(buf));
                    info_row("", "%s km polar", buf);
                    info_row("", "%s (mean)", buf2);
                } else {
                    format_grouped(a, buf, sizeof(buf));
                    info_row("Radius", "%s km", buf);
                    info_row("", "%s", buf2);
                }
            }
            if (body.kind == BodyKind::BlackHole) {
                format_distance(a, buf, sizeof(buf));
                info_row("Horizon", "%s", buf);
                info_row("Spin a/M", "%.3g", body.spin);
                if (body.spin > 0.0) {
                    format_distance(2.0 * gm / (kSpeedOfLightKmS * kSpeedOfLightKmS), buf, sizeof(buf));
                    info_row("Schwarzschild r", "%s", buf);
                }
            }
            if (gm > 0.0) {
                const char* label = body.kind == BodyKind::Barycenter ? "Total mass" : "Mass";
                if (gm >= 1e-3 * kSunGmKm3S2) {
                    info_row(label, "%.4g M_sun", gm / kSunGmKm3S2);
                } else {
                    info_row(label, "%.4g M_earth", gm / kEarthGmKm3S2);
                }
                const double kg = gm / kGravitationalConstantKm3KgS2;
                info_row("", "%.4g kg", kg);
                if (has_surface) {
                    const double volume_cm3 = 4.0 / 3.0 * kPi * a * b * c * 1e15;
                    info_row("Density", "%.3g g/cm^3", kg * 1e3 / volume_cm3);
                    info_row("Surface gravity", "%.3g m/s^2", gm / (mean_radius * mean_radius) * 1e3);
                }
            }
            if (has_surface && body.pm_rate_deg_per_day != 0.0 && !body.pulsar.enabled) {
                format_duration(360.0 / std::fabs(body.pm_rate_deg_per_day) * kSecondsPerDay, buf, sizeof(buf));
                info_row("Rotation", "%s%s", buf, body.pm_rate_deg_per_day < 0.0 ? " (retrograde)" : "");
            }
            if (has_surface && body.free_precession.enabled) {
                const Body::FreePrecession& fp = body.free_precession;
                format_duration(360.0 / fp.precession_deg_per_day * kSecondsPerDay, buf, sizeof(buf));
                info_row("Precession", "%s, at %.1f deg", buf, fp.nutation_deg);
                format_duration(360.0 / fp.spin_deg_per_day * kSecondsPerDay, buf, sizeof(buf));
                info_row("Spin (long axis)", "%s", buf);
            }
            if (body.kind == BodyKind::Star) {
                info_row("Temperature", "%.0f K", body.temperature_k);
                info_row("Luminosity", "%.3g L_sun", body.luminosity_solar);
            }
            ImGui::EndTable();
        }
    }

    // --- Orbit (osculating two-body orbit around the center) ---
    if (orbiting) {
        std::snprintf(buf, sizeof(buf), "Orbit around %s", host_name);
        ImGui::SeparatorText(buf);
        if (begin_info_table()) {
            format_distance(glm::length(relative.position), buf, sizeof(buf));
            info_row("Distance", "%s", buf);
            info_row("Speed", "%.3g km/s", glm::length(relative.velocity));
            OrbitElements el;
            if (osculating_elements(relative, center.gm, &el)) {
                info_row("Eccentricity", "%.3g%s", el.eccentricity, el.period_s > 0.0 ? "" : " (unbound)");
                format_distance(el.periapsis_km, buf, sizeof(buf));
                info_row("Periapsis", "%s", buf);
                if (el.period_s > 0.0) {
                    format_distance(el.apoapsis_km, buf, sizeof(buf));
                    info_row("Apoapsis", "%s", buf);
                    format_duration(el.period_s, buf, sizeof(buf));
                    info_row("Period", "%s", buf);
                }
            }
            ImGui::EndTable();
        }
    }

    // --- Camera ---
    if (body.visible) {
        ImGui::SeparatorText("From camera");
        if (begin_info_table()) {
            const double d = glm::length(body.world_position - m_view.camera.position);
            format_distance(d, buf, sizeof(buf));
            info_row("Distance", "%s", buf);
            format_duration(d / kSpeedOfLightKmS, buf, sizeof(buf));
            info_row("Light time", "%s", buf);
            ImGui::EndTable();
        }
    }

    ImGui::End();
}

void App::build_labels()
{
    const Scene& scene = m_sim.scene();
    const ImGuiIO& io = ImGui::GetIO();
    const float font_size = ImGui::GetFontSize();
    m_labels.update(scene, m_view, m_sim.camera().target(), glm::vec2(io.DisplaySize.x, io.DisplaySize.y),
                    font_size,
                    [](const std::string& text) {
                        const ImVec2 size = ImGui::CalcTextSize(text.c_str());
                        return glm::vec2(size.x, size.y);
                    },
                    io.DeltaTime);

    ImDrawList* draw = ImGui::GetBackgroundDrawList();

    // Lagrange points and other markers: a small diamond and a name.
    for (const LabelLayout::SceneMarker& marker : m_labels.markers()) {
        const ImVec2 s(marker.position.x, marker.position.y);
        const float r = 4.0f;
        draw->AddQuad(ImVec2(s.x, s.y - r), ImVec2(s.x + r, s.y), ImVec2(s.x, s.y + r), ImVec2(s.x - r, s.y),
                      IM_COL32(180, 200, 255, 170), 1.2f);
        draw->AddText(ImVec2(s.x + r + 3.0f, s.y - font_size * 0.5f), IM_COL32(180, 200, 255, 170),
                      scene.markers[static_cast<size_t>(marker.marker)].name.c_str());
    }

    const float dot_radius = LabelLayout::kMarkerRadius;
    for (const LabelLayout::BodyMark& mark : m_labels.bodies()) {
        const Body& body = scene.bodies[static_cast<size_t>(mark.body)];
        const ImVec2 center(mark.position.x, mark.position.y);
        if (mark.dot == LabelLayout::Dot::Hollow) {
            draw->AddCircle(center, dot_radius + 1.0f, to_imgui_color(body.orbit_color, 0.95f * mark.fade), 0, 1.5f);
        } else if (mark.dot == LabelLayout::Dot::Filled) {
            draw->AddCircleFilled(center, dot_radius,
                                  to_imgui_color(glm::mix(body.color, glm::vec3(1.0f), 0.3f), 0.95f * mark.fade));
        }
        if (mark.text_alpha > 0.0f) {
            draw->AddText(ImVec2(mark.text_position.x, mark.text_position.y),
                          IM_COL32(220, 225, 235, static_cast<int>(190.0f * mark.text_alpha)), body.name.c_str());
        }
    }

    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left) && !io.WantCaptureMouse) {
        const int picked = m_labels.pick(glm::vec2(io.MousePos.x, io.MousePos.y), kLabelPickRadius);
        if (picked >= 0) {
            m_sim.set_focus(picked);
        }
    }
}

} // namespace astraxis
