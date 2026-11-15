#include "app/hr_diagram_overlay.hpp"

#include "core/color.hpp"
#include "scene/hr_diagram.hpp"
#include "scene/scene.hpp"

#include <imgui.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <vector>

namespace astraxis {

namespace {

// The diagram's span, the same in every scene so that stars can be compared
// between them: hot on the left, as is customary.
constexpr double kLogTLeft = 4.85;  // 70,800 K
constexpr double kLogTRight = 3.35; // 2,240 K
constexpr double kLogLTop = 7.0;
constexpr double kLogLBottom = -5.0;

ImU32 gray(int value, float alpha)
{
    return IM_COL32(value, value + 4, value + 10, static_cast<int>(255.0f * std::clamp(alpha, 0.0f, 1.0f)));
}

// A star's color for the panel: its black body, sRGB encoded.
ImU32 star_color(double kelvin, float alpha)
{
    const glm::dvec3 c = blackbody_linear_srgb(kelvin);
    const auto encode = [](double v) {
        v = std::clamp(v, 0.0, 1.0);
        v = v <= 0.0031308 ? 12.92 * v : 1.055 * std::pow(v, 1.0 / 2.4) - 0.055;
        return static_cast<int>(std::lround(255.0 * v));
    };
    return IM_COL32(encode(c.r), encode(c.g), encode(c.b), static_cast<int>(255.0f * alpha));
}

const HrBand& ms_band()
{
    static const HrBand band = main_sequence_band();
    return band;
}

const HrBand& wd_band()
{
    static const HrBand band = white_dwarf_band();
    return band;
}

// The band's lower edge (its first half: the zero-age main sequence, or the
// white dwarfs' smaller radius) at `log_t`, by straight lines between samples.
double lower_edge_at(const HrBand& band, size_t count, double log_t)
{
    double best = -FLT_MAX;
    for (size_t i = 1; i < count; ++i) {
        const HrPoint& a = band.outline[i - 1];
        const HrPoint& b = band.outline[i];
        if ((log_t - a.log_t) * (log_t - b.log_t) <= 0.0 && a.log_t != b.log_t) {
            best = std::max(best, a.log_l + (log_t - a.log_t) / (b.log_t - a.log_t) * (b.log_l - a.log_l));
        }
    }
    return best;
}

} // namespace

bool shows_hr_diagram(const Body& body)
{
    return body.kind == BodyKind::Star && !body.neutron_star && body.hr_diagram && body.temperature_k > 0.0 &&
           body.luminosity_solar > 0.0;
}

void draw_hr_diagram(const Scene& scene, int focus, float width)
{
    const float em = ImGui::GetFontSize();
    const float small = 0.85f * em;
    const float exponent = 0.65f * em;
    ImFont* font = ImGui::GetFont();
    ImDrawList* draw = ImGui::GetWindowDrawList();

    const ImVec2 origin = ImGui::GetCursorScreenPos();
    const float height = 0.8f * width;
    ImGui::Dummy(ImVec2(width, height));

    // Room for the luminosity labels on the left and the temperatures below.
    const float left = origin.x + 2.2f * em;
    const float right = origin.x + width - 0.3f * em;
    const float top = origin.y + 0.5f * em;
    const float bottom = origin.y + height - small - 0.4f * em;
    const auto to_x = [&](double log_t) {
        return left + static_cast<float>((kLogTLeft - log_t) / (kLogTLeft - kLogTRight)) * (right - left);
    };
    const auto to_y = [&](double log_l) {
        return top + static_cast<float>((kLogLTop - log_l) / (kLogLTop - kLogLBottom)) * (bottom - top);
    };
    const auto from_x = [&](float x) {
        return kLogTLeft - static_cast<double>((x - left) / (right - left)) * (kLogTLeft - kLogTRight);
    };
    const auto to_screen = [&](HrPoint p) { return ImVec2(to_x(p.log_t), to_y(p.log_l)); };

    char text[32];
    // Luminosity: a faint line every 4 dex, labelled as powers of ten.
    for (int e = -4; e <= 4; e += 4) {
        const float y = to_y(e);
        draw->AddLine(ImVec2(left, y), ImVec2(right, y), gray(150, 0.12f), 1.0f);
        if (e == 0) {
            const float w = font->CalcTextSizeA(small, FLT_MAX, 0.0f, "1").x;
            draw->AddText(font, small, ImVec2(left - 0.4f * em - w, y - 0.5f * small), gray(170, 0.7f), "1");
            continue;
        }
        std::snprintf(text, sizeof(text), "%d", e);
        const float we = font->CalcTextSizeA(exponent, FLT_MAX, 0.0f, text).x;
        const float wb = font->CalcTextSizeA(small, FLT_MAX, 0.0f, "10").x;
        const float x = left - 0.4f * em - we - wb;
        draw->AddText(font, small, ImVec2(x, y - 0.5f * small), gray(170, 0.7f), "10");
        draw->AddText(font, exponent, ImVec2(x + wb, y - 0.5f * small - 0.3f * exponent), gray(170, 0.7f), text);
    }
    draw->AddText(font, exponent, ImVec2(origin.x, origin.y), gray(170, 0.6f), "L_sun");
    // Temperature: a faint line and a label at each.
    for (const double kelvin : {30000.0, 10000.0, 3000.0}) {
        const float x = to_x(std::log10(kelvin));
        draw->AddLine(ImVec2(x, top), ImVec2(x, bottom), gray(150, 0.12f), 1.0f);
        std::snprintf(text, sizeof(text), "%.0f %03.0f K", std::floor(kelvin / 1000.0), std::fmod(kelvin, 1000.0));
        const float w = font->CalcTextSizeA(small, FLT_MAX, 0.0f, text).x;
        const float tx = std::clamp(x - 0.5f * w, origin.x, origin.x + width - w);
        draw->AddText(font, small, ImVec2(tx, bottom + 0.3f * em), gray(170, 0.7f), text);
    }

    draw->PushClipRect(ImVec2(left, top), ImVec2(right, bottom), true);

    // The bands, each named just below its lower edge: the main sequence's name
    // ending at 10,000 K, beside the A stars; the white dwarfs' starting at the
    // left, where nothing else is.
    std::vector<ImVec2> points;
    const auto draw_band = [&](const HrBand& band, float fill, const char* name, float label_x, bool right_aligned) {
        points.clear();
        for (const HrPoint& p : band.outline) {
            points.push_back(to_screen(p));
        }
        draw->AddConcavePolyFilled(points.data(), static_cast<int>(points.size()), gray(150, fill));
        draw->AddPolyline(points.data(), static_cast<int>(points.size()), gray(170, 2.0f * fill), ImDrawFlags_Closed,
                          1.0f);
        const float w = font->CalcTextSizeA(small, FLT_MAX, 0.0f, name).x;
        const float x = right_aligned ? label_x - w : label_x;
        // The edge falls to the right, so clear it at the label's right end.
        const double edge = lower_edge_at(band, band.outline.size() / 2, from_x(x + w));
        draw->AddText(font, small, ImVec2(x, to_y(edge) + 0.2f * em), gray(170, 0.55f), name);
    };
    draw_band(wd_band(), 0.09f, "white dwarfs", left + 0.4f * em, false);
    draw_band(ms_band(), 0.16f, "main sequence", to_x(4.0), true);

    // The Sun's place, for reference.
    const ImVec2 sun = to_screen({std::log10(5772.0), 0.0});
    draw->AddCircle(sun, 0.32f * em, gray(200, 0.6f), 0, 1.0f);

    // The scene's stars, the focused one last.
    const float dot = 0.17f * em;
    for (int i = 0; i < static_cast<int>(scene.bodies.size()); ++i) {
        const Body& b = scene.bodies[static_cast<size_t>(i)];
        if (i != focus && shows_hr_diagram(b)) {
            draw->AddCircleFilled(to_screen({std::log10(b.temperature_k), std::log10(b.luminosity_solar)}), dot,
                                  star_color(b.temperature_k, 0.75f));
        }
    }
    if (focus >= 0 && focus < static_cast<int>(scene.bodies.size())) {
        const Body& b = scene.bodies[static_cast<size_t>(focus)];
        if (shows_hr_diagram(b)) {
            const ImVec2 p = to_screen({std::log10(b.temperature_k), std::log10(b.luminosity_solar)});
            draw->AddCircle(p, 0.5f * em, gray(230, 0.85f), 0, 1.2f);
            // A pulsating star's cycle, over the ring: it is small at this scale.
            if (b.pulsation && b.pulsation->samples.size() > 1) {
                points.clear();
                for (const Pulsation::Sample& s : b.pulsation->samples) {
                    points.push_back(to_screen({std::log10(s.temperature_k), std::log10(s.luminosity_solar)}));
                }
                draw->AddPolyline(points.data(), static_cast<int>(points.size()), star_color(b.temperature_k, 0.9f),
                                  ImDrawFlags_Closed, 1.0f);
            }
            draw->AddCircleFilled(p, 0.22f * em, star_color(b.temperature_k, 1.0f));
        }
    }

    draw->PopClipRect();
}

} // namespace astraxis
