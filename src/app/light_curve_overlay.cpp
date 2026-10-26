#include "app/light_curve_overlay.hpp"

#include "core/time.hpp"
#include "scene/light_curve.hpp"
#include "scene/scene.hpp"

#include <imgui.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace astraxis {

namespace {

// In font sizes (em).
constexpr float kWidth = 20.0f;
constexpr float kHeight = 8.0f;
constexpr float kPadding = 0.7f;
constexpr float kTextSize = 1.0f;

ImU32 gray(int value, float alpha)
{
    return IM_COL32(value, value + 4, value + 10, static_cast<int>(255.0f * std::clamp(alpha, 0.0f, 1.0f)));
}

bool overlaps(const ScreenRect& a, const ScreenRect& b)
{
    return a.min_x < b.max_x && b.min_x < a.max_x && a.min_y < b.max_y && b.min_y < a.max_y;
}

// A recorded light curve (Scene::observed_light_curve) over its whole span:
// apparent magnitude, brighter up, the years along the bottom, the points
// joined where they are close enough, and the clock's moment marked with its
// year and magnitude (none outside the record or across a gap).
void draw_observed_light_curve(const CaptionFonts& fonts, const ObservedLightCurve& lc, double t_tdb,
                               const ScreenRect& rect, float alpha)
{
    if (alpha <= 0.0f) {
        return;
    }
    ImFont* font = fonts.book ? fonts.book : ImGui::GetFont();
    const float em = ImGui::GetFontSize();
    const float text_size = kTextSize * em;
    const float small = 0.85f * em;
    const float pad = kPadding * em;
    ImDrawList* draw = ImGui::GetBackgroundDrawList();
    draw->AddRectFilled(ImVec2(rect.min_x, rect.min_y), ImVec2(rect.max_x, rect.max_y),
                        IM_COL32(6, 8, 12, static_cast<int>(140.0f * alpha)), 0.4f * em);

    // Heading, and the moment's year and magnitude on the right.
    char text[96];
    std::snprintf(text, sizeof(text), "Seen from Earth, %s band%s%s", lc.band.c_str(), lc.source.empty() ? "" : " - ",
                  lc.source.c_str());
    draw->AddText(font, text_size, ImVec2(rect.min_x + pad, rect.min_y + 0.6f * pad), gray(200, 0.85f * alpha), text);
    const double now_year = 2000.0 + t_tdb / (kDaysPerJulianYear * kSecondsPerDay);
    const double now_mag = observed_magnitude(lc, t_tdb);
    if (std::isnan(now_mag)) {
        std::snprintf(text, sizeof(text), "%.0f", std::floor(now_year));
    } else {
        std::snprintf(text, sizeof(text), "%.0f  %s %.1f", std::floor(now_year), lc.band.c_str(), now_mag);
    }
    const float readout_width = font->CalcTextSizeA(text_size, FLT_MAX, 0.0f, text).x;
    draw->AddText(font, text_size, ImVec2(rect.max_x - pad - readout_width, rect.min_y + 0.6f * pad),
                  IM_COL32(255, 222, 170, static_cast<int>(230.0f * alpha)), text);

    // The magnitude range of the points shown, with a margin, in whole magnitudes.
    const std::vector<double>& t = lc.t_tdb;
    float lo_mag = FLT_MAX;
    float hi_mag = -FLT_MAX;
    for (size_t i = 0; i < t.size(); ++i) {
        if (t[i] >= lc.from_tdb && t[i] <= lc.to_tdb) {
            lo_mag = std::min(lo_mag, lc.magnitude[i]);
            hi_mag = std::max(hi_mag, lc.magnitude[i]);
        }
    }
    if (lo_mag > hi_mag) {
        return;
    }
    constexpr double kMarginMag = 0.3;
    const double top_mag = std::floor(static_cast<double>(lo_mag) - kMarginMag);
    const double bottom_mag = std::max(std::ceil(static_cast<double>(hi_mag) + kMarginMag), top_mag + 1.0);

    const float left = rect.min_x + pad + 1.6f * em;
    const float right = rect.max_x - pad;
    const float top = rect.min_y + 0.6f * pad + text_size + 0.6f * em;
    const float bottom = rect.max_y - pad - small - 0.3f * em;
    const auto to_y = [&](double mag) {
        return top + static_cast<float>((mag - top_mag) / (bottom_mag - top_mag)) * (bottom - top);
    };
    const auto to_x = [&](double time) {
        return left + static_cast<float>((time - lc.from_tdb) / (lc.to_tdb - lc.from_tdb)) * (right - left);
    };

    const double step = bottom_mag - top_mag > 6.0 ? 2.0 : 1.0;
    for (double m = std::ceil(top_mag / step) * step; m <= bottom_mag + 1e-6; m += step) {
        const float y = to_y(m);
        draw->AddLine(ImVec2(left, y), ImVec2(right, y), gray(150, 0.18f * alpha), 1.0f);
        std::snprintf(text, sizeof(text), "%.0f", std::abs(m) < 1e-9 ? 0.0 : m); // no "-0"
        const float w = font->CalcTextSizeA(small, FLT_MAX, 0.0f, text).x;
        draw->AddText(font, small, ImVec2(left - 0.4f * em - w, y - 0.5f * small), gray(170, 0.6f * alpha), text);
    }
    // Years: a tick every 10, labels every 50.
    const double year_s = kDaysPerJulianYear * kSecondsPerDay;
    const double first_year = 2000.0 + lc.from_tdb / year_s;
    const double last_year = 2000.0 + lc.to_tdb / year_s;
    for (double y = std::ceil(first_year / 10.0) * 10.0; y <= last_year; y += 10.0) {
        const float tx = to_x((y - 2000.0) * year_s);
        const bool label = std::fmod(y, 50.0) == 0.0;
        draw->AddLine(ImVec2(tx, bottom), ImVec2(tx, bottom + (label ? 0.35f : 0.2f) * em),
                      gray(170, 0.45f * alpha), 1.0f);
        if (label) {
            std::snprintf(text, sizeof(text), "%.0f", y);
            const float w = font->CalcTextSizeA(small, FLT_MAX, 0.0f, text).x;
            draw->AddText(font, small, ImVec2(tx - 0.5f * w, bottom + 0.35f * em), gray(170, 0.6f * alpha), text);
        }
    }

    // The record: runs of points joined, each point a dot.
    const ImU32 curve_color = IM_COL32(255, 222, 170, static_cast<int>(200.0f * alpha));
    const float line_width = 1.2f * ImGui::GetStyle().FontScaleDpi;
    std::vector<ImVec2> run;
    for (size_t i = 0; i < t.size(); ++i) {
        if (t[i] < lc.from_tdb || t[i] > lc.to_tdb) {
            continue;
        }
        const ImVec2 p(to_x(t[i]), to_y(lc.magnitude[i]));
        if (!run.empty() && t[i] - t[i - 1] > lc.gap_s) {
            draw->AddPolyline(run.data(), static_cast<int>(run.size()), curve_color, 0, line_width);
            run.clear();
        }
        run.push_back(p);
        draw->AddCircleFilled(p, 0.12f * em, curve_color);
    }
    if (run.size() > 1) {
        draw->AddPolyline(run.data(), static_cast<int>(run.size()), curve_color, 0, line_width);
    }

    // The clock's moment.
    if (t_tdb >= lc.from_tdb && t_tdb <= lc.to_tdb) {
        const float nx = to_x(t_tdb);
        draw->AddLine(ImVec2(nx, top), ImVec2(nx, bottom), gray(220, 0.3f * alpha), 1.0f);
        if (!std::isnan(now_mag)) {
            draw->AddCircleFilled(ImVec2(nx, to_y(now_mag)), 0.25f * em,
                                  IM_COL32(255, 240, 215, static_cast<int>(255.0f * alpha)));
        }
    }
}

} // namespace

ImVec2 light_curve_plot_size()
{
    const float em = ImGui::GetFontSize();
    return ImVec2(kWidth * em, kHeight * em);
}

ScreenRect light_curve_plot_rect(const ImVec2& display, float bottom, float margin, float caption_width,
                                 std::span<const ScreenRect> avoid)
{
    const ImVec2 size = light_curve_plot_size();
    const float left = display.x - margin - size.x;
    ScreenRect rect{left, bottom - size.y, display.x - margin, bottom};
    bool clear = left >= 0.5f * (display.x + caption_width) + margin;
    for (const ScreenRect& r : avoid) {
        clear = clear && !overlaps(rect, r);
    }
    if (!clear) {
        rect = {left, margin, display.x - margin, margin + size.y};
    }
    return rect;
}

void draw_light_curve(const CaptionFonts& fonts, const Scene& scene, const LightCurveSamples& samples, double t_tdb,
                      double rate, const ScreenRect& rect, float alpha)
{
    if (scene.observed_light_curve.enabled()) {
        draw_observed_light_curve(fonts, scene.observed_light_curve, t_tdb, rect, alpha);
        return;
    }
    const LightCurve& lc = scene.light_curve;
    alpha *= static_cast<float>(light_curve_visibility(lc, rate));
    if (!lc.enabled() || samples.empty() || alpha <= 0.0f) {
        return;
    }
    const std::vector<float>& mags = samples.magnitudes();
    const int count = static_cast<int>(mags.size());
    const double t0 = samples.time(0);
    const double t1 = samples.time(count - 1);

    ImFont* font = fonts.book ? fonts.book : ImGui::GetFont();
    const float em = ImGui::GetFontSize();
    const float text_size = kTextSize * em;
    const float pad = kPadding * em;
    ImDrawList* draw = ImGui::GetBackgroundDrawList();
    draw->AddRectFilled(ImVec2(rect.min_x, rect.min_y), ImVec2(rect.max_x, rect.max_y),
                        IM_COL32(6, 8, 12, static_cast<int>(140.0f * alpha)), 0.4f * em);

    // Heading: where it is seen from and in which band; the dip now on the right.
    const std::string& observer = scene.bodies[static_cast<size_t>(lc.observer)].name;
    char text[96];
    std::snprintf(text, sizeof(text), "Seen from %s, %s band", observer == "Sun" ? "Earth" : observer.c_str(),
                  lc.band.c_str());
    const ImU32 text_color = gray(200, 0.85f * alpha);
    draw->AddText(font, text_size, ImVec2(rect.min_x + pad, rect.min_y + 0.6f * pad), text_color, text);

    const double x = std::clamp((t_tdb - t0) / (t1 - t0), 0.0, 1.0) * (count - 1);
    const int i0 = std::min(static_cast<int>(x), count - 2);
    const double now_mag = mags[static_cast<size_t>(i0)] +
                           (mags[static_cast<size_t>(i0 + 1)] - mags[static_cast<size_t>(i0)]) * (x - i0);
    std::snprintf(text, sizeof(text), fonts.book ? "\xce\x94m %.2f" : "dm %.2f", now_mag);
    const float readout_width = font->CalcTextSizeA(text_size, FLT_MAX, 0.0f, text).x;
    draw->AddText(font, text_size, ImVec2(rect.max_x - pad - readout_width, rect.min_y + 0.6f * pad),
                  IM_COL32(255, 222, 170, static_cast<int>(230.0f * alpha)), text);

    // The plot: 0 (nothing hidden) at the top, depth_mag at the bottom.
    const float left = rect.min_x + pad + 1.6f * em;
    const float right = rect.max_x - pad;
    const float top = rect.min_y + 0.6f * pad + text_size + 0.6f * em;
    const float bottom = rect.max_y - pad - 0.4f * em;
    const auto to_y = [&](double mag) {
        return top + static_cast<float>(std::clamp(mag / lc.depth_mag, 0.0, 1.05)) * (bottom - top);
    };
    const auto to_x = [&](double t) { return left + static_cast<float>((t - t0) / (t1 - t0)) * (right - left); };

    const double step = lc.depth_mag <= 0.3 ? 0.1 : lc.depth_mag <= 1.5 ? 0.5 : 1.0;
    const float small = 0.85f * em;
    for (double m = 0.0; m <= lc.depth_mag + 1e-6; m += step) {
        const float y = to_y(m);
        draw->AddLine(ImVec2(left, y), ImVec2(right, y), gray(150, 0.18f * alpha), 1.0f);
        std::snprintf(text, sizeof(text), step < 0.5 ? "%.1f" : m == std::floor(m) ? "%.0f" : "%.1f", m);
        const float w = font->CalcTextSizeA(small, FLT_MAX, 0.0f, text).x;
        draw->AddText(font, small, ImVec2(left - 0.4f * em - w, y - 0.5f * small), gray(170, 0.6f * alpha), text);
    }
    // A tick for every day (TDB), if they are not too dense.
    if (t1 - t0 <= 40.0 * kSecondsPerDay) {
        for (double d = std::ceil(t0 / kSecondsPerDay); d * kSecondsPerDay <= t1; d += 1.0) {
            const float tx = to_x(d * kSecondsPerDay);
            draw->AddLine(ImVec2(tx, bottom), ImVec2(tx, bottom + 0.35f * em), gray(170, 0.45f * alpha), 1.0f);
        }
    }

    // The curve, and the clock's moment on it.
    std::vector<ImVec2> points(static_cast<size_t>(count));
    for (int i = 0; i < count; ++i) {
        points[static_cast<size_t>(i)] = ImVec2(to_x(samples.time(i)), to_y(mags[static_cast<size_t>(i)]));
    }
    draw->AddPolyline(points.data(), count, IM_COL32(255, 222, 170, static_cast<int>(220.0f * alpha)), 0,
                      1.4f * ImGui::GetStyle().FontScaleDpi);
    const float nx = to_x(t_tdb);
    draw->AddLine(ImVec2(nx, top), ImVec2(nx, bottom), gray(220, 0.3f * alpha), 1.0f);
    draw->AddCircleFilled(ImVec2(nx, to_y(now_mag)), 0.25f * em, IM_COL32(255, 240, 215, static_cast<int>(255.0f * alpha)));
}

} // namespace astraxis
