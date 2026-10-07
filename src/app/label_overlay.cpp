#include "app/label_overlay.hpp"

#include "scene/label_layout.hpp"
#include "scene/scene.hpp"
#include "scene/simulation.hpp"

#include <glm/glm.hpp>
#include <imgui.h>

#include <algorithm>
#include <cfloat>
#include <system_error>
#include <utility>
#include <vector>

namespace astraxis {

namespace {

// Title card timing, real seconds after the scene is loaded: a short wait for
// the first frames, then fade in, hold, fade out.
constexpr double kTitleDelay = 0.5;
constexpr double kTitleFadeIn = 1.5;
constexpr double kTitleHold = 2.5;
constexpr double kTitleFadeOut = 2.0;
constexpr double kTitleEnd = kTitleDelay + kTitleFadeIn + kTitleHold + kTitleFadeOut;

constexpr float kTitleHeight = 0.06f;        // cap-to-descender size, fraction of the output height
constexpr float kTitleTrackingStart = 0.20f; // letter spacing in ems; widens slowly while shown
constexpr float kTitleTrackingEnd = 0.28f;
constexpr float kTitleCenterY = 0.70f;       // in the lower third, clear of the focus at the centre
constexpr float kTitleMaxWidth = 0.8f;       // fraction of the output width

ImU32 to_imgui_color(const glm::vec3& c, float alpha)
{
    return ImGui::ColorConvertFloat4ToU32(ImVec4(c.r, c.g, c.b, alpha));
}

double smoothstep(double edge0, double edge1, double x)
{
    const double t = std::clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

} // namespace

void draw_labels(LabelLayout& labels, const Scene& scene, const OutputView& view, int focus)
{
    const ImGuiIO& io = ImGui::GetIO();
    const float font_size = ImGui::GetFontSize();
    labels.update(scene, view, focus, glm::vec2(io.DisplaySize.x, io.DisplaySize.y), font_size,
                  [](const std::string& text) {
                      const ImVec2 size = ImGui::CalcTextSize(text.c_str());
                      return glm::vec2(size.x, size.y);
                  },
                  io.DeltaTime);

    ImDrawList* draw = ImGui::GetBackgroundDrawList();

    // Lagrange points and other markers: a small diamond and a name.
    for (const LabelLayout::SceneMarker& marker : labels.markers()) {
        const ImVec2 s(marker.position.x, marker.position.y);
        const float r = 4.0f;
        draw->AddQuad(ImVec2(s.x, s.y - r), ImVec2(s.x + r, s.y), ImVec2(s.x, s.y + r), ImVec2(s.x - r, s.y),
                      IM_COL32(180, 200, 255, 170), 1.2f);
        draw->AddText(ImVec2(s.x + r + 3.0f, s.y - font_size * 0.5f), IM_COL32(180, 200, 255, 170),
                      scene.markers[static_cast<size_t>(marker.marker)].name.c_str());
    }

    const float dot_radius = LabelLayout::kMarkerRadius;
    for (const LabelLayout::BodyMark& mark : labels.bodies()) {
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
}

ImFont* load_fonts(const std::filesystem::path& asset_dir)
{
    ImFontAtlas* atlas = ImGui::GetIO().Fonts;
    atlas->AddFontDefault();
    const std::filesystem::path file = asset_dir / "fonts" / "Jost-300-Light.ttf";
    std::error_code ec;
    if (!std::filesystem::is_regular_file(file, ec)) {
        return nullptr; // ImGui asserts on a missing file
    }
    // ImGui opens it as UTF-8. The atlas is dynamic (ImGui 1.92): glyphs are
    // baked at whatever size they are drawn.
    const std::u8string path = file.u8string();
    return atlas->AddFontFromFileTTF(reinterpret_cast<const char*>(path.c_str()), 48.0f);
}

bool scene_title_visible(double age)
{
    return age < kTitleEnd;
}

void draw_scene_title(ImFont* font, const std::string& title, double age)
{
    if (title.empty() || !scene_title_visible(age)) {
        return;
    }
    const double fade_in = smoothstep(kTitleDelay, kTitleDelay + kTitleFadeIn, age);
    const double fade_out = 1.0 - smoothstep(kTitleEnd - kTitleFadeOut, kTitleEnd, age);
    const float alpha = static_cast<float>(std::min(fade_in, fade_out));
    if (alpha <= 0.0f) {
        return;
    }

    // Upper case for ASCII letters only: the bytes of multi-byte UTF-8
    // sequences are left alone (and a Greek letter in a star's name stays
    // lower case, as it should).
    std::string text = title;
    for (char& c : text) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    // The title's code points as byte ranges: a UTF-8 sequence runs up to the
    // next byte that is not a continuation byte (10xxxxxx).
    std::vector<std::pair<const char*, const char*>> glyphs;
    for (size_t i = 0; i < text.size();) {
        size_t next = i + 1;
        while (next < text.size() && (static_cast<unsigned char>(text[next]) & 0xC0) == 0x80) {
            ++next;
        }
        glyphs.emplace_back(text.data() + i, text.data() + next);
        i = next;
    }

    const ImGuiIO& io = ImGui::GetIO();
    if (!font) {
        font = ImGui::GetFont();
    }
    const float size = kTitleHeight * io.DisplaySize.y;
    const float progress = static_cast<float>(std::clamp((age - kTitleDelay) / (kTitleEnd - kTitleDelay), 0.0, 1.0));
    const float tracking = size * (kTitleTrackingStart + (kTitleTrackingEnd - kTitleTrackingStart) * progress);

    // ImGui has no letter spacing: lay out one code point at a time. The last
    // one's spacing is not counted, so the title stays centred.
    float width = 0.0f;
    for (const auto& [begin, end] : glyphs) {
        width += font->CalcTextSizeA(size, FLT_MAX, 0.0f, begin, end).x;
    }
    width += tracking * static_cast<float>(glyphs.size() - 1);
    // A long title in a narrow output: shrink to fit.
    const float fit = std::min(1.0f, kTitleMaxWidth * io.DisplaySize.x / std::max(width, 1.0f));

    const float scaled = size * fit;
    const float left = 0.5f * (io.DisplaySize.x - width * fit);
    const float top = kTitleCenterY * io.DisplaySize.y - 0.5f * scaled;
    ImDrawList* draw = ImGui::GetBackgroundDrawList();

    // A faint dark halo (the text drawn around its place) keeps the thin
    // strokes readable over bright orbits and stars.
    const float r = std::max(1.0f, 0.04f * scaled);
    const ImVec2 halo[] = {{-r, 0.0f}, {r, 0.0f}, {0.0f, -r}, {0.0f, r},
                           {-0.7f * r, -0.7f * r}, {0.7f * r, -0.7f * r}, {-0.7f * r, 0.7f * r}, {0.7f * r, 0.7f * r}};
    const ImU32 shadow = IM_COL32(0, 0, 0, static_cast<int>(70.0f * alpha));
    const ImU32 color = IM_COL32(236, 240, 246, static_cast<int>(235.0f * alpha));
    for (int pass = 0; pass < 2; ++pass) {
        float x = left;
        for (const auto& [begin, end] : glyphs) {
            if (pass == 0) {
                for (const ImVec2& d : halo) {
                    draw->AddText(font, scaled, ImVec2(x + d.x, top + d.y), shadow, begin, end);
                }
            } else {
                draw->AddText(font, scaled, ImVec2(x, top), color, begin, end);
            }
            x += font->CalcTextSizeA(scaled, FLT_MAX, 0.0f, begin, end).x + tracking * fit;
        }
    }
}

} // namespace astraxis
