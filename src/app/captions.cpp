#include "app/captions.hpp"

#include "scene/scene.hpp"

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

// Event caption timing, real seconds after the jump: the camera flies for 2 s
// (Simulation::jump_to_event); the caption stays for a leisurely reading time.
constexpr double kCaptionDelay = 1.5;
constexpr double kCaptionFadeIn = 0.8;
constexpr double kCaptionFadeOut = 1.5;
constexpr double kCaptionWordsPerSecond = 3.0;
constexpr double kCaptionExtraSeconds = 4.0;
constexpr double kCaptionMinSeconds = 8.0;
constexpr double kCaptionMaxSeconds = 45.0;

// Caption sizes, relative to the UI font size (which includes the DPI scale).
constexpr float kHeadingSize = 1.9f;
constexpr float kSubheadingSize = 1.1f;
constexpr float kBodySize = 1.3f;
constexpr float kLineSpacing = 1.35f; // body line advance, in body sizes
constexpr float kPadding = 1.2f;

constexpr ImU32 kAccent = IM_COL32(200, 170, 110, 255); // the rule and the subheading

using Span = std::pair<const char*, const char*>;

double smoothstep(double edge0, double edge1, double x)
{
    const double t = std::clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

ImFont* add_font(const std::filesystem::path& file)
{
    std::error_code ec;
    if (!std::filesystem::is_regular_file(file, ec)) {
        return nullptr; // ImGui asserts on a missing file
    }
    // ImGui opens it as UTF-8. The atlas is dynamic (ImGui 1.92): glyphs are
    // baked at whatever size they are drawn.
    const std::u8string path = file.u8string();
    return ImGui::GetIO().Fonts->AddFontFromFileTTF(reinterpret_cast<const char*>(path.c_str()), 24.0f);
}

bool is_continuation(char c)
{
    return (static_cast<unsigned char>(c) & 0xC0) == 0x80; // 10xxxxxx
}

// The code points of UTF-8 text as byte ranges.
std::vector<Span> code_points(const std::string& text)
{
    std::vector<Span> out;
    for (size_t i = 0; i < text.size();) {
        size_t next = i + 1;
        while (next < text.size() && is_continuation(text[next])) {
            ++next;
        }
        out.emplace_back(text.data() + i, text.data() + next);
        i = next;
    }
    return out;
}

// Text broken into lines no wider than `width` at `size`; a newline breaks a line.
std::vector<Span> wrap_lines(ImFont* font, float size, const std::string& text, float width)
{
    std::vector<Span> lines;
    const char* s = text.data();
    const char* const end = s + text.size();
    while (s < end) {
        const char* e = font->CalcWordWrapPosition(size, s, end, width);
        if (e == s && *s != '\n') { // a glyph wider than the line: take it anyway
            ++e;
            while (e < end && is_continuation(*e)) {
                ++e;
            }
        }
        const char* trimmed = e;
        while (trimmed > s && (trimmed[-1] == ' ' || trimmed[-1] == '\n')) {
            --trimmed;
        }
        lines.emplace_back(s, trimmed);
        s = e;
        if (s < end && *s == '\n') {
            ++s;
        }
        while (s < end && *s == ' ') {
            ++s;
        }
    }
    return lines;
}

// The text drawn around its place: a faint dark halo that keeps thin strokes
// readable over bright orbits and stars.
void add_halo(ImDrawList* draw, ImFont* font, float size, const ImVec2& pos, float r, ImU32 color, const char* begin,
              const char* end)
{
    const ImVec2 offsets[] = {{-r, 0.0f}, {r, 0.0f}, {0.0f, -r}, {0.0f, r},
                              {-0.7f * r, -0.7f * r}, {0.7f * r, -0.7f * r}, {-0.7f * r, 0.7f * r}, {0.7f * r, 0.7f * r}};
    for (const ImVec2& d : offsets) {
        draw->AddText(font, size, ImVec2(pos.x + d.x, pos.y + d.y), color, begin, end);
    }
}

ImU32 with_alpha(ImU32 color, float alpha)
{
    const auto a = static_cast<ImU32>(255.0f * std::clamp(alpha, 0.0f, 1.0f));
    return (color & ~IM_COL32_A_MASK) | (a << IM_COL32_A_SHIFT);
}

} // namespace

CaptionFonts load_fonts(const std::filesystem::path& asset_dir)
{
    ImGui::GetIO().Fonts->AddFontDefault();
    CaptionFonts fonts;
    fonts.light = add_font(asset_dir / "fonts" / "Jost-300-Light.ttf");
    fonts.book = add_font(asset_dir / "fonts" / "Jost-400-Book.ttf");
    return fonts;
}

bool scene_title_visible(double age)
{
    return age < kTitleEnd;
}

void draw_scene_title(const CaptionFonts& fonts, const std::string& title, double age)
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
    const std::vector<Span> glyphs = code_points(text);

    const ImGuiIO& io = ImGui::GetIO();
    ImFont* font = fonts.light ? fonts.light : ImGui::GetFont();
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

    const float r = std::max(1.0f, 0.04f * scaled);
    const ImU32 shadow = IM_COL32(0, 0, 0, static_cast<int>(70.0f * alpha));
    const ImU32 color = IM_COL32(236, 240, 246, static_cast<int>(235.0f * alpha));
    for (int pass = 0; pass < 2; ++pass) {
        float x = left;
        for (const auto& [begin, end] : glyphs) {
            if (pass == 0) {
                add_halo(draw, font, scaled, ImVec2(x, top), r, shadow, begin, end);
            } else {
                draw->AddText(font, scaled, ImVec2(x, top), color, begin, end);
            }
            x += font->CalcTextSizeA(scaled, FLT_MAX, 0.0f, begin, end).x + tracking * fit;
        }
    }
}

CaptionTimes caption_times(const std::string& caption)
{
    size_t words = 0;
    bool in_word = false;
    for (const char c : caption) {
        const bool space = c == ' ' || c == '\n';
        if (!space && !in_word) {
            ++words;
        }
        in_word = !space;
    }
    const double reading = std::clamp(static_cast<double>(words) / kCaptionWordsPerSecond + kCaptionExtraSeconds,
                                      kCaptionMinSeconds, kCaptionMaxSeconds);
    CaptionTimes t;
    t.fade_in = kCaptionDelay;
    t.shown = t.fade_in + kCaptionFadeIn;
    t.fade_out = t.shown + reading;
    t.gone = t.fade_out + kCaptionFadeOut;
    return t;
}

void draw_event_caption(const CaptionFonts& fonts, const SceneEvent& event, double age, const ImVec2& bottom_center,
                        float width)
{
    if (event.caption.empty()) {
        return;
    }
    const CaptionTimes t = caption_times(event.caption);
    const float alpha = static_cast<float>(
        std::min(smoothstep(t.fade_in, t.shown, age), 1.0 - smoothstep(t.fade_out, t.gone, age)));
    if (alpha <= 0.0f) {
        return;
    }

    // "Voyager 2 at Neptune (CA 1989-08-25 04:03, 33,661 km)": the parenthesis
    // goes on a line of its own under the heading.
    std::string heading = event.name;
    std::string subheading;
    const size_t open = heading.rfind(" (");
    if (open != std::string::npos && heading.back() == ')') {
        subheading = heading.substr(open + 2, heading.size() - open - 3);
        heading.resize(open);
    }

    ImFont* light = fonts.light ? fonts.light : ImGui::GetFont();
    ImFont* book = fonts.book ? fonts.book : ImGui::GetFont();
    const float em = ImGui::GetFontSize();
    const float heading_size = kHeadingSize * em;
    const float sub_size = kSubheadingSize * em;
    const float body_size = kBodySize * em;
    const float body_line = kLineSpacing * body_size;
    const float pad = kPadding * em;
    const float inner = width - 2.0f * pad;

    const std::vector<Span> heading_lines = wrap_lines(light, heading_size, heading, inner);
    const std::vector<Span> body_lines = wrap_lines(book, body_size, event.caption, inner);
    float height = 2.0f * pad + static_cast<float>(heading_lines.size()) * heading_size;
    if (!subheading.empty()) {
        height += 0.3f * em + sub_size;
    }
    height += 0.8f * em + static_cast<float>(body_lines.size()) * body_line;

    ImDrawList* draw = ImGui::GetBackgroundDrawList();
    const ImVec2 min(bottom_center.x - 0.5f * width, bottom_center.y - height);
    const ImVec2 max(bottom_center.x + 0.5f * width, bottom_center.y);
    draw->AddRectFilled(min, max, IM_COL32(6, 8, 12, static_cast<int>(170.0f * alpha)), 0.5f * em);
    // A thin rule down the left edge, as on a museum label.
    draw->AddRectFilled(ImVec2(min.x + 0.5f * pad, min.y + pad), ImVec2(min.x + 0.5f * pad + 0.12f * em, max.y - pad),
                        with_alpha(kAccent, 0.8f * alpha));

    const ImU32 text = IM_COL32(236, 240, 246, static_cast<int>(240.0f * alpha));
    float y = min.y + pad;
    for (const auto& [begin, end] : heading_lines) {
        draw->AddText(light, heading_size, ImVec2(min.x + pad, y), text, begin, end);
        y += heading_size;
    }
    if (!subheading.empty()) {
        y += 0.3f * em;
        draw->AddText(book, sub_size, ImVec2(min.x + pad, y), with_alpha(kAccent, 0.9f * alpha), subheading.c_str());
        y += sub_size;
    }
    y += 0.8f * em;
    for (const auto& [begin, end] : body_lines) {
        draw->AddText(book, body_size, ImVec2(min.x + pad, y), text, begin, end);
        y += body_line;
    }
}

} // namespace astraxis
