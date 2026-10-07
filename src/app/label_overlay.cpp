#include "app/label_overlay.hpp"

#include "scene/label_layout.hpp"
#include "scene/scene.hpp"
#include "scene/simulation.hpp"

#include <glm/glm.hpp>
#include <imgui.h>

namespace astraxis {

namespace {

ImU32 to_imgui_color(const glm::vec3& c, float alpha)
{
    return ImGui::ColorConvertFloat4ToU32(ImVec4(c.r, c.g, c.b, alpha));
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

} // namespace astraxis
