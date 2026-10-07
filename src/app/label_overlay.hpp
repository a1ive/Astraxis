#pragma once

#include <filesystem>
#include <string>

struct ImFont;

namespace astraxis {

class LabelLayout;
class Scene;
struct OutputView;

// Lays out the labels of one output (sized as the current ImGui frame's
// DisplaySize) and draws them into the background draw list: scene markers,
// body dots and body names. Shared by all hosts.
void draw_labels(LabelLayout& labels, const Scene& scene, const OutputView& view, int focus);

// Adds the default font and then the title font (assets/fonts) to the current
// ImGui context, so the UI keeps the default. Returns the title font, or
// nullptr if it could not be read (the title then uses the default).
ImFont* load_fonts(const std::filesystem::path& asset_dir);

// The scene's name as a film's title card: large and letter-spaced, fading in
// and out over the first seconds after the scene is loaded (`age`, real
// seconds; Simulation::scene_age). Draws into the background draw list.
void draw_scene_title(ImFont* font, const std::string& title, double age);
// False once the title has faded out.
bool scene_title_visible(double age);

} // namespace astraxis
