#pragma once

struct ImFont;

namespace astraxis {

class LabelLayout;
class Scene;
struct OutputView;

// Lays out the labels of one output (sized as the current ImGui frame's
// DisplaySize) and draws them into the background draw list: scene markers,
// body dots and body names, with `alpha` (the picture's fade). Shared by all hosts.
// Below a spacecraft's name: its burn's delta-v during an event with a ghost
// (Scene::burn_delta_v), in `readout_font` (which has a Greek delta) if given.
void draw_labels(LabelLayout& labels, const Scene& scene, const OutputView& view, int focus, float alpha = 1.0f,
                 ImFont* readout_font = nullptr);

} // namespace astraxis
