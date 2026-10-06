#pragma once

struct ImDrawList;

namespace astraxis {

class LabelLayout;
class Scene;
struct OutputView;

// Lays out the labels of one output (sized as the current ImGui frame's
// DisplaySize) and draws them into the background draw list: scene markers,
// body dots and body names. Shared by all hosts.
void draw_labels(LabelLayout& labels, const Scene& scene, const OutputView& view, int focus);

} // namespace astraxis
