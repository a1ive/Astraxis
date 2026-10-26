#pragma once

#include "app/captions.hpp"

#include <span>

struct ImVec2;

namespace astraxis {

class Scene;
class LightCurveSamples;

// Size of the light curve plot, in pixels (for the current ImGui frame's font size).
ImVec2 light_curve_plot_size();

// Where the plot goes on an output of `display` size: in the corner at
// `bottom_right` (minus `margin`), unless it would cover the caption, centred
// at the bottom and `caption_width` wide, or one of the `avoid` rectangles;
// then in the top right corner.
ScreenRect light_curve_plot_rect(const ImVec2& display, float bottom, float margin, float caption_width,
                                 std::span<const ScreenRect> avoid = {});

// A recorded light curve (Scene::observed_light_curve) takes its place: the
// whole record, apparent magnitudes, the clock's moment marked.
//
// The scene's light curve (Scene::light_curve) as a small plot: magnitudes
// below the unhidden brightness (brighter up) over the samples' span, with
// the clock's moment marked. Draws into the background draw list, with
// `alpha` (the picture's fade); nothing if the scene has no light curve. At
// `rate` (simulated seconds per second; 0 when paused) the curve would be
// replaced within a few seconds, so it fades out when the clock runs faster.
void draw_light_curve(const CaptionFonts& fonts, const Scene& scene, const LightCurveSamples& samples, double t_tdb,
                      double rate, const ScreenRect& rect, float alpha);

} // namespace astraxis
