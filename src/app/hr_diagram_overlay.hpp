#pragma once

namespace astraxis {

class Scene;
struct Body;

// Whether the info panel shows an H-R diagram for `body`: a star with its own
// temperature (not a neutron star, nor one left off with hr_diagram = false).
bool shows_hr_diagram(const Body& body);

// A Hertzsprung-Russell diagram at the cursor of the current ImGui window,
// `width` pixels wide: the main-sequence and white-dwarf bands
// (scene/hr_diagram.hpp), the Sun's place for reference, the scene's stars in
// their colors, and the star `focus` marked (with its cycle, if it pulsates).
void draw_hr_diagram(const Scene& scene, int focus, float width);

} // namespace astraxis
