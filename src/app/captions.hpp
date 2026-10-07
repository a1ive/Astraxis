#pragma once

#include <filesystem>
#include <string>

struct ImFont;
struct ImVec2;

namespace astraxis {

struct SceneEvent;

// The fonts of the title card and the event captions (assets/fonts). Either
// may be nullptr if its file could not be read; the default font is used then.
struct CaptionFonts {
    ImFont* light = nullptr; // titles and headings
    ImFont* book = nullptr;  // caption text
};

// Adds the default font and then the caption fonts to the current ImGui
// context, so the UI keeps the default.
CaptionFonts load_fonts(const std::filesystem::path& asset_dir);

// The scene's name as a film's title card: large and letter-spaced, fading in
// and out over the first seconds after the scene is loaded (`age`, real
// seconds; Simulation::scene_age). Draws into the background draw list.
void draw_scene_title(const CaptionFonts& fonts, const std::string& title, double age);
// False once the title has faded out.
bool scene_title_visible(double age);

// When an event's caption is shown, in real seconds after the jump: it waits
// for the camera's flight, fades in, stays for a reading time that grows with
// its length, and fades out.
struct CaptionTimes {
    double fade_in = 0.0;  // starts fading in
    double shown = 0.0;    // fully shown
    double fade_out = 0.0; // starts fading out
    double gone = 0.0;
};
CaptionTimes caption_times(const std::string& caption);

// A museum placard for an event with a caption: the event's name as a heading
// (a trailing parenthesis, such as the date, set apart below it) and the
// caption, on a dark card `width` wide whose bottom edge is centred on
// `bottom_center`. `age` is in real seconds since the jump. Draws into the
// background draw list; nothing if the event has no caption.
void draw_event_caption(const CaptionFonts& fonts, const SceneEvent& event, double age, const ImVec2& bottom_center,
                        float width);

} // namespace astraxis
