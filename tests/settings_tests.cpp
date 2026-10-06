// Tests for the settings file (config.toml) and the scene list.

#include "test_util.hpp"

#include "app/settings.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace astraxis;

namespace {

// config.toml: [view] shared, mode sections override it, saving keeps the rest.
void test_settings()
{
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "astraxis_settings_test";
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
    std::filesystem::create_directories(dir);
    const std::filesystem::path path = dir / "config.toml";

    Settings s;
    s.scene = "jupiter";
    check(load_settings(path, SettingsSection::Window, s) && s.scene == "solar_system" && s.view.orbits,
          "missing config gives the defaults");

    {
        std::ofstream file(path, std::ios::binary);
        file << "# by hand\n"
                "[view]\nexposure = 2\nlabels = false\nstars = 9.0\nunknown = 1\n"
                "[wallpaper]\nlabels = true\nscene = \"sgr_a\"\nfps = 30\n"
                "[window]\nscene = \"jupiter\"\nfullscreen = true\nexposure = 0.5\n";
    }
    Settings window;
    Settings wallpaper;
    check(load_settings(path, SettingsSection::Window, window), "config loads");
    check(load_settings(path, SettingsSection::Wallpaper, wallpaper), "config loads (wallpaper)");
    check(window.scene == "jupiter" && window.fullscreen && !window.labels, "[window] keys over [view]");
    check(window.view.post.exposure == 0.5f, "[window] overrides [view] exposure", window.view.post.exposure);
    check(window.view.star_brightness == 2.0f, "out-of-range value is clamped", window.view.star_brightness);
    check(wallpaper.scene == "sgr_a" && wallpaper.fps == 30 && wallpaper.labels, "[wallpaper] keys");
    check(wallpaper.view.post.exposure == 2.0f, "[wallpaper] inherits [view]", wallpaper.view.post.exposure);
    check(!wallpaper.fullscreen, "[wallpaper] ignores window-only keys");

    // Saving [window] writes the view keys to [view] and drops its overrides there.
    window.view.post.exposure = 1.25f;
    window.view.post.bloom_strength = 0.04f;
    window.view.orbits = false;
    window.display = 2;
    window.display_name = "Test Monitor";
    check(save_settings(path, SettingsSection::Window, window), "config saves");
    Settings again;
    Settings wallpaper_again;
    load_settings(path, SettingsSection::Window, again);
    load_settings(path, SettingsSection::Wallpaper, wallpaper_again);
    check(again.view.post.exposure == 1.25f && !again.view.orbits && again.display == 2 &&
              again.display_name == "Test Monitor" && again.scene == "jupiter" && again.fullscreen,
          "saved [window] settings read back");
    check(wallpaper_again.scene == "sgr_a" && wallpaper_again.fps == 30 && wallpaper_again.labels &&
              wallpaper_again.view.post.exposure == 1.25f,
          "saving [window] keeps [wallpaper]");
    {
        std::ifstream file(path, std::ios::binary);
        const std::string text((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        check(text.find("unknown = 1") != std::string::npos, "saving keeps unknown keys");
        check(text.find("bloom = 0.04\n") != std::string::npos, "floats are written tidily");
    }

    // A file that does not parse is neither read nor overwritten.
    {
        std::ofstream file(path, std::ios::binary);
        file << "[view\nexposure = 3\n";
    }
    Settings broken;
    std::string error;
    check(!load_settings(path, SettingsSection::Window, broken, &error) && !error.empty() &&
              broken.view.post.exposure == 1.0f,
          "unparsable config is reported, defaults used");
    check(!save_settings(path, SettingsSection::Window, window), "unparsable config is not overwritten");

    // The settings dialog: view keys go to the section, only where they differ from [view].
    {
        std::ofstream file(path, std::ios::binary);
        file << "[view]\nexposure = 2\nlabels = false\n[screensaver]\nstars = 0.5\norbits = false\n";
    }
    Settings saver;
    load_settings(path, SettingsSection::Screensaver, saver);
    saver.view.post.exposure = 2.0f;     // as in [view]: no override
    saver.view.star_brightness = 1.0f;   // back to [view] (the default): the override goes
    saver.view.post.bloom_strength = 0.1f; // differs: an override
    saver.view.orbits = false;            // still differs: stays
    saver.scene = "saturn";
    check(save_settings(path, SettingsSection::Screensaver, saver, nullptr, ViewScope::Overrides),
          "dialog settings save");
    {
        std::ifstream file(path, std::ios::binary);
        const std::string text((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        const size_t at = text.find("\n[screensaver]"); // not the header comment's mention
        const size_t end = at == std::string::npos ? at : text.find("\n[", at + 1);
        const std::string section = at == std::string::npos ? "" : text.substr(at, end - at);
        check(section.find("exposure") == std::string::npos && section.find("stars") == std::string::npos &&
                  section.find("labels") == std::string::npos,
              "no overrides where the dialog matches [view]");
        check(section.find("bloom = 0.1") != std::string::npos && section.find("orbits = false") != std::string::npos,
              "overrides where the dialog differs from [view]");
        check(text.find("exposure = 2") != std::string::npos, "[view] untouched by the dialog");
    }
    Settings saver_again;
    load_settings(path, SettingsSection::Screensaver, saver_again);
    check(saver_again.scene == "saturn" && saver_again.view.post.bloom_strength == 0.1f &&
              saver_again.view.post.exposure == 2.0f && !saver_again.labels && !saver_again.view.orbits &&
              saver_again.view.star_brightness == 1.0f,
          "dialog settings read back");

    // Power keys: [wallpaper] / [screensaver] only; unknown battery names are ignored.
    {
        std::ofstream file(path, std::ios::binary);
        file << "[wallpaper]\nrender_scale = 0.1\nbattery = \"pause\"\npause_covered = false\n"
                "[screensaver]\nbattery = \"sometimes\"\n[window]\nrender_scale = 0.5\n";
    }
    Settings wall;
    Settings saver2;
    Settings win;
    load_settings(path, SettingsSection::Wallpaper, wall);
    load_settings(path, SettingsSection::Screensaver, saver2);
    load_settings(path, SettingsSection::Window, win);
    check(wall.render_scale == 0.25f && wall.battery == BatteryPolicy::Pause && !wall.pause_covered,
          "power keys read (render scale clamped)");
    check(saver2.battery == BatteryPolicy::Limit, "unknown battery policy keeps the default");
    check(win.render_scale == 1.0f, "[window] ignores render_scale");
    wall.battery = BatteryPolicy::Run;
    wall.render_scale = 0.5f;
    save_settings(path, SettingsSection::Wallpaper, wall, nullptr, ViewScope::Overrides);
    Settings wall_again;
    load_settings(path, SettingsSection::Wallpaper, wall_again);
    check(wall_again.battery == BatteryPolicy::Run && wall_again.render_scale == 0.5f && !wall_again.pause_covered,
          "power keys saved");

    // Scene names for the dialog.
    const std::vector<SceneEntry> scenes = list_scenes(ASTRAXIS_ASSET_DIR);
    const auto earth_moon = std::find_if(scenes.begin(), scenes.end(),
                                         [](const SceneEntry& s) { return s.stem == "earth_moon"; });
    check(scenes.size() > 10 && earth_moon != scenes.end() && earth_moon->name == "Earth-Moon",
          "scene list with names");

    std::filesystem::remove_all(dir, ec);
}

} // namespace

void run_settings_tests()
{
    test_settings();
}
