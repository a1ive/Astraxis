#include "app/settings.hpp"

#include <toml++/toml.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <system_error>

namespace astraxis {

namespace {

constexpr const char* kHeader =
    "# Astraxis settings, rewritten by the program (comments are not kept).\n"
    "# [view] is shared by all modes; [window], [wallpaper] and [screensaver]\n"
    "# hold each mode's own keys and may override any [view] key.\n\n";

const char* section_name(SettingsSection section)
{
    switch (section) {
    case SettingsSection::Window:
        return "window";
    case SettingsSection::Wallpaper:
        return "wallpaper";
    case SettingsSection::Screensaver:
        return "screensaver";
    }
    return "window";
}

// Calls `f(key, field[, min, max])` for every view key. `S` is Settings or
// const Settings, so that one list serves reading and writing.
template <typename S, typename F>
void visit_view_keys(S& s, F&& f)
{
    f("orbits", s.view.orbits);
    f("labels", s.labels);
    f("belts", s.view.belts);
    f("atmospheres", s.view.atmospheres);
    f("plumes", s.view.plumes);
    f("comets", s.view.comets);
    // Ranges as on the control panel's sliders.
    f("stars", s.view.star_brightness, 0.0f, 2.0f);
    f("line_width", s.view.line_width, 0.5f, 4.0f);
    f("exposure", s.view.post.exposure, 0.1f, 8.0f);
    f("bloom", s.view.post.bloom_strength, 0.0f, 0.2f);
}

template <typename S, typename F>
void visit_mode_keys(S& s, SettingsSection section, F&& f)
{
    f("scene", s.scene);
    f("display", s.display, 0, 64);
    f("display_name", s.display_name);
    f("fps", s.fps, 0, 1000);
    if (section == SettingsSection::Window) {
        f("fullscreen", s.fullscreen);
        f("auto_tour", s.auto_tour);
        f("info", s.info);
    }
}

// Reads the keys present in a table; values of the wrong type are ignored.
struct Reader {
    const toml::table& table;

    void operator()(const char* key, bool& value) const
    {
        if (const auto v = table[key].value<bool>()) {
            value = *v;
        }
    }
    void operator()(const char* key, float& value, float lo, float hi) const
    {
        if (const auto v = table[key].value<double>()) {
            value = std::clamp(static_cast<float>(*v), lo, hi);
        }
    }
    void operator()(const char* key, int& value, int lo, int hi) const
    {
        if (const auto v = table[key].value<int64_t>()) {
            value = static_cast<int>(std::clamp<int64_t>(*v, lo, hi));
        }
    }
    void operator()(const char* key, std::string& value) const
    {
        if (const auto v = table[key].value<std::string>()) {
            value = *v;
        }
    }
};

struct Writer {
    toml::table& table;

    void operator()(const char* key, bool value) const { table.insert_or_assign(key, value); }
    void operator()(const char* key, float value, float, float) const
    {
        // Through 6 significant digits, so that 0.04f is written as 0.04.
        char text[32];
        std::snprintf(text, sizeof(text), "%.6g", static_cast<double>(value));
        table.insert_or_assign(key, std::strtod(text, nullptr));
    }
    void operator()(const char* key, int value, int, int) const
    {
        table.insert_or_assign(key, static_cast<int64_t>(value));
    }
    void operator()(const char* key, const std::string& value) const
    {
        if (value.empty()) {
            table.erase(key);
        } else {
            table.insert_or_assign(key, value);
        }
    }
};

struct Eraser {
    toml::table& table;

    template <typename... Args>
    void operator()(const char* key, Args&&...) const
    {
        table.erase(key);
    }
};

bool parse_file(const std::filesystem::path& path, toml::table& out, std::string* error)
{
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        if (error) {
            *error = "cannot open " + path.string();
        }
        return false;
    }
    std::ostringstream text;
    text << file.rdbuf();
    const std::string source = path.string();
    try {
        out = toml::parse(text.str(), source);
        return true;
    } catch (const toml::parse_error& e) {
        if (error) {
            std::ostringstream msg;
            msg << source << ":" << e.source().begin.line << ": " << e.description();
            *error = msg.str();
        }
    }
    return false;
}

toml::table& child_table(toml::table& root, const char* key)
{
    if (!root[key].is_table()) {
        root.insert_or_assign(key, toml::table{});
    }
    return *root[key].as_table();
}

} // namespace

bool load_settings(const std::filesystem::path& path, SettingsSection section, Settings& out, std::string* error)
{
    out = Settings{};
    std::error_code ec;
    if (!std::filesystem::exists(path, ec)) {
        return true;
    }
    toml::table root;
    if (!parse_file(path, root, error)) {
        return false;
    }
    if (const toml::table* view = root["view"].as_table()) {
        visit_view_keys(out, Reader{*view});
    }
    if (const toml::table* mode = root[section_name(section)].as_table()) {
        visit_view_keys(out, Reader{*mode});
        visit_mode_keys(out, section, Reader{*mode});
    }
    return true;
}

bool save_settings(const std::filesystem::path& path, SettingsSection section, const Settings& settings,
                   std::string* error)
{
    toml::table root;
    std::error_code ec;
    if (std::filesystem::exists(path, ec) && !parse_file(path, root, error)) {
        return false;
    }

    toml::table& view = child_table(root, "view");
    toml::table& mode = child_table(root, section_name(section));
    visit_view_keys(settings, Writer{view});
    visit_view_keys(settings, Eraser{mode});
    visit_mode_keys(settings, section, Writer{mode});

    // Write a temporary file and rename it over the old one, so that a failed
    // write never leaves a truncated config.
    std::filesystem::path temp = path;
    temp += ".tmp";
    {
        std::ofstream file(temp, std::ios::binary | std::ios::trunc);
        file << kHeader << root << '\n';
        if (!file) {
            if (error) {
                *error = "cannot write " + temp.string();
            }
            return false;
        }
    }
    std::filesystem::rename(temp, path, ec);
    if (ec) {
        if (error) {
            *error = "cannot replace " + path.string() + ": " + ec.message();
        }
        std::filesystem::remove(temp, ec);
        return false;
    }
    return true;
}

} // namespace astraxis
