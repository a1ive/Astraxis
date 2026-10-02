#pragma once

#include "scene/scene.hpp"

#include <filesystem>
#include <string>
#include <string_view>

namespace astraxis {

// Loads a scene description (TOML, see assets/scenes/jupiter.toml for the format).
// On failure returns false and sets `error`; `out` is left in an unspecified state.
bool load_scene_file(const std::filesystem::path& path, Scene& out, std::string* error);
// `asset_root` resolves relative asset paths (e.g. ephemeris files).
bool load_scene_string(std::string_view text, std::string_view source_name, Scene& out, std::string* error,
                       const std::filesystem::path& asset_root = {});

} // namespace astraxis
