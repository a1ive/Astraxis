#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// The scene files in <asset_dir>/scenes, sorted by name.
std::vector<std::filesystem::path> list_scene_files(const std::filesystem::path& asset_dir);

// Index of the scene whose file stem is `name`, or -1.
int find_scene(const std::vector<std::filesystem::path>& files, const std::string& name);

} // namespace astraxis
