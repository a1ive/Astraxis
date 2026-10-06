#include "app/scene_list.hpp"

#include <algorithm>
#include <system_error>

namespace astraxis {

std::vector<std::filesystem::path> list_scene_files(const std::filesystem::path& asset_dir)
{
    std::vector<std::filesystem::path> files;
    std::error_code ec;
    for (const auto& entry : std::filesystem::directory_iterator(asset_dir / "scenes", ec)) {
        if (entry.path().extension() == ".toml") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

int find_scene(const std::vector<std::filesystem::path>& files, const std::string& name)
{
    for (size_t i = 0; i < files.size(); ++i) {
        if (files[i].stem() == name) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

} // namespace astraxis
