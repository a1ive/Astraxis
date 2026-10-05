#include "platform/paths.hpp"

#include <SDL3/SDL_filesystem.h>

#include <string>
#include <system_error>

namespace astraxis {

std::filesystem::path path_from_utf8(const char* s)
{
    const std::string str(s);
    return std::filesystem::path(std::u8string(str.begin(), str.end()));
}

std::filesystem::path asset_directory()
{
    std::error_code ec;
    if (const char* base = SDL_GetBasePath()) {
        const std::filesystem::path candidate = path_from_utf8(base) / "assets";
        if (std::filesystem::is_directory(candidate, ec)) {
            return candidate;
        }
    }
    return std::filesystem::current_path(ec) / "assets";
}

std::filesystem::path config_path()
{
    if (const char* base = SDL_GetBasePath()) {
        return path_from_utf8(base) / "config.toml";
    }
    std::error_code ec;
    return std::filesystem::current_path(ec) / "config.toml";
}

} // namespace astraxis
