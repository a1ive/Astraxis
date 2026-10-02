#include "platform/paths.hpp"

#include <SDL3/SDL_filesystem.h>

#include <string>
#include <system_error>

namespace astraxis {

namespace {

std::filesystem::path from_utf8(const char* s)
{
    const std::string str(s);
    return std::filesystem::path(std::u8string(str.begin(), str.end()));
}

} // namespace

std::filesystem::path asset_directory()
{
    std::error_code ec;
    if (const char* base = SDL_GetBasePath()) {
        const std::filesystem::path candidate = from_utf8(base) / "assets";
        if (std::filesystem::is_directory(candidate, ec)) {
            return candidate;
        }
    }
    return std::filesystem::current_path(ec) / "assets";
}

} // namespace astraxis
