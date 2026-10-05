#pragma once

#include <filesystem>

namespace astraxis {

// A path from a UTF-8 string (SDL strings and, through SDL_main, argv).
std::filesystem::path path_from_utf8(const char* s);

// Directory holding runtime assets: "<exe dir>/assets" (copied there by the
// build), falling back to "<working dir>/assets".
std::filesystem::path asset_directory();

// The settings file: "<exe dir>/config.toml" (portable; not saved when that
// directory is read-only), falling back to "<working dir>/config.toml".
std::filesystem::path config_path();

} // namespace astraxis
