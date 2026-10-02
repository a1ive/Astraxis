#pragma once

#include <filesystem>

namespace astraxis {

// Directory holding runtime assets: "<exe dir>/assets" (copied there by the
// build), falling back to "<working dir>/assets".
std::filesystem::path asset_directory();

} // namespace astraxis
