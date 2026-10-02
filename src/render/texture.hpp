#pragma once

#include <SDL3/SDL_gpu.h>

#include <cstdint>
#include <filesystem>
#include <string>

namespace astraxis {

// Loads a JPEG/PNG as an sRGB RGBA8 texture with a full mip chain (sampling
// returns linear values). Returns nullptr and sets `error` on failure;
// `width` (optional) receives the image width.
SDL_GPUTexture* load_texture_srgb(SDL_GPUDevice* device, const std::filesystem::path& path, std::string* error,
                                  uint32_t* width = nullptr);

// RGBA8 texture without mipmaps; `srgb` selects the sRGB format (sampling
// decodes to linear), otherwise values are used as they are.
SDL_GPUTexture* create_texture_rgba8(SDL_GPUDevice* device, const uint8_t* rgba, uint32_t width, uint32_t height,
                                     bool srgb);

// 1x1 texture of the given sRGB color.
SDL_GPUTexture* create_solid_texture(SDL_GPUDevice* device, uint8_t r, uint8_t g, uint8_t b, uint8_t a);

} // namespace astraxis
