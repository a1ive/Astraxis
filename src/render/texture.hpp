#pragma once

#include <SDL3/SDL_gpu.h>
#include <glm/vec2.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// Loads a JPEG/PNG as an sRGB RGBA8 texture with a full mip chain (sampling
// returns linear values). Returns nullptr and sets `error` on failure;
// `width` (optional) receives the image width.
SDL_GPUTexture* load_texture_srgb(SDL_GPUDevice* device, const std::filesystem::path& path, std::string* error,
                                  uint32_t* width = nullptr);

// Loads a JPEG/PNG whose channels are plain numbers, not colors (a relief map): RGBA8 with a
// full mip chain, sampled as stored.
SDL_GPUTexture* load_texture_linear(SDL_GPUDevice* device, const std::filesystem::path& path, std::string* error);

// RGBA8 texture with optional mipmaps; `srgb` selects the sRGB format (sampling
// decodes to linear), otherwise values are used as they are.
SDL_GPUTexture* create_texture_rgba8(SDL_GPUDevice* device, const uint8_t* rgba, uint32_t width, uint32_t height,
                                     bool srgb, bool mipmaps = false);

// One-row R16G16_FLOAT texture of a radial profile (e.g. a ring's optical
// depth), with a mip chain of pairwise averages so that distant rings do not
// shimmer. Linear values, no sRGB decoding.
SDL_GPUTexture* create_profile_texture(SDL_GPUDevice* device, const std::vector<glm::vec2>& samples);

// R16_FLOAT texture of a table of values, without mips, row after row.
SDL_GPUTexture* create_float_texture(SDL_GPUDevice* device, const float* values, uint32_t width, uint32_t height);
// Two channels per texel (R16G16_FLOAT), interleaved.
SDL_GPUTexture* create_float2_texture(SDL_GPUDevice* device, const float* values, uint32_t width, uint32_t height);

// 1x1 texture of the given sRGB color.
SDL_GPUTexture* create_solid_texture(SDL_GPUDevice* device, uint8_t r, uint8_t g, uint8_t b, uint8_t a);

} // namespace astraxis
