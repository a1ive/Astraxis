#pragma once

#include <SDL3/SDL_gpu.h>
#include <glm/vec2.hpp>

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace astraxis {

// A JPEG/PNG decoded to RGBA8 on the CPU, not yet a texture: decoding is the
// slow part of a texture load and may run on any thread, the upload may not.
struct DecodedImage {
    struct Free {
        void operator()(uint8_t* pixels) const;
    };
    std::unique_ptr<uint8_t, Free> pixels; // null if the file could not be read or decoded
    uint32_t width = 0;
    uint32_t height = 0;
    std::string source; // the file, for messages
    std::string error;  // why `pixels` is null
};

DecodedImage decode_image(const std::filesystem::path& path);

// An RGBA8 texture of a decoded image with a full mip chain: sRGB for colors
// (sampling returns linear values), plain for numbers such as a relief map.
// Returns nullptr and sets `error` if the image is empty or the upload fails.
SDL_GPUTexture* upload_image(SDL_GPUDevice* device, const DecodedImage& image, bool srgb, std::string* error);

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
