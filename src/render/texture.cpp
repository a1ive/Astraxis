#include "render/texture.hpp"

#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#define STBI_ONLY_PNG
#define STBI_NO_STDIO
#include <stb_image.h>

#include <SDL3/SDL_log.h>

#include <algorithm>
#include <bit>
#include <cstring>
#include <fstream>
#include <vector>

namespace astraxis {

namespace {

SDL_GPUTexture* create_and_upload(SDL_GPUDevice* device, const uint8_t* rgba, uint32_t width, uint32_t height,
                                  bool mipmaps, bool srgb = true)
{
    const uint32_t levels = mipmaps ? static_cast<uint32_t>(std::bit_width(std::max(width, height))) : 1u;

    SDL_GPUTextureCreateInfo info = {};
    info.type = SDL_GPU_TEXTURETYPE_2D;
    info.format = srgb ? SDL_GPU_TEXTUREFORMAT_R8G8B8A8_UNORM_SRGB : SDL_GPU_TEXTUREFORMAT_R8G8B8A8_UNORM;
    info.usage = SDL_GPU_TEXTUREUSAGE_SAMPLER | (levels > 1 ? SDL_GPU_TEXTUREUSAGE_COLOR_TARGET : 0);
    info.width = width;
    info.height = height;
    info.layer_count_or_depth = 1;
    info.num_levels = levels;
    info.sample_count = SDL_GPU_SAMPLECOUNT_1;
    SDL_GPUTexture* texture = SDL_CreateGPUTexture(device, &info);
    if (!texture) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_CreateGPUTexture failed: %s", SDL_GetError());
        return nullptr;
    }

    const uint32_t size = width * height * 4;
    SDL_GPUTransferBufferCreateInfo transfer_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD, .size = size};
    SDL_GPUTransferBuffer* transfer = SDL_CreateGPUTransferBuffer(device, &transfer_info);
    if (!transfer) {
        SDL_ReleaseGPUTexture(device, texture);
        return nullptr;
    }
    void* mapped = SDL_MapGPUTransferBuffer(device, transfer, false);
    std::memcpy(mapped, rgba, size);
    SDL_UnmapGPUTransferBuffer(device, transfer);

    SDL_GPUCommandBuffer* cmd = SDL_AcquireGPUCommandBuffer(device);
    SDL_GPUCopyPass* copy = SDL_BeginGPUCopyPass(cmd);
    SDL_GPUTextureTransferInfo src = {};
    src.transfer_buffer = transfer;
    src.pixels_per_row = width;
    src.rows_per_layer = height;
    SDL_GPUTextureRegion dst = {};
    dst.texture = texture;
    dst.w = width;
    dst.h = height;
    dst.d = 1;
    SDL_UploadToGPUTexture(copy, &src, &dst, false);
    SDL_EndGPUCopyPass(copy);
    if (levels > 1) {
        SDL_GenerateMipmapsForGPUTexture(cmd, texture);
    }
    SDL_SubmitGPUCommandBuffer(cmd);
    SDL_ReleaseGPUTransferBuffer(device, transfer);
    return texture;
}

} // namespace

SDL_GPUTexture* load_texture_srgb(SDL_GPUDevice* device, const std::filesystem::path& path, std::string* error,
                                  uint32_t* width)
{
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        if (error) {
            *error = "cannot open " + path.string();
        }
        return nullptr;
    }
    const std::vector<char> bytes((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    int w = 0;
    int h = 0;
    int channels = 0;
    stbi_uc* pixels = stbi_load_from_memory(reinterpret_cast<const stbi_uc*>(bytes.data()),
                                            static_cast<int>(bytes.size()), &w, &h, &channels, 4);
    if (!pixels) {
        if (error) {
            *error = path.string() + ": " + stbi_failure_reason();
        }
        return nullptr;
    }

    SDL_GPUTexture* texture =
        create_and_upload(device, pixels, static_cast<uint32_t>(w), static_cast<uint32_t>(h), true);
    stbi_image_free(pixels);
    if (!texture && error) {
        *error = path.string() + ": GPU upload failed";
    }
    if (texture && width) {
        *width = static_cast<uint32_t>(w);
    }
    return texture;
}

SDL_GPUTexture* create_texture_rgba8(SDL_GPUDevice* device, const uint8_t* rgba, uint32_t width, uint32_t height,
                                     bool srgb)
{
    return create_and_upload(device, rgba, width, height, false, srgb);
}

SDL_GPUTexture* create_solid_texture(SDL_GPUDevice* device, uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    const uint8_t rgba[4] = {r, g, b, a};
    return create_and_upload(device, rgba, 1, 1, false);
}

} // namespace astraxis
