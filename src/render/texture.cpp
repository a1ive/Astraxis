#include "render/texture.hpp"

#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#define STBI_ONLY_PNG
#define STBI_NO_STDIO
#include <stb_image.h>

#include <SDL3/SDL_log.h>
#include <glm/packing.hpp>

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
    SDL_GPUTransferBufferCreateInfo transfer_info = {
        .usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD, .size = size, .props = 0};
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

SDL_GPUTexture* create_profile_texture(SDL_GPUDevice* device, const std::vector<glm::vec2>& samples)
{
    if (samples.empty()) {
        return nullptr;
    }
    // Mip levels: pairwise averages down to a single texel.
    std::vector<std::vector<glm::vec2>> levels{samples};
    while (levels.back().size() > 1) {
        const std::vector<glm::vec2>& prev = levels.back();
        std::vector<glm::vec2> next((prev.size() + 1) / 2);
        for (size_t i = 0; i < next.size(); ++i) {
            const glm::vec2 b = 2 * i + 1 < prev.size() ? prev[2 * i + 1] : prev[2 * i];
            next[i] = 0.5f * (prev[2 * i] + b);
        }
        levels.push_back(std::move(next));
    }

    SDL_GPUTextureCreateInfo info = {};
    info.type = SDL_GPU_TEXTURETYPE_2D;
    info.format = SDL_GPU_TEXTUREFORMAT_R16G16_FLOAT;
    info.usage = SDL_GPU_TEXTUREUSAGE_SAMPLER;
    info.width = static_cast<uint32_t>(samples.size());
    info.height = 1;
    info.layer_count_or_depth = 1;
    info.num_levels = static_cast<uint32_t>(levels.size());
    info.sample_count = SDL_GPU_SAMPLECOUNT_1;
    SDL_GPUTexture* texture = SDL_CreateGPUTexture(device, &info);
    if (!texture) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_CreateGPUTexture failed: %s", SDL_GetError());
        return nullptr;
    }

    // Each level at a 512-byte aligned offset (the strictest backend requirement).
    constexpr uint32_t kAlign = 512;
    std::vector<uint32_t> offsets;
    uint32_t size = 0;
    for (const auto& level : levels) {
        offsets.push_back(size);
        size += (static_cast<uint32_t>(level.size()) * 4 + kAlign - 1) / kAlign * kAlign;
    }
    SDL_GPUTransferBufferCreateInfo transfer_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD, .size = size, .props = 0};
    SDL_GPUTransferBuffer* transfer = SDL_CreateGPUTransferBuffer(device, &transfer_info);
    if (!transfer) {
        SDL_ReleaseGPUTexture(device, texture);
        return nullptr;
    }
    auto* mapped = static_cast<uint8_t*>(SDL_MapGPUTransferBuffer(device, transfer, false));
    for (size_t l = 0; l < levels.size(); ++l) {
        auto* texels = reinterpret_cast<uint32_t*>(mapped + offsets[l]);
        for (size_t i = 0; i < levels[l].size(); ++i) {
            texels[i] = glm::packHalf2x16(levels[l][i]);
        }
    }
    SDL_UnmapGPUTransferBuffer(device, transfer);

    SDL_GPUCommandBuffer* cmd = SDL_AcquireGPUCommandBuffer(device);
    SDL_GPUCopyPass* copy = SDL_BeginGPUCopyPass(cmd);
    for (size_t l = 0; l < levels.size(); ++l) {
        SDL_GPUTextureTransferInfo src = {};
        src.transfer_buffer = transfer;
        src.offset = offsets[l];
        src.pixels_per_row = static_cast<uint32_t>(levels[l].size());
        src.rows_per_layer = 1;
        SDL_GPUTextureRegion dst = {};
        dst.texture = texture;
        dst.mip_level = static_cast<uint32_t>(l);
        dst.w = static_cast<uint32_t>(levels[l].size());
        dst.h = 1;
        dst.d = 1;
        SDL_UploadToGPUTexture(copy, &src, &dst, false);
    }
    SDL_EndGPUCopyPass(copy);
    SDL_SubmitGPUCommandBuffer(cmd);
    SDL_ReleaseGPUTransferBuffer(device, transfer);
    return texture;
}

SDL_GPUTexture* create_solid_texture(SDL_GPUDevice* device, uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    const uint8_t rgba[4] = {r, g, b, a};
    return create_and_upload(device, rgba, 1, 1, false);
}

} // namespace astraxis
