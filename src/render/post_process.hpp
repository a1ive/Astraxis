#pragma once

#include <SDL3/SDL_gpu.h>

#include <cstdint>
#include <vector>

namespace astraxis {

struct PostSettings {
    float exposure = 1.0f;
    float bloom_strength = 0.04f;
    float bloom_threshold = 1.5f; // linear HDR; only brighter light blooms
};

// HDR -> bloom (downsample/upsample chain) -> ACES tonemap -> swapchain.
class PostProcess {
public:
    bool init(SDL_GPUDevice* device, SDL_GPUTextureFormat hdr_format, SDL_GPUTextureFormat output_format);
    void shutdown();

    // Records all post passes; `output` is fully overwritten.
    void run(SDL_GPUCommandBuffer* cmd, SDL_GPUTexture* hdr, uint32_t width, uint32_t height,
             SDL_GPUTexture* output, const PostSettings& settings);

private:
    struct Level {
        SDL_GPUTexture* texture = nullptr;
        uint32_t width = 0;
        uint32_t height = 0;
    };

    bool ensure_chain(uint32_t width, uint32_t height);
    void release_chain();
    void fullscreen_pass(SDL_GPUCommandBuffer* cmd, SDL_GPUTexture* target, SDL_GPULoadOp load,
                         SDL_GPUGraphicsPipeline* pipeline, SDL_GPUTexture* source, const void* uniforms,
                         uint32_t uniforms_size);

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUTextureFormat m_hdr_format = SDL_GPU_TEXTUREFORMAT_INVALID;
    SDL_GPUGraphicsPipeline* m_downsample = nullptr;
    SDL_GPUGraphicsPipeline* m_upsample = nullptr;
    SDL_GPUGraphicsPipeline* m_composite = nullptr;
    SDL_GPUSampler* m_linear_clamp = nullptr;

    std::vector<Level> m_levels;
    uint32_t m_chain_width = 0;
    uint32_t m_chain_height = 0;
};

} // namespace astraxis
