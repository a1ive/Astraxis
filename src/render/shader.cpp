#include "render/shader.hpp"

#include <SDL3/SDL_log.h>

namespace astraxis {

SDL_GPUShader* create_shader(SDL_GPUDevice* device, const ShaderDesc& desc)
{
    if (!(SDL_GetGPUShaderFormats(device) & SDL_GPU_SHADERFORMAT_DXIL)) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "GPU device does not accept DXIL shaders");
        return nullptr;
    }

    SDL_GPUShaderCreateInfo info = {};
    info.code = desc.dxil.data();
    info.code_size = desc.dxil.size();
    info.entrypoint = "main";
    info.format = SDL_GPU_SHADERFORMAT_DXIL;
    info.stage = desc.stage;
    info.num_samplers = desc.num_samplers;
    info.num_storage_textures = desc.num_storage_textures;
    info.num_storage_buffers = desc.num_storage_buffers;
    info.num_uniform_buffers = desc.num_uniform_buffers;

    SDL_GPUShader* shader = SDL_CreateGPUShader(device, &info);
    if (!shader) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_CreateGPUShader failed: %s", SDL_GetError());
    }
    return shader;
}

} // namespace astraxis
