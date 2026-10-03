#pragma once

#include <SDL3/SDL_gpu.h>

#include <cstdint>
#include <span>

namespace astraxis {

// Embedded binaries come from <shaders/<name>.<stage>.h> (k<Name><Stage>Dxil / Spirv).
struct ShaderDesc {
    std::span<const unsigned char> dxil;  // Direct3D 12
    std::span<const unsigned char> spirv; // Vulkan
    SDL_GPUShaderStage stage = SDL_GPU_SHADERSTAGE_VERTEX;
    uint32_t num_samplers = 0;
    uint32_t num_storage_textures = 0;
    uint32_t num_storage_buffers = 0;
    uint32_t num_uniform_buffers = 0;
};

// Creates a shader (entry point `main`) in the format the device accepts,
// preferring DXIL. Returns nullptr on failure.
SDL_GPUShader* create_shader(SDL_GPUDevice* device, const ShaderDesc& desc);

} // namespace astraxis
