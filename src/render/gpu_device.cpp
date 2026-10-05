#include "render/gpu_device.hpp"

#include <SDL3/SDL_log.h>

namespace astraxis {

bool GpuDevice::init(bool debug)
{
    // DXIL / SPIR-V for our shaders (Direct3D 12 / Vulkan); DXBC because the
    // ImGui backend ships DXBC for D3D12. SDL tries Direct3D 12 before Vulkan;
    // SDL_GPU_DRIVER=vulkan forces Vulkan (e.g. to test the SPIR-V path on Windows).
    m_device = SDL_CreateGPUDevice(SDL_GPU_SHADERFORMAT_DXIL | SDL_GPU_SHADERFORMAT_DXBC | SDL_GPU_SHADERFORMAT_SPIRV,
                                   debug, nullptr);
    if (!m_device) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "SDL_CreateGPUDevice failed: %s", SDL_GetError());
        return false;
    }

    m_scene_format.color = SDL_GPU_TEXTUREFORMAT_R16G16B16A16_FLOAT;
    m_scene_format.depth = SDL_GPUTextureSupportsFormat(m_device, SDL_GPU_TEXTUREFORMAT_D32_FLOAT,
                                                        SDL_GPU_TEXTURETYPE_2D,
                                                        SDL_GPU_TEXTUREUSAGE_DEPTH_STENCIL_TARGET)
                               ? SDL_GPU_TEXTUREFORMAT_D32_FLOAT
                               : SDL_GPU_TEXTUREFORMAT_D24_UNORM;
    const bool msaa4 = SDL_GPUTextureSupportsSampleCount(m_device, m_scene_format.color, SDL_GPU_SAMPLECOUNT_4) &&
                       SDL_GPUTextureSupportsSampleCount(m_device, m_scene_format.depth, SDL_GPU_SAMPLECOUNT_4);
    m_scene_format.samples = msaa4 ? SDL_GPU_SAMPLECOUNT_4 : SDL_GPU_SAMPLECOUNT_1;

    SDL_LogInfo(SDL_LOG_CATEGORY_GPU, "GPU driver: %s, MSAA %s", driver_name(), msaa4 ? "4x" : "off");
    return true;
}

void GpuDevice::shutdown()
{
    if (m_device) {
        SDL_WaitForGPUIdle(m_device);
        SDL_DestroyGPUDevice(m_device);
        m_device = nullptr;
    }
}

const char* GpuDevice::driver_name() const
{
    return m_device ? SDL_GetGPUDeviceDriver(m_device) : "none";
}

} // namespace astraxis
