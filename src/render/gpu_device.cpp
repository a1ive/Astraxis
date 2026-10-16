#include "render/gpu_device.hpp"

#include <SDL3/SDL_log.h>

namespace astraxis {

namespace {

SDL_GPUSampleCount sample_count_of(int samples)
{
    switch (samples) {
    case 2:
        return SDL_GPU_SAMPLECOUNT_2;
    case 4:
        return SDL_GPU_SAMPLECOUNT_4;
    default:
        return SDL_GPU_SAMPLECOUNT_1;
    }
}

} // namespace

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
    m_max_msaa = supports_msaa(4) ? 4 : supports_msaa(2) ? 2 : 1;
    set_msaa(4);

    const char* name = SDL_GetStringProperty(SDL_GetGPUDeviceProperties(m_device), SDL_PROP_GPU_DEVICE_NAME_STRING, "");
    SDL_LogInfo(SDL_LOG_CATEGORY_GPU, "GPU: %s (%s), MSAA up to %dx", name, driver_name(), m_max_msaa);
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

bool GpuDevice::set_msaa(int samples)
{
    const int count = samples >= 4 && supports_msaa(4) ? 4 : samples >= 2 && supports_msaa(2) ? 2 : 1;
    const SDL_GPUSampleCount sample_count = sample_count_of(count);
    if (sample_count == m_scene_format.samples) {
        return false;
    }
    m_scene_format.samples = sample_count;
    return true;
}

int GpuDevice::msaa() const
{
    return 1 << static_cast<int>(m_scene_format.samples);
}

bool GpuDevice::supports_msaa(int samples) const
{
    const SDL_GPUSampleCount count = sample_count_of(samples);
    return SDL_GPUTextureSupportsSampleCount(m_device, m_scene_format.color, count) &&
           SDL_GPUTextureSupportsSampleCount(m_device, m_scene_format.depth, count);
}

const char* GpuDevice::driver_name() const
{
    return m_device ? SDL_GetGPUDeviceDriver(m_device) : "none";
}

} // namespace astraxis
