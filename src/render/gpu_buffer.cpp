#include "render/gpu_buffer.hpp"

#include <SDL3/SDL_log.h>

#include <cstring>

namespace astraxis {

SDL_GPUBuffer* create_static_buffer(SDL_GPUDevice* device, SDL_GPUBufferUsageFlags usage, const void* data,
                                    uint32_t size)
{
    SDL_GPUBufferCreateInfo buffer_info = {.usage = usage, .size = size};
    SDL_GPUBuffer* buffer = SDL_CreateGPUBuffer(device, &buffer_info);

    SDL_GPUTransferBufferCreateInfo transfer_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD, .size = size};
    SDL_GPUTransferBuffer* transfer = SDL_CreateGPUTransferBuffer(device, &transfer_info);
    if (!buffer || !transfer) {
        SDL_LogError(SDL_LOG_CATEGORY_GPU, "Static buffer creation failed: %s", SDL_GetError());
        SDL_ReleaseGPUBuffer(device, buffer);
        SDL_ReleaseGPUTransferBuffer(device, transfer);
        return nullptr;
    }

    void* mapped = SDL_MapGPUTransferBuffer(device, transfer, false);
    std::memcpy(mapped, data, size);
    SDL_UnmapGPUTransferBuffer(device, transfer);

    SDL_GPUCommandBuffer* cmd = SDL_AcquireGPUCommandBuffer(device);
    SDL_GPUCopyPass* copy = SDL_BeginGPUCopyPass(cmd);
    SDL_GPUTransferBufferLocation src = {.transfer_buffer = transfer, .offset = 0};
    SDL_GPUBufferRegion dst = {.buffer = buffer, .offset = 0, .size = size};
    SDL_UploadToGPUBuffer(copy, &src, &dst, false);
    SDL_EndGPUCopyPass(copy);
    SDL_SubmitGPUCommandBuffer(cmd);
    SDL_ReleaseGPUTransferBuffer(device, transfer);
    return buffer;
}

void StreamBuffer::shutdown()
{
    if (m_device) {
        SDL_ReleaseGPUBuffer(m_device, m_buffer);
        SDL_ReleaseGPUTransferBuffer(m_device, m_transfer);
    }
    m_buffer = nullptr;
    m_transfer = nullptr;
    m_capacity = 0;
}

bool StreamBuffer::upload(SDL_GPUCommandBuffer* cmd, const void* data, uint32_t size)
{
    if (size == 0) {
        return true;
    }
    if (size > m_capacity) {
        // Released resources stay alive until the GPU is done with them.
        shutdown();
        uint32_t capacity = 4096;
        while (capacity < size) {
            capacity *= 2;
        }
        SDL_GPUBufferCreateInfo buffer_info = {.usage = m_usage, .size = capacity};
        m_buffer = SDL_CreateGPUBuffer(m_device, &buffer_info);
        SDL_GPUTransferBufferCreateInfo transfer_info = {.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD,
                                                         .size = capacity};
        m_transfer = SDL_CreateGPUTransferBuffer(m_device, &transfer_info);
        if (!m_buffer || !m_transfer) {
            SDL_LogError(SDL_LOG_CATEGORY_GPU, "Stream buffer creation failed: %s", SDL_GetError());
            shutdown();
            return false;
        }
        m_capacity = capacity;
    }

    void* mapped = SDL_MapGPUTransferBuffer(m_device, m_transfer, true);
    std::memcpy(mapped, data, size);
    SDL_UnmapGPUTransferBuffer(m_device, m_transfer);

    SDL_GPUCopyPass* copy = SDL_BeginGPUCopyPass(cmd);
    SDL_GPUTransferBufferLocation src = {.transfer_buffer = m_transfer, .offset = 0};
    SDL_GPUBufferRegion dst = {.buffer = m_buffer, .offset = 0, .size = size};
    SDL_UploadToGPUBuffer(copy, &src, &dst, true);
    SDL_EndGPUCopyPass(copy);
    return true;
}

} // namespace astraxis
