#pragma once

#include <SDL3/SDL_gpu.h>

#include <cstdint>

namespace astraxis {

// Creates a GPU buffer and uploads `data` into it (blocking submit).
SDL_GPUBuffer* create_static_buffer(SDL_GPUDevice* device, SDL_GPUBufferUsageFlags usage, const void* data,
                                    uint32_t size);

// GPU buffer rewritten every frame through a transfer buffer.
class StreamBuffer {
public:
    void init(SDL_GPUDevice* device, SDL_GPUBufferUsageFlags usage)
    {
        m_device = device;
        m_usage = usage;
    }
    void shutdown();

    // Records a copy pass; must be called outside any render pass. Grows as needed.
    bool upload(SDL_GPUCommandBuffer* cmd, const void* data, uint32_t size);

    SDL_GPUBuffer* buffer() const { return m_buffer; }

private:
    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUBufferUsageFlags m_usage = 0;
    SDL_GPUBuffer* m_buffer = nullptr;
    SDL_GPUTransferBuffer* m_transfer = nullptr;
    uint32_t m_capacity = 0;
};

} // namespace astraxis
