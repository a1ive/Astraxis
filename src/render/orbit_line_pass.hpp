#pragma once

#include "render/gpu_buffer.hpp"
#include "render/renderer.hpp"
#include "scene/camera.hpp"

#include <SDL3/SDL_gpu.h>
#include <glm/vec4.hpp>

#include <vector>

namespace astraxis {

// Draws polylines (orbits) with constant screen-space width.
// Usage per frame: begin() -> add_line()... -> upload() (outside render passes) -> draw().
class OrbitLinePass {
public:
    bool init(SDL_GPUDevice* device, const SceneTargetFormat& format);
    void shutdown();

    void begin();
    // `points`: xyz = camera-relative position, w = fade (0 head .. 1 tail).
    void add_line(const std::vector<glm::vec4>& points, const glm::vec4& color);
    void upload(SDL_GPUCommandBuffer* cmd);
    void draw(SDL_GPUCommandBuffer* cmd, SDL_GPURenderPass* pass, const CameraView& view, float width_px) const;

private:
    struct Line {
        uint32_t first = 0;
        uint32_t count = 0;
        glm::vec4 color{1.0f};
    };

    SDL_GPUDevice* m_device = nullptr;
    SDL_GPUGraphicsPipeline* m_pipeline = nullptr;
    StreamBuffer m_points_buffer;
    std::vector<glm::vec4> m_points;
    std::vector<Line> m_lines;
};

} // namespace astraxis
