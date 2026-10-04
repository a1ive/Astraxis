#include "scene/shape_model.hpp"

#include <glm/geometric.hpp>

#include <array>
#include <cstring>
#include <fstream>
#include <map>

namespace astraxis {

static_assert(sizeof(glm::vec3) == 3 * sizeof(float), "positions are read as packed floats");

// File format (tools/shapes/axmesh.py), little endian:
//   char[8] "AXMESH2\0"; uint32 vertex count; uint32 triangle count; uint32 flags;
//   float32 x, y, z (km) per vertex; float32 relative albedo per vertex;
//   float32 map u per vertex (if flags bit 0);
//   uint32 vertex indices, 3 per triangle (counter-clockwise seen from outside).
bool load_shape_model(const std::filesystem::path& path, ShapeModel& out, std::string* error)
{
    auto fail = [&](const std::string& what) {
        if (error) {
            *error = path.string() + ": " + what;
        }
        return false;
    };
    std::ifstream in(path, std::ios::binary);
    char magic[8];
    uint32_t counts[3];
    if (!in || !in.read(magic, 8) || std::memcmp(magic, "AXMESH2", 8) != 0 ||
        !in.read(reinterpret_cast<char*>(counts), sizeof(counts))) {
        return fail("not a shape model");
    }
    if (counts[0] < 4 || counts[1] < 4) {
        return fail("empty mesh");
    }
    out.positions.resize(counts[0]);
    out.albedo.resize(counts[0]);
    out.map_u.resize((counts[2] & 1u) ? counts[0] : 0u);
    out.indices.resize(static_cast<size_t>(counts[1]) * 3);
    if (!in.read(reinterpret_cast<char*>(out.positions.data()),
                 static_cast<std::streamsize>(out.positions.size() * sizeof(glm::vec3))) ||
        !in.read(reinterpret_cast<char*>(out.albedo.data()),
                 static_cast<std::streamsize>(out.albedo.size() * sizeof(float))) ||
        !in.read(reinterpret_cast<char*>(out.map_u.data()),
                 static_cast<std::streamsize>(out.map_u.size() * sizeof(float))) ||
        !in.read(reinterpret_cast<char*>(out.indices.data()),
                 static_cast<std::streamsize>(out.indices.size() * sizeof(uint32_t)))) {
        return fail("truncated");
    }
    for (const uint32_t i : out.indices) {
        if (i >= counts[0]) {
            return fail("vertex index out of range");
        }
    }

    // Vertices at the same position (a map seam, a pole row) share one normal, or
    // the shading would show the seam. The cross product's length is twice the face
    // area, which weights the average.
    std::map<std::array<float, 3>, uint32_t> first;
    std::vector<uint32_t> shared(out.positions.size());
    for (uint32_t i = 0; i < out.positions.size(); ++i) {
        const glm::vec3& p = out.positions[i];
        shared[i] = first.try_emplace({p.x, p.y, p.z}, i).first->second;
    }
    std::vector<glm::vec3> sum(out.positions.size(), glm::vec3(0.0f));
    for (size_t k = 0; k < out.indices.size(); k += 3) {
        const uint32_t a = out.indices[k], b = out.indices[k + 1], c = out.indices[k + 2];
        const glm::vec3 n = glm::cross(out.positions[b] - out.positions[a], out.positions[c] - out.positions[a]);
        sum[shared[a]] += n;
        sum[shared[b]] += n;
        sum[shared[c]] += n;
    }
    out.normals.resize(out.positions.size());
    for (size_t i = 0; i < out.positions.size(); ++i) {
        const glm::vec3& n = sum[shared[i]];
        const float len = glm::length(n);
        out.normals[i] = len > 0.0f ? n / len : glm::vec3(0.0f, 0.0f, 1.0f);
    }
    return true;
}

} // namespace astraxis
