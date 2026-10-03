#include "scene/shape_model.hpp"

#include <glm/geometric.hpp>

#include <cstring>
#include <fstream>

namespace astraxis {

static_assert(sizeof(glm::vec3) == 3 * sizeof(float), "positions are read as packed floats");

// File format, little endian:
//   char[8] "AXMESH1\0"; uint32 vertex count; uint32 triangle count;
//   float32 x, y, z (km) per vertex; float32 relative albedo per vertex;
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
    uint32_t counts[2];
    if (!in || !in.read(magic, 8) || std::memcmp(magic, "AXMESH1", 8) != 0 ||
        !in.read(reinterpret_cast<char*>(counts), sizeof(counts))) {
        return fail("not a shape model");
    }
    if (counts[0] < 4 || counts[1] < 4) {
        return fail("empty mesh");
    }
    out.positions.resize(counts[0]);
    out.albedo.resize(counts[0]);
    out.indices.resize(static_cast<size_t>(counts[1]) * 3);
    if (!in.read(reinterpret_cast<char*>(out.positions.data()),
                 static_cast<std::streamsize>(out.positions.size() * sizeof(glm::vec3))) ||
        !in.read(reinterpret_cast<char*>(out.albedo.data()),
                 static_cast<std::streamsize>(out.albedo.size() * sizeof(float))) ||
        !in.read(reinterpret_cast<char*>(out.indices.data()),
                 static_cast<std::streamsize>(out.indices.size() * sizeof(uint32_t)))) {
        return fail("truncated");
    }
    for (const uint32_t i : out.indices) {
        if (i >= counts[0]) {
            return fail("vertex index out of range");
        }
    }

    // The cross product's length is twice the face area, which weights the average.
    out.normals.assign(out.positions.size(), glm::vec3(0.0f));
    for (size_t k = 0; k < out.indices.size(); k += 3) {
        const uint32_t a = out.indices[k], b = out.indices[k + 1], c = out.indices[k + 2];
        const glm::vec3 n = glm::cross(out.positions[b] - out.positions[a], out.positions[c] - out.positions[a]);
        out.normals[a] += n;
        out.normals[b] += n;
        out.normals[c] += n;
    }
    for (glm::vec3& n : out.normals) {
        const float len = glm::length(n);
        n = len > 0.0f ? n / len : glm::vec3(0.0f, 0.0f, 1.0f);
    }
    return true;
}

} // namespace astraxis
