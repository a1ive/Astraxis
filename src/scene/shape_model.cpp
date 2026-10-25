#include "scene/shape_model.hpp"

#include <glm/geometric.hpp>

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <map>
#include <numeric>

namespace astraxis {

static_assert(sizeof(glm::vec3) == 3 * sizeof(float), "positions are read as packed floats");

// File format (tools/shapes/axmesh.py), little endian:
//   char[8] "AXMESH2\0"; uint32 vertex count; uint32 triangle count; uint32 flags
//   (bit 0: map coordinates; bit 1: hard edges, see below);
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
    // the shading would show the seam. A mesh with hard edges (flags bit 1: a
    // built shape such as a spacecraft, whose faces have their own vertices) is
    // smooth only where triangles share a vertex index. The cross product's length
    // is twice the face area, which weights the average.
    const bool hard_edges = (counts[2] & 2u) != 0;
    std::map<std::array<float, 3>, uint32_t> first;
    std::vector<uint32_t> shared(out.positions.size());
    for (uint32_t i = 0; i < out.positions.size(); ++i) {
        const glm::vec3& p = out.positions[i];
        shared[i] = hard_edges ? i : first.try_emplace({p.x, p.y, p.z}, i).first->second;
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

std::vector<ShapeOcclusionNode> build_shape_occlusion(const ShapeModel& shape)
{
    static_assert(sizeof(ShapeOcclusionNode) == 48);
    const size_t count = shape.indices.size() / 3;
    if (count == 0) {
        return {};
    }
    std::vector<size_t> triangles(count);
    std::iota(triangles.begin(), triangles.end(), size_t{0});
    std::vector<ShapeOcclusionNode> nodes;
    nodes.reserve(2 * count - 1);
    auto vertex = [&](size_t triangle, size_t corner) {
        return shape.positions[shape.indices[3 * triangle + corner]];
    };
    auto build = [&](auto&& self, size_t begin, size_t end) -> void {
        const size_t index = nodes.size();
        nodes.emplace_back();
        if (end - begin == 1) {
            const size_t t = triangles[begin];
            nodes[index].a = vertex(t, 0);
            nodes[index].b = vertex(t, 1);
            nodes[index].c = vertex(t, 2);
            nodes[index].leaf = 1;
        } else {
            glm::vec3 lo = vertex(triangles[begin], 0), hi = lo;
            for (size_t i = begin; i < end; ++i) {
                for (size_t k = 0; k < 3; ++k) {
                    lo = glm::min(lo, vertex(triangles[i], k));
                    hi = glm::max(hi, vertex(triangles[i], k));
                }
            }
            // Pad the slabs for float ray transforms and shared-edge rounding.
            const glm::vec3 margin = 1e-6f * (glm::max(glm::abs(lo), glm::abs(hi)) + 1.0f);
            nodes[index].a = lo - margin;
            nodes[index].b = hi + margin;
            const glm::vec3 extent = hi - lo;
            const int axis = extent.x >= extent.y && extent.x >= extent.z ? 0 : (extent.y >= extent.z ? 1 : 2);
            const size_t middle = begin + (end - begin) / 2;
            std::nth_element(triangles.begin() + static_cast<std::ptrdiff_t>(begin),
                             triangles.begin() + static_cast<std::ptrdiff_t>(middle),
                             triangles.begin() + static_cast<std::ptrdiff_t>(end), [&](size_t a, size_t b) {
                return vertex(a, 0)[axis] + vertex(a, 1)[axis] + vertex(a, 2)[axis] <
                       vertex(b, 0)[axis] + vertex(b, 1)[axis] + vertex(b, 2)[axis];
            });
            self(self, begin, middle);
            self(self, middle, end);
        }
        nodes[index].escape = static_cast<uint32_t>(nodes.size());
    };
    build(build, 0, count);
    return nodes;
}

} // namespace astraxis
