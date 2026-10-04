#pragma once

#include <glm/vec3.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// An irregular body's surface as a triangle mesh in its body-fixed frame
// (x toward the prime meridian, z along the spin pole), drawn in place of the
// ellipsoid.
struct ShapeModel {
    std::vector<glm::vec3> positions; // km
    std::vector<glm::vec3> normals;   // unit, area-weighted average of the adjacent faces
    std::vector<float> albedo;        // relative to the mean; scales the body color
    // Optional (empty if none): east longitude / 360 deg of each vertex, for a
    // surface map in planetocentric coordinates; continuous across each triangle.
    std::vector<float> map_u;
    std::vector<uint32_t> indices;    // 3 per triangle, counter-clockwise seen from outside
};

// Loads a mesh written by the tools in tools/shapes/ (axmesh.py).
bool load_shape_model(const std::filesystem::path& path, ShapeModel& out, std::string* error);

} // namespace astraxis
