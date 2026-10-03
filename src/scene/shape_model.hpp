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
    std::vector<uint32_t> indices;    // 3 per triangle, counter-clockwise seen from outside
};

// Loads a mesh written by tools/shapes/make_arrokoth.py.
bool load_shape_model(const std::filesystem::path& path, ShapeModel& out, std::string* error);

} // namespace astraxis
