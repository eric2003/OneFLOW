#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace ONEFLOW {

enum class StlFormat
{
    Auto,
    ASCII,
    Binary
};

struct StlTriangle
{
    std::array<double, 3> normal{};
    std::array<std::array<double, 3>, 3> vertices{};
    std::uint16_t attribute{0};
};

struct StlMesh
{
    // ASCII solid name or the opaque 80-byte binary header.
    std::string header;
    std::vector<StlTriangle> triangles;
};

// Recover the single closed XY boundary represented by a planar STL surface.
// STL triangulation edges are removed by edge-incidence counting; only edges
// used by one triangle are returned as the outer loop.
struct StlPlanarBoundary
{
    std::vector<std::array<double, 3>> vertices;
};

StlPlanarBoundary ExtractStlPlanarBoundary(
    const StlMesh& mesh,
    double weldTolerance = 1.0e-9);

class StlMeshReader
{
public:
    static StlMesh Read(
        const std::string& fileName,
        StlFormat format = StlFormat::Auto);
};

class StlMeshWriter
{
public:
    static void WriteASCII(const std::string& fileName, const StlMesh& mesh);
    static void WriteBinary(const std::string& fileName, const StlMesh& mesh);
};

} // namespace ONEFLOW
