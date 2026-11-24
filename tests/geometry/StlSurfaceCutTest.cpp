#include "StlSurfaceCut.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>

namespace {

ONEFLOW::StlMesh MakeOutwardTetrahedron()
{
    using Point = std::array<double, 3>;
    const Point v0{0.0, 0.0, 0.0};
    const Point v1{1.0, 0.0, 0.0};
    const Point v2{0.0, 1.0, 0.0};
    const Point v3{0.0, 0.0, 1.0};
    ONEFLOW::StlMesh mesh;
    mesh.header = "tetrahedron";
    for (const auto& vertices : {
        std::array<Point, 3>{v1, v2, v3},
        std::array<Point, 3>{v0, v1, v3},
        std::array<Point, 3>{v0, v3, v2},
        std::array<Point, 3>{v0, v2, v1}})
    {
        ONEFLOW::StlTriangle triangle;
        triangle.vertices = vertices;
        mesh.triangles.push_back(triangle);
    }
    return mesh;
}

} // namespace

TEST(StlSurfaceCutTest, ReportsClosedOutwardInputTopology)
{
    const auto report = ONEFLOW::AnalyzeStlSurface(MakeOutwardTetrahedron());
    EXPECT_EQ(report.connectedComponentCount, 1u);
    EXPECT_EQ(report.boundaryEdgeCount, 0u);
    EXPECT_EQ(report.nonManifoldEdgeCount, 0u);
    EXPECT_TRUE(report.HasConsistentOutwardOrientation());
    EXPECT_NEAR(report.signedVolume, 1.0 / 6.0, 1.0e-12);
}

TEST(StlSurfaceCutTest, ClipsAndCapsBetweenTwoPlanes)
{
    const auto source = MakeOutwardTetrahedron();
    ONEFLOW::StlSurfaceCutOptions options;
    options.inlet.point = {0.2, 0.0, 0.0};
    options.inlet.normal = {1.0, 0.0, 0.0};
    options.inlet.keepPositive = true;
    options.outlet.point = {0.6, 0.0, 0.0};
    options.outlet.normal = {1.0, 0.0, 0.0};
    options.outlet.keepPositive = false;

    const auto result = ONEFLOW::CutStlSurfaceBetweenPlanes(source, options);
    EXPECT_FALSE(result.wall.triangles.empty());
    EXPECT_EQ(result.inletLoop.size(), 3u);
    EXPECT_EQ(result.outletLoop.size(), 3u);
    EXPECT_EQ(result.inletCap.triangles.size(), 1u);
    EXPECT_EQ(result.outletCap.triangles.size(), 1u);
    EXPECT_TRUE(result.topology.IsClosedTwoManifold());
    EXPECT_EQ(result.topology.connectedComponentCount, 1u);
    EXPECT_GT(result.topology.signedVolume, 0.0);
}

TEST(StlSurfaceCutTest, RejectsZeroAreaTriangles)
{
    auto source = MakeOutwardTetrahedron();
    source.triangles.front().vertices[2] = {
        0.5 * (source.triangles.front().vertices[0][0] +
               source.triangles.front().vertices[1][0]),
        0.5 * (source.triangles.front().vertices[0][1] +
               source.triangles.front().vertices[1][1]),
        0.5 * (source.triangles.front().vertices[0][2] +
               source.triangles.front().vertices[1][2])};
    EXPECT_THROW(ONEFLOW::AnalyzeStlSurface(source), std::invalid_argument);
}

TEST(StlSurfaceCutTest, RejectsOpenInputShell)
{
    auto source = MakeOutwardTetrahedron();
    source.triangles.pop_back();
    ONEFLOW::StlSurfaceCutOptions options;
    options.inlet.point = {0.2, 0.0, 0.0};
    options.outlet.point = {0.6, 0.0, 0.0};
    EXPECT_THROW(ONEFLOW::CutStlSurfaceBetweenPlanes(source, options), std::invalid_argument);
}

TEST(StlSurfaceCutTest, CutsBlindCylinderFixtureAndMeasuresObliqueEllipse)
{
#ifdef ONEFLOW_SOURCE_DIR
    const std::filesystem::path sourcePath =
        std::filesystem::path(ONEFLOW_SOURCE_DIR) /
        "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python/input/model.stl";
    const auto source = ONEFLOW::StlMeshReader::Read(sourcePath.string());
    ONEFLOW::StlSurfaceCutOptions options;
    options.inlet.point = {-1.2, 0.0, 0.0};
    options.inlet.normal = {0.8660254037844386, 0.5, 0.0};
    options.inlet.keepPositive = true;
    options.outlet.point = {1.5, 0.0, 0.0};
    options.outlet.normal = {1.0, 0.0, 0.0};
    options.outlet.keepPositive = false;
    options.tolerance = 1.0e-8;

    const auto result = ONEFLOW::CutStlSurfaceBetweenPlanes(source, options);
    EXPECT_GT(result.inletLoop.size(), 8u);
    EXPECT_GT(result.outletLoop.size(), 8u);
    EXPECT_TRUE(result.topology.IsClosedTwoManifold());
    EXPECT_EQ(result.topology.connectedComponentCount, 1u);
    EXPECT_GT(result.topology.signedVolume, 0.0);

    const auto extent = [](const auto& loop, const std::array<double, 3>& axis)
    {
        double minimum = std::numeric_limits<double>::max();
        double maximum = std::numeric_limits<double>::lowest();
        for (const auto& point : loop)
        {
            const double projection =
                point[0] * axis[0] + point[1] * axis[1] + point[2] * axis[2];
            minimum = std::min(minimum, projection);
            maximum = std::max(maximum, projection);
        }
        return 0.5 * (maximum - minimum);
    };
    const std::array<double, 3> inletMajor{-0.5, 0.8660254037844386, 0.0};
    const std::array<double, 3> inletMinor{0.0, 0.0, 1.0};
    const double ellipseRatio =
        extent(result.inletLoop, inletMajor) / extent(result.inletLoop, inletMinor);
    EXPECT_NEAR(ellipseRatio, 1.1547005383792517, 0.05);

    const std::array<double, 3> outletY{0.0, 1.0, 0.0};
    const std::array<double, 3> outletZ{0.0, 0.0, 1.0};
    const double circleRatio =
        extent(result.outletLoop, outletY) / extent(result.outletLoop, outletZ);
    EXPECT_NEAR(circleRatio, 1.0, 0.05);
#endif
}
