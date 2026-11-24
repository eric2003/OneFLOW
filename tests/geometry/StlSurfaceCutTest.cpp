#include "StlSurfaceCut.h"

#include <gtest/gtest.h>
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

TEST(StlSurfaceCutTest, RejectsOpenInputShell)
{
    auto source = MakeOutwardTetrahedron();
    source.triangles.pop_back();
    ONEFLOW::StlSurfaceCutOptions options;
    options.inlet.point = {0.2, 0.0, 0.0};
    options.outlet.point = {0.6, 0.0, 0.0};
    EXPECT_THROW(ONEFLOW::CutStlSurfaceBetweenPlanes(source, options), std::invalid_argument);
}
