#include <gtest/gtest.h>
#include "GmshMeshReader.h"
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <iterator>
#include <sstream>
namespace {
const std::string mesh=R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
4 4 10 40
0 1 0 1
10
0 0 0
0 2 0 1
20
1 0 0
0 3 0 1
30
1 1 0
0 4 0 1
40
0 1 0
$EndNodes
$Elements
5 6 1 6
1 1 1 1
1 10 20
1 2 1 1
2 20 30
1 3 1 1
3 30 40
1 4 1 1
4 40 10
2 1 2 2
5 10 20 30
6 10 30 40
$EndElements
)MSH";
}
TEST(GmshMeshReaderTest, ReadsCoordinatesEntityAndOrderedConnectivity) {
 auto m=ONEFLOW::ReadGmsh41AsciiMesh(mesh);ASSERT_EQ(m.nodes.size(),4u);ASSERT_EQ(m.elements.size(),6u);
 EXPECT_EQ(m.nodes[0].tag,10);EXPECT_DOUBLE_EQ(m.nodes[2].coordinates[1],1.0);
 EXPECT_EQ(m.elements[4].tag,5);EXPECT_EQ(m.elements[4].type,2);EXPECT_EQ(m.elements[4].entityDimension,2);
 EXPECT_EQ(m.elements[4].nodeTags,(std::vector<long long>{10,20,30}));
}
TEST(GmshMeshReaderTest, RejectsDuplicateNodeTags) {
 auto bad=mesh;auto p=bad.find("20\n1 0 0");ASSERT_NE(p,std::string::npos);bad.replace(p,8,"10\n1 0 0");
 EXPECT_THROW(ONEFLOW::ReadGmsh41AsciiMesh(bad),std::runtime_error);
}
TEST(GmshMeshReaderTest, RejectsUnknownNodeReferences) {
 auto bad=mesh;auto p=bad.find("5 10 20 30");ASSERT_NE(p,std::string::npos);bad.replace(p,10,"5 10 20 999");
 EXPECT_THROW(ONEFLOW::ReadGmsh41AsciiMesh(bad),std::runtime_error);
}
TEST(GmshMeshReaderTest, RejectsUnsupportedElementTypes) {
 auto bad=mesh;auto p=bad.find("2 1 2 2");ASSERT_NE(p,std::string::npos);bad.replace(p,7,"2 1 99 2");
 EXPECT_THROW(ONEFLOW::ReadGmsh41AsciiMesh(bad),std::runtime_error);
}


TEST(GmshMeshReaderTest, ReadsBundledSquareMshFileAndBuildsDiagnosticReport) {
    const std::string path =
        std::string(ONEFLOW_SOURCE_DIR) + "/grid/conversion/gmsh_to_cgns/grid/square.msh";
    std::ifstream input(path);
    ASSERT_TRUE(input.good()) << path;
    const std::string content{
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    };
    const auto mesh = ONEFLOW::ReadGmsh41AsciiMesh(content);
    ASSERT_EQ(mesh.nodes.size(), 4u);
    ASSERT_EQ(mesh.elements.size(), 6u);
    ASSERT_EQ(mesh.nodeTagToIndex.size(), 4u);
    EXPECT_EQ(mesh.nodeTagToIndex.at(10), 0u);
    EXPECT_EQ(mesh.nodeTagToIndex.at(40), 3u);
    EXPECT_EQ(mesh.nodes.at(mesh.nodeTagToIndex.at(30)).tag, 30);

    const auto report = ONEFLOW::BuildGmshMeshReadReport(mesh);
    EXPECT_EQ(report.nodeCount, 4u);
    EXPECT_EQ(report.elementCount, 6u);
    EXPECT_EQ(report.elementTypeCounts.at(1), 4u);
    EXPECT_EQ(report.elementTypeCounts.at(2), 2u);
    EXPECT_EQ(report.triangleCount, 2u);
    EXPECT_EQ(report.degenerateTriangleCount, 0u);
    EXPECT_EQ(report.positiveXYTriangleCount, 2u);
    EXPECT_DOUBLE_EQ(report.minimumCoordinates[0], 0.0);
    EXPECT_DOUBLE_EQ(report.maximumCoordinates[0], 1.0);
    EXPECT_DOUBLE_EQ(report.minimumCoordinates[1], 0.0);
    EXPECT_DOUBLE_EQ(report.maximumCoordinates[1], 1.0);
    EXPECT_DOUBLE_EQ(report.minimumTriangleArea, 0.5);
    EXPECT_DOUBLE_EQ(report.maximumTriangleArea, 0.5);

    const std::string json = ONEFLOW::GmshMeshReadReportToJson(report);
    EXPECT_NE(json.find("\"status\": \"READ_SUCCESS\""), std::string::npos);
    EXPECT_NE(json.find("\"physical_groups\": \"PARSED\""), std::string::npos);
    EXPECT_NE(json.find("\"physical_group_count\": 5"), std::string::npos);
    ASSERT_EQ(mesh.physicalNames.size(), 5u);
    EXPECT_EQ(mesh.physicalNames.at({1, 1}), "inflow");
    EXPECT_EQ(mesh.physicalNames.at({1, 2}), "outflow");
    EXPECT_EQ(mesh.physicalNames.at({1, 3}), "symmetry");
    EXPECT_EQ(mesh.physicalNames.at({1, 4}), "farfield");
    EXPECT_EQ(mesh.physicalNames.at({2, 5}), "fluid_domain");
    ASSERT_EQ(mesh.entityPhysicalTags.size(), 9u);
    EXPECT_EQ(mesh.entityPhysicalTags.at({1, 1}), (std::vector<long long>{3}));
    EXPECT_EQ(mesh.entityPhysicalTags.at({1, 2}), (std::vector<long long>{2}));
    EXPECT_EQ(mesh.entityPhysicalTags.at({1, 3}), (std::vector<long long>{4}));
    EXPECT_EQ(mesh.entityPhysicalTags.at({1, 4}), (std::vector<long long>{1}));
    EXPECT_EQ(mesh.entityPhysicalTags.at({2, 1}), (std::vector<long long>{5}));
    const std::string markdown = GmshMeshReportToMarkdown(mesh, "square.cgns");
    EXPECT_NE(markdown.find("## Overall Summary"), std::string::npos);
    EXPECT_NE(markdown.find("## Geometric Entities"), std::string::npos);
    EXPECT_NE(markdown.find("inflow"), std::string::npos);
    EXPECT_NE(markdown.find("4 nodes, 6 elements"), std::string::npos);
}

TEST(GmshMeshReaderTest, ReadsTwoQuadRouteFixtureAndPreservesInternalInterfaceSemantics) {
    const std::string path =
        std::string(ONEFLOW_SOURCE_DIR) + "/grid/conversion/gmsh_to_cgns/grid/two_quads.msh";
    std::ifstream input(path);
    ASSERT_TRUE(input.good()) << path;
    const std::string content{
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    };

    const auto mesh = ONEFLOW::ReadGmsh41AsciiMesh(content);
    ASSERT_EQ(mesh.nodes.size(), 6u);
    ASSERT_EQ(mesh.elements.size(), 9u);
    ASSERT_EQ(mesh.physicalNames.size(), 5u);
    EXPECT_EQ(mesh.physicalNames.at({1, 1}), "inflow");
    EXPECT_EQ(mesh.physicalNames.at({1, 2}), "outflow");
    EXPECT_EQ(mesh.physicalNames.at({1, 3}), "symmetry");
    EXPECT_EQ(mesh.physicalNames.at({1, 4}), "farfield");
    EXPECT_EQ(mesh.physicalNames.at({2, 5}), "fluid_domain");

    const auto report = ONEFLOW::BuildGmshMeshReadReport(mesh);
    EXPECT_EQ(report.nodeCount, 6u);
    EXPECT_EQ(report.elementCount, 9u);
    EXPECT_EQ(report.elementTypeCounts.at(1), 7u);
    EXPECT_EQ(report.elementTypeCounts.at(3), 2u);

    ASSERT_EQ(mesh.entityPhysicalTags.size(), 8u);
    EXPECT_EQ(mesh.entityPhysicalTags.at({1, 1}), (std::vector<long long>{1}));
    EXPECT_EQ(mesh.entityPhysicalTags.at({1, 2}), (std::vector<long long>{3}));
    EXPECT_EQ(mesh.entityPhysicalTags.count({1, 3}), 0u);
    EXPECT_EQ(mesh.entityPhysicalTags.at({1, 4}), (std::vector<long long>{3}));
    EXPECT_EQ(mesh.entityPhysicalTags.at({1, 5}), (std::vector<long long>{2}));
    EXPECT_EQ(mesh.entityPhysicalTags.at({1, 6}), (std::vector<long long>{4}));
    EXPECT_EQ(mesh.entityPhysicalTags.at({1, 7}), (std::vector<long long>{4}));
    EXPECT_EQ(mesh.entityPhysicalTags.at({2, 1}), (std::vector<long long>{5}));
    EXPECT_EQ(mesh.entityPhysicalTags.at({2, 2}), (std::vector<long long>{5}));

    std::size_t quadCells = 0;
    std::size_t boundaryLines = 0;
    for (const auto& element : mesh.elements) {
        if (element.entityDimension == 2 && element.type == 3) {
            ++quadCells;
            EXPECT_EQ(element.nodeTags.size(), 4u);
        }
        if (element.entityDimension == 1 && element.type == 1) {
            ++boundaryLines;
            EXPECT_EQ(element.nodeTags.size(), 2u);
        }
    }
    EXPECT_EQ(quadCells, 2u);
    EXPECT_EQ(boundaryLines, 7u);
}
