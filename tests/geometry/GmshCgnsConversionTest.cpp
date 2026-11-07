#include <gtest/gtest.h>
#include "GridConversion.h"
#include "GridTypes.h"
#include "HXCgns.h"
#include <chrono>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#ifdef ENABLE_CGNS
using namespace ONEFLOW;

namespace {
std::filesystem::path MakeTemporaryCaseDirectory() {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto path = std::filesystem::temp_directory_path() /
        ("oneflow-gmsh-cgns-" + std::to_string(stamp));
    std::filesystem::create_directories(path);
    return path;
}
}

TEST(GmshCgnsConversionTest, WritesAndIndependentlyReadsTriangleZone) {
    const std::filesystem::path source =
        std::filesystem::path(ONEFLOW_SOURCE_DIR) /
        "grid" / "gmsh_to_cgns" / "grid" / "square.msh";
    ASSERT_TRUE(std::filesystem::is_regular_file(source));

    const std::filesystem::path caseDirectory = MakeTemporaryCaseDirectory();
    const std::filesystem::path input = caseDirectory / "square.msh";
    const std::filesystem::path resultsDirectory = caseDirectory / "results";
    std::filesystem::create_directories(resultsDirectory);
    const std::filesystem::path output = resultsDirectory / "square.cgns";
    const std::filesystem::path report =
        resultsDirectory / "square.gmsh-read-report.json";
    const std::filesystem::path markdownReport =
        resultsDirectory / "square.gmsh-report.md";
    std::filesystem::copy_file(source, input);

    ONEFLOW::GridConfig config;
    config.sourceType = ONEFLOW::GridFileType::Gmsh;
    config.targetType = ONEFLOW::GridFileType::CGNS;
    config.sourceFile = input.filename().string();
    config.targetFile = (std::filesystem::path("results") / output.filename()).string();
    ONEFLOW::ConvertGrid(config, caseDirectory.string());

    ASSERT_TRUE(std::filesystem::is_regular_file(output));
    ASSERT_TRUE(std::filesystem::is_regular_file(report));
    ASSERT_TRUE(std::filesystem::is_regular_file(markdownReport));
    {
        std::ifstream reportInput(markdownReport);
        const std::string markdown{std::istreambuf_iterator<char>(reportInput),
                                   std::istreambuf_iterator<char>()};
        EXPECT_NE(markdown.find("## Physical Groups"), std::string::npos);
        EXPECT_NE(markdown.find("fluid_domain"), std::string::npos);
        EXPECT_NE(markdown.find("BAR_2"), std::string::npos);
        EXPECT_NE(markdown.find("CGNS writer completed successfully"), std::string::npos);
        EXPECT_NE(markdown.find("| inflow | BCInflow |"), std::string::npos);
        EXPECT_NE(markdown.find("| outflow | BCOutflow |"), std::string::npos);
        EXPECT_NE(markdown.find("| symmetry | BCSymmetryPlane |"), std::string::npos);
        EXPECT_NE(markdown.find("| farfield | BCFarfield |"), std::string::npos);
    }

    int fileId = -1;
    ASSERT_EQ(cg_open(output.string().c_str(), CG_MODE_READ, &fileId), CG_OK);
    int baseCount = 0;
    ASSERT_EQ(cg_nbases(fileId, &baseCount), CG_OK);
    ASSERT_EQ(baseCount, 1);

    char baseName[33] = {};
    int cellDimension = 0;
    int physicalDimension = 0;
    ASSERT_EQ(cg_base_read(fileId, 1, baseName, &cellDimension, &physicalDimension), CG_OK);
    EXPECT_EQ(cellDimension, 2);
    EXPECT_EQ(physicalDimension, 3);

    int zoneCount = 0;
    ASSERT_EQ(cg_nzones(fileId, 1, &zoneCount), CG_OK);
    ASSERT_EQ(zoneCount, 1);
    char zoneName[33] = {};
    cgsize_t zoneSize[9] = {};
    ASSERT_EQ(cg_zone_read(fileId, 1, 1, zoneName, zoneSize), CG_OK);
    EXPECT_EQ(zoneSize[0], 4);
    EXPECT_EQ(zoneSize[1], 2);

    int sectionCount = 0;
    ASSERT_EQ(cg_nsections(fileId, 1, 1, &sectionCount), CG_OK);
    ASSERT_EQ(sectionCount, 2);
    char sectionName[33] = {};
    ElementType_t elementType = TRI_3;
    cgsize_t start = 0;
    cgsize_t end = 0;
    int boundaryCount = 0;
    int parentFlag = 0;
    ASSERT_EQ(cg_section_read(fileId, 1, 1, 1, sectionName, &elementType,
                              &start, &end, &boundaryCount, &parentFlag), CG_OK);
    EXPECT_EQ(elementType, TRI_3);
    EXPECT_EQ(start, 1);
    EXPECT_EQ(end, 2);

    cgsize_t connectivity[6] = {};
    ASSERT_EQ(cg_elements_read(fileId, 1, 1, 1, connectivity, nullptr), CG_OK);
    const std::vector<cgsize_t> expected{1, 2, 3, 1, 3, 4};
    EXPECT_EQ(std::vector<cgsize_t>(connectivity, connectivity + 6), expected);

    cgsize_t rangeMinimum[3] = {1, 0, 0};
    cgsize_t rangeMaximum[3] = {4, 0, 0};
    double x[4] = {};
    double y[4] = {};
    ASSERT_EQ(cg_coord_read(fileId, 1, 1, "CoordinateX", RealDouble,
                            rangeMinimum, rangeMaximum, x), CG_OK);
    ASSERT_EQ(cg_coord_read(fileId, 1, 1, "CoordinateY", RealDouble,
                            rangeMinimum, rangeMaximum, y), CG_OK);
    EXPECT_DOUBLE_EQ(x[0], 0.0);
    EXPECT_DOUBLE_EQ(x[1], 1.0);
    EXPECT_DOUBLE_EQ(x[2], 1.0);
    EXPECT_DOUBLE_EQ(x[3], 0.0);
    EXPECT_DOUBLE_EQ(y[0], 0.0);
    EXPECT_DOUBLE_EQ(y[1], 0.0);
    EXPECT_DOUBLE_EQ(y[2], 1.0);
    EXPECT_DOUBLE_EQ(y[3], 1.0);

    char boundarySectionName[33] = {};
    ElementType_t boundaryType = NODE;
    cgsize_t boundaryStart = 0;
    cgsize_t boundaryEnd = 0;
    ASSERT_EQ(cg_section_read(fileId, 1, 1, 2, boundarySectionName,
                              &boundaryType, &boundaryStart, &boundaryEnd,
                              &boundaryCount, &parentFlag), CG_OK);
    EXPECT_EQ(boundaryType, BAR_2);
    EXPECT_EQ(boundaryStart, 3);
    EXPECT_EQ(boundaryEnd, 6);

    int bocoCount = 0;
    ASSERT_EQ(cg_nbocos(fileId, 1, 1, &bocoCount), CG_OK);
    ASSERT_EQ(bocoCount, 4);
    cgsize_t boundaryConnectivity[8] = {};
    ASSERT_EQ(cg_elements_read(fileId, 1, 1, 2, boundaryConnectivity, nullptr), CG_OK);

    struct ExpectedBoundary
    {
        const char * name;
        BCType_t type;
        cgsize_t nodeA;
        cgsize_t nodeB;
    };
    const std::vector<ExpectedBoundary> expectedBoundaries{
        {"inflow", BCInflow, 4, 1},
        {"outflow", BCOutflow, 2, 3},
        {"symmetry", BCSymmetryPlane, 1, 2},
        {"farfield", BCFarfield, 3, 4}
    };
    std::vector<cgsize_t> observedBoundaryElementIds;
    for (int bocoIndex = 1; bocoIndex <= bocoCount; ++bocoIndex)
    {
        char bocoName[33] = {};
        BCType_t bocoType = BCTypeNull;
        PointSetType_t pointSetType = PointSetTypeNull;
        cgsize_t pointCount = 0;
        int normalIndex[3] = {};
        cgsize_t normalListSize = 0;
        DataType_t normalDataType = DataTypeNull;
        int dataSetCount = 0;
        ASSERT_EQ(cg_boco_info(fileId, 1, 1, bocoIndex, bocoName, &bocoType,
                               &pointSetType, &pointCount, normalIndex,
                               &normalListSize, &normalDataType, &dataSetCount), CG_OK);
        ASSERT_EQ(pointSetType, PointList);
        ASSERT_EQ(pointCount, 1);
        cgsize_t boundaryElementId = 0;
        ASSERT_EQ(cg_boco_read(fileId, 1, 1, bocoIndex, &boundaryElementId, nullptr), CG_OK);
        observedBoundaryElementIds.push_back(boundaryElementId);

        const auto expected = std::find_if(
            expectedBoundaries.begin(), expectedBoundaries.end(),
            [bocoName](const ExpectedBoundary & item) {
                return std::string(item.name) == bocoName;
            });
        ASSERT_NE(expected, expectedBoundaries.end()) << "Unexpected BC_t name: " << bocoName;
        EXPECT_EQ(bocoType, expected->type);
        ASSERT_GE(boundaryElementId, 3);
        ASSERT_LE(boundaryElementId, 6);
        const std::size_t connectivityOffset =
            static_cast<std::size_t>(boundaryElementId - 3) * 2;
        const cgsize_t node0 = boundaryConnectivity[connectivityOffset];
        const cgsize_t node1 = boundaryConnectivity[connectivityOffset + 1];
        EXPECT_EQ(std::minmax(node0, node1),
                  std::minmax(expected->nodeA, expected->nodeB))
            << "Boundary group does not reference its expected geometric edge: " << bocoName;
    }
    std::sort(observedBoundaryElementIds.begin(), observedBoundaryElementIds.end());
    EXPECT_EQ(observedBoundaryElementIds, (std::vector<cgsize_t>{3, 4, 5, 6}));
    EXPECT_EQ(cg_close(fileId), CG_OK);

    std::filesystem::remove_all(caseDirectory);
}
#else
TEST(GmshCgnsConversionTest, CGNSDisabled) {
    GTEST_SKIP() << "CGNS is not enabled in this build.";
}
#endif
