#include <gtest/gtest.h>
#include "GridConversion.h"
#include "GridTypes.h"
#include "HXCgns.h"
#include <chrono>
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
    const std::filesystem::path output = caseDirectory / "square.cgns";
    const std::filesystem::path report =
        caseDirectory / "square.gmsh-read-report.json";
    const std::filesystem::path markdownReport =
        caseDirectory / "square.gmsh-report.md";
    std::filesystem::copy_file(source, input);

    ONEFLOW::GridConfig config;
    config.sourceType = ONEFLOW::GridFileType::Gmsh;
    config.targetType = ONEFLOW::GridFileType::CGNS;
    config.sourceFile = input.filename().string();
    config.targetFile = output.filename().string();
    ONEFLOW::ConvertGrid(config, caseDirectory.string());

    ASSERT_TRUE(std::filesystem::is_regular_file(output));
    ASSERT_TRUE(std::filesystem::is_regular_file(report));
    ASSERT_TRUE(std::filesystem::is_regular_file(markdownReport));
    {
        std::ifstream reportInput(markdownReport);
        const std::string markdown{std::istreambuf_iterator<char>(reportInput),
                                   std::istreambuf_iterator<char>()};
        EXPECT_NE(markdown.find("## Physical Groups"), std::string::npos);
        EXPECT_NE(markdown.find("outer boundary"), std::string::npos);
        EXPECT_NE(markdown.find("BAR_2"), std::string::npos);
        EXPECT_NE(markdown.find("CGNS writer completed successfully"), std::string::npos);
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
    ASSERT_EQ(bocoCount, 1);
    char bocoName[33] = {};
    BCType_t bocoType = BCTypeNull;
    PointSetType_t pointSetType = PointSetTypeNull;
    cgsize_t pointCount = 0;
    cgsize_t normalIndex[3] = {};
    cgsize_t normalListSize = 0;
    DataType_t normalDataType = DataTypeNull;
    int dataSetCount = 0;
    ASSERT_EQ(cg_boco_info(fileId, 1, 1, 1, bocoName, &bocoType,
                           &pointSetType, &pointCount, normalIndex,
                           &normalListSize, &normalDataType, &dataSetCount), CG_OK);
    EXPECT_STREQ(bocoName, "outer boundary");
    EXPECT_EQ(bocoType, BCTypeUserDefined);
    EXPECT_EQ(pointSetType, PointList);
    EXPECT_EQ(pointCount, 4);
    cgsize_t boundaryElements[4] = {};
    ASSERT_EQ(cg_boco_read(fileId, 1, 1, 1, boundaryElements, nullptr), CG_OK);
    EXPECT_EQ(std::vector<cgsize_t>(boundaryElements, boundaryElements + 4),
              (std::vector<cgsize_t>{3, 4, 5, 6}));
    EXPECT_EQ(cg_close(fileId), CG_OK);

    std::filesystem::remove_all(caseDirectory);
}
#else
TEST(GmshCgnsConversionTest, CGNSDisabled) {
    GTEST_SKIP() << "CGNS is not enabled in this build.";
}
#endif
