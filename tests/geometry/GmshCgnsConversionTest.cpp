#include <gtest/gtest.h>
#include "GridConversion.h"
#include "GridTypes.h"
#include "HXCgns.h"
#include "CgnsFactory.h"
#include "CgnsBase.h"
#include "GridElem.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "CellMesh.h"
#include "NodeMesh.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "Dimension.h"
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

class ScopedDimensionRestore
{
public:
    ScopedDimensionRestore() : previousDimension( Dim::GetDimension() ) {}
    ~ScopedDimensionRestore() { Dim::SetDimension( previousDimension ); }

private:
    int previousDimension;
};
}

TEST(CgnsBaseDimensionTest, UsesActiveTopologyDimensionAndThreeCoordinateComponents)
{
    ScopedDimensionRestore restoreDimension;
    ONEFLOW::Dim::SetDimension(ONEFLOW::TWO_D);

    ONEFLOW::CgnsBase base;
    base.SetDefaultCgnsBaseBasicInfo();

    EXPECT_EQ(base.celldim, ONEFLOW::TWO_D);
    EXPECT_EQ(base.phydim, ONEFLOW::THREE_D);
}

TEST(GmshCgnsConversionTest, WritesAndIndependentlyReadsTriangleZone) {
    const std::filesystem::path source =
        std::filesystem::path(ONEFLOW_SOURCE_DIR) /
        "grid" / "conversion" / "gmsh_to_cgns" / "grid" / "square.msh";
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
        GridLocation_t boundaryLocation = GridLocationNull;
        ASSERT_EQ(cg_boco_gridlocation_read(fileId, 1, 1, bocoIndex,
                                            &boundaryLocation), CG_OK);
        EXPECT_EQ(boundaryLocation, EdgeCenter);
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

TEST(GmshCgnsConversionTest, ConvertsGmshThroughCgnsIntoOneFlowGrid) {
    const std::filesystem::path source =
        std::filesystem::path(ONEFLOW_SOURCE_DIR) /
        "grid" / "conversion" / "gmsh_to_cgns" / "grid" / "square.msh";
    ASSERT_TRUE(std::filesystem::is_regular_file(source));

    const std::filesystem::path caseDirectory = MakeTemporaryCaseDirectory();
    const std::filesystem::path input = caseDirectory / "square.msh";
    std::filesystem::copy_file(source, input);

    ONEFLOW::GridConfig config;
    config.sourceType = ONEFLOW::GridFileType::Gmsh;
    config.targetType = ONEFLOW::GridFileType::OneFLOW;
    config.sourceFile = input.filename().string();
    config.topology = ONEFLOW::GridTopology::Unstructured;
    const std::filesystem::path resultsDirectory = caseDirectory / "results";
    std::filesystem::create_directories(resultsDirectory);
    const std::filesystem::path nativeGrid =
        resultsDirectory / "square.oneflow.grid";
    config.targetFile = nativeGrid.string();

    ASSERT_NO_THROW(ONEFLOW::ConvertGrid(config, caseDirectory.string()));

    const std::filesystem::path interchange =
        resultsDirectory / "square.gmsh-to-oneflow.cgns";
    ASSERT_TRUE(std::filesystem::is_regular_file(interchange));
    EXPECT_GT(std::filesystem::file_size(interchange), 0u);
    ASSERT_TRUE(std::filesystem::is_regular_file(nativeGrid));
    EXPECT_GT(std::filesystem::file_size(nativeGrid), 0u);

    // The CGNS interchange must survive the OneFLOW output stage unchanged.
    int fileId = -1;
    ASSERT_EQ(cg_open(interchange.string().c_str(), CG_MODE_READ, &fileId), CG_OK);
    int baseCount = 0;
    ASSERT_EQ(cg_nbases(fileId, &baseCount), CG_OK);
    ASSERT_EQ(baseCount, 1);
    int zoneCount = 0;
    ASSERT_EQ(cg_nzones(fileId, 1, &zoneCount), CG_OK);
    ASSERT_EQ(zoneCount, 1);
    int boundaryCount = 0;
    ASSERT_EQ(cg_nbocos(fileId, 1, 1, &boundaryCount), CG_OK);
    EXPECT_EQ(boundaryCount, 4);
    EXPECT_EQ(cg_close(fileId), CG_OK);

    // Re-run the same CGNS reader boundary and inspect the actual OneFLOW
    // unstructured-grid model, not just the existence of output files.
    ONEFLOW::GridConfig readConfig = config;
    readConfig.sourceType = ONEFLOW::GridFileType::CGNS;
    readConfig.sourceFile = std::filesystem::relative(
        interchange, caseDirectory).generic_string();
    readConfig.targetType = ONEFLOW::GridFileType::OneFLOW;

    ONEFLOW::CgnsFactory cgnsFactory;
    ASSERT_NO_THROW(cgnsFactory.ReadCgnsGrid(readConfig, caseDirectory.string()));
    ASSERT_NO_THROW(cgnsFactory.ProcessCgnsBases());
    ASSERT_NE(cgnsFactory.gridElemGroup, nullptr);
    // Deliberately leave a stale 3D process setting: the CGNS Base cell
    // dimension must take precedence when the 2D unstructured grid is built.
    ScopedDimensionRestore restoreDimension;
    ONEFLOW::Dim::SetDimension(ONEFLOW::THREE_D);
    ONEFLOW::Grids grids =
        cgnsFactory.gridElemGroup->GenerateLocalOneFlowGrids(readConfig);
    EXPECT_EQ(ONEFLOW::Dim::GetDimension(), ONEFLOW::TWO_D);
    ASSERT_EQ(grids.size(), 1u);

    const auto * grid = dynamic_cast<const ONEFLOW::UnsGrid *>(grids.front().get());
    ASSERT_NE(grid, nullptr);
    EXPECT_EQ(grid->nNodes, 4);
    EXPECT_EQ(grid->nCells, 2);
    EXPECT_EQ(grid->nBFaces, 4);
    EXPECT_EQ(grid->GetFaceTopo().GetBoundaryFaceCount(), 4);
    EXPECT_EQ(grid->GetFaceTopo().GetFaceCount(), 5);
    ASSERT_NE(grid->nodeMesh, nullptr);
    ASSERT_EQ(grid->nodeMesh->xN.size(), 4u);
    ASSERT_EQ(grid->nodeMesh->yN.size(), 4u);
    EXPECT_DOUBLE_EQ(grid->nodeMesh->xN[0], 0.0);
    EXPECT_DOUBLE_EQ(grid->nodeMesh->xN[1], 1.0);
    EXPECT_DOUBLE_EQ(grid->nodeMesh->xN[2], 1.0);
    EXPECT_DOUBLE_EQ(grid->nodeMesh->xN[3], 0.0);
    EXPECT_DOUBLE_EQ(grid->nodeMesh->yN[0], 0.0);
    EXPECT_DOUBLE_EQ(grid->nodeMesh->yN[1], 0.0);
    EXPECT_DOUBLE_EQ(grid->nodeMesh->yN[2], 1.0);
    EXPECT_DOUBLE_EQ(grid->nodeMesh->yN[3], 1.0);
    const auto & bcRecord = grid->GetFaceTopo().GetBcRecord();
    ASSERT_EQ(bcRecord.bcType.size(), 4u);
    ASSERT_EQ(bcRecord.bcNameId.size(), 4u);

    struct ExpectedOneFlowBoundary
    {
        const char * name;
        int bcType;
        int nodeA;
        int nodeB;
    };
    const std::vector<ExpectedOneFlowBoundary> expectedOneFlowBoundaries{
        { "inflow", BC::INFLOW, 3, 0 },
        { "outflow", BC::OUTFLOW, 1, 2 },
        { "symmetry", BC::SYMMETRY, 0, 1 },
        { "farfield", BC::FARFIELD, 2, 3 }
    };
    const auto & boundaryFaces = grid->GetFaceTopo().GetFaces();
    for ( const auto & expected : expectedOneFlowBoundaries )
    {
        const int expectedNameId = RegionNameMap::FindRegionId( expected.name );
        ASSERT_GE( expectedNameId, 0 )
            << "Missing registered OneFLOW boundary name: " << expected.name;

        int matchedFace = -1;
        int matchCount = 0;
        for ( int iFace = 0; iFace < grid->nBFaces; ++ iFace )
        {
            if ( bcRecord.bcType[ iFace ] == expected.bcType &&
                 bcRecord.bcNameId[ iFace ] == expectedNameId )
            {
                matchedFace = iFace;
                ++ matchCount;
            }
        }
        ASSERT_EQ( matchCount, 1 )
            << "Expected exactly one boundary face for " << expected.name;
        ASSERT_LT( matchedFace, static_cast<int>( boundaryFaces.size() ) );
        ASSERT_EQ( boundaryFaces[ matchedFace ].size(), 2u );

        std::vector<int> actualNodes = boundaryFaces[ matchedFace ];
        std::sort( actualNodes.begin(), actualNodes.end() );
        std::vector<int> expectedNodes{ expected.nodeA, expected.nodeB };
        std::sort( expectedNodes.begin(), expectedNodes.end() );
        EXPECT_EQ( actualNodes, expectedNodes )
            << "Boundary name/type was attached to the wrong geometric edge: "
            << expected.name;
    }
    EXPECT_EQ(grid->GetCellMesh().GetCellTopo().eTypes.size(), 2u);

    grids.clear();
    std::filesystem::remove_all(caseDirectory);
}


TEST(GmshCgnsConversionTest, ConvertsSingleHexahedronWithFaceCenterBoundaries) {
    const std::filesystem::path source =
        std::filesystem::path(ONEFLOW_SOURCE_DIR) /
        "grid" / "conversion" / "gmsh_to_cgns" / "grid" / "cube_single_hex.msh";
    ASSERT_TRUE(std::filesystem::is_regular_file(source));

    const std::filesystem::path caseDirectory = MakeTemporaryCaseDirectory();
    const std::filesystem::path input = caseDirectory / "cube_single_hex.msh";
    const std::filesystem::path resultsDirectory = caseDirectory / "results";
    std::filesystem::create_directories(resultsDirectory);
    std::filesystem::copy_file(source, input);

    ONEFLOW::GridConfig config;
    config.sourceType = ONEFLOW::GridFileType::Gmsh;
    config.targetType = ONEFLOW::GridFileType::CGNS;
    config.sourceFile = input.filename().string();
    config.targetFile = (std::filesystem::path("results") / "cube_single_hex.cgns").string();
    ONEFLOW::ConvertGrid(config, caseDirectory.string());

    const std::filesystem::path output = resultsDirectory / "cube_single_hex.cgns";
    ASSERT_TRUE(std::filesystem::is_regular_file(output));

    int fileId = -1;
    ASSERT_EQ(cg_open(output.string().c_str(), CG_MODE_READ, &fileId), CG_OK);
    int cellDimension = 0;
    int physicalDimension = 0;
    char baseName[33] = {};
    ASSERT_EQ(cg_base_read(fileId, 1, baseName, &cellDimension, &physicalDimension), CG_OK);
    EXPECT_EQ(cellDimension, 3);
    EXPECT_EQ(physicalDimension, 3);

    char zoneName[33] = {};
    cgsize_t zoneSize[9] = {};
    ASSERT_EQ(cg_zone_read(fileId, 1, 1, zoneName, zoneSize), CG_OK);
    EXPECT_EQ(zoneSize[0], 8);
    EXPECT_EQ(zoneSize[1], 1);
    EXPECT_EQ(zoneSize[2], 0);

    int sectionCount = 0;
    ASSERT_EQ(cg_nsections(fileId, 1, 1, &sectionCount), CG_OK);
    ASSERT_EQ(sectionCount, 2);
    char sectionName[33] = {};
    ElementType_t elementType = ElementTypeNull;
    cgsize_t start = 0;
    cgsize_t end = 0;
    int boundaryCount = 0;
    int parentFlag = 0;
    ASSERT_EQ(cg_section_read(fileId, 1, 1, 1, sectionName, &elementType,
                              &start, &end, &boundaryCount, &parentFlag), CG_OK);
    EXPECT_EQ(elementType, HEXA_8);
    EXPECT_EQ(start, 1);
    EXPECT_EQ(end, 1);
    ASSERT_EQ(cg_section_read(fileId, 1, 1, 2, sectionName, &elementType,
                              &start, &end, &boundaryCount, &parentFlag), CG_OK);
    EXPECT_EQ(elementType, QUAD_4);
    EXPECT_EQ(start, 2);
    EXPECT_EQ(end, 7);

    int bocoCount = 0;
    ASSERT_EQ(cg_nbocos(fileId, 1, 1, &bocoCount), CG_OK);
    ASSERT_EQ(bocoCount, 6);
    struct ExpectedCubeBoundary
    {
        const char * name;
        BCType_t type;
        cgsize_t elementId;
    };
    const std::vector<ExpectedCubeBoundary> expectedBoundaries{
        { "xmin_inflow", BCInflow, 2 },
        { "xmax_outflow", BCOutflow, 3 },
        { "zmin_symmetry", BCSymmetryPlane, 4 },
        { "ymax_farfield", BCFarfield, 5 },
        { "ymin_wall", BCWall, 6 },
        { "zmax_symmetry", BCSymmetryPlane, 7 }
    };
    std::vector<cgsize_t> observedBoundaryElementIds;
    int wallFaceCount = 0;
    for (int bocoIndex = 1; bocoIndex <= bocoCount; ++bocoIndex) {
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
        EXPECT_EQ(pointSetType, PointList);
        GridLocation_t boundaryLocation = GridLocationNull;
        ASSERT_EQ(cg_boco_gridlocation_read(fileId, 1, 1, bocoIndex,
                                            &boundaryLocation), CG_OK);
        EXPECT_EQ(boundaryLocation, FaceCenter);
        ASSERT_EQ(pointCount, 1) << "Each cube face must remain an independent BC_t patch";

        const auto expected = std::find_if(
            expectedBoundaries.begin(), expectedBoundaries.end(),
            [bocoName](const ExpectedCubeBoundary & item) {
                return std::string(item.name) == bocoName;
            });
        ASSERT_NE(expected, expectedBoundaries.end())
            << "Unexpected cube BC_t name: " << bocoName;
        EXPECT_EQ(bocoType, expected->type);
        std::vector<cgsize_t> ids(static_cast<std::size_t>(pointCount));
        ASSERT_EQ(cg_boco_read(fileId, 1, 1, bocoIndex, ids.data(), nullptr), CG_OK);
        EXPECT_EQ(ids.front(), expected->elementId);
        observedBoundaryElementIds.insert(observedBoundaryElementIds.end(),
                                          ids.begin(), ids.end());
        if (std::string(bocoName) == "ymin_wall")
            wallFaceCount += static_cast<int>(pointCount);
    }
    EXPECT_EQ(wallFaceCount, 1);
    std::sort(observedBoundaryElementIds.begin(), observedBoundaryElementIds.end());
    EXPECT_EQ(observedBoundaryElementIds,
              (std::vector<cgsize_t>{2, 3, 4, 5, 6, 7}));
    EXPECT_EQ(cg_close(fileId), CG_OK);

    ONEFLOW::GridConfig readConfig = config;
    readConfig.sourceType = ONEFLOW::GridFileType::CGNS;
    readConfig.sourceFile = "results/cube_single_hex.cgns";
    readConfig.targetType = ONEFLOW::GridFileType::OneFLOW;
    ONEFLOW::CgnsFactory cgnsFactory;
    ASSERT_NO_THROW(cgnsFactory.ReadCgnsGrid(readConfig, caseDirectory.string()));
    ASSERT_NO_THROW(cgnsFactory.ProcessCgnsBases());
    ASSERT_NE(cgnsFactory.gridElemGroup, nullptr);

    ScopedDimensionRestore restoreDimension;
    ONEFLOW::Dim::SetDimension(ONEFLOW::TWO_D);
    ONEFLOW::Grids grids =
        cgnsFactory.gridElemGroup->GenerateLocalOneFlowGrids(readConfig);
    EXPECT_EQ(ONEFLOW::Dim::GetDimension(), ONEFLOW::THREE_D);
    ASSERT_EQ(grids.size(), 1u);
    const auto * grid = dynamic_cast<const ONEFLOW::UnsGrid *>(grids.front().get());
    ASSERT_NE(grid, nullptr);
    EXPECT_EQ(grid->nNodes, 8);
    EXPECT_EQ(grid->nCells, 1);
    EXPECT_EQ(grid->nBFaces, 6);
    EXPECT_EQ(grid->GetFaceTopo().GetBoundaryFaceCount(), 6);
    EXPECT_EQ(grid->GetFaceTopo().GetFaceCount(), 6);
    ASSERT_NE(grid->nodeMesh, nullptr);
    EXPECT_EQ(grid->nodeMesh->xN.size(), 8u);
    EXPECT_EQ(grid->nodeMesh->yN.size(), 8u);
    EXPECT_EQ(grid->nodeMesh->zN.size(), 8u);
    EXPECT_EQ(grid->GetCellMesh().GetCellTopo().eTypes.size(), 1u);

    const auto& cubeBcRecord = grid->GetFaceTopo().GetBcRecord();
    ASSERT_EQ(cubeBcRecord.bcType.size(), 6u);
    ASSERT_EQ(cubeBcRecord.bcNameId.size(), 6u);
    const std::vector<std::pair<std::string, int>> expectedOneFlowBoundaries{
        {"xmin_inflow", BC::INFLOW},
        {"xmax_outflow", BC::OUTFLOW},
        {"ymin_wall", BC::SOLID_SURFACE},
        {"ymax_farfield", BC::FARFIELD},
        {"zmin_symmetry", BC::SYMMETRY},
        {"zmax_symmetry", BC::SYMMETRY},
    };
    for (const auto& expected : expectedOneFlowBoundaries)
    {
        const int expectedNameId = RegionNameMap::FindRegionId(expected.first);
        ASSERT_GE(expectedNameId, 0)
            << "Missing OneFLOW boundary patch name: " << expected.first;
        int matchCount = 0;
        for (std::size_t face = 0; face < cubeBcRecord.bcType.size(); ++face)
        {
            if (cubeBcRecord.bcNameId[face] == expectedNameId &&
                cubeBcRecord.bcType[face] == expected.second)
                ++matchCount;
        }
        EXPECT_EQ(matchCount, 1)
            << "Expected exactly one OneFLOW boundary face for " << expected.first;
    }

    grids.clear();
    std::filesystem::remove_all(caseDirectory);
}

#else
TEST(GmshCgnsConversionTest, CGNSDisabled) {
    GTEST_SKIP() << "CGNS is not enabled in this build.";
}
#endif
