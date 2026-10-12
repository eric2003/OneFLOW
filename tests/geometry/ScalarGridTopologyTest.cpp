#include <gtest/gtest.h>
#include <stdexcept>

#include "ScalarGrid.h"
#include "Boundary.h"
#include "ScalarIFace.h"
#include "DataBook.h"
#include "DataBaseIO.h"

TEST(ScalarGridTopologyTest, CalcC2CRejectsInvalidInternalCellReferences)
{
    ONEFLOW::ScalarGrid grid;
    grid.eTypes.AddData(1);
    grid.faces.Resize(1);
    grid.lc.AddData(0);
    grid.rc.AddData(1);

    ONEFLOW::EList adjacency;
    EXPECT_THROW(grid.CalcC2C(adjacency), std::runtime_error);
    EXPECT_EQ(adjacency.GetNElements(), 0u);
}

TEST(ScalarGridTopologyTest, CalcC2CRejectsSelfConnectedInternalFace)
{
    ONEFLOW::ScalarGrid grid;
    grid.eTypes.AddData(1);
    grid.faces.Resize(1);
    grid.lc.AddData(0);
    grid.rc.AddData(0);

    ONEFLOW::EList adjacency;
    EXPECT_THROW(grid.CalcC2C(adjacency), std::runtime_error);
    EXPECT_EQ(adjacency.GetNElements(), 0u);
}

TEST(ScalarGridTopologyTest, CalcC2CBuildsSymmetricInternalAdjacency)
{
    ONEFLOW::ScalarGrid grid;
    grid.eTypes.AddData(1);
    grid.eTypes.AddData(1);
    grid.faces.Resize(1);
    grid.lc.AddData(0);
    grid.rc.AddData(1);

    ONEFLOW::EList adjacency;
    grid.CalcC2C(adjacency);

    ASSERT_EQ(adjacency.GetNElements(), 2u);
    ASSERT_EQ(adjacency[0].size(), 1u);
    ASSERT_EQ(adjacency[1].size(), 1u);
    EXPECT_EQ(adjacency[0][0], 1);
    EXPECT_EQ(adjacency[1][0], 0);
}

TEST(ScalarGridTopologyTest, CalcInterfaceToBcFaceRejectsInconsistentFaceArrays)
{
    ONEFLOW::ScalarGrid grid;
    grid.faces.Resize(1);
    grid.lc.AddData(0);
    grid.rc.AddData(1);
    grid.bcTypes.AddData(ONEFLOW::BC::INTERFACE);
    grid.scalarIFace->interface_to_bcface = {42};

    EXPECT_THROW(grid.CalcInterfaceToBcFace(), std::runtime_error);
    ASSERT_EQ(grid.scalarIFace->interface_to_bcface.size(), 1u);
    EXPECT_EQ(grid.scalarIFace->interface_to_bcface[0], 42);
}

TEST(ScalarGridTopologyTest, CalcInterfaceToBcFaceRejectsMissingInterfaceTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.faces.Resize(1);
    grid.lc.AddData(0);
    grid.rc.AddData(1);
    grid.fBcTypes.AddData(ONEFLOW::BC::INTERFACE);
    grid.bcTypes.AddData(ONEFLOW::BC::INTERFACE);
    grid.scalarIFace.reset();

    EXPECT_THROW(grid.CalcInterfaceToBcFace(), std::logic_error);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsUnmatchedBoundaryFaceWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(3);
    grid.yn.Resize(3);
    grid.zn.Resize(3);
    grid.elements.AddElem(std::vector<int>{0, 1});
    grid.elements.AddElem(std::vector<int>{1, 2});
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.faces.Resize(1);
    grid.lc.AddData(23);

    auto boundaryCondition = std::make_unique< ONEFLOW::ScalarBcco >();
    boundaryCondition->AddBcPoint(0);
    grid.scalarBccos->AddBcco(std::move(boundaryCondition));

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 23);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsBoundaryFaceMatchingMultipleConditions)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.yn.Resize(2);
    grid.zn.Resize(2);
    grid.elements.AddElem(std::vector<int>{0, 1});
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.faces.Resize(1);
    grid.lc.AddData(19);

    auto firstBoundary = std::make_unique< ONEFLOW::ScalarBcco >();
    firstBoundary->AddBcPoint(0);
    firstBoundary->AddBcPoint(1);
    grid.scalarBccos->AddBcco(std::move(firstBoundary));

    auto overlappingBoundary = std::make_unique< ONEFLOW::ScalarBcco >();
    overlappingBoundary->AddBcPoint(0);
    grid.scalarBccos->AddBcco(std::move(overlappingBoundary));

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 19);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsNullBoundaryConditionWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.yn.Resize(2);
    grid.zn.Resize(2);
    grid.faces.Resize(1);
    grid.lc.AddData(37);
    grid.scalarBccos->bccos.push_back(nullptr);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 37);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsInvalidBoundaryNodeWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.yn.Resize(2);
    grid.zn.Resize(2);
    grid.faces.Resize(1);
    grid.lc.AddData(41);

    auto boundaryCondition = std::make_unique< ONEFLOW::ScalarBcco >();
    boundaryCondition->AddBcPoint(2);
    grid.scalarBccos->AddBcco(std::move(boundaryCondition));

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 41);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsMismatchedCoordinateArraysWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.yn.Resize(1);
    grid.zn.Resize(2);
    grid.faces.Resize(1);
    grid.lc.AddData(31);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 31);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsMismatchedCellArraysWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.faces.Resize(1);
    grid.lc.AddData(17);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 17);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsOutOfRangeNodeBeforeResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.elements.AddElem(std::vector<int>{0, 2});
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.faces.Resize(1);
    grid.lc.AddData(17);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 17);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsNonManifoldFaceWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(4);
    grid.elements.AddElem(std::vector<int>{0, 1});
    grid.elements.AddElem(std::vector<int>{0, 2});
    grid.elements.AddElem(std::vector<int>{0, 3});
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.faces.Resize(1);
    grid.lc.AddData(23);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 23);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsDuplicateCellNodesWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.elements.AddElem(std::vector<int>{0, 0});
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.faces.Resize(1);
    grid.lc.AddData(29);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 29);
}


namespace
{
void SetScalarFaceAdjacency( ONEFLOW::ScalarGrid & grid )
{
    grid.nCells = 2;
    grid.nFaces = 3;
    grid.nBFaces = 2;
    grid.lc.AddData( 0 );
    grid.lc.AddData( 1 );
    grid.lc.AddData( 0 );
    grid.rc.AddData( -1 );
    grid.rc.AddData( -1 );
    grid.rc.AddData( 1 );
}
}

TEST( ScalarGridTopologyTest, NormalizeBcAssignsGhostCellsAndPreservesInternalAdjacency )
{
    ONEFLOW::ScalarGrid grid;
    SetScalarFaceAdjacency( grid );

    grid.NormalizeBc();

    ASSERT_EQ( grid.rc.GetNElements(), 3u );
    EXPECT_EQ( grid.rc[ 0 ], 2 );
    EXPECT_EQ( grid.rc[ 1 ], 3 );
    EXPECT_EQ( grid.rc[ 2 ], 1 );
}

TEST( ScalarGridTopologyTest, NormalizeBcRejectsSelfConnectedInternalFaceBeforeMutation )
{
    ONEFLOW::ScalarGrid grid;
    SetScalarFaceAdjacency( grid );
    grid.rc[ 2 ] = 0;

    EXPECT_THROW( grid.NormalizeBc(), std::logic_error );

    EXPECT_EQ( grid.rc[ 0 ], -1 );
    EXPECT_EQ( grid.rc[ 1 ], -1 );
    EXPECT_EQ( grid.rc[ 2 ], 0 );
}

TEST( ScalarGridTopologyTest, NormalizeBcRejectsInvalidLeftCellBeforeMutation )
{
    ONEFLOW::ScalarGrid grid;
    SetScalarFaceAdjacency( grid );
    grid.lc[ 1 ] = 2;

    EXPECT_THROW( grid.NormalizeBc(), std::out_of_range );

    EXPECT_EQ( grid.rc[ 0 ], -1 );
    EXPECT_EQ( grid.rc[ 1 ], -1 );
}

TEST( ScalarGridTopologyTest, NormalizeBcRejectsBoundaryFaceWithUnexpectedRightCell )
{
    ONEFLOW::ScalarGrid grid;
    SetScalarFaceAdjacency( grid );
    grid.rc[ 0 ] = 1;

    EXPECT_THROW( grid.NormalizeBc(), std::logic_error );

    EXPECT_EQ( grid.rc[ 0 ], 1 );
    EXPECT_EQ( grid.rc[ 1 ], -1 );
}

TEST( ScalarGridTopologyTest, NormalizeBcRejectsConnectivitySizeMismatch )
{
    ONEFLOW::ScalarGrid grid;
    SetScalarFaceAdjacency( grid );
    grid.rc.data.pop_back();

    EXPECT_THROW( grid.NormalizeBc(), std::logic_error );

    EXPECT_EQ( grid.rc.GetNElements(), 2u );
    EXPECT_EQ( grid.rc[ 0 ], -1 );
    EXPECT_EQ( grid.rc[ 1 ], -1 );
}

TEST( ScalarGridTopologyTest, WriteGridFaceTopologyRejectsInconsistentParallelArrays )
{
    ONEFLOW::ScalarGrid grid;
    grid.nNodes = 2;
    grid.nCells = 1;
    grid.nFaces = 1;
    grid.xn.AddData( 0.0 );
    grid.xn.AddData( 1.0 );
    grid.yn.AddData( 0.0 );
    grid.yn.AddData( 0.0 );
    grid.zn.AddData( 0.0 );
    grid.zn.AddData( 0.0 );
    grid.faces.AddElem( std::vector< int >{ 0 } );
    grid.lc.AddData( 0 );
    grid.rc.AddData( -1 );
    // fTypes intentionally remains empty to violate the face-record invariant.

    ONEFLOW::DataBook databook;
    EXPECT_THROW( grid.WriteGridFaceTopology( &databook ), std::logic_error );
}

TEST( ScalarGridTopologyTest, WriteGridFaceTopologyRejectsDuplicateNodeIndices )
{
    ONEFLOW::ScalarGrid grid;
    grid.nNodes = 2;
    grid.nCells = 1;
    grid.nFaces = 1;
    grid.xn.AddData( 0.0 );
    grid.xn.AddData( 1.0 );
    grid.yn.AddData( 0.0 );
    grid.yn.AddData( 0.0 );
    grid.zn.AddData( 0.0 );
    grid.zn.AddData( 0.0 );
    grid.eTypes.AddData( ONEFLOW::BAR_2 );
    grid.faces.AddElem( std::vector< int >{ 0, 0 } );
    grid.fTypes.AddData( ONEFLOW::NODE );
    grid.lc.AddData( 0 );
    grid.rc.AddData( -1 );

    ONEFLOW::DataBook databook;
    EXPECT_THROW( grid.WriteGridFaceTopology( &databook ), std::logic_error );
}

TEST( ScalarGridTopologyTest, ReadGridFaceTopologyRejectsEmptyFace )
{
    ONEFLOW::ScalarGrid grid;
    grid.nNodes = 2;
    grid.nCells = 1;
    grid.nFaces = 1;

    ONEFLOW::DataBook databook;
    ONEFLOW::IntField faceTypes{ ONEFLOW::NODE };
    ONEFLOW::IntField faceNodeCounts{ 0 };
    ONEFLOW::HXWrite( &databook, faceTypes );
    ONEFLOW::HXWrite( &databook, faceNodeCounts );
    databook.MoveToBegin();

    EXPECT_THROW( grid.ReadGridFaceTopology( &databook ), std::runtime_error );
}

TEST( ScalarGridTopologyTest, ReadGridFaceTopologyRejectsDuplicateNodeIndices )
{
    ONEFLOW::ScalarGrid grid;
    grid.nNodes = 2;
    grid.nCells = 1;
    grid.nFaces = 1;

    ONEFLOW::DataBook databook;
    ONEFLOW::IntField faceTypes{ ONEFLOW::BAR_2 };
    ONEFLOW::IntField faceNodeCounts{ 2 };
    ONEFLOW::IntField faceNodeIds{ 0, 0 };
    ONEFLOW::HXWrite( &databook, faceTypes );
    ONEFLOW::HXWrite( &databook, faceNodeCounts );
    ONEFLOW::HXWrite( &databook, faceNodeIds );
    databook.MoveToBegin();

    EXPECT_THROW( grid.ReadGridFaceTopology( &databook ), std::runtime_error );
    EXPECT_EQ( grid.faces.GetNElements(), 0u );
    EXPECT_EQ( grid.lc.GetNElements(), 0u );
    EXPECT_EQ( grid.rc.GetNElements(), 0u );
    EXPECT_EQ( grid.fTypes.GetNElements(), 0u );
}

TEST( ScalarGridTopologyTest, ReadGridReplacesMeshAndKeepsCountsConsistent )
{
    ONEFLOW::ScalarGrid grid;
    grid.nNodes = 1;
    grid.nCells = 1;
    grid.nFaces = 1;
    grid.nBFaces = 1;
    grid.xn.AddData( 99.0 );
    grid.yn.AddData( 0.0 );
    grid.zn.AddData( 0.0 );
    grid.xfc.AddData( 99.0 );
    grid.eTypes.AddData( ONEFLOW::BAR_2 );

    ONEFLOW::DataBook databook;
    const int nodeCount = 2;
    const int faceCount = 2;
    const int cellCount = 1;
    const std::vector< ONEFLOW::Real > x{ 0.0, 1.0 };
    const std::vector< ONEFLOW::Real > y{ 0.0, 0.0 };
    const std::vector< ONEFLOW::Real > z{ 0.0, 0.0 };
    const int volumeBcType = -1;
    const std::vector< int > elementTypes{ ONEFLOW::BAR_2 };
    const std::vector< int > faceTypes{ ONEFLOW::NODE, ONEFLOW::NODE };
    const std::vector< int > faceNodeCounts{ 1, 1 };
    const std::vector< int > faceNodeIds{ 0, 1 };
    const std::vector< int > leftCells{ 0, 0 };
    const std::vector< int > rightCells{ -1, -1 };
    const int boundaryFaceCount = 2;
    const std::vector< int > boundaryTypes{ 0, 0 };
    const std::vector< int > boundaryNameIds{ 0, 0 };
    const int interfaceCount = 0;

    ONEFLOW::HXWrite( &databook, nodeCount );
    ONEFLOW::HXWrite( &databook, faceCount );
    ONEFLOW::HXWrite( &databook, cellCount );
    ONEFLOW::HXWrite( &databook, x );
    ONEFLOW::HXWrite( &databook, y );
    ONEFLOW::HXWrite( &databook, z );
    ONEFLOW::HXWrite( &databook, volumeBcType );
    ONEFLOW::HXWrite( &databook, elementTypes );
    ONEFLOW::HXWrite( &databook, faceTypes );
    ONEFLOW::HXWrite( &databook, faceNodeCounts );
    ONEFLOW::HXWrite( &databook, faceNodeIds );
    ONEFLOW::HXWrite( &databook, leftCells );
    ONEFLOW::HXWrite( &databook, rightCells );
    ONEFLOW::HXWrite( &databook, boundaryFaceCount );
    ONEFLOW::HXWrite( &databook, boundaryTypes );
    ONEFLOW::HXWrite( &databook, boundaryNameIds );
    ONEFLOW::HXWrite( &databook, interfaceCount );
    databook.MoveToBegin();

    grid.ReadGrid( &databook );

    EXPECT_EQ( grid.nNodes, grid.GetNNodes() );
    EXPECT_EQ( grid.nCells, grid.GetNCells() );
    EXPECT_EQ( grid.nFaces, grid.GetNFaces() );
    EXPECT_EQ( grid.nBFaces, grid.GetNBFaces() );
    ASSERT_EQ( grid.xn.GetNElements(), 2u );
    EXPECT_DOUBLE_EQ( grid.xn[ 0 ], 0.0 );
    EXPECT_DOUBLE_EQ( grid.xn[ 1 ], 1.0 );
    EXPECT_EQ( grid.xfc.GetNElements(), 0u );
}

TEST( ScalarGridTopologyTest, ReadGridPreservesExistingMeshWhenPayloadIsTruncated )
{
	ONEFLOW::ScalarGrid grid;
	grid.nNodes = 1;
	grid.nCells = 1;
	grid.nFaces = 1;
	grid.nBFaces = 1;
	grid.xn.AddData( 12.5 );
	grid.yn.AddData( 23.5 );
	grid.zn.AddData( 34.5 );
	grid.eTypes.AddData( ONEFLOW::BAR_2 );
	grid.faces.AddElem( std::vector< int >{ 0 } );
	grid.bcTypes.AddData( 77 );

	ONEFLOW::DataBook databook;
	const int nodeCount = 2;
	const int faceCount = 1;
	const int cellCount = 1;
	ONEFLOW::HXWrite( &databook, nodeCount );
	ONEFLOW::HXWrite( &databook, faceCount );
	ONEFLOW::HXWrite( &databook, cellCount );
	databook.MoveToBegin();

	EXPECT_THROW( grid.ReadGrid( &databook ), std::exception );

	EXPECT_EQ( grid.nNodes, 1 );
	EXPECT_EQ( grid.nCells, 1 );
	EXPECT_EQ( grid.nFaces, 1 );
	EXPECT_EQ( grid.nBFaces, 1 );
	ASSERT_EQ( grid.xn.GetNElements(), 1u );
	EXPECT_DOUBLE_EQ( grid.xn[ 0 ], 12.5 );
	ASSERT_EQ( grid.eTypes.GetNElements(), 1u );
	EXPECT_EQ( grid.eTypes[ 0 ], ONEFLOW::BAR_2 );
	ASSERT_EQ( grid.bcTypes.GetNElements(), 1u );
	EXPECT_EQ( grid.bcTypes[ 0 ], 77 );
}

TEST( ScalarGridTopologyTest, ReadBoundaryTopologyRejectsInterfaceCountMismatch )
{
    ONEFLOW::ScalarGrid grid;
    grid.nFaces = 1;
    grid.nBFaces = 0;
    grid.bcTypes.AddData( 77 );
    grid.bcNameIds.AddData( 88 );
    grid.scalarIFace->zones.push_back( 123 );

    ONEFLOW::DataBook databook;
    const int nBoundaryFaces = 1;
    ONEFLOW::IntField boundaryTypes{ ONEFLOW::BC::INTERFACE };
    ONEFLOW::IntField boundaryNameIds{ 0 };
    const int nInterfaces = 0;
    ONEFLOW::HXWrite( &databook, nBoundaryFaces );
    ONEFLOW::HXWrite( &databook, boundaryTypes );
    ONEFLOW::HXWrite( &databook, boundaryNameIds );
    ONEFLOW::HXWrite( &databook, nInterfaces );
    databook.MoveToBegin();

    EXPECT_THROW( grid.ReadBoundaryTopology( &databook ), std::runtime_error );
    EXPECT_EQ( grid.nBFaces, 0 );
    ASSERT_EQ( grid.bcTypes.GetNElements(), 1u );
    EXPECT_EQ( grid.bcTypes[ 0 ], 77 );
    ASSERT_EQ( grid.bcNameIds.GetNElements(), 1u );
    EXPECT_EQ( grid.bcNameIds[ 0 ], 88 );
    ASSERT_EQ( grid.scalarIFace->zones.size(), 1u );
    EXPECT_EQ( grid.scalarIFace->zones[ 0 ], 123 );
}
