#include <gtest/gtest.h>
#include <stdexcept>
#include <limits>
#include <cmath>

#include "ScalarGrid.h"
#include "Boundary.h"
#include "ScalarIFace.h"
#include "DataBook.h"
#include "DataBaseIO.h"
#include "Constant.h"
#include "ElementHome.h"

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

TEST( ScalarGridTopologyTest, WriteGridFaceTopologyRejectsSelfConnectedFace )
{
    ONEFLOW::ScalarGrid grid;
    grid.nNodes = 1;
    grid.nCells = 1;
    grid.nFaces = 1;
    grid.xn.AddData( 0.0 );
    grid.yn.AddData( 0.0 );
    grid.zn.AddData( 0.0 );
    grid.eTypes.AddData( ONEFLOW::BAR_2 );
    grid.faces.AddElem( std::vector< int >{ 0 } );
    grid.fTypes.AddData( ONEFLOW::NODE );
    grid.lc.AddData( 0 );
    grid.rc.AddData( 0 );

    ONEFLOW::DataBook databook;
    EXPECT_THROW( grid.WriteGridFaceTopology( &databook ), std::logic_error );
}

TEST( ScalarGridTopologyTest, FaceTopologyRoundTripReconstructsBoundaryGhostCells )
{
    ONEFLOW::ScalarGrid source;
    source.nNodes = 2;
    source.nCells = 1;
    source.nFaces = 2;
    source.nBFaces = 2;
    source.xn.AddData( 0.0 );
    source.xn.AddData( 1.0 );
    source.yn.AddData( 0.0 );
    source.yn.AddData( 0.0 );
    source.zn.AddData( 0.0 );
    source.zn.AddData( 0.0 );
    source.eTypes.AddData( ONEFLOW::BAR_2 );
    source.faces.AddElem( std::vector< int >{ 0 } );
    source.faces.AddElem( std::vector< int >{ 1 } );
    source.fTypes.AddData( ONEFLOW::NODE );
    source.fTypes.AddData( ONEFLOW::NODE );
    source.lc.AddData( 0 );
    source.lc.AddData( 0 );
    // These right-cell IDs are derived by NormalizeBc, not persisted topology.
    source.rc.AddData( 1 );
    source.rc.AddData( 2 );
    source.bcTypes.AddData( 0 );
    source.bcTypes.AddData( 0 );

    ONEFLOW::DataBook databook;
    source.WriteGridFaceTopology( &databook );
    databook.MoveToBegin();

    ONEFLOW::ScalarGrid restored;
    restored.nNodes = 2;
    restored.nCells = 1;
    restored.nFaces = 2;
    restored.ReadGridFaceTopology( &databook );

    ASSERT_EQ( restored.rc.GetNElements(), 2u );
    EXPECT_EQ( restored.rc[ 0 ], -1 );
    EXPECT_EQ( restored.rc[ 1 ], -1 );

    restored.nBFaces = 2;
    restored.NormalizeBc();
    EXPECT_EQ( restored.rc[ 0 ], 1 );
    EXPECT_EQ( restored.rc[ 1 ], 2 );
}

TEST( ScalarGridTopologyTest, ReadGridFaceTopologyRejectsSelfConnectedFace )
{
    ONEFLOW::ScalarGrid grid;
    grid.nNodes = 1;
    grid.nCells = 1;
    grid.nFaces = 1;

    ONEFLOW::DataBook databook;
    const ONEFLOW::IntField faceTypes{ ONEFLOW::NODE };
    const ONEFLOW::IntField faceNodeCounts{ 1 };
    const ONEFLOW::IntField faceNodeIds{ 0 };
    const ONEFLOW::IntField leftCells{ 0 };
    const ONEFLOW::IntField rightCells{ 0 };
    ONEFLOW::HXWrite( &databook, faceTypes );
    ONEFLOW::HXWrite( &databook, faceNodeCounts );
    ONEFLOW::HXWrite( &databook, faceNodeIds );
    ONEFLOW::HXWrite( &databook, leftCells );
    ONEFLOW::HXWrite( &databook, rightCells );
    databook.MoveToBegin();

    EXPECT_THROW( grid.ReadGridFaceTopology( &databook ), std::runtime_error );
    EXPECT_EQ( grid.faces.GetNElements(), 0u );
    EXPECT_EQ( grid.lc.GetNElements(), 0u );
    EXPECT_EQ( grid.rc.GetNElements(), 0u );
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


TEST( ScalarGridTopologyTest, CalcTopologyResetsInterfaceTopology )
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.AddData( 0.0 );
    grid.xn.AddData( 1.0 );
    grid.yn.AddData( 0.0 );
    grid.yn.AddData( 0.0 );
    grid.zn.AddData( 0.0 );
    grid.zn.AddData( 0.0 );
    grid.elements.AddElem( std::vector< int >{ 0, 1 } );
    grid.eTypes.AddData( ONEFLOW::BAR_2 );

    auto leftBoundary = std::make_unique< ONEFLOW::ScalarBcco >();
    leftBoundary->bcType = ONEFLOW::BCOutflow;
    leftBoundary->AddBcPoint( 0 );
    grid.scalarBccos->AddBcco( std::move( leftBoundary ) );

    auto rightBoundary = std::make_unique< ONEFLOW::ScalarBcco >();
    rightBoundary->bcType = ONEFLOW::BCOutflow;
    rightBoundary->AddBcPoint( 1 );
    grid.scalarBccos->AddBcco( std::move( rightBoundary ) );

    grid.scalarIFace->iglobalfaces = { 17 };
    grid.scalarIFace->zones = { 2 };
    grid.scalarIFace->cells = { 4 };
    grid.scalarIFace->target_interfaces = { 6 };
    grid.scalarIFace->global_to_local_interfaces.emplace( 17, 0 );
    grid.scalarIFace->local_to_global_interfaces.emplace( 0, 17 );
    grid.scalarIFace->interface_to_bcface = { 1 };
    ONEFLOW::ScalarIFaceIJ neighbor;
    neighbor.zonej = 2;
    neighbor.ifaces = { 0 };
    neighbor.recv_ifaces = { 1 };
    grid.scalarIFace->data.push_back( neighbor );

    // Geometry metrics are derived from the face/cell topology and must not
    // survive a successful topology rebuild.
    grid.xfc.AddData( 10.0 );
    grid.yfc.AddData( 11.0 );
    grid.zfc.AddData( 12.0 );
    grid.xfn.AddData( 13.0 );
    grid.yfn.AddData( 14.0 );
    grid.zfn.AddData( 15.0 );
    grid.area.AddData( 16.0 );
    grid.xcc.AddData( 17.0 );
    grid.ycc.AddData( 18.0 );
    grid.zcc.AddData( 19.0 );
    grid.vol.AddData( 20.0 );
    grid.nTCells = 21;

    grid.CalcTopology();

    EXPECT_TRUE( grid.xfc.data.empty() );
    EXPECT_TRUE( grid.yfc.data.empty() );
    EXPECT_TRUE( grid.zfc.data.empty() );
    EXPECT_TRUE( grid.xfn.data.empty() );
    EXPECT_TRUE( grid.yfn.data.empty() );
    EXPECT_TRUE( grid.zfn.data.empty() );
    EXPECT_TRUE( grid.area.data.empty() );
    EXPECT_TRUE( grid.xcc.data.empty() );
    EXPECT_TRUE( grid.ycc.data.empty() );
    EXPECT_TRUE( grid.zcc.data.empty() );
    EXPECT_TRUE( grid.vol.data.empty() );
    EXPECT_EQ( grid.nTCells, 0 );

    EXPECT_TRUE( grid.scalarIFace->iglobalfaces.empty() );
    EXPECT_TRUE( grid.scalarIFace->zones.empty() );
    EXPECT_TRUE( grid.scalarIFace->cells.empty() );
    EXPECT_TRUE( grid.scalarIFace->target_interfaces.empty() );
    EXPECT_TRUE( grid.scalarIFace->global_to_local_interfaces.empty() );
    EXPECT_TRUE( grid.scalarIFace->local_to_global_interfaces.empty() );
    EXPECT_TRUE( grid.scalarIFace->interface_to_bcface.empty() );
    EXPECT_TRUE( grid.scalarIFace->data.empty() );
}


TEST( ScalarGridTopologyTest, PushElementInvalidatesDerivedTopologyAndGeometry )
{
	ONEFLOW::ScalarGrid grid;
	grid.xn.Resize( 2 );
	grid.yn.Resize( 2 );
	grid.zn.Resize( 2 );
	grid.faces.AddElem( std::vector< int >{ 0 } );
	grid.fTypes.AddData( ONEFLOW::NODE );
	grid.lc.AddData( 0 );
	grid.rc.AddData( -1 );
	grid.bcTypes.AddData( ONEFLOW::BCOutflow );
	grid.xfc.AddData( 10.0 );
	grid.xcc.AddData( 20.0 );
	grid.vol.AddData( 30.0 );
	grid.nFaces = 1;
	grid.nBFaces = 1;
	grid.nTCells = 1;

	grid.PushElement( 0, 1, ONEFLOW::BAR_2 );

	ASSERT_EQ( grid.elements.GetNElements(), 1u );
	ASSERT_EQ( grid.eTypes.GetNElements(), 1u );
	EXPECT_TRUE( grid.faces.data.empty() );
	EXPECT_TRUE( grid.fTypes.data.empty() );
	EXPECT_TRUE( grid.lc.data.empty() );
	EXPECT_TRUE( grid.rc.data.empty() );
	EXPECT_TRUE( grid.bcTypes.data.empty() );
	EXPECT_TRUE( grid.xfc.data.empty() );
	EXPECT_TRUE( grid.xcc.data.empty() );
	EXPECT_TRUE( grid.vol.data.empty() );
	EXPECT_EQ( grid.nFaces, 0 );
	EXPECT_EQ( grid.nBFaces, 0 );
	EXPECT_EQ( grid.nTCells, 0 );
}

TEST( ScalarGridTopologyTest, PushBoundaryFaceInvalidatesDerivedTopologyAndGeometry )
{
	ONEFLOW::ScalarGrid grid;
	grid.faces.AddElem( std::vector< int >{ 0 } );
	grid.fTypes.AddData( ONEFLOW::NODE );
	grid.lc.AddData( 0 );
	grid.rc.AddData( -1 );
	grid.bcTypes.AddData( ONEFLOW::BCOutflow );
	grid.xfc.AddData( 10.0 );
	grid.nFaces = 1;
	grid.nBFaces = 1;

	grid.PushBoundaryFace( 0, ONEFLOW::NODE );

	ASSERT_EQ( grid.boundaryElements.GetNElements(), 1u );
	ASSERT_EQ( grid.bcETypes.GetNElements(), 1u );
	EXPECT_TRUE( grid.faces.data.empty() );
	EXPECT_TRUE( grid.lc.data.empty() );
	EXPECT_TRUE( grid.rc.data.empty() );
	EXPECT_TRUE( grid.bcTypes.data.empty() );
	EXPECT_TRUE( grid.xfc.data.empty() );
	EXPECT_EQ( grid.nFaces, 0 );
	EXPECT_EQ( grid.nBFaces, 0 );
}

TEST( ScalarGridTopologyTest, CreateNodesInvalidatesTopologyAndGeometry )
{
	ONEFLOW::ScalarGrid grid;
	grid.xn.AddData( 0.0 );
	grid.yn.AddData( 0.0 );
	grid.zn.AddData( 0.0 );
	grid.faces.AddElem( std::vector< int >{ 0 } );
	grid.fTypes.AddData( ONEFLOW::NODE );
	grid.lc.AddData( 0 );
	grid.rc.AddData( -1 );
	grid.bcTypes.AddData( ONEFLOW::BCOutflow );
	grid.xfc.AddData( 10.0 );
	grid.xcc.AddData( 20.0 );
	grid.vol.AddData( 30.0 );
	grid.nFaces = 1;
	grid.nBFaces = 1;
	grid.nTCells = 1;

	grid.CreateNodes( 2 );

	EXPECT_EQ( grid.xn.GetNElements(), 2u );
	EXPECT_EQ( grid.yn.GetNElements(), 2u );
	EXPECT_EQ( grid.zn.GetNElements(), 2u );
	EXPECT_TRUE( grid.faces.data.empty() );
	EXPECT_TRUE( grid.fTypes.data.empty() );
	EXPECT_TRUE( grid.lc.data.empty() );
	EXPECT_TRUE( grid.rc.data.empty() );
	EXPECT_TRUE( grid.bcTypes.data.empty() );
	EXPECT_TRUE( grid.xfc.data.empty() );
	EXPECT_TRUE( grid.xcc.data.empty() );
	EXPECT_TRUE( grid.vol.data.empty() );
	EXPECT_EQ( grid.nFaces, 0 );
	EXPECT_EQ( grid.nBFaces, 0 );
	EXPECT_EQ( grid.nTCells, 0 );
}

TEST( ScalarGridTopologyTest, CreateNodesRejectsNegativeCountWithoutInvalidatingTopology )
{
	ONEFLOW::ScalarGrid grid;
	grid.faces.AddElem( std::vector< int >{ 0 } );
	grid.lc.AddData( 0 );
	grid.xfc.AddData( 10.0 );
	grid.nFaces = 1;

	EXPECT_THROW( grid.CreateNodes( -1 ), std::invalid_argument );

	ASSERT_EQ( grid.faces.GetNElements(), 1u );
	ASSERT_EQ( grid.lc.GetNElements(), 1u );
	ASSERT_EQ( grid.xfc.GetNElements(), 1u );
	EXPECT_EQ( grid.nFaces, 1 );
}

TEST( ScalarGridTopologyTest, NormalizeRejectsOutOfRangeAdjacencyBeforeMutation )
{
	ONEFLOW::ScalarGrid grid;
	grid.eTypes.AddData( ONEFLOW::BAR_2 );
	grid.faces.AddElem( std::vector< int >{ 0, 1 } );
	grid.lc.AddData( 0 );
	grid.rc.AddData( 2 );

	EXPECT_THROW( grid.Normalize(), std::runtime_error );
	ASSERT_EQ( grid.faces.GetNElements(), 1u );
	EXPECT_EQ( grid.faces[ 0 ], ( std::vector< int >{ 0, 1 } ) );
	EXPECT_EQ( grid.lc[ 0 ], 0 );
	EXPECT_EQ( grid.rc[ 0 ], 2 );
}

TEST( ScalarGridTopologyTest, NormalizeRejectsMissingFaceConnectivityBeforeMutation )
{
	ONEFLOW::ScalarGrid grid;
	grid.eTypes.AddData( ONEFLOW::BAR_2 );
	grid.bcTypes.AddData( 0 );
	grid.lc.AddData( ONEFLOW::INVALID_INDEX );
	grid.rc.AddData( 0 );

	EXPECT_THROW( grid.Normalize(), std::runtime_error );
	EXPECT_EQ( grid.lc[ 0 ], ONEFLOW::INVALID_INDEX );
	EXPECT_EQ( grid.rc[ 0 ], 0 );
}


TEST( ScalarGridTopologyTest, PushElementRejectsOutOfRangeCgnsNodeWithoutInvalidatingDerivedState )
{
	ONEFLOW::ScalarGrid grid;
	grid.xn.Resize( 2 );
	grid.yn.Resize( 2 );
	grid.zn.Resize( 2 );
	grid.elements.AddElem( std::vector< int >{ 0, 1 } );
	grid.eTypes.AddData( ONEFLOW::BAR_2 );
	grid.faces.AddElem( std::vector< int >{ 0 } );
	grid.fTypes.AddData( ONEFLOW::NODE );
	grid.lc.AddData( 0 );
	grid.rc.AddData( -1 );
	grid.bcTypes.AddData( ONEFLOW::BCOutflow );
	grid.xfc.AddData( 10.0 );
	grid.nFaces = 1;
	grid.nBFaces = 1;

	ONEFLOW::CgIntField invalidConnectivity{ 1, 3 };
	EXPECT_THROW( grid.PushElement( invalidConnectivity, ONEFLOW::BAR_2 ),
		std::invalid_argument );

	ASSERT_EQ( grid.elements.GetNElements(), 1u );
	EXPECT_EQ( grid.elements[ 0 ], ( std::vector< int >{ 0, 1 } ) );
	ASSERT_EQ( grid.eTypes.GetNElements(), 1u );
	EXPECT_EQ( grid.eTypes[ 0 ], ONEFLOW::BAR_2 );
	ASSERT_EQ( grid.faces.GetNElements(), 1u );
	EXPECT_EQ( grid.faces[ 0 ], ( std::vector< int >{ 0 } ) );
	ASSERT_EQ( grid.lc.GetNElements(), 1u );
	EXPECT_EQ( grid.lc[ 0 ], 0 );
	ASSERT_EQ( grid.rc.GetNElements(), 1u );
	EXPECT_EQ( grid.rc[ 0 ], -1 );
	ASSERT_EQ( grid.xfc.GetNElements(), 1u );
	EXPECT_DOUBLE_EQ( grid.xfc[ 0 ], 10.0 );
	EXPECT_EQ( grid.nFaces, 1 );
	EXPECT_EQ( grid.nBFaces, 1 );
}


TEST( ScalarGridTopologyTest, PushElementRejectsInvalidNodeWithoutInvalidatingDerivedState )
{
	ONEFLOW::ScalarGrid grid;
	grid.xn.Resize( 2 );
	grid.yn.Resize( 2 );
	grid.zn.Resize( 2 );
	grid.elements.AddElem( std::vector< int >{ 0, 1 } );
	grid.eTypes.AddData( ONEFLOW::BAR_2 );
	grid.faces.AddElem( std::vector< int >{ 0 } );
	grid.fTypes.AddData( ONEFLOW::NODE );
	grid.lc.AddData( 0 );
	grid.rc.AddData( -1 );
	grid.bcTypes.AddData( ONEFLOW::BCOutflow );
	grid.xfc.AddData( 10.0 );
	grid.nFaces = 1;
	grid.nBFaces = 1;

	EXPECT_THROW( grid.PushElement( 0, 2, ONEFLOW::BAR_2 ), std::out_of_range );
	EXPECT_THROW( grid.PushElement( 0, 0, ONEFLOW::BAR_2 ), std::invalid_argument );
	EXPECT_THROW( grid.PushElement( 0, 1, -1 ), std::invalid_argument );

	ASSERT_EQ( grid.elements.GetNElements(), 1u );
	EXPECT_EQ( grid.elements[ 0 ], ( std::vector< int >{ 0, 1 } ) );
	ASSERT_EQ( grid.eTypes.GetNElements(), 1u );
	EXPECT_EQ( grid.eTypes[ 0 ], ONEFLOW::BAR_2 );
	ASSERT_EQ( grid.faces.GetNElements(), 1u );
	EXPECT_EQ( grid.faces[ 0 ], ( std::vector< int >{ 0 } ) );
	ASSERT_EQ( grid.lc.GetNElements(), 1u );
	EXPECT_EQ( grid.lc[ 0 ], 0 );
	ASSERT_EQ( grid.xfc.GetNElements(), 1u );
	EXPECT_DOUBLE_EQ( grid.xfc[ 0 ], 10.0 );
	EXPECT_EQ( grid.nFaces, 1 );
	EXPECT_EQ( grid.nBFaces, 1 );
}

TEST( ScalarGridTopologyTest, GenerateGridRejectsInvalidDomainWithoutClearingExistingMesh )
{
	ONEFLOW::ScalarGrid grid;
	grid.xn.AddData( 12.0 );
	grid.yn.AddData( 0.0 );
	grid.zn.AddData( 0.0 );
	grid.elements.AddElem( std::vector< int >{ 0, 0 } );
	grid.eTypes.AddData( ONEFLOW::BAR_2 );
	grid.faces.AddElem( std::vector< int >{ 0 } );
	grid.fTypes.AddData( ONEFLOW::NODE );
	grid.lc.AddData( 0 );
	grid.rc.AddData( -1 );
	grid.nFaces = 1;
	grid.nBFaces = 1;

	EXPECT_THROW( grid.GenerateGrid( 3, 2.0, 2.0 ), std::invalid_argument );
	EXPECT_THROW( grid.GenerateGrid( 3, 3.0, 2.0 ), std::invalid_argument );

	ASSERT_EQ( grid.xn.GetNElements(), 1u );
	EXPECT_DOUBLE_EQ( grid.xn[ 0 ], 12.0 );
	ASSERT_EQ( grid.elements.GetNElements(), 1u );
	EXPECT_EQ( grid.elements[ 0 ], ( std::vector< int >{ 0, 0 } ) );
	ASSERT_EQ( grid.faces.GetNElements(), 1u );
	EXPECT_EQ( grid.faces[ 0 ], ( std::vector< int >{ 0 } ) );
	EXPECT_EQ( grid.nFaces, 1 );
	EXPECT_EQ( grid.nBFaces, 1 );
}

TEST( ScalarGridTopologyTest, GenerateGridRejectsNonFiniteDomainWithoutClearingExistingMesh )
{
	ONEFLOW::ScalarGrid grid;
	grid.xn.AddData( 7.0 );
	grid.yn.AddData( 0.0 );
	grid.zn.AddData( 0.0 );

	EXPECT_THROW( grid.GenerateGrid( 3, 0.0,
		std::numeric_limits< ONEFLOW::Real >::infinity() ), std::invalid_argument );

	ASSERT_EQ( grid.xn.GetNElements(), 1u );
	EXPECT_DOUBLE_EQ( grid.xn[ 0 ], 7.0 );
}

TEST( ScalarGridTopologyTest, GenerateGridRejectsUnrepresentableSpacingWithoutClearingExistingMesh )
{
	ONEFLOW::ScalarGrid grid;
	grid.xn.AddData( 42.0 );
	grid.yn.AddData( 0.0 );
	grid.zn.AddData( 0.0 );
	grid.elements.AddElem( std::vector< int >{ 0, 0 } );
	grid.eTypes.AddData( ONEFLOW::BAR_2 );

	const ONEFLOW::Real xmin = 1.0;
	const ONEFLOW::Real xmax = std::nextafter( xmin, 2.0 );
	EXPECT_THROW( grid.GenerateGrid( 3, xmin, xmax ), std::invalid_argument );

	ASSERT_EQ( grid.xn.GetNElements(), 1u );
	EXPECT_DOUBLE_EQ( grid.xn[ 0 ], 42.0 );
	ASSERT_EQ( grid.elements.GetNElements(), 1u );
	EXPECT_EQ( grid.elements[ 0 ], ( std::vector< int >{ 0, 0 } ) );
}
 
TEST( ScalarGridTopologyTest, GenerateGridBuildsOneDimensionalUnstructuredMesh )
{
	ONEFLOW::ScalarGrid grid;
	grid.GenerateGrid( 4, 0.0, 3.0 );

	ASSERT_EQ( grid.GetNNodes(), 4 );
	ASSERT_EQ( grid.GetNCells(), 3 );
	EXPECT_EQ( grid.xn[ 0 ], 0.0 );
	EXPECT_EQ( grid.xn[ 1 ], 1.0 );
	EXPECT_EQ( grid.xn[ 2 ], 2.0 );
	EXPECT_EQ( grid.xn[ 3 ], 3.0 );
	EXPECT_EQ( grid.elements.GetNElements(), 3u );
	EXPECT_EQ( grid.elements[ 0 ], ( std::vector< int >{ 0, 1 } ) );
	EXPECT_EQ( grid.elements[ 1 ], ( std::vector< int >{ 1, 2 } ) );
	EXPECT_EQ( grid.elements[ 2 ], ( std::vector< int >{ 2, 3 } ) );
	EXPECT_EQ( grid.eTypes[ 0 ], ONEFLOW::BAR_2 );
	EXPECT_EQ( grid.eTypes[ 1 ], ONEFLOW::BAR_2 );
	EXPECT_EQ( grid.eTypes[ 2 ], ONEFLOW::BAR_2 );

	// Source-mesh generation does not eagerly construct derived face topology.
	EXPECT_EQ( grid.GetNFaces(), 0 );
	EXPECT_EQ( grid.GetNBFaces(), 0 );
	EXPECT_TRUE( grid.xfc.data.empty() );
	ASSERT_EQ( grid.scalarBccos->bccos.size(), 2u );
	EXPECT_EQ( grid.scalarBccos->bccos[ 0 ]->vertexList, ( std::vector< int >{ 0 } ) );
	EXPECT_EQ( grid.scalarBccos->bccos[ 1 ]->vertexList, ( std::vector< int >{ 3 } ) );
}


TEST( UnitElementInitializationTest, ReinitializationReplacesPreviousTopology )
{
    ONEFLOW::UnitElement element;

    element.Initialize( ONEFLOW::BAR_3 );
    EXPECT_EQ( element.GetElementType(), ONEFLOW::BAR_3 );
    EXPECT_EQ( element.GetElementFaceNumber(), 2 );
    EXPECT_EQ( element.GetChildElementNumbers(), 2 );

    element.Initialize( ONEFLOW::TRI_3 );
    EXPECT_EQ( element.GetElementType(), ONEFLOW::TRI_3 );
    EXPECT_EQ( element.GetElementFaceNumber(), 3 );
    EXPECT_EQ( element.GetChildElementNumbers(), 0 );
}

TEST( UnitElementInitializationTest, UnsupportedTypePreservesExistingTopology )
{
    ONEFLOW::UnitElement element;
    element.Initialize( ONEFLOW::BAR_3 );

    EXPECT_THROW( element.Initialize( -1 ), std::invalid_argument );
    EXPECT_EQ( element.GetElementType(), ONEFLOW::BAR_3 );
    EXPECT_EQ( element.GetElementFaceNumber(), 2 );
    EXPECT_EQ( element.GetChildElementNumbers(), 2 );

    EXPECT_THROW( element.Initialize( 9999 ), std::invalid_argument );
    EXPECT_EQ( element.GetElementType(), ONEFLOW::BAR_3 );
    EXPECT_EQ( element.GetElementFaceNumber(), 2 );
    EXPECT_EQ( element.GetChildElementNumbers(), 2 );
}

TEST( ElementDimensionTest, ReportsIntrinsicDimensionIndependentOfMeshDimension )
{
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::NODE ), 0 );

    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::BAR_2 ), 1 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::BAR_3 ), 1 );

    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::TRI_3 ), 2 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::TRI_6 ), 2 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::QUAD_4 ), 2 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::QUAD_8 ), 2 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::QUAD_9 ), 2 );

    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::TETRA_4 ), 3 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::TETRA_10 ), 3 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::PYRA_5 ), 3 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::PYRA_14 ), 3 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::PENTA_6 ), 3 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::PENTA_15 ), 3 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::PENTA_18 ), 3 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::HEXA_8 ), 3 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::HEXA_20 ), 3 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::HEXA_27 ), 3 );
}

TEST( ElementDimensionTest, RejectsInvalidOrUnsupportedElementTypes )
{
    EXPECT_EQ( ONEFLOW::GetElementDimension( -1 ), -1 );
    // Use the geometry registry count instead of duplicating the CGNS element-count macro.
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::ElementHome::GetUnitElementSlotCount() ), -1 );
}

TEST( ElementDimensionTest, DimensionLookupUsesRegistryStorageAsSourceOfTruth )
{
    ONEFLOW::ElementHome::Free();

    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::BAR_2 ), -1 );

    ONEFLOW::ElementHome::Initialize();
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::BAR_2 ), 1 );
}

TEST( ElementDimensionTest, RegistryRejectsOutOfRangeElementTypes )
{
    EXPECT_THROW( ONEFLOW::ElementHome::GetUnitElementDefinition( -1 ), std::out_of_range );
    EXPECT_THROW(
        ONEFLOW::ElementHome::GetUnitElementDefinition( ONEFLOW::ElementHome::GetUnitElementSlotCount() ),
        std::out_of_range );
    EXPECT_THROW(
        ONEFLOW::GetElementNodeNumbers( ONEFLOW::ElementHome::GetUnitElementSlotCount() ),
        std::out_of_range );

    EXPECT_FALSE( ONEFLOW::IsBasicVolumeElementType( -1 ) );
    EXPECT_FALSE( ONEFLOW::IsBasicVolumeElementType(
        ONEFLOW::ElementHome::GetUnitElementSlotCount() ) );
}

TEST( ElementDimensionTest, RegistryDistinguishesUnsupportedSlotsFromOutOfRangeIndices )
{
    bool foundUnsupportedSlot = false;
    for ( int elementType = 0;
          elementType < ONEFLOW::ElementHome::GetUnitElementSlotCount();
          ++ elementType )
    {
        if ( ONEFLOW::ElementHome::HasUnitElement( elementType ) ) continue;

        foundUnsupportedSlot = true;
        EXPECT_EQ( ONEFLOW::GetElementNodeNumbers( elementType ), 0 );
        EXPECT_EQ( ONEFLOW::GetElementDimension( elementType ), -1 );
        EXPECT_FALSE( ONEFLOW::IsBasicVolumeElementType( elementType ) );
        EXPECT_THROW(
            ONEFLOW::ElementHome::GetUnitElementDefinition( elementType ),
            std::out_of_range );
        break;
    }

    EXPECT_TRUE( foundUnsupportedSlot );
}

TEST( ElementDimensionTest, RegistryFreeClearsStorageAndCanBeReinitialized )
{
    const int expectedSlotCount = ONEFLOW::ElementHome::GetUnitElementSlotCount();
    ASSERT_GT( expectedSlotCount, 0 );

    ONEFLOW::ElementHome::Free();

    EXPECT_EQ( ONEFLOW::ElementHome::GetUnitElementSlotCount(), 0 );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::BAR_2 ), -1 );
    EXPECT_FALSE( ONEFLOW::IsBasicVolumeElementType( ONEFLOW::BAR_2 ) );
    EXPECT_THROW(
        ONEFLOW::ElementHome::GetUnitElementDefinition( ONEFLOW::BAR_2 ),
        std::out_of_range );

    ONEFLOW::ElementHome::Initialize();

    EXPECT_EQ( ONEFLOW::ElementHome::GetUnitElementSlotCount(), expectedSlotCount );
    EXPECT_EQ( ONEFLOW::GetElementDimension( ONEFLOW::BAR_2 ), 1 );
}

TEST( ElementDimensionTest, RegistryInitializationCanBeRepeated )
{
    const int expectedSlotCount = ONEFLOW::ElementHome::GetUnitElementSlotCount();
    ASSERT_GT( expectedSlotCount, 0 );

    ONEFLOW::ElementHome::Initialize();
    EXPECT_EQ( ONEFLOW::ElementHome::GetUnitElementSlotCount(), expectedSlotCount );
    EXPECT_EQ(
        ONEFLOW::ElementHome::GetUnitElementDefinition( ONEFLOW::BAR_2 ).GetElementType(),
        ONEFLOW::BAR_2 );

    ONEFLOW::ElementHome::Initialize();
    EXPECT_EQ( ONEFLOW::ElementHome::GetUnitElementSlotCount(), expectedSlotCount );
}
