#include <gtest/gtest.h>
#include <stdexcept>

#include "GridSummaryAdapter.h"
#include "UnsGridStorageValidation.h"
#include "UnsGrid.h"
#include "StrGrid.h"
#include "Dimension.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "InterFace.h"

TEST( GridSummaryAdapterTest, DefaultGridHasDeterministicEmptyMetadata )
{
    const ONEFLOW::UnsGrid grid;

    EXPECT_FALSE( grid.HasMeshStorage() );

    const ONEFLOW::GridSummary summary = ONEFLOW::MakeGridSummary( grid );

    EXPECT_EQ( summary.dimension, ONEFLOW::THREE_D );
    EXPECT_EQ( summary.nodeCount, 0 );
    EXPECT_EQ( summary.faceCount, 0 );
    EXPECT_EQ( summary.cellCount, 0 );
    EXPECT_EQ( summary.boundaryFaceCount, 0 );
    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( summary ) );
    EXPECT_NO_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ) );
}

TEST( GridSummaryAdapterTest, CopiesMetadataWhenCoordinateStorageMatches )
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    grid.nodeMesh->CreateNodes( 12 );
    grid.dimension = 2;
    grid.nNodes = 12;
    grid.nFaces = 17;
    grid.nCells = 6;
    grid.nBFaces = 12;

    const ONEFLOW::GridSummary summary = ONEFLOW::MakeGridSummary( grid );

    EXPECT_EQ( summary.dimension, 2 );
    EXPECT_EQ( summary.nodeCount, 12 );
    EXPECT_EQ( summary.faceCount, 17 );
    EXPECT_EQ( summary.cellCount, 6 );
    EXPECT_EQ( summary.boundaryFaceCount, 12 );
    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( summary ) );
    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );
}

TEST( GridSummaryAdapterTest, RejectsInterfaceCountThatDiffersFromMappingStorage )
{
    ONEFLOW::UnsGrid grid;
    grid.Init();

    grid.nFaces = 1;
    grid.nBFaces = 1;
    grid.nIFaces = 1;

    ONEFLOW::FaceTopo & faceTopo = grid.GetFaceTopo();
    faceTopo.GetFaceTypes().resize( 1 );
    faceTopo.GetFaces().resize( 1 );
    faceTopo.GetLeftCells().resize( 1 );
    faceTopo.GetRightCells().resize( 1 );
    faceTopo.SetNBFaces( 1 );
    faceTopo.GetBcRecord().bcType.resize( 1 );
    faceTopo.GetBcRecord().bcNameId.resize( 1 );

    grid.interFace->Set( 1, &grid );
    EXPECT_NO_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ) );

    grid.interFace->i2b[ 0 ] = 1;
    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );
    grid.interFace->i2b[ 0 ] = 0;

    grid.interFace->localCellId.clear();
    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );
    grid.interFace->localCellId.resize( 1 );

    grid.interFace->i2b.clear();
    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );
}

TEST( GridSummaryAdapterTest, RejectsNegativeMeshCountsBeforeSizeConversion )
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    grid.nFaces = -1;

    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );
}

TEST( GridSummaryAdapterTest, RejectsPositiveNodeCountWithoutCoordinateStorage )
{
    ONEFLOW::UnsGrid grid;
    grid.nNodes = 1;

    EXPECT_THROW( ONEFLOW::MakeGridSummary( grid ), std::invalid_argument );
}

TEST( GridSummaryAdapterTest, RejectsNodeCountThatDiffersFromCoordinateStorage )
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    grid.nodeMesh->CreateNodes( 3 );
    grid.nNodes = 4;

    EXPECT_THROW( ONEFLOW::MakeGridSummary( grid ), std::invalid_argument );
    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );
}

TEST( GridSummaryAdapterTest, RejectsFaceNodeIdsOutsideCoordinateStorage )
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    grid.nodeMesh->CreateNodes( 2 );
    grid.nNodes = 2;
    grid.nFaces = 1;

    ONEFLOW::FaceTopo & faceTopo = grid.GetFaceTopo();
    faceTopo.GetFaceTypes().resize( 1 );
    faceTopo.GetFaces().resize( 1 );
    faceTopo.GetLeftCells().resize( 1 );
    faceTopo.GetRightCells().resize( 1 );
    faceTopo.GetFaces()[ 0 ] = { 2 };

    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );

    faceTopo.GetFaces()[ 0 ][ 0 ] = -1;
    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );
}

TEST( GridSummaryAdapterTest, RejectsCellNodeIdsOutsideCoordinateStorage )
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    grid.nodeMesh->CreateNodes( 3 );
    grid.nNodes = 3;
    grid.nCells = 1;

    ONEFLOW::CellTopo & cellTopo = grid.GetCellMesh().GetCellTopo();
    cellTopo.eTypes.resize( 1 );
    cellTopo.blank.resize( 1 );
    cellTopo.elements.resize( 1 );
    cellTopo.elements[ 0 ] = { 0, 1, 3 };

    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );

    cellTopo.elements[ 0 ] = { 0, -1, 2 };
    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );

    cellTopo.elements[ 0 ] = { 0, 1, 2 };
    EXPECT_NO_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ) );
}

TEST( GridSummaryAdapterTest, RejectsInconsistentMetadata )
{
    ONEFLOW::UnsGrid grid;
    grid.dimension = 2;
    grid.nNodes = 4;
    grid.nFaces = 3;
    grid.nCells = 1;
    grid.nBFaces = 4;

    EXPECT_THROW( ONEFLOW::MakeGridSummary( grid ), std::invalid_argument );
    EXPECT_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ), std::logic_error );

    grid.nBFaces = 2;
    grid.nNodes = -1;
    EXPECT_THROW( ONEFLOW::MakeGridSummary( grid ), std::invalid_argument );
}

TEST( GridSummaryAdapterTest, FreeClearsGridCountsAndAllowsReinitialization )
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    EXPECT_TRUE( grid.HasMeshStorage() );
    grid.nNodes = 12;
    grid.nFaces = 17;
    grid.nCells = 6;
    grid.nBFaces = 12;
    grid.nIFaces = 2;

    ONEFLOW::Grid & baseGrid = grid;
    baseGrid.Free();
    EXPECT_FALSE( grid.HasMeshStorage() );

    const ONEFLOW::GridSummary releasedSummary = ONEFLOW::MakeGridSummary( grid );
    EXPECT_EQ( releasedSummary.nodeCount, 0 );
    EXPECT_EQ( releasedSummary.faceCount, 0 );
    EXPECT_EQ( releasedSummary.cellCount, 0 );
    EXPECT_EQ( releasedSummary.boundaryFaceCount, 0 );
    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( releasedSummary ) );
    EXPECT_NO_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ) );

    EXPECT_NO_THROW( grid.Init() );
    EXPECT_TRUE( grid.HasMeshStorage() );
    const ONEFLOW::GridSummary reinitializedSummary = ONEFLOW::MakeGridSummary( grid );
    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( reinitializedSummary ) );
    EXPECT_NO_THROW( ONEFLOW::ValidateUnsGridStorageCounts( grid ) );
}

TEST( GridSummaryAdapterTest, StructuredGridCleanupCanBeCalledThroughBaseInterface )
{
    ONEFLOW::StrGrid grid;
    grid.Init();
    grid.ni = 4;
    grid.nj = 3;
    grid.nk = 1;
    grid.nNodes = 12;

    ONEFLOW::Grid & baseGrid = grid;
    baseGrid.Free();

    EXPECT_EQ( grid.ni, 0 );
    EXPECT_EQ( grid.nj, 0 );
    EXPECT_EQ( grid.nk, 0 );
    EXPECT_EQ( grid.nNodes, 0 );
    EXPECT_FALSE( static_cast< bool >( grid.bcRegionGroup ) );

    EXPECT_NO_THROW( grid.Init() );
    EXPECT_TRUE( static_cast< bool >( grid.bcRegionGroup ) );
}


TEST( GridSummaryAdapterTest, GenerateI2BRequiresMatchingInterfaceStorage )
{
    ONEFLOW::BcRecord boundaryRecord;
    boundaryRecord.bcType = {
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::SOLID_SURFACE,
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::INTERFACE
    };
    boundaryRecord.bcNameId.resize( 4 );

    ONEFLOW::InterFace interFace;
    interFace.Set( 2 );

    EXPECT_THROW( boundaryRecord.GenerateI2B( interFace ), std::logic_error );

    interFace.Set( 3 );
    EXPECT_NO_THROW( boundaryRecord.GenerateI2B( interFace ) );
    ASSERT_EQ( interFace.i2b.size(), 3 );
    EXPECT_EQ( interFace.i2b[ 0 ], 0 );
    EXPECT_EQ( interFace.i2b[ 1 ], 2 );
    EXPECT_EQ( interFace.i2b[ 2 ], 3 );
    EXPECT_NO_THROW( boundaryRecord.ValidateI2B( interFace ) );

    interFace.i2b[ 1 ] = 1;
    EXPECT_THROW( boundaryRecord.ValidateI2B( interFace ), std::logic_error );
    interFace.i2b[ 1 ] = 2;

    interFace.Set( 2 );
    EXPECT_THROW( boundaryRecord.ValidateI2B( interFace ), std::logic_error );
}
