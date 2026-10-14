#include <gtest/gtest.h>

#include "GridSummaryAdapter.h"
#include "UnsGrid.h"
#include "StrGrid.h"

TEST( GridSummaryAdapterTest, DefaultGridHasDeterministicEmptyMetadata )
{
    const ONEFLOW::UnsGrid grid;

    const ONEFLOW::GridSummary summary = ONEFLOW::MakeGridSummary( grid );

    EXPECT_EQ( summary.dimension, ONEFLOW::THREE_D );
    EXPECT_EQ( summary.nodeCount, 0 );
    EXPECT_EQ( summary.faceCount, 0 );
    EXPECT_EQ( summary.cellCount, 0 );
    EXPECT_EQ( summary.boundaryFaceCount, 0 );
    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( summary ) );
}

TEST( GridSummaryAdapterTest, CopiesCurrentGridCountsWithoutOwningStorage )
{
    ONEFLOW::UnsGrid grid;
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
}

TEST( GridSummaryAdapterTest, FreeClearsGridCountsAndAllowsReinitialization )
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    grid.nNodes = 12;
    grid.nFaces = 17;
    grid.nCells = 6;
    grid.nBFaces = 12;
    grid.nIFaces = 2;

    ONEFLOW::Grid & baseGrid = grid;
    baseGrid.Free();

    const ONEFLOW::GridSummary releasedSummary = ONEFLOW::MakeGridSummary( grid );
    EXPECT_EQ( releasedSummary.nodeCount, 0 );
    EXPECT_EQ( releasedSummary.faceCount, 0 );
    EXPECT_EQ( releasedSummary.cellCount, 0 );
    EXPECT_EQ( releasedSummary.boundaryFaceCount, 0 );
    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( releasedSummary ) );

    EXPECT_NO_THROW( grid.Init() );
    const ONEFLOW::GridSummary reinitializedSummary = ONEFLOW::MakeGridSummary( grid );
    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( reinitializedSummary ) );
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
    EXPECT_EQ( grid.bcRegionGroup, nullptr );

    EXPECT_NO_THROW( grid.Init() );
    EXPECT_NE( grid.bcRegionGroup, nullptr );
}
