#include <gtest/gtest.h>
#include <stdexcept>

#include "UnsGrid.h"
#include "FaceTopo.h"

namespace
{
void SetFaceCellConnectivity( ONEFLOW::UnsGrid & grid )
{
    grid.nCells = 2;
    grid.nFaces = 3;
    grid.nBFaces = 2;

    auto & faceTopo = grid.GetFaceTopo();
    faceTopo.GetLeftCells().assign( { 0, 1, 0 } );
    faceTopo.GetRightCells().assign( { -1, -1, 1 } );
}
}

TEST(UnsGridTopologyTest, BoundaryFacesReceiveGhostCellsAndInternalAdjacencyIsPreserved)
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    SetFaceCellConnectivity( grid );

    grid.NormalizeBc();

    const auto & faceTopo = grid.GetFaceTopo();
    const auto & leftCells = faceTopo.GetLeftCells();
    const auto & rightCells = faceTopo.GetRightCells();

    ASSERT_EQ( leftCells.size(), 3u );
    ASSERT_EQ( rightCells.size(), 3u );
    EXPECT_EQ( leftCells[ 0 ], 0 );
    EXPECT_EQ( leftCells[ 1 ], 1 );
    EXPECT_EQ( leftCells[ 2 ], 0 );
    EXPECT_EQ( rightCells[ 0 ], 2 );
    EXPECT_EQ( rightCells[ 1 ], 3 );
    EXPECT_EQ( rightCells[ 2 ], 1 );
}

TEST(UnsGridTopologyTest, RejectsInvalidInternalRightCellBeforeChangingBoundaryCells)
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    SetFaceCellConnectivity( grid );
    grid.GetFaceTopo().GetRightCells()[ 2 ] = 2;

    EXPECT_THROW( grid.NormalizeBc(), std::out_of_range );

    const auto & rightCells = grid.GetFaceTopo().GetRightCells();
    EXPECT_EQ( rightCells[ 0 ], -1 );
    EXPECT_EQ( rightCells[ 1 ], -1 );
    EXPECT_EQ( rightCells[ 2 ], 2 );
}

TEST(UnsGridTopologyTest, RejectsInvalidLeftCellBeforeChangingBoundaryCells)
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    SetFaceCellConnectivity( grid );
    grid.GetFaceTopo().GetLeftCells()[ 1 ] = 2;

    EXPECT_THROW( grid.NormalizeBc(), std::out_of_range );

    const auto & rightCells = grid.GetFaceTopo().GetRightCells();
    EXPECT_EQ( rightCells[ 0 ], -1 );
    EXPECT_EQ( rightCells[ 1 ], -1 );
}

TEST(UnsGridTopologyTest, RejectsConnectivityArrayLengthMismatch)
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    SetFaceCellConnectivity( grid );
    grid.GetFaceTopo().GetRightCells().pop_back();

    EXPECT_THROW( grid.NormalizeBc(), std::logic_error );
}
