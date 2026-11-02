#include <gtest/gtest.h>

#include "ULimiterRange.h"

using namespace ONEFLOW;

namespace
{
    bool IsInterfaceType( int bcType )
    {
        return bcType == 7;
    }
}

TEST( LimiterRange, IncludesInterfaceAndInteriorNeighborsAsOffsets )
{
    const RealField q = { 10.0, 20.0, 100.0 };
    const IntField leftCells = { 0, 0 };
    const IntField rightCells = { 2, 1 };
    const IntField boundaryTypes = { 7, 3 };
    RealField minOffsets( 3, 0.0 );
    RealField maxOffsets( 3, 0.0 );
    Real minOffsetRange = 0.0;
    Real maxOffsetRange = 0.0;

    BuildCellNeighborOffsetBounds(
        q, leftCells, rightCells, boundaryTypes, 2, 1, 2,
        IsInterfaceType, minOffsets, maxOffsets, minOffsetRange, maxOffsetRange );

    EXPECT_DOUBLE_EQ( minOffsets[ 0 ], 0.0 );
    EXPECT_DOUBLE_EQ( maxOffsets[ 0 ], 90.0 );
    EXPECT_DOUBLE_EQ( minOffsets[ 1 ], -10.0 );
    EXPECT_DOUBLE_EQ( maxOffsets[ 1 ], 0.0 );

    // Ghost entries remain untouched; only physical-cell offsets enter the reduction.
    EXPECT_DOUBLE_EQ( minOffsets[ 2 ], 0.0 );
    EXPECT_DOUBLE_EQ( maxOffsets[ 2 ], 0.0 );
    EXPECT_DOUBLE_EQ( minOffsetRange, -10.0 );
    EXPECT_DOUBLE_EQ( maxOffsetRange, 90.0 );
}

TEST( LimiterRange, IgnoresNonInterfaceBoundaryFaces )
{
    const RealField q = { 10.0, 20.0, 100.0 };
    const IntField leftCells = { 0, 0 };
    const IntField rightCells = { 2, 1 };
    const IntField boundaryTypes = { 3, 3 };
    RealField minOffsets( 3, 0.0 );
    RealField maxOffsets( 3, 0.0 );
    Real minOffsetRange = 0.0;
    Real maxOffsetRange = 0.0;

    BuildCellNeighborOffsetBounds(
        q, leftCells, rightCells, boundaryTypes, 2, 1, 2,
        IsInterfaceType, minOffsets, maxOffsets, minOffsetRange, maxOffsetRange );

    EXPECT_DOUBLE_EQ( minOffsets[ 0 ], 0.0 );
    EXPECT_DOUBLE_EQ( maxOffsets[ 0 ], 10.0 );
    EXPECT_DOUBLE_EQ( minOffsets[ 1 ], -10.0 );
    EXPECT_DOUBLE_EQ( maxOffsets[ 1 ], 0.0 );
    EXPECT_DOUBLE_EQ( minOffsetRange, -10.0 );
    EXPECT_DOUBLE_EQ( maxOffsetRange, 10.0 );
}
