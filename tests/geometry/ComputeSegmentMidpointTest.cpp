// tests/geometry/ComputeSegmentMidpointTest.cpp
// Locks pure R^3 segment midpoint used by line interface splits.

#include <gtest/gtest.h>

#include "FaceSearch.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( ComputeSegmentMidpointTest, MidpointOfUnitSegment )
{
    Real mx = 0, my = 0, mz = 0;
    ComputeSegmentMidpoint( 0, 0, 0, 2, 0, 0, mx, my, mz );
    EXPECT_DOUBLE_EQ( mx, 1.0 );
    EXPECT_DOUBLE_EQ( my, 0.0 );
    EXPECT_DOUBLE_EQ( mz, 0.0 );
}

TEST( ComputeSegmentMidpointTest, MidpointIsSymmetric )
{
    Real a0, a1, a2, b0, b1, b2;
    ComputeSegmentMidpoint( 1, 2, 3, 4, 5, 6, a0, a1, a2 );
    ComputeSegmentMidpoint( 4, 5, 6, 1, 2, 3, b0, b1, b2 );
    EXPECT_DOUBLE_EQ( a0, b0 );
    EXPECT_DOUBLE_EQ( a1, b1 );
    EXPECT_DOUBLE_EQ( a2, b2 );
    EXPECT_DOUBLE_EQ( a0, 2.5 );
    EXPECT_DOUBLE_EQ( a1, 3.5 );
    EXPECT_DOUBLE_EQ( a2, 4.5 );
}

TEST( ComputeSegmentMidpointTest, DegenerateSegmentIsThePoint )
{
    Real mx, my, mz;
    ComputeSegmentMidpoint( 7, 8, 9, 7, 8, 9, mx, my, mz );
    EXPECT_DOUBLE_EQ( mx, 7.0 );
    EXPECT_DOUBLE_EQ( my, 8.0 );
    EXPECT_DOUBLE_EQ( mz, 9.0 );
}
