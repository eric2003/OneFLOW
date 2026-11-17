// tests/geometry/BindCellToFaceSidesTest.cpp
// Locks left-then-right cell binding for an existing shared face.

#include <gtest/gtest.h>

#include "FaceSolver.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( BindCellToFaceSidesTest, FillsEmptyLeft )
{
    int left = INVALID_INDEX;
    int right = INVALID_INDEX;
    EXPECT_EQ( BindCellToFaceSides( 7, left, right ), 0 );
    EXPECT_EQ( left, 7 );
    EXPECT_EQ( right, INVALID_INDEX );
}

TEST( BindCellToFaceSidesTest, FillsRightWhenLeftHeldByOtherCell )
{
    int left = 3;
    int right = INVALID_INDEX;
    EXPECT_EQ( BindCellToFaceSides( 8, left, right ), 1 );
    EXPECT_EQ( left, 3 );
    EXPECT_EQ( right, 8 );
}

TEST( BindCellToFaceSidesTest, SameCellOnLeftDoesNotWriteRight )
{
    int left = 4;
    int right = INVALID_INDEX;
    EXPECT_EQ( BindCellToFaceSides( 4, left, right ), -1 );
    EXPECT_EQ( left, 4 );
    EXPECT_EQ( right, INVALID_INDEX );
}

TEST( BindCellToFaceSidesTest, FullyBoundIsNoOp )
{
    int left = 1;
    int right = 2;
    EXPECT_EQ( BindCellToFaceSides( 9, left, right ), -1 );
    EXPECT_EQ( left, 1 );
    EXPECT_EQ( right, 2 );
}
