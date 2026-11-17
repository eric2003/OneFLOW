// tests/geometry/RebuildFaceCellSidesTest.cpp
// Locks left/right cell rebuild when interface boundary faces split.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceTopo.h"
#include "Boundary.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( RebuildFaceCellSidesTest, NoInterfaceBoundaryUnchangedPattern )
{
    // 2 boundary + 1 interior; nCells = 100
    const IntField bcType{ BC::SOLID_SURFACE, BC::INFLOW };
    const IntField left{ 1, 2, 3 };
    const IntField right{ -1, -1, 4 };
    IntField outL, outR;

    RebuildFaceCellSidesAfterInterfaceSplit(
        bcType, left, right, /*nBFaces=*/2, /*nFaces=*/3,
        IntField{}, /*nCells=*/100, outL, outR );

    EXPECT_EQ( outL, ( IntField{ 1, 2, 3 } ) );
    EXPECT_EQ( outR, ( IntField{ 100, 101, 4 } ) );
}

TEST( RebuildFaceCellSidesTest, SplitInterfaceRepeatsLeftAndGhostRights )
{
    const IntField bcType{ BC::INTERFACE, BC::SOLID_SURFACE };
    const IntField left{ 7, 8, 9 };
    const IntField right{ -1, -1, 10 };
    const IntField childCounts{ 3 };
    IntField outL, outR;

    RebuildFaceCellSidesAfterInterfaceSplit(
        bcType, left, right, 2, 3, childCounts, /*nCells=*/50, outL, outR );

    // boundary: 3 copies of left 7, then left 8; rights 50,51,52,53; interior 9/10
    EXPECT_EQ( outL, ( IntField{ 7, 7, 7, 8, 9 } ) );
    EXPECT_EQ( outR, ( IntField{ 50, 51, 52, 53, 10 } ) );
}

TEST( RebuildFaceCellSidesTest, UnsplitInterfaceEmitsOnce )
{
    const IntField bcType{ BC::INTERFACE };
    const IntField left{ 3 };
    const IntField right{ -1 };
    const IntField childCounts{ 0 };
    IntField outL, outR;

    RebuildFaceCellSidesAfterInterfaceSplit(
        bcType, left, right, 1, 1, childCounts, 20, outL, outR );

    EXPECT_EQ( outL, ( IntField{ 3 } ) );
    EXPECT_EQ( outR, ( IntField{ 20 } ) );
}

TEST( RebuildFaceCellSidesTest, RejectsChildCountShortfall )
{
    const IntField bcType{ BC::INTERFACE };
    const IntField left{ 1 };
    const IntField right{ -1 };
    IntField outL, outR;
    EXPECT_THROW(
        RebuildFaceCellSidesAfterInterfaceSplit(
            bcType, left, right, 1, 1, IntField{}, 10, outL, outR ),
        std::out_of_range );
}
