// tests/geometry/FaceSolverStageTest.cpp
// Locks StageNewFace / SetStagedFaceBc invariants for the polygon-style path
// without a full CGNS scan. Product sources under codes/ are unchanged.

#include <gtest/gtest.h>
#include <stdexcept>

#include "Constant.h"
#include "FaceSolver.h"
#include "FaceTopo.h"

using namespace ONEFLOW;

TEST( FaceSolverStageTest, StageNewFaceGrowsFaceFlagsWithTopology )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, /*faceType=*/5, /*leftCellId=*/0 );
    solver.StageNewFace( IntField{ 4, 5, 6, 7 }, /*faceType=*/7, INVALID_INDEX );

    ASSERT_TRUE( solver.StagedBcMatchesFaceCount() );

    const FaceTopo & topo = solver.GetFaceTopo();
    EXPECT_EQ( topo.GetFaces().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( topo.GetFaceFlags().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( topo.GetLeftCells().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( topo.GetRightCells().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( solver.faceBcTypes.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( solver.faceBcNameIds.size(), static_cast< std::size_t >( 2 ) );

    EXPECT_EQ( topo.GetFaceFlags()[ 0 ], 0 );
    EXPECT_EQ( topo.GetFaceFlags()[ 1 ], 0 );
    EXPECT_EQ( topo.GetLeftCells()[ 1 ], INVALID_INDEX );
    EXPECT_EQ( solver.faceBcTypes[ 0 ], INVALID_INDEX );
    EXPECT_EQ( solver.faceBcNameIds[ 1 ], INVALID_INDEX );
}

// Polygon path: stage with INVALID left cell, then claim ownership and set BC
// defaults only through SetStagedFaceBc (no direct faceBcTypes[id] writes).
TEST( FaceSolverStageTest, PolygonStyleClaimUsesSetStagedFaceBc )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3, 4 }, /*faceType=*/7, INVALID_INDEX );
    ASSERT_TRUE( solver.StagedBcMatchesFaceCount() );

    const int faceId = 0;
    FaceTopo & topo = solver.GetFaceTopo();
    topo.GetFaceFlags()[ faceId ] = 1;
    topo.GetLeftCells()[ faceId ] = /*polyhedron cell*/ 9;
    topo.GetRightCells()[ faceId ] = INVALID_INDEX;
    solver.SetStagedFaceBc( faceId, INVALID_INDEX, INVALID_INDEX );

    EXPECT_EQ( topo.GetLeftCells()[ faceId ], 9 );
    EXPECT_EQ( solver.faceBcTypes[ faceId ], INVALID_INDEX );
    EXPECT_EQ( solver.faceBcNameIds[ faceId ], INVALID_INDEX );
    EXPECT_TRUE( solver.StagedBcMatchesFaceCount() );
}

TEST( FaceSolverStageTest, SetStagedFaceBcRejectsOutOfRange )
{
    FaceSolver solver;
    EXPECT_THROW( solver.SetStagedFaceBc( 0, 1, 1 ), std::out_of_range );
}
