// tests/geometry/FaceSolverStageTest.cpp
// Locks StageNewFace / ClaimStagedFaceAsLeft / AttachStagedFaceAsRight
// invariants for the polygon-style path without a full CGNS scan.

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

    EXPECT_EQ( topo.GetFaceFlags()[ 0 ], FaceSolver::kFaceUnclaimed );
    EXPECT_EQ( topo.GetFaceFlags()[ 1 ], FaceSolver::kFaceUnclaimed );
    EXPECT_EQ( topo.GetLeftCells()[ 1 ], INVALID_INDEX );
    EXPECT_EQ( solver.faceBcTypes[ 0 ], INVALID_INDEX );
    EXPECT_EQ( solver.faceBcNameIds[ 1 ], INVALID_INDEX );
}

TEST( FaceSolverStageTest, ClaimStagedFaceAsLeftAssignsOwnership )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3, 4 }, /*faceType=*/7, INVALID_INDEX );
    ASSERT_TRUE( solver.StagedBcMatchesFaceCount() );

    solver.ClaimStagedFaceAsLeft( /*faceId=*/0, /*cellId=*/9 );

    const FaceTopo & topo = solver.GetFaceTopo();
    EXPECT_EQ( topo.GetFaceFlags()[ 0 ], FaceSolver::kFaceClaimed );
    EXPECT_EQ( topo.GetLeftCells()[ 0 ], 9 );
    EXPECT_EQ( topo.GetRightCells()[ 0 ], INVALID_INDEX );
    EXPECT_EQ( solver.faceBcTypes[ 0 ], INVALID_INDEX );
    EXPECT_EQ( solver.faceBcNameIds[ 0 ], INVALID_INDEX );
    EXPECT_TRUE( solver.StagedBcMatchesFaceCount() );
}

TEST( FaceSolverStageTest, AttachStagedFaceAsRightSetsNeighbor )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, /*faceType=*/5, INVALID_INDEX );
    solver.ClaimStagedFaceAsLeft( 0, /*cellId=*/3 );
    solver.AttachStagedFaceAsRight( 0, /*cellId=*/8 );

    const FaceTopo & topo = solver.GetFaceTopo();
    EXPECT_EQ( topo.GetLeftCells()[ 0 ], 3 );
    EXPECT_EQ( topo.GetRightCells()[ 0 ], 8 );
    EXPECT_EQ( topo.GetFaceFlags()[ 0 ], FaceSolver::kFaceClaimed );
}

TEST( FaceSolverStageTest, ClaimStagedFaceAsLeftRejectsDoubleClaim )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, /*faceType=*/5, INVALID_INDEX );
    solver.ClaimStagedFaceAsLeft( 0, 1 );
    EXPECT_THROW( solver.ClaimStagedFaceAsLeft( 0, 2 ), std::logic_error );
}

TEST( FaceSolverStageTest, AttachStagedFaceAsRightRejectsUnclaimed )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, /*faceType=*/5, INVALID_INDEX );
    EXPECT_THROW( solver.AttachStagedFaceAsRight( 0, 1 ), std::logic_error );
}

TEST( FaceSolverStageTest, SetStagedFaceBcRejectsOutOfRange )
{
    FaceSolver solver;
    EXPECT_THROW( solver.SetStagedFaceBc( 0, 1, 1 ), std::out_of_range );
}
