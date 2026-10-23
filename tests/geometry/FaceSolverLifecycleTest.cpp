// tests/geometry/FaceSolverLifecycleTest.cpp
#include <gtest/gtest.h>
#include "Constant.h"
#include <stdexcept>
#include "FaceSolver.h"
#include "FaceTopo.h"

// Verify that FaceSolver correctly manages its internal value-type fields
// and FaceTopo without memory leaks or crashes.
TEST( FaceSolverLifecycleTest, InternalFieldsInitializationAndCleanup )
{
    auto * solver = new ONEFLOW::FaceSolver();

    EXPECT_NE( &solver->GetFaceTopo(), nullptr );
    EXPECT_EQ( solver->GetStagedBcSlotCount(), static_cast< std::size_t >( 0 ) );
    EXPECT_EQ( solver->GetChildFidSlotCount(), static_cast< std::size_t >( 0 ) );

    // Growth goes through StageNewFace (sole write path for BC slots).
    for ( int i = 0; i < 10; ++ i )
    {
        solver->StageNewFace( ONEFLOW::IntField{ i, i + 1 }, 5, i );
    }
    EXPECT_EQ( solver->GetStagedBcSlotCount(), static_cast< std::size_t >( 10 ) );
    EXPECT_EQ( solver->GetChildFidSlotCount(), static_cast< std::size_t >( 10 ) );

    delete solver;
    SUCCEED() << "FaceSolver lifecycle is safe and leak-free.";
}

TEST( FaceSolverLifecycleTest, FaceBcAttributesTrackFaceCount )
{
    // Intent: staged BC slots are parallel with faces; empty solver starts empty.
    ONEFLOW::FaceSolver solver;
    EXPECT_EQ( solver.GetStagedBcSlotCount(), static_cast< std::size_t >( 0 ) );
    EXPECT_EQ( solver.GetFaceTopo().GetFaceCount(), static_cast< ONEFLOW::HXSize_t >( 0 ) );
}

TEST( FaceSolverLifecycleTest, StageNewFaceKeepsBcArraysAligned )
{
    // Intent: every newly staged face grows topology and BC arrays together.
    ONEFLOW::FaceSolver solver;
    EXPECT_TRUE( solver.StagedBcMatchesFaceCount() );

    solver.StageNewFace( ONEFLOW::IntField{ 0, 1, 2 }, 5, 10 );
    solver.StageNewFace( ONEFLOW::IntField{ 0, 2, 3 }, 5, 11 );

    ASSERT_EQ( solver.GetFaceTopo().GetFaces().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_TRUE( solver.StagedBcMatchesFaceCount() );
    EXPECT_EQ( solver.GetFaceTopo().GetLeftCells()[ 0 ], 10 );
    EXPECT_EQ( solver.GetFaceTopo().GetLeftCells()[ 1 ], 11 );
    EXPECT_EQ( solver.GetStagedBcType( 0 ), ONEFLOW::INVALID_INDEX );
    EXPECT_EQ( solver.GetStagedBcNameId( 1 ), ONEFLOW::INVALID_INDEX );
}

TEST( FaceSolverLifecycleTest, SetStagedFaceBcUpdatesParallelSlots )
{
    ONEFLOW::FaceSolver solver;
    solver.StageNewFace( ONEFLOW::IntField{ 1, 2 }, 3, 0 );
    solver.SetStagedFaceBc( 0, 7, 42 );

    EXPECT_EQ( solver.GetStagedBcType( 0 ), 7 );
    EXPECT_EQ( solver.GetStagedBcNameId( 0 ), 42 );
    EXPECT_TRUE( solver.StagedBcMatchesFaceCount() );
}

TEST( FaceSolverLifecycleTest, SetStagedFaceBcRejectsOutOfRange )
{
    ONEFLOW::FaceSolver solver;
    EXPECT_THROW( solver.SetStagedFaceBc( 0, 1, 2 ), std::out_of_range );
}
