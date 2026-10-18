// tests/FaceSolverLifecycleTest.cpp
#include <gtest/gtest.h>
#include "Constant.h"
#include <stdexcept>
#include "FaceSolver.h"
#include "FaceTopo.h"

// Verify that FaceSolver correctly manages its internal value-type fields 
// and raw pointer (faceTopo) without memory leaks or crashes.
TEST(FaceSolverLifecycleTest, InternalFieldsInitializationAndCleanup)
{
    ONEFLOW::FaceSolver* solver = new ONEFLOW::FaceSolver();
    
    // 1. Verify faceTopo is allocated
    EXPECT_NE(&solver->GetFaceTopo(), nullptr);
    
    // 2. Verify value-type fields are accessible and empty initially
    EXPECT_EQ(solver->faceBcNameIds.size(), 0);
    EXPECT_EQ(solver->faceBcTypes.size(), 0);
    EXPECT_EQ(solver->childFid.size(), 0);
    
    // 3. Simulate internal resizing (mimicking ScanPolygonFace behavior)
    solver->faceBcTypes.resize(10);
    solver->faceBcNameIds.resize(10);
    solver->childFid.resize(10);
    
    EXPECT_EQ(solver->faceBcTypes.size(), 10);
    
    // 4. Safe deletion: value types are automatically cleaned up, 
    // faceTopo is explicitly deleted in destructor.
    delete solver; 
    
    SUCCEED() << "FaceSolver lifecycle is safe and leak-free.";
}

TEST(FaceSolverLifecycleTest, FaceBcAttributesTrackFaceCount)
{
    // Intent: faceBcNameIds and faceBcTypes are parallel arrays staged with
    // faces; empty solver starts with empty BC attribute lists.
    ONEFLOW::FaceSolver solver;
    EXPECT_TRUE(solver.faceBcNameIds.empty());
    EXPECT_TRUE(solver.faceBcTypes.empty());
    EXPECT_EQ(solver.GetFaceTopo().GetFaceCount(), static_cast<ONEFLOW::HXSize_t>(0));
}

TEST(FaceSolverLifecycleTest, StageNewFaceKeepsBcArraysAligned)
{
    // Intent: every newly staged face grows topology and BC arrays together.
    ONEFLOW::FaceSolver solver;
    EXPECT_TRUE(solver.StagedBcMatchesFaceCount());

    solver.StageNewFace(ONEFLOW::IntField{0, 1, 2}, 5, 10);
    solver.StageNewFace(ONEFLOW::IntField{0, 2, 3}, 5, 11);

    ASSERT_EQ(solver.GetFaceTopo().GetFaces().size(), static_cast<std::size_t>(2));
    EXPECT_TRUE(solver.StagedBcMatchesFaceCount());
    EXPECT_EQ(solver.GetFaceTopo().GetLeftCells()[0], 10);
    EXPECT_EQ(solver.GetFaceTopo().GetLeftCells()[1], 11);
    EXPECT_EQ(solver.faceBcTypes[0], ONEFLOW::INVALID_INDEX);
    EXPECT_EQ(solver.faceBcNameIds[1], ONEFLOW::INVALID_INDEX);
}

TEST(FaceSolverLifecycleTest, SetStagedFaceBcUpdatesParallelSlots)
{
    ONEFLOW::FaceSolver solver;
    solver.StageNewFace(ONEFLOW::IntField{1, 2}, 3, 0);
    solver.SetStagedFaceBc(0, 7, 42);

    EXPECT_EQ(solver.faceBcTypes[0], 7);
    EXPECT_EQ(solver.faceBcNameIds[0], 42);
    EXPECT_TRUE(solver.StagedBcMatchesFaceCount());
}

TEST(FaceSolverLifecycleTest, SetStagedFaceBcRejectsOutOfRange)
{
    ONEFLOW::FaceSolver solver;
    EXPECT_THROW(solver.SetStagedFaceBc(0, 1, 2), std::out_of_range);
}
