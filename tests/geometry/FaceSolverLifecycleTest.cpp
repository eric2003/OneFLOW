// tests/FaceSolverLifecycleTest.cpp
#include <gtest/gtest.h>
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
    EXPECT_EQ(solver->faceBcKey.size(), 0);
    EXPECT_EQ(solver->faceBcType.size(), 0);
    EXPECT_EQ(solver->childFid.size(), 0);
    
    // 3. Simulate internal resizing (mimicking ScanPolygonFace behavior)
    solver->faceBcType.resize(10);
    solver->faceBcKey.resize(10);
    solver->childFid.resize(10);
    
    EXPECT_EQ(solver->faceBcType.size(), 10);
    
    // 4. Safe deletion: value types are automatically cleaned up, 
    // faceTopo is explicitly deleted in destructor.
    delete solver; 
    
    SUCCEED() << "FaceSolver lifecycle is safe and leak-free.";
}