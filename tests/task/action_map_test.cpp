// tests/test_action_map.cpp
#include <gtest/gtest.h>
#include "ActionMap.h"

// Test fixture for ActionMap lifecycle management
class ActionMapTest : public ::testing::Test 
{
protected:
    void SetUp() override 
    {
        // Initialize the singleton state before each test
        ONEFLOW::ActionMap::Init();
    }

    void TearDown() override 
    {
        // Free resources after each test to prevent memory leaks
        ONEFLOW::ActionMap::Free();
    }
};

TEST_F(ActionMapTest, RegisterAndGetActionId) 
{
    // Register action using underlying implementation
    ONEFLOW::ActionMap::imp->Register("ComputeFlux");
    ONEFLOW::ActionMap::imp->Register("UpdateBoundary");

    // Verify IDs match registration order
    EXPECT_EQ(ONEFLOW::ActionMap::GetActionId("ComputeFlux"), 0);
    EXPECT_EQ(ONEFLOW::ActionMap::GetActionId("UpdateBoundary"), 1);

    // Verify non-existent action returns -1
    EXPECT_EQ(ONEFLOW::ActionMap::GetActionId("NonExistent"), -1);
}

TEST_F(ActionMapTest, GetActionNameById) 
{
    ONEFLOW::ActionMap::imp->Register("SolveNavierStokes");

    EXPECT_EQ(ONEFLOW::ActionMap::GetActionName(0), "SolveNavierStokes");
    EXPECT_EQ(ONEFLOW::ActionMap::GetActionName(999), "");
}

// ActionMapImp::Register/GetActionId/GetActionName is pure logic:
// takes a string, returns an id, with no MPI and no file I/O (ReadFile can be
// skipped for now). This is the best starting point in the whole codebase
// for writing the first unit test.
TEST( ActionMapImpTest, RegisterAssignsSequentialIds )
{
    ONEFLOW::ActionMapImp imp;
    imp.Register( "solve" );
    imp.Register( "output" );

    EXPECT_EQ( imp.GetActionId( "solve" ), 0 );
    EXPECT_EQ( imp.GetActionId( "output" ), 1 );
    EXPECT_EQ( imp.GetActionName( 0 ), "solve" );
}

TEST( ActionMapImpTest, UnknownNameReturnsNegativeOne )
{
    ONEFLOW::ActionMapImp imp;
    EXPECT_EQ( imp.GetActionId( "not_registered" ), -1 );
}