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
    // CHANGED: ActionMap::imp is now private; use the public Register()
    // entry point instead (added in this refactor).
    ONEFLOW::ActionMap::Register("ComputeFlux");
    ONEFLOW::ActionMap::Register("UpdateBoundary");

    EXPECT_EQ(ONEFLOW::ActionMap::GetActionId("ComputeFlux"), 0);
    EXPECT_EQ(ONEFLOW::ActionMap::GetActionId("UpdateBoundary"), 1);
    EXPECT_EQ(ONEFLOW::ActionMap::GetActionId("NonExistent"), -1);
}

TEST_F(ActionMapTest, GetActionNameById)
{
    // CHANGED: same as above.
    ONEFLOW::ActionMap::Register("SolveNavierStokes");

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

// New tests covering behavior introduced by this refactor:
// Register() rejecting empty names, GetActionName() with a vector-backed
// out-of-range id, and Clear() resetting state correctly.

TEST( ActionMapImpTest, RegisterIgnoresEmptyName )
{
    ONEFLOW::ActionMapImp imp;
    imp.Register( "" );

    EXPECT_EQ( imp.GetActionId( "" ), -1 );
    EXPECT_EQ( imp.GetActionName( 0 ), "" ); // id 0 was never assigned
}

TEST( ActionMapImpTest, GetActionNameWithNegativeIdReturnsEmptyString )
{
    // Regression-style test for the vector-backed idToName: a negative
    // index must never be used to index into the vector directly.
    ONEFLOW::ActionMapImp imp;
    imp.Register( "solve" );

    EXPECT_EQ( imp.GetActionName( -1 ), "" );
}

TEST( ActionMapImpTest, ClearResetsAllState )
{
    ONEFLOW::ActionMapImp imp;
    imp.Register( "solve" );
    imp.Register( "output" );

    imp.Clear();

    EXPECT_EQ( imp.GetActionId( "solve" ), -1 );
    EXPECT_EQ( imp.GetActionName( 0 ), "" );

    // Registering again after Clear() should restart ids from 0.
    imp.Register( "restart" );
    EXPECT_EQ( imp.GetActionId( "restart" ), 0 );
}

TEST_F( ActionMapTest, FreeThenInitAgainStartsFromCleanState )
{
    // Exercises the new Init()/Free() semantics (Clear()-based instead of
    // delete/new-based) through the static ActionMap facade.
    ONEFLOW::ActionMap::Register( "solve" );
    ONEFLOW::ActionMap::Free();
    ONEFLOW::ActionMap::Init();

    EXPECT_EQ( ONEFLOW::ActionMap::GetActionId( "solve" ), -1 );

    ONEFLOW::ActionMap::Register( "restart" );
    EXPECT_EQ( ONEFLOW::ActionMap::GetActionId( "restart" ), 0 );
}