#include <gtest/gtest.h>
#include "SolverMap.h"
#include "SolverNameList.h"
#include "GridState.h"
#include <string>

using namespace ONEFLOW;

// Index-map + SelectSolverNames seams only - no SafeClone, no MPI, no file I/O.

class SolverMapIndexTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        SolverMap::ClearIndexMaps();
        SolverNameClass::Reset();
    }

    void TearDown() override
    {
        SolverMap::ClearIndexMaps();
        SolverNameClass::Reset();
    }
};

TEST_F( SolverMapIndexTest, AddSolverInfo_BidirectionalLookup )
{
    // Mimic two solvers registered in order (as BuildSolversInBucket would).
    SolverMap::AddSolverInfo( /*solverType=*/10, /*solverIndex=*/0 );
    SolverMap::AddSolverInfo( /*solverType=*/20, /*solverIndex=*/1 );

    EXPECT_EQ( SolverMap::GetSolverIndexBySolverType( 10 ), 0 );
    EXPECT_EQ( SolverMap::GetSolverIndexBySolverType( 20 ), 1 );
    EXPECT_EQ( SolverMap::GetSolverTypeBySolverIndex( 0 ), 10 );
    EXPECT_EQ( SolverMap::GetSolverTypeBySolverIndex( 1 ), 20 );

    ASSERT_EQ( SolverMap::solverTypes.size(), 2u );
    EXPECT_EQ( SolverMap::solverTypes[ 0 ], 10 );
    EXPECT_EQ( SolverMap::solverTypes[ 1 ], 20 );
}

TEST_F( SolverMapIndexTest, AddSolverInfo_DuplicateTypeKeepsFirstIndex )
{
    SolverMap::AddSolverInfo( 10, 0 );
    SolverMap::AddSolverInfo( 10, 5 );  // same type, different index - first wins

    EXPECT_EQ( SolverMap::GetSolverIndexBySolverType( 10 ), 0 );
    EXPECT_EQ( SolverMap::solverTypes.size(), 1u );
}

TEST_F( SolverMapIndexTest, ClearIndexMaps_EmptiesAll )
{
    SolverMap::AddSolverInfo( 1, 0 );
    SolverMap::AddSolverInfo( 2, 1 );
    SolverMap::ClearIndexMaps();

    EXPECT_TRUE( SolverMap::solverTypes.empty() );
    EXPECT_TRUE( SolverMap::solverTypeToIndex.empty() );
    EXPECT_TRUE( SolverMap::solverIndexToType.empty() );
}

TEST_F( SolverMapIndexTest, SelectSolverNames_InjectedOverridesGlobal )
{
    StringField injected;
    injected.push_back( "UNsSolver" );
    injected.push_back( "UTurbSolver" );

    // Global list deliberately different / empty after Reset
    const StringField & selected =
        SolverMap::SelectSolverNames( ONEFLOW::UMESH, &injected );

    ASSERT_EQ( selected.size(), 2u );
    EXPECT_EQ( selected[ 0 ], "UNsSolver" );
    EXPECT_EQ( selected[ 1 ], "UTurbSolver" );
    // Must be the same object (reference to injected), not a copy of global
    EXPECT_EQ( &selected, &injected );
}

TEST_F( SolverMapIndexTest, SelectSolverNames_NullUsesSolverNameClass )
{
    StringField base;
    base.push_back( "NsSolver" );
    SolverNameClass::LoadFromBaseNames( base );

    const StringField & selected =
        SolverMap::SelectSolverNames( ONEFLOW::UMESH, nullptr );

    ASSERT_EQ( selected.size(), 1u );
    EXPECT_EQ( selected[ 0 ], "UNsSolver" );
    EXPECT_EQ( &selected, &SolverNameClass::GetSolverNames( ONEFLOW::UMESH ) );
}
