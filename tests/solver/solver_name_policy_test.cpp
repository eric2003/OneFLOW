#include <gtest/gtest.h>
#include "SolverNamePolicy.h"
#include "SolverNameList.h"
#include "SolverMap.h"
#include <string>
#include <vector>

using namespace ONEFLOW;

TEST( SolverNamePolicy, Plateuns2dslau2_NsBecomesUNs )
{
    EXPECT_EQ( MakeUnstructuredSolverName( "NsSolver" ), "UNsSolver" );
    EXPECT_EQ( MakeStructuredSolverName( "NsSolver" ), "SNsSolver" );
}

TEST( SolverNamePolicy, ExpandPreservesOrder )
{
    std::vector<std::string> base = { "NsSolver", "TurbSolver" };
    auto uns = ExpandSolverNamesForUnstructured( base );
    ASSERT_EQ( uns.size(), 2u );
    EXPECT_EQ( uns[ 0 ], "UNsSolver" );
    EXPECT_EQ( uns[ 1 ], "UTurbSolver" );
}

// ---------------------------------------------------------------------------
// Production-path contract (same steps as SolverNameClass::ReadSolverNames
// after parsing script/solver.txt). No file I/O; mirrors expand + dual lists.
// ---------------------------------------------------------------------------
TEST( SolverNamePolicy, ReadSolverNamesPath_DualListsMatchSafeCloneNames )
{
    // Fake base names as if read from script/solver.txt
    const std::vector<std::string> baseNames = {
        "NsSolver",
        "TurbSolver"
    };

    // Same calls production uses after ReadSolverNames(StringField&)
    const auto uns = ExpandSolverNamesForUnstructured( baseNames );
    const auto str = ExpandSolverNamesForStructured( baseNames );

    ASSERT_EQ( uns.size(), baseNames.size() );
    ASSERT_EQ( str.size(), baseNames.size() );

    // Names handed to Solver::SafeClone for UMESH / SMESH
    EXPECT_EQ( uns[ 0 ], "UNsSolver" );
    EXPECT_EQ( uns[ 1 ], "UTurbSolver" );
    EXPECT_EQ( str[ 0 ], "SNsSolver" );
    EXPECT_EQ( str[ 1 ], "STurbSolver" );

    // Pairwise: each base expands to one U* and one S* at the same index
    for ( size_t i = 0; i < baseNames.size(); ++ i )
    {
        EXPECT_EQ( uns[ i ], MakeUnstructuredSolverName( baseNames[ i ] ) );
        EXPECT_EQ( str[ i ], MakeStructuredSolverName( baseNames[ i ] ) );
    }
}

TEST( SolverNamePolicy, ReadSolverNamesPath_EmptyInputYieldsEmptyLists )
{
    const std::vector<std::string> baseNames;
    const auto uns = ExpandSolverNamesForUnstructured( baseNames );
    const auto str = ExpandSolverNamesForStructured( baseNames );
    EXPECT_TRUE( uns.empty() );
    EXPECT_TRUE( str.empty() );
}

