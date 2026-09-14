#include <gtest/gtest.h>
#include "SolverNamePolicy.h"

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