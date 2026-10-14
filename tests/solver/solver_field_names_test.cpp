#include <gtest/gtest.h>
#include "NsFieldNames.h"
#include "INsFieldNames.h"
#include "TurbFieldNames.h"
#include "SolverInfo.h"
#include <string>

using namespace ONEFLOW;

// ---------------------------------------------------------------------------
// Per-solver tables live in ns/ins/turb modules (happy ignorance).
// ---------------------------------------------------------------------------

TEST( SolverFieldNames, NsMatchesAllocUnsteadyTxt )
{
    const auto u = NsSolverUnsteadyFieldNames();
    ASSERT_EQ( u.flow.size(), 3u );
    EXPECT_EQ( u.flow[ 0 ], "q" );
    EXPECT_EQ( u.flow[ 1 ], "q1" );
    EXPECT_EQ( u.flow[ 2 ], "q2" );
    EXPECT_EQ( u.residual[ 0 ], "res" );
    EXPECT_EQ( u.residual[ 1 ], "res1" );
}

TEST( SolverFieldNames, TurbUsesPrefixNotBareQ )
{
    const auto u = TurbSolverUnsteadyFieldNames();
    const auto L = TurbSolverFieldLogisticsNames();
    EXPECT_EQ( u.flow[ 0 ], "turbq" );
    EXPECT_EQ( L.residualName, "turbres" );
    EXPECT_EQ( L.limiterName, "turblimiter" );
    EXPECT_TRUE( L.bcSolutionName.empty() );
}

TEST( SolverFieldNames, INsSharesNsStyleNames )
{
    const auto ns = NsSolverUnsteadyFieldNames();
    const auto ins = INsSolverUnsteadyFieldNames();
    EXPECT_EQ( ns.flow[ 0 ], ins.flow[ 0 ] );
    EXPECT_EQ( ns.residual[ 0 ], ins.residual[ 0 ] );
}

TEST( SolverFieldNames, NsAndTurbSolutionNamesDiffer )
{
    EXPECT_NE(
        NsSolverFieldLogisticsNames().solutionGradNames[ 0 ],
        TurbSolverFieldLogisticsNames().solutionGradNames[ 0 ] );
}

TEST( SolverFieldNames, Apply_NsModuleTablesOntoSolverInfo )
{
    SolverInfo info;
    ApplyFieldNameContract(
        &info,
        NsSolverFieldLogisticsNames(),
        NsSolverUnsteadyFieldNames() );
    EXPECT_EQ( info.residualName, "res" );
    EXPECT_EQ( info.stageSolutionNames[ 0 ], "q1" );
}
