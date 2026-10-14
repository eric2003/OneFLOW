#include <gtest/gtest.h>
#include "NsFieldNames.h"
#include "INsFieldNames.h"
#include "TurbFieldNames.h"
#include "SolverInfo.h"
#include <string>

using namespace ONEFLOW;

TEST( SolverFieldNames, NsDualTimeNamesHaveExpectedOrder )
{
    const auto names = NsSolverUnsteadyFieldNames();

    ASSERT_EQ( names.flow.size(), 3u );
    EXPECT_EQ( names.flow[ 0 ], "q" );
    EXPECT_EQ( names.flow[ 1 ], "q1" );
    EXPECT_EQ( names.flow[ 2 ], "q2" );

    ASSERT_EQ( names.residual.size(), 3u );
    EXPECT_EQ( names.residual[ 0 ], "res" );
    EXPECT_EQ( names.residual[ 1 ], "res1" );
    EXPECT_EQ( names.residual[ 2 ], "res2" );
}

TEST( SolverFieldNames, TurbUsesPrefixedDualTimeNames )
{
    const auto names = TurbSolverUnsteadyFieldNames();

    ASSERT_EQ( names.flow.size(), 3u );
    EXPECT_EQ( names.flow[ 0 ], "turbq" );
    EXPECT_EQ( names.flow[ 1 ], "turbq1" );
    EXPECT_EQ( names.flow[ 2 ], "turbq2" );

    ASSERT_EQ( names.residual.size(), 3u );
    EXPECT_EQ( names.residual[ 0 ], "turbres" );
    EXPECT_EQ( names.residual[ 1 ], "turbres1" );
    EXPECT_EQ( names.residual[ 2 ], "turbres2" );
}

TEST( SolverFieldNames, TurbLogisticsUsePrefixedNames )
{
    const auto names = TurbSolverFieldLogisticsNames();

    EXPECT_EQ(
        names.solutionGradNames,
        ( StringField{ "turbq", "turbdqdx", "turbdqdy", "turbdqdz" } ) );
    EXPECT_EQ( names.implicitFieldNames, ( StringField{ "turbq", "turbdq" } ) );
    EXPECT_EQ(
        names.bcGradFieldNames,
        ( StringField{ "turbbcdqdx", "turbbcdqdy", "turbbcdqdz" } ) );
    EXPECT_EQ( names.residualName, "turbres" );
    EXPECT_EQ( names.limiterName, "turblimiter" );
    EXPECT_TRUE( names.bcSolutionName.empty() );
}

TEST( SolverFieldNames, INsAndNsHaveEquivalentLogisticsNames )
{
    const auto ns = NsSolverFieldLogisticsNames();
    const auto ins = INsSolverFieldLogisticsNames();

    EXPECT_EQ( ns.solutionGradNames, ins.solutionGradNames );
    EXPECT_EQ( ns.implicitFieldNames, ins.implicitFieldNames );
    EXPECT_EQ( ns.bcGradFieldNames, ins.bcGradFieldNames );
    EXPECT_EQ( ns.residualName, ins.residualName );
    EXPECT_EQ( ns.limiterName, ins.limiterName );
    EXPECT_EQ( ns.bcSolutionName, ins.bcSolutionName );
}

TEST( SolverFieldNames, INsAndNsHaveEquivalentDualTimeNames )
{
    const auto ns = NsSolverUnsteadyFieldNames();
    const auto ins = INsSolverUnsteadyFieldNames();

    EXPECT_EQ( ns.flow, ins.flow );
    EXPECT_EQ( ns.residual, ins.residual );
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
    info.resFileName = "keep_me.dat";

    ApplyFieldNameContract(
        &info,
        NsSolverFieldLogisticsNames(),
        NsSolverUnsteadyFieldNames() );

    EXPECT_EQ( info.residualName, "res" );
    EXPECT_EQ( info.solutionGradNames, ( StringField{ "q", "dqdx", "dqdy", "dqdz" } ) );
    EXPECT_EQ( info.implicitFieldNames, ( StringField{ "q", "dq" } ) );
    EXPECT_EQ( info.bcGradFieldNames, ( StringField{ "bcdqdx", "bcdqdy", "bcdqdz" } ) );
    EXPECT_EQ( info.stageSolutionNames, ( StringField{ "q1", "q2" } ) );
    EXPECT_EQ( info.stageResidualNames, ( StringField{ "res1", "res2" } ) );
    EXPECT_EQ( info.limiterName, "limiter" );
    EXPECT_EQ( info.bcSolutionName, "bc_q" );
    EXPECT_EQ( info.resFileName, "keep_me.dat" );
}
