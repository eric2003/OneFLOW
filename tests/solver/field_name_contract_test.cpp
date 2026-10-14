#include <gtest/gtest.h>
#include "FieldNameContract.h"
#include "SolverInfo.h"
#include <string>
#include <vector>

using namespace ONEFLOW;

namespace
{
    std::vector<std::string> ToVector( const StringField & sf )
    {
        std::vector<std::string> out;
        out.reserve( sf.size() );
        for ( std::size_t i = 0; i < sf.size(); ++ i )
        {
            out.push_back( sf[ i ] );
        }
        return out;
    }
}

// ---------------------------------------------------------------------------
// Dual-time order: Current, Previous, Old (flow + residual).
// ---------------------------------------------------------------------------

TEST( FieldNameContract, MakeUnsteady_OrdersCurrentPreviousOld )
{
    StringField stageQ;
    stageQ.push_back( "q1" );
    stageQ.push_back( "q2" );
    StringField stageR;
    stageR.push_back( "res1" );
    stageR.push_back( "res2" );

    const UnsteadyFieldNames names =
        MakeUnsteadyFieldNames( "q", stageQ, "res", stageR );

    ASSERT_EQ( names.flow.size(), 3u );
    EXPECT_EQ( names.flow[ 0 ], "q" );
    EXPECT_EQ( names.flow[ 1 ], "q1" );
    EXPECT_EQ( names.flow[ 2 ], "q2" );

    ASSERT_EQ( names.residual.size(), 3u );
    EXPECT_EQ( names.residual[ 0 ], "res" );
    EXPECT_EQ( names.residual[ 1 ], "res1" );
    EXPECT_EQ( names.residual[ 2 ], "res2" );
}

TEST( FieldNameContract, NsMatchesSystemAllocUnsteadyTxt )
{
    const auto names = ExpectedNsUnsteadyFieldNames();
    EXPECT_EQ( ToVector( names.flow ), ( std::vector<std::string>{ "q", "q1", "q2" } ) );
    EXPECT_EQ( ToVector( names.residual ), ( std::vector<std::string>{ "res", "res1", "res2" } ) );
}

TEST( FieldNameContract, INsSharesNsDualTimeNames )
{
    const auto ns = ExpectedNsUnsteadyFieldNames();
    const auto ins = ExpectedINsUnsteadyFieldNames();
    EXPECT_EQ( ToVector( ns.flow ), ToVector( ins.flow ) );
    EXPECT_EQ( ToVector( ns.residual ), ToVector( ins.residual ) );
}

TEST( FieldNameContract, TurbMatchesSystemAllocUnsteadyTxt )
{
    const auto names = ExpectedTurbUnsteadyFieldNames();
    EXPECT_EQ( ToVector( names.flow ), ( std::vector<std::string>{ "turbq", "turbq1", "turbq2" } ) );
    EXPECT_EQ( ToVector( names.residual ), ( std::vector<std::string>{ "turbres", "turbres1", "turbres2" } ) );
}

TEST( FieldNameContract, TurbUsesPrefixStrategyNotBareQ )
{
    const auto turb = ExpectedTurbUnsteadyFieldNames();
    for ( std::size_t i = 0; i < turb.flow.size(); ++ i )
    {
        EXPECT_NE( turb.flow[ i ].find( "turb" ), std::string::npos );
    }
    for ( std::size_t i = 0; i < turb.residual.size(); ++ i )
    {
        EXPECT_NE( turb.residual[ i ].find( "turb" ), std::string::npos );
    }
}

// ---------------------------------------------------------------------------
// Logistics name table (solutionGrad / implicit / bcGrad / residual / limiter).
// Mirrors NsSolver / INsSolver / TurbSolver StaticInit.
// ---------------------------------------------------------------------------

TEST( FieldNameContract, NsLogistics_SolutionGradOrder )
{
    const auto n = ExpectedNsSolverLogisticsNames();
    ASSERT_EQ( n.solutionGradNames.size(), 4u );
    EXPECT_EQ( n.solutionGradNames[ 0 ], "q" );
    EXPECT_EQ( n.solutionGradNames[ 1 ], "dqdx" );
    EXPECT_EQ( n.solutionGradNames[ 2 ], "dqdy" );
    EXPECT_EQ( n.solutionGradNames[ 3 ], "dqdz" );
}

TEST( FieldNameContract, NsLogistics_ImplicitIsSolutionAndDelta )
{
    const auto n = ExpectedNsSolverLogisticsNames();
    ASSERT_EQ( n.implicitFieldNames.size(), 2u );
    EXPECT_EQ( n.implicitFieldNames[ 0 ], "q" );
    EXPECT_EQ( n.implicitFieldNames[ 1 ], "dq" );
    EXPECT_EQ( n.residualName, "res" );
    EXPECT_EQ( n.limiterName, "limiter" );
    EXPECT_EQ( n.bcSolutionName, "bc_q" );
}

TEST( FieldNameContract, NsLogistics_BcGradNames )
{
    const auto n = ExpectedNsSolverLogisticsNames();
    ASSERT_EQ( n.bcGradFieldNames.size(), 3u );
    EXPECT_EQ( n.bcGradFieldNames[ 0 ], "bcdqdx" );
    EXPECT_EQ( n.bcGradFieldNames[ 1 ], "bcdqdy" );
    EXPECT_EQ( n.bcGradFieldNames[ 2 ], "bcdqdz" );
}

TEST( FieldNameContract, TurbLogistics_PrefixedAndEmptyBcSolution )
{
    const auto n = ExpectedTurbSolverLogisticsNames();
    EXPECT_EQ( n.solutionGradNames[ 0 ], "turbq" );
    EXPECT_EQ( n.solutionGradNames[ 1 ], "turbdqdx" );
    EXPECT_EQ( n.implicitFieldNames[ 1 ], "turbdq" );
    EXPECT_EQ( n.residualName, "turbres" );
    EXPECT_EQ( n.limiterName, "turblimiter" );
    EXPECT_TRUE( n.bcSolutionName.empty() );
    EXPECT_EQ( n.bcGradFieldNames[ 0 ], "turbbcdqdx" );
}

TEST( FieldNameContract, INsLogisticsMatchesNs )
{
    const auto ns = ExpectedNsSolverLogisticsNames();
    const auto ins = ExpectedINsSolverLogisticsNames();
    EXPECT_EQ( ToVector( ns.solutionGradNames ), ToVector( ins.solutionGradNames ) );
    EXPECT_EQ( ToVector( ns.implicitFieldNames ), ToVector( ins.implicitFieldNames ) );
    EXPECT_EQ( ToVector( ns.bcGradFieldNames ), ToVector( ins.bcGradFieldNames ) );
    EXPECT_EQ( ns.residualName, ins.residualName );
    EXPECT_EQ( ns.limiterName, ins.limiterName );
    EXPECT_EQ( ns.bcSolutionName, ins.bcSolutionName );
}

TEST( FieldNameContract, TurbAndNsDoNotShareBareSolutionName )
{
    const auto ns = ExpectedNsSolverLogisticsNames();
    const auto turb = ExpectedTurbSolverLogisticsNames();
    EXPECT_NE( ns.solutionGradNames[ 0 ], turb.solutionGradNames[ 0 ] );
    EXPECT_NE( ns.residualName, turb.residualName );
}

TEST( FieldNameContract, Apply_NsFillsSolverInfoSlots )
{
    SolverInfo info;
    ApplyFieldNameContract(
        &info,
        ExpectedNsSolverLogisticsNames(),
        ExpectedNsUnsteadyFieldNames() );

    EXPECT_EQ( info.residualName, "res" );
    ASSERT_EQ( info.solutionGradNames.size(), 4u );
    EXPECT_EQ( info.solutionGradNames[ 0 ], "q" );
    EXPECT_EQ( info.implicitFieldNames[ 1 ], "dq" );
    ASSERT_EQ( info.stageSolutionNames.size(), 2u );
    EXPECT_EQ( info.stageSolutionNames[ 0 ], "q1" );
    EXPECT_EQ( info.stageResidualNames[ 0 ], "res1" );
    EXPECT_EQ( info.limiterName, "limiter" );
    EXPECT_EQ( info.bcSolutionName, "bc_q" );
}

TEST( FieldNameContract, Apply_TurbUsesPrefixedNames )
{
    SolverInfo info;
    ApplyFieldNameContract(
        &info,
        ExpectedTurbSolverLogisticsNames(),
        ExpectedTurbUnsteadyFieldNames() );

    EXPECT_EQ( info.residualName, "turbres" );
    EXPECT_EQ( info.solutionGradNames[ 0 ], "turbq" );
    EXPECT_EQ( info.stageSolutionNames[ 0 ], "turbq1" );
    EXPECT_EQ( info.stageResidualNames[ 1 ], "turbres2" );
    EXPECT_TRUE( info.bcSolutionName.empty() );
}

TEST( FieldNameContract, Apply_DoesNotTouchResFileName )
{
    SolverInfo info;
    info.resFileName = "keep_me.dat";
    ApplyFieldNameContract(
        &info,
        ExpectedNsSolverLogisticsNames(),
        ExpectedNsUnsteadyFieldNames() );
    EXPECT_EQ( info.resFileName, "keep_me.dat" );
}
