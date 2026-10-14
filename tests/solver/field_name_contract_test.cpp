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
// Generic machinery only (no solver module knowledge).
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

TEST( FieldNameContract, Apply_FillsSlotsFromCallerTables )
{
    // Caller supplies tables; FieldNameContract does not know "NS" vs "Turb".
    const SolverLogisticsNames logistics = MakeSolverLogisticsNames(
        "q", "dqdx", "dqdy", "dqdz", "dq", "res",
        "bcdqdx", "bcdqdy", "bcdqdz",
        "limiter", "bc_q" );
    StringField stageQ;
    stageQ.push_back( "q1" );
    stageQ.push_back( "q2" );
    StringField stageR;
    stageR.push_back( "res1" );
    stageR.push_back( "res2" );
    const UnsteadyFieldNames unsteady =
        MakeUnsteadyFieldNames( "q", stageQ, "res", stageR );

    SolverInfo info;
    info.resFileName = "keep_me.dat";
    ApplyFieldNameContract( &info, logistics, unsteady );

    EXPECT_EQ( info.residualName, "res" );
    EXPECT_EQ( info.solutionGradNames[ 0 ], "q" );
    EXPECT_EQ( info.implicitFieldNames[ 1 ], "dq" );
    EXPECT_EQ( info.stageSolutionNames[ 0 ], "q1" );
    EXPECT_EQ( info.stageResidualNames[ 1 ], "res2" );
    EXPECT_EQ( info.resFileName, "keep_me.dat" );
}

TEST( FieldNameContract, Apply_TurbStyleTablesStayPrefixed )
{
    const SolverLogisticsNames logistics = MakeSolverLogisticsNames(
        "turbq", "turbdqdx", "turbdqdy", "turbdqdz", "turbdq", "turbres",
        "turbbcdqdx", "turbbcdqdy", "turbbcdqdz",
        "turblimiter", "" );
    StringField stageQ;
    stageQ.push_back( "turbq1" );
    stageQ.push_back( "turbq2" );
    StringField stageR;
    stageR.push_back( "turbres1" );
    stageR.push_back( "turbres2" );
    const UnsteadyFieldNames unsteady =
        MakeUnsteadyFieldNames( "turbq", stageQ, "turbres", stageR );

    SolverInfo info;
    ApplyFieldNameContract( &info, logistics, unsteady );

    EXPECT_EQ( info.residualName, "turbres" );
    EXPECT_NE( info.solutionGradNames[ 0 ].find( "turb" ), std::string::npos );
    EXPECT_TRUE( info.bcSolutionName.empty() );
}
