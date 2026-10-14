#include <gtest/gtest.h>
#include "FieldNameContract.h"
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
// Dual-time order contract: Current, Previous, Old for both flow and residual.
// Mirrors system/ns|ins|turb/alloc/unsteady.txt and SolverInfo StaticInit.
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
    const auto flow = ToVector( names.flow );
    const auto res = ToVector( names.residual );

    EXPECT_EQ( flow, ( std::vector<std::string>{ "q", "q1", "q2" } ) );
    EXPECT_EQ( res, ( std::vector<std::string>{ "res", "res1", "res2" } ) );
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
    const auto flow = ToVector( names.flow );
    const auto res = ToVector( names.residual );

    EXPECT_EQ( flow, ( std::vector<std::string>{ "turbq", "turbq1", "turbq2" } ) );
    EXPECT_EQ( res, ( std::vector<std::string>{ "turbres", "turbres1", "turbres2" } ) );
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
