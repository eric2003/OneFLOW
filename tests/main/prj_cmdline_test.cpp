#include "Prj.h"
#include <gtest/gtest.h>
#include <stdexcept>

using ONEFLOW::Prj;
using ONEFLOW::CmdLineOptions;

TEST( PrjParseCmdLineArgs, NormalRun )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs( { "OneFlow.exe", "r", "test/plateuns2dslau2/" } );
    EXPECT_FALSE( opt.debug );
    EXPECT_EQ( opt.prjName, "test/plateuns2dslau2/" );
}

TEST( PrjParseCmdLineArgs, DebugFlagRecognized )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs( { "OneFlow.exe", "d", "test/plateuns2dslau2/" } );
    EXPECT_TRUE( opt.debug );
    EXPECT_EQ( opt.prjName, "test/plateuns2dslau2/" );
}

TEST( PrjParseCmdLineArgs, NonDFlagIsNotDebug )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs( { "OneFlow.exe", "x", "test/plateuns2dslau2/" } );
    EXPECT_FALSE( opt.debug );
}

TEST( PrjParseCmdLineArgs, TooFewArgsThrows )
{
    EXPECT_THROW( Prj::ParseCmdLineArgs( { "OneFlow.exe", "d" } ), std::invalid_argument );
}

TEST( PrjParseCmdLineArgs, EmptyArgsThrows )
{
    EXPECT_THROW( Prj::ParseCmdLineArgs( {} ), std::invalid_argument );
}

TEST( PrjParseCmdLineArgs, ExtraArgsBeyondThirdAreIgnored )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs(
        { "OneFlow.exe", "d", "test/plateuns2dslau2/", "extra_ignored" } );
    EXPECT_TRUE( opt.debug );
    EXPECT_EQ( opt.prjName, "test/plateuns2dslau2/" );
}