#include "CliParser.h"
#include "SimulationCommandLine.h"

#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>
TEST( CliParserTest, ParsesRequiredAndRepeatedPositionals )
{
    ONEFLOW::CliParser parser( "OneFLOW" );
    parser.AddPositional( "mode", "Run mode" );
    parser.AddPositional(
        "caseDir", "Case directory", 1, ONEFLOW::CliParser::Unlimited );

    const ONEFLOW::ParsedCommandLine parsed = parser.Parse(
        { "r", "caseA", "caseA", "caseB" } );

    EXPECT_EQ( parsed.GetOne( "mode" ), "r" );
    const std::vector<std::string> & caseDirs = parsed.GetAll( "caseDir" );
    ASSERT_EQ( caseDirs.size(), 3u );
    EXPECT_EQ( caseDirs[ 0 ], "caseA" );
    EXPECT_EQ( caseDirs[ 1 ], "caseA" );
    EXPECT_EQ( caseDirs[ 2 ], "caseB" );
}

TEST( CliParserTest, OptionalPositionalCanBeAbsentOrPresent )
{
    ONEFLOW::CliParser parser( "OneFLOW" );
    parser.AddPositional( "mode", "Run mode", 0, 1 );

    const ONEFLOW::ParsedCommandLine absent = parser.Parse( {} );
    EXPECT_FALSE( absent.Has( "mode" ) );

    const ONEFLOW::ParsedCommandLine present = parser.Parse( { "d" } );
    EXPECT_TRUE( present.Has( "mode" ) );
    EXPECT_EQ( present.GetOne( "mode" ), "d" );
}

TEST( CliParserTest, RejectsMissingRequiredPositional )
{
    ONEFLOW::CliParser parser( "OneFLOW" );
    parser.AddPositional( "mode", "Run mode" );
    parser.AddPositional( "caseDir", "Case directory" );

    EXPECT_THROW( parser.Parse( { "r" } ), std::invalid_argument );
}

TEST( CliParserTest, RejectsUnexpectedExtraPositional )
{
    ONEFLOW::CliParser parser( "OneFLOW" );
    parser.AddPositional( "mode", "Run mode" );

    EXPECT_THROW(
        parser.Parse( { "r", "unexpected" } ), std::invalid_argument );
}

TEST( CliParserTest, RepeatablePositionalMustBeLast )
{
    ONEFLOW::CliParser parser( "OneFLOW" );
    parser.AddPositional(
        "caseDir", "Case directory", 1, ONEFLOW::CliParser::Unlimited );

    EXPECT_THROW(
        parser.AddPositional( "mode", "Run mode" ),
        std::invalid_argument );
}

TEST( CliParserTest, BuildsUsageFromArgumentDefinitions )
{
    ONEFLOW::CliParser parser( "OneFLOW" );
    parser.AddPositional( "mode", "Run mode" );
    parser.AddPositional(
        "caseDir", "Case directory", 1, ONEFLOW::CliParser::Unlimited );

    const std::string usage = parser.Usage();
    EXPECT_NE( usage.find( "Usage: OneFLOW <mode> <caseDir>..." ),
        std::string::npos );
    EXPECT_NE( usage.find( "Case directory" ), std::string::npos );
}

TEST( SimulationInvocationTest, OneArgumentSelectsLightweightTest )
{
    EXPECT_EQ(
        ONEFLOW::ClassifySimulationInvocation( { "OneFlow.exe" } ),
        ONEFLOW::SimulationInvocationKind::LightweightTest );
}

TEST( SimulationInvocationTest, OneAdditionalArgumentIsRejected )
{
    EXPECT_EQ(
        ONEFLOW::ClassifySimulationInvocation( { "OneFlow.exe", "d" } ),
        ONEFLOW::SimulationInvocationKind::InvalidArgumentCount );
}

TEST( SimulationInvocationTest, MultipleCaseArgumentsSelectFullSimulation )
{
    EXPECT_EQ(
        ONEFLOW::ClassifySimulationInvocation(
            { "OneFlow.exe", "r", "caseA", "caseB" } ),
        ONEFLOW::SimulationInvocationKind::FullSimulation );
}

TEST( SimulationCommandLineOptionsTest, ParsesDebugModeAndRepeatedCases )
{
    const ONEFLOW::SimulationCommandLineOptions options =
        ONEFLOW::ParseSimulationCommandLineOptions(
            { "OneFlow.exe", "d", "caseA", "caseA", "caseB" } );

    EXPECT_TRUE( options.debug );
    ASSERT_EQ( options.caseDirectories.size(), 3u );
    EXPECT_EQ( options.caseDirectories[ 0 ], "caseA" );
    EXPECT_EQ( options.caseDirectories[ 1 ], "caseA" );
    EXPECT_EQ( options.caseDirectories[ 2 ], "caseB" );
}

TEST( SimulationCommandLineOptionsTest, RunModeRemainsNonDebug )
{
    const ONEFLOW::SimulationCommandLineOptions options =
        ONEFLOW::ParseSimulationCommandLineOptions(
            { "OneFlow.exe", "r", "caseA" } );

    EXPECT_FALSE( options.debug );
    ASSERT_EQ( options.caseDirectories.size(), 1u );
    EXPECT_EQ( options.caseDirectories.front(), "caseA" );
}

