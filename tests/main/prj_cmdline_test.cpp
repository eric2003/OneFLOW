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

namespace
{

    std::filesystem::path RemoveTrailingSeparator(
        const std::filesystem::path & path )
    {
        std::string value = path.string();

        while ( value.size() > 1 &&
            ( value.back() == '/' || value.back() == '\\' ) )
        {
            value.pop_back();
        }

        return std::filesystem::path( value );
    }

}

TEST( PrjSetPrjBaseDir, RelativePath )
{
    Prj::current_dir = std::filesystem::current_path().string();

    Prj::SetPrjBaseDir( "plate" );

    const std::filesystem::path actual =
        RemoveTrailingSeparator( Prj::prjBaseDir );

    const std::filesystem::path expected =
        std::filesystem::current_path() / "plate";

    EXPECT_EQ(
        actual.lexically_normal(),
        expected.lexically_normal() );

    EXPECT_TRUE(
        !Prj::prjBaseDir.empty() &&
        ( Prj::prjBaseDir.back() == '/' ||
            Prj::prjBaseDir.back() == '\\' ) );
}


TEST( PrjSetPrjBaseDir, NestedRelativePath )
{
    Prj::current_dir = std::filesystem::current_path().string();

    Prj::SetPrjBaseDir( "cases/plate" );

    const std::filesystem::path actual =
        RemoveTrailingSeparator( Prj::prjBaseDir );

    const std::filesystem::path expected =
        std::filesystem::current_path() / "cases" / "plate";

    EXPECT_EQ(
        actual.lexically_normal(),
        expected.lexically_normal() );

    EXPECT_TRUE(
        !Prj::prjBaseDir.empty() &&
        ( Prj::prjBaseDir.back() == '/' ||
            Prj::prjBaseDir.back() == '\\' ) );
}


TEST( PrjSetPrjBaseDir, AbsolutePath )
{
    Prj::current_dir = std::filesystem::current_path().string();

    const std::filesystem::path absoluteCase =
        std::filesystem::temp_directory_path()
        / "OneFLOW"
        / "2026"
        / "plate";

    Prj::SetPrjBaseDir( absoluteCase.string() );

    const std::filesystem::path actual =
        RemoveTrailingSeparator( Prj::prjBaseDir );

    EXPECT_EQ(
        actual.lexically_normal(),
        absoluteCase.lexically_normal() );

    EXPECT_TRUE(
        !Prj::prjBaseDir.empty() &&
        ( Prj::prjBaseDir.back() == '/' ||
            Prj::prjBaseDir.back() == '\\' ) );
}


TEST( PrjSetPrjBaseDir, AbsolutePathIsNotPrefixedByCurrentDirectory )
{
    Prj::current_dir = std::filesystem::current_path().string();

    const std::filesystem::path absoluteCase =
        std::filesystem::temp_directory_path()
        / "OneFLOW"
        / "2026"
        / "plate";

    Prj::SetPrjBaseDir( absoluteCase.string() );

    const std::filesystem::path actual =
        RemoveTrailingSeparator( Prj::prjBaseDir );

    EXPECT_EQ(
        actual.lexically_normal(),
        absoluteCase.lexically_normal() );
}