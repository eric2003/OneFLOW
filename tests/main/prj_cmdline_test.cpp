#include "Prj.h"
#include "FileUtils.h"
#include <gtest/gtest.h>
#include <stdexcept>

using ONEFLOW::Prj;
using ONEFLOW::CmdLineOptions;

TEST( PrjParseCmdLineArgs, NormalRun )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs(
        { "OneFlow.exe", "r", "test/plateuns2dslau2/" } );

    EXPECT_FALSE( opt.debug );
    EXPECT_EQ( opt.caseDir, "test/plateuns2dslau2/" );
}

TEST( PrjParseCmdLineArgs, DebugFlagRecognized )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs(
        { "OneFlow.exe", "d", "test/plateuns2dslau2/" } );

    EXPECT_TRUE( opt.debug );
    EXPECT_EQ( opt.caseDir, "test/plateuns2dslau2/" );
}

TEST( PrjParseCmdLineArgs, NonDFlagIsNotDebug )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs(
        { "OneFlow.exe", "x", "test/plateuns2dslau2/" } );

    EXPECT_FALSE( opt.debug );
}

TEST( PrjParseCmdLineArgs, TooFewArgsThrows )
{
    EXPECT_THROW(
        Prj::ParseCmdLineArgs( { "OneFlow.exe", "d" } ),
        std::invalid_argument );
}

TEST( PrjParseCmdLineArgs, EmptyArgsThrows )
{
    EXPECT_THROW(
        Prj::ParseCmdLineArgs( {} ),
        std::invalid_argument );
}
TEST( PrjParseCmdLineArgs, ExtraArgsBeyondThirdAreIgnored )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs(
        {
            "OneFlow.exe",
            "d",
            "test/plateuns2dslau2/",
            "extra_ignored"
        } );

    EXPECT_TRUE( opt.debug );
    EXPECT_EQ( opt.caseDir, "test/plateuns2dslau2/" );
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


TEST( PrjCasePath, LeadingSlashIsCaseRelative )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjCasePathTest"
        / "case";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    const std::filesystem::path relative =
        Prj::GetPrjFileName( "grid/test.dat" );

    const std::filesystem::path leadingSlash =
        Prj::GetPrjFileName( "/grid/test.dat" );

    EXPECT_EQ(
        relative.lexically_normal(),
        leadingSlash.lexically_normal() );
}

TEST( PrjCasePath, OpenPrjFileUsesCaseRelativePath )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjOpenFileTest"
        / "case";

    std::filesystem::remove_all( caseDir );

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    std::fstream file;

    Prj::OpenPrjFile(
        file,
        "/grid/test.dat",
        std::ios_base::out );

    ASSERT_TRUE( file.is_open() );

    file << "OneFLOW";
    Prj::CloseFile( file );

    const std::filesystem::path expected =
        caseDir / "grid" / "test.dat";

    EXPECT_TRUE( std::filesystem::is_regular_file( expected ) );

    std::filesystem::remove_all( caseDir );
}

TEST( PrjCasePath, MakePrjDirUsesCaseRelativePath )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjMakeDirTest"
        / "case";

    std::filesystem::remove_all( caseDir );

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    Prj::MakePrjDir( "/output/data" );

    const std::filesystem::path expected =
        caseDir / "output" / "data";

    EXPECT_TRUE( std::filesystem::is_directory( expected ) );

    std::filesystem::remove_all( caseDir );
}

TEST( PrjSystemPath, GetSystemFileNameJoinsWithSystemRoot )
{
    const std::string savedRoot = Prj::system_root;

    const std::filesystem::path systemRoot =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjSystemPathTest"
        / "system";

    Prj::system_root = systemRoot.string() + "/";

    const std::filesystem::path expected = systemRoot / "action" / "actionFileList.txt";

    EXPECT_EQ(
        std::filesystem::path( Prj::GetSystemFileName( "action/actionFileList.txt" ) )
        .lexically_normal(),
        expected.lexically_normal() );

    // A leading slash should behave the same as a system-root-relative path,
    // matching the equivalent leading-slash handling in GetPrjFileName.
    EXPECT_EQ(
        std::filesystem::path( Prj::GetSystemFileName( "/action/actionFileList.txt" ) )
        .lexically_normal(),
        expected.lexically_normal() );

    Prj::system_root = savedRoot;
}

TEST( PrjPath, GetDirName )
{
    EXPECT_EQ( Prj::GetDirName( "grid/test.dat" ), "grid" );
    EXPECT_EQ( Prj::GetDirName( "/grid/test.dat" ), "/grid" );
    EXPECT_EQ( Prj::GetDirName( "test.dat" ), "" );
    EXPECT_EQ( Prj::GetDirName( "grid\\test.dat" ), "grid" );
}

TEST( PrjPathUtils, SlashDetection )
{
    EXPECT_FALSE( ONEFLOW::EndWithSlash( "" ) );

    EXPECT_TRUE( ONEFLOW::EndWithSlash( "/" ) );
    EXPECT_TRUE( ONEFLOW::EndWithSlash( "\\" ) );

    EXPECT_TRUE( ONEFLOW::EndWithForwardSlash( "grid/" ) );
    EXPECT_TRUE( ONEFLOW::EndWithBackwardSlash( "grid\\" ) );
    EXPECT_TRUE( ONEFLOW::StartWithForwardSlash( "/grid" ) );
    EXPECT_FALSE( ONEFLOW::StartWithForwardSlash( "\\grid" ) );
}

TEST( PrjPathUtils, SlashRemoval )
{
    EXPECT_EQ( ONEFLOW::RemoveFirstSlash( "/grid/test.dat" ),
        "grid/test.dat" );

    EXPECT_EQ( ONEFLOW::RemoveFirstSlash( "\\grid\\test.dat" ),
        "grid\\test.dat" );

    EXPECT_EQ( ONEFLOW::RemoveFirstSlash( "grid/test.dat" ),
        "grid/test.dat" );

    EXPECT_EQ( ONEFLOW::RemoveFirstSlash( "" ), "" );

    EXPECT_EQ( ONEFLOW::RemoveEndSlash( "grid/" ), "grid" );
    EXPECT_EQ( ONEFLOW::RemoveEndSlash( "grid\\" ), "grid" );
    EXPECT_EQ( ONEFLOW::RemoveEndSlash( "grid" ), "grid" );
    EXPECT_EQ( ONEFLOW::RemoveEndSlash( "" ), "" );
}