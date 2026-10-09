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
TEST( PrjParseCmdLineArgs, AdditionalCaseDirectoriesArePreserved )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs(
        {
            "OneFlow.exe",
            "d",
            "test/caseA/",
            "test/caseB/"
        } );

    EXPECT_TRUE( opt.debug );
    EXPECT_EQ( opt.caseDir, "test/caseA/" );
    ASSERT_EQ( opt.caseDirs.size(), 2u );
    EXPECT_EQ( opt.caseDirs[ 0 ], "test/caseA/" );
    EXPECT_EQ( opt.caseDirs[ 1 ], "test/caseB/" );

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


TEST( PrjCasePath, RootedPathIsNotCaseRelative )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjCasePathTest"
        / "case";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    const std::filesystem::path rootedPath =
        std::filesystem::path( "/grid/test.dat" );

    const std::filesystem::path actual =
        Prj::GetPrjFileName( rootedPath.string() );

    // A rooted path is not case-relative. On POSIX it is absolute; on
    // Windows it is root-relative and keeps the current drive.
    const std::filesystem::path expected =
        ( caseDir / rootedPath ).lexically_normal();

    EXPECT_EQ(
        actual.lexically_normal(),
        expected );
}

TEST( PrjCasePath, AbsolutePathIsNotCaseRelative )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjAbsolutePathTest"
        / "case";

    const std::filesystem::path absolutePath =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjAbsolutePathTest"
        / "external"
        / "grid.dat";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    const std::filesystem::path actual =
        Prj::GetPrjFileName( absolutePath.string() );

    EXPECT_EQ(
        actual.lexically_normal(),
        absolutePath.lexically_normal() );
}

TEST( PrjCasePath, RelativeCasePathIsResolvedFromCurrentCase )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjRelativeCasePathTest"
        / "consumer";

    const std::filesystem::path sourceCaseDir =
        caseDir.parent_path() / "producer";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    const std::filesystem::path actual =
        Prj::GetCaseFileName(
            "../producer",
            "grid/test.dat" );

    const std::filesystem::path expected =
        sourceCaseDir / "grid" / "test.dat";

    EXPECT_EQ(
        std::filesystem::path( actual ).lexically_normal(),
        expected.lexically_normal() );
}

TEST( PrjCasePath, OpenPrjFileUsesRelativePath )
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
        "grid/test.dat",
        std::ios_base::out );

    ASSERT_TRUE( file.is_open() );

    file << "OneFLOW";
    Prj::CloseFile( file );

    const std::filesystem::path expected =
        caseDir / "grid" / "test.dat";

    EXPECT_TRUE( std::filesystem::is_regular_file( expected ) );

    std::filesystem::remove_all( caseDir );
}

TEST( PrjCasePath, OpenPrjFileInCaseRoot )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjOpenRootFileTest"
        / "case";

    std::filesystem::remove_all( caseDir );
    std::filesystem::create_directories( caseDir );

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    std::fstream file;

    Prj::OpenPrjFile(
        file,
        "test.dat",
        std::ios_base::out );

    ASSERT_TRUE( file.is_open() );

    file << "OneFLOW";
    Prj::CloseFile( file );

    EXPECT_TRUE( std::filesystem::is_regular_file( caseDir / "test.dat" ) );

    std::filesystem::remove_all( caseDir );
}

TEST( PrjCasePath, OpenCaseFileInCaseRoot )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjOpenCaseRootFileTest"
        / "case";

    std::filesystem::remove_all( caseDir );
    std::filesystem::create_directories( caseDir );

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    std::fstream file;

    Prj::OpenCaseFile(
        file,
        caseDir.string(),
        "test.dat",
        std::ios_base::out );

    ASSERT_TRUE( file.is_open() );

    file << "OneFLOW";
    Prj::CloseFile( file );

    EXPECT_TRUE( std::filesystem::is_regular_file( caseDir / "test.dat" ) );

    std::filesystem::remove_all( caseDir );
}

TEST( PrjCasePath, OpenPrjFileUsesAbsolutePath )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjOpenAbsoluteFileTest"
        / "case";

    const std::filesystem::path absoluteFile =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjOpenAbsoluteFileTest"
        / "external"
        / "test.dat";

    std::filesystem::remove_all( caseDir );
    std::filesystem::remove_all( absoluteFile.parent_path() );

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    std::fstream file;

    Prj::OpenPrjFile(
        file,
        absoluteFile.string(),
        std::ios_base::out );

    ASSERT_TRUE( file.is_open() );

    file << "OneFLOW";
    Prj::CloseFile( file );

    EXPECT_TRUE( std::filesystem::is_regular_file( absoluteFile ) );
    EXPECT_FALSE(
        std::filesystem::is_regular_file(
            caseDir / "external" / "test.dat" ) );

    std::filesystem::remove_all( caseDir );
    std::filesystem::remove_all( absoluteFile.parent_path() );
}

TEST( PrjCasePath, MakePrjDirFailsWhenPathCannotBeCreated )
{
    const std::filesystem::path testRoot =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjMakeDirFailureTest";
    const std::filesystem::path blockingFile = testRoot / "not_a_directory";

    std::filesystem::remove_all( testRoot );
    std::filesystem::create_directories( testRoot );
    {
        std::ofstream file( blockingFile );
        ASSERT_TRUE( file.is_open() );
        file << "block";
    }

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( testRoot.string() );

    EXPECT_THROW(
        Prj::MakePrjDir( "not_a_directory/child" ),
        std::runtime_error );

    std::filesystem::remove_all( testRoot );
}

TEST( PrjCasePath, MakePrjDirUsesRelativePath )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjMakeDirTest"
        / "case";

    std::filesystem::remove_all( caseDir );

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    Prj::MakePrjDir( "output/data" );

    const std::filesystem::path expected =
        caseDir / "output" / "data";

    EXPECT_TRUE( std::filesystem::is_directory( expected ) );

    std::filesystem::remove_all( caseDir );
}

TEST( PrjCasePath, MakePrjDirUsesAbsolutePath )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjMakeAbsoluteDirTest"
        / "case";

    const std::filesystem::path absoluteDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjMakeAbsoluteDirTest"
        / "external"
        / "output"
        / "data";

    std::filesystem::remove_all( caseDir );
    std::filesystem::remove_all( absoluteDir.parent_path().parent_path() );

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    Prj::MakePrjDir( absoluteDir.string() );

    EXPECT_TRUE( std::filesystem::is_directory( absoluteDir ) );
    EXPECT_FALSE(
        std::filesystem::is_directory(
            caseDir / "external" / "output" / "data" ) );

    std::filesystem::remove_all( caseDir );
    std::filesystem::remove_all( absoluteDir.parent_path().parent_path() );
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

    const std::filesystem::path absoluteFile =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjSystemPathTest"
        / "external"
        / "actionFileList.txt";

    EXPECT_EQ(
        std::filesystem::path(
            Prj::GetSystemFileName( absoluteFile.string() ) )
        .lexically_normal(),
        absoluteFile.lexically_normal() );

    const std::filesystem::path rootedPath =
        std::filesystem::path( "/action/actionFileList.txt" );

    // A rooted path is not system-root-relative. On POSIX it is absolute;
    // on Windows it is root-relative and keeps the current drive.
    const std::filesystem::path rootedExpected =
        ( systemRoot / rootedPath ).lexically_normal();

    EXPECT_EQ(
        std::filesystem::path(
            Prj::GetSystemFileName( rootedPath.string() ) )
        .lexically_normal(),
        rootedExpected );

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

TEST( PrjSetPrjBaseDir, RelativePathRequiresInitializedCurrentDirectory )
{
    const std::string savedCurrentDir = Prj::current_dir;
    Prj::current_dir.clear();

    EXPECT_THROW(
        Prj::ResolveCaseDir( "plate" ),
        std::runtime_error );

    Prj::current_dir = savedCurrentDir;
}
