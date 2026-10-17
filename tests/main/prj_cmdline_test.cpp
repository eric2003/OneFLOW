#include "Prj.h"
#include "FileUtils.h"
#include <gtest/gtest.h>
#include <stdexcept>

using ONEFLOW::Prj;
using ONEFLOW::CmdLineOptions;

TEST( PrjParseCmdLineArgs, NormalRun )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs(
        { "OneFlow.exe", "r", "test/cfd/plateuns2dslau2/" } );

    EXPECT_FALSE( opt.debug );
    EXPECT_EQ( opt.caseDir, "test/cfd/plateuns2dslau2/" );
}

TEST( PrjParseCmdLineArgs, DebugFlagRecognized )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs(
        { "OneFlow.exe", "d", "test/cfd/plateuns2dslau2/" } );

    EXPECT_TRUE( opt.debug );
    EXPECT_EQ( opt.caseDir, "test/cfd/plateuns2dslau2/" );
}

TEST( PrjParseCmdLineArgs, NonDFlagIsNotDebug )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs(
        { "OneFlow.exe", "x", "test/cfd/plateuns2dslau2/" } );

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

TEST( FileUtilsCreateDirectory, ExistingDirectoryIsSuccess )
{
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path()
        / "OneFLOW_HXCreateDirectoryExistingTest";

    std::filesystem::remove_all( directory );
    ASSERT_TRUE( std::filesystem::create_directories( directory ) );

    EXPECT_TRUE( ONEFLOW::HX_CreateDirectory( directory.string() ));

    std::filesystem::remove_all( directory );
}

TEST( FileUtilsCreateDirectory, ExistingFileIsFailure )
{
    const std::filesystem::path root =
        std::filesystem::temp_directory_path()
        / "OneFLOW_HXCreateDirectoryFileTest";
    const std::filesystem::path filePath = root / "not_a_directory";

    std::filesystem::remove_all( root );
    std::filesystem::create_directories( root );
    {
        std::ofstream file( filePath );
        ASSERT_TRUE( file.is_open() );
        file << "block";
    }

    EXPECT_FALSE( ONEFLOW::HX_CreateDirectory( filePath.string() ));

    std::filesystem::remove_all( root );
}



TEST( PrjIOPaths, InputAndOutputDefaultToCurrentCase )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIOPathsDefaultTest"
        / "consumer";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );
    Prj::SetIOBaseDirs( "", "" );

    EXPECT_EQ(
        std::filesystem::path( Prj::GetInputFileName( "grid/input.ofl" ) ).lexically_normal(),
        ( caseDir / "grid" / "input.ofl" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetOutputFileName( "results/flow.dat" ) ).lexically_normal(),
        ( caseDir / "results" / "flow.dat" ).lexically_normal() );
}

TEST( PrjIOPaths, InputAndOutputCanUseDifferentRoots )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIOPathsDifferentRootsTest"
        / "consumer";
    const std::filesystem::path producerDir = caseDir.parent_path() / "producer";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );
    Prj::SetIOBaseDirs( "../producer", "case" );

    EXPECT_EQ(
        std::filesystem::path( Prj::GetInputFileName( "results/grid.ofl" ) ).lexically_normal(),
        ( producerDir / "results" / "grid.ofl" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetOutputFileName( "results/flow.dat" ) ).lexically_normal(),
        ( caseDir / "results" / "flow.dat" ).lexically_normal() );
}

TEST( PrjIOPaths, AbsoluteFilesBypassConfiguredRoots )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIOPathsAbsoluteTest"
        / "consumer";
    const std::filesystem::path externalFile =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIOPathsAbsoluteTest"
        / "shared"
        / "grid.ofl";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );
    Prj::SetIOBaseDirs( "case", "case" );

    EXPECT_EQ(
        std::filesystem::path( Prj::GetInputFileName( externalFile.string() ) ).lexically_normal(),
        externalFile.lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetOutputFileName( externalFile.string() ) ).lexically_normal(),
        externalFile.lexically_normal() );
}

TEST( PrjIOPaths, ProjectRootMustBeExplicitlyConfigured )
{
    const std::filesystem::path projectDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIOPathsProjectTest"
        / "project";
    const std::filesystem::path caseDir = projectDir / "cases" / "consumer";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );
    Prj::SetProjectBaseDir( projectDir.string() );
    Prj::SetIOBaseDirs( "project", "project" );

    EXPECT_EQ(
        std::filesystem::path( Prj::GetInputFileName( "shared/grid.ofl" ) ).lexically_normal(),
        ( projectDir / "shared" / "grid.ofl" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetOutputFileName( "results/flow.dat" ) ).lexically_normal(),
        ( projectDir / "results" / "flow.dat" ).lexically_normal() );
}


TEST( PrjIOPaths, OpenPrjFileUsesModeSpecificRoots )
{
    const std::filesystem::path root =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIOPathsOpenModeTest";
    const std::filesystem::path caseDir = root / "consumer";
    const std::filesystem::path producerDir = root / "producer";
    std::filesystem::remove_all( root );
    std::filesystem::create_directories( producerDir / "results" );

    {
        std::ofstream source( producerDir / "results" / "grid.ofl" );
        source << "source-grid";
    }

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );
    Prj::SetIOBaseDirs( "../producer", "case" );

    std::fstream input;
    Prj::OpenPrjFile( input, "results/grid.ofl", std::ios_base::in );
    std::string value;
    input >> value;
    Prj::CloseFile( input );
    EXPECT_EQ( value, "source-grid" );

    std::fstream output;
    Prj::OpenPrjFile( output, "results/flow.dat", std::ios_base::out );
    output << "flow-result";
    Prj::CloseFile( output );

    EXPECT_TRUE( std::filesystem::is_regular_file( caseDir / "results" / "flow.dat" ) );
    EXPECT_FALSE( std::filesystem::exists( producerDir / "results" / "flow.dat" ) );

    std::filesystem::remove_all( root );
}


TEST( PrjIOPaths, DomainRootsInheritGenericRoots )
{
    const std::filesystem::path root =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIODomainInheritanceTest";
    const std::filesystem::path caseDir = root / "consumer";
    const std::filesystem::path genericInput = root / "shared-input";
    const std::filesystem::path genericOutput = root / "shared-output";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );
    Prj::SetIOBaseDirs( genericInput.string(), genericOutput.string() );
    Prj::SetDomainIOBaseDirs( "", "", "", "" );

    EXPECT_EQ(
        std::filesystem::path( Prj::GetMeshInputFileName( "grid.ofl" ) ).lexically_normal(),
        ( genericInput / "grid.ofl" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetFlowInputFileName( "restart.dat" ) ).lexically_normal(),
        ( genericInput / "restart.dat" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetMeshOutputFileName( "grid.ofl" ) ).lexically_normal(),
        ( genericOutput / "grid.ofl" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetFlowOutputFileName( "flow.dat" ) ).lexically_normal(),
        ( genericOutput / "flow.dat" ).lexically_normal() );
}

TEST( PrjIOPaths, DomainRootsCanBeConfiguredIndependently )
{
    const std::filesystem::path root =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIODomainRootsTest";
    const std::filesystem::path caseDir = root / "cases" / "consumer";
    const std::filesystem::path projectDir = root / "project";
    const std::filesystem::path meshInput = root / "mesh-source";
    const std::filesystem::path flowInput = root / "flow-source";
    const std::filesystem::path meshOutput = root / "mesh-output";
    const std::filesystem::path flowOutput = root / "flow-output";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );
    Prj::SetProjectBaseDir( projectDir.string() );
    Prj::SetIOBaseDirs( "case", "case" );
    Prj::SetDomainIOBaseDirs(
        meshInput.string(), flowInput.string(),
        meshOutput.string(), flowOutput.string() );

    EXPECT_EQ(
        std::filesystem::path( Prj::GetMeshInputFileName( "grid.ofl" ) ).lexically_normal(),
        ( meshInput / "grid.ofl" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetFlowInputFileName( "restart.dat" ) ).lexically_normal(),
        ( flowInput / "restart.dat" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetMeshOutputFileName( "grid.ofl" ) ).lexically_normal(),
        ( meshOutput / "grid.ofl" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetFlowOutputFileName( "flow.dat" ) ).lexically_normal(),
        ( flowOutput / "flow.dat" ).lexically_normal() );
}

TEST( PrjIOPaths, DiagnosticsReportEffectiveDomainRoots )
{
    const std::filesystem::path root =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIODiagnosticsTest";
    const std::filesystem::path caseDir = root / "case";
    const std::filesystem::path meshInput = root / "mesh-input";
    const std::filesystem::path flowInput = root / "flow-input";
    const std::filesystem::path meshOutput = root / "mesh-output";
    const std::filesystem::path flowOutput = root / "flow-output";
    const std::filesystem::path report = root / "reports" / "io_paths.txt";

    std::filesystem::remove_all( root );
    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );
    Prj::SetIOBaseDirs( "case", "case" );
    Prj::SetDomainIOBaseDirs(
        meshInput.string(), flowInput.string(),
        meshOutput.string(), flowOutput.string() );

    Prj::WriteIOPathInfo( report.string() );

    std::ifstream file( report );
    ASSERT_TRUE( file.is_open() );
    std::ostringstream contents;
    contents << file.rdbuf();
    const std::string text = contents.str();
    file.close();

    const auto formattedDirectory = []( std::filesystem::path path )
    {
        path = path.lexically_normal();
        path.make_preferred();
        std::string value = path.string();
        std::filesystem::path rootPath = path.root_path();
        rootPath.make_preferred();
        const std::string rootName = rootPath.string();
        while ( value.size() > rootName.size()
            && ! value.empty()
            && ( value.back() == '/' || value.back() == '\\' ) )
        {
            value.pop_back();
        }
        return value;
    };

    EXPECT_NE(
        text.find( "case_root=" + formattedDirectory( caseDir ) ),
        std::string::npos );
    EXPECT_NE(
        text.find( "mesh_input_root=" + formattedDirectory( meshInput ) ),
        std::string::npos );
    EXPECT_NE(
        text.find( "flow_input_root=" + formattedDirectory( flowInput ) ),
        std::string::npos );
    EXPECT_NE(
        text.find( "mesh_output_root=" + formattedDirectory( meshOutput ) ),
        std::string::npos );
    EXPECT_NE(
        text.find( "flow_output_root=" + formattedDirectory( flowOutput ) ),
        std::string::npos );

    EXPECT_EQ( text.find( "mesh_input_root=" + formattedDirectory( meshInput ) + "\\" ),
        std::string::npos );

    std::filesystem::remove_all( root );
}

TEST( PrjIOPaths, DomainRootsResetWhenBindingNextCase )
{
    const std::filesystem::path root =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIODomainCaseResetTest";
    const std::filesystem::path firstCase = root / "first";
    const std::filesystem::path secondCase = root / "second";
    const std::filesystem::path externalMesh = root / "mesh";
    const std::filesystem::path externalFlow = root / "flow";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( firstCase.string() );
    Prj::SetProjectBaseDir( root.string() );
    Prj::SetIOBaseDirs( "project", "project" );
    Prj::SetDomainIOBaseDirs(
        externalMesh.string(), externalFlow.string(),
        externalMesh.string(), externalFlow.string() );

    Prj::SetPrjBaseDir( secondCase.string() );

    EXPECT_TRUE( Prj::projectBaseDir.empty() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetMeshInputFileName( "grid.ofl" ) ).lexically_normal(),
        ( secondCase / "grid.ofl" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetFlowInputFileName( "restart.dat" ) ).lexically_normal(),
        ( secondCase / "restart.dat" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetMeshOutputFileName( "grid.ofl" ) ).lexically_normal(),
        ( secondCase / "grid.ofl" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetFlowOutputFileName( "flow.dat" ) ).lexically_normal(),
        ( secondCase / "flow.dat" ).lexically_normal() );
}

TEST( PrjIOPaths, DomainRootCanInheritProjectRoot )
{
    const std::filesystem::path root =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIODomainProjectTest";
    const std::filesystem::path caseDir = root / "cases" / "consumer";
    const std::filesystem::path projectDir = root / "project";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );
    Prj::SetProjectBaseDir( projectDir.string() );
    Prj::SetIOBaseDirs( "case", "case" );
    Prj::SetDomainIOBaseDirs( "project", "", "", "" );

    EXPECT_EQ(
        std::filesystem::path( Prj::GetMeshInputFileName( "grid.ofl" ) ).lexically_normal(),
        ( projectDir / "grid.ofl" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( Prj::GetFlowInputFileName( "restart.dat" ) ).lexically_normal(),
        ( caseDir / "restart.dat" ).lexically_normal() );
}

TEST( PrjIOPaths, ProjectRootMustBeDeclared )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "OneFLOW_PrjIOPathsUndeclaredProjectTest"
        / "consumer";

    Prj::current_dir = std::filesystem::current_path().string();
    Prj::SetPrjBaseDir( caseDir.string() );

    EXPECT_THROW(
        Prj::SetIOBaseDirs( "project", "case" ),
        std::logic_error );
}
