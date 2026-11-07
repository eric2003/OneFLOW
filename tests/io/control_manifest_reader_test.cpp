#include "gtest/gtest.h"
#include "ParaFile.h"
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <string>

using namespace ONEFLOW;

namespace {

class TemporaryScriptDirectory
{
public:
    TemporaryScriptDirectory()
    {
        static unsigned int nextId = 0;
        path_ = std::filesystem::temp_directory_path()
            / ( "oneflow_script_discovery_" + std::to_string( nextId++ ) );
        std::filesystem::remove_all( path_ );
        std::filesystem::create_directories( path_ / "script" );
    }

    ~TemporaryScriptDirectory()
    {
        std::error_code error;
        std::filesystem::remove_all( path_, error );
    }

    std::filesystem::path ScriptPath( const std::string & name ) const
    {
        return path_ / "script" / name;
    }

    std::string CasePath() const { return path_.string(); }

    void Write( const std::string & name, const std::string & content ) const
    {
        std::ofstream output( ScriptPath( name ) );
        ASSERT_TRUE( output.good() );
        output << content;
        ASSERT_TRUE( output.good() );
    }

private:
    std::filesystem::path path_;
};

} // namespace

TEST( ControlManifestReader, ReadsRealFilesInFilenameOrderAndChecksManifest )
{
    TemporaryScriptDirectory testCase;
    testCase.Write( "control.txt", "\"second.txt\";\n\"first.txt\";\n" );
    testCase.Write( "first.txt", "int first = 1;\n" );
    testCase.Write( "second.txt", "int second = 2;\n" );
    testCase.Write( "extra.txt", "int extra = 3;\n" );
    testCase.Write( "solver.txt", "NsSolver\n" );

    const ScriptFileDiscovery discovery =
        DiscoverCaseScriptFiles( testCase.CasePath() );

    ASSERT_EQ( discovery.mode, ScriptFileSourceMode::Manifest );
    EXPECT_FALSE( discovery.hasFatalErrors );
    ASSERT_EQ( discovery.filePaths.size(), 3u );
    EXPECT_EQ(
        std::filesystem::path( discovery.filePaths[ 0 ] ).filename().string(),
        "extra.txt" );
    EXPECT_EQ(
        std::filesystem::path( discovery.filePaths[ 1 ] ).filename().string(),
        "first.txt" );
    EXPECT_EQ(
        std::filesystem::path( discovery.filePaths[ 2 ] ).filename().string(),
        "second.txt" );
    EXPECT_TRUE( std::any_of( discovery.diagnostics.begin(), discovery.diagnostics.end(),
        []( const std::string & item ) {
            return item.find( "solver.txt" ) != std::string::npos;
        } ) );
    EXPECT_TRUE( std::any_of( discovery.diagnostics.begin(), discovery.diagnostics.end(),
        []( const std::string & item ) {
            return item.find( "extra.txt" ) != std::string::npos;
        } ) );
    EXPECT_TRUE( std::any_of( discovery.diagnostics.begin(), discovery.diagnostics.end(),
        []( const std::string & item ) {
            return item.find( "filename order" ) != std::string::npos;
        } ) );
}

TEST( ControlManifestReader, ReportsMissingAndDuplicateManifestEntries )
{
    TemporaryScriptDirectory testCase;
    testCase.Write(
        "control.txt",
        "\"present.txt\";\n\"missing.txt\";\n\"present.txt\";\n" );
    testCase.Write( "present.txt", "int present = 1;\n" );

    const ScriptFileDiscovery discovery =
        DiscoverCaseScriptFiles( testCase.CasePath() );

    ASSERT_EQ( discovery.filePaths.size(), 1u );
    EXPECT_EQ(
        std::filesystem::path( discovery.filePaths[ 0 ] ).filename().string(),
        "present.txt" );
    EXPECT_FALSE( discovery.hasFatalErrors );
    EXPECT_TRUE( std::any_of( discovery.diagnostics.begin(), discovery.diagnostics.end(),
        []( const std::string & item ) {
            return item.find( "missing.txt" ) != std::string::npos;
        } ) );
}

TEST( ControlManifestReader, ReadsParameterStyleControlAsAnOrdinaryFile )
{
    TemporaryScriptDirectory testCase;
    testCase.Write(
        "control.txt",
        "string simutask = \"Grid\";\nint dimension = 2;\n" );
    testCase.Write( "grid.txt", "int gridScale = 1;\n" );

    const ScriptFileDiscovery discovery =
        DiscoverCaseScriptFiles( testCase.CasePath() );

    EXPECT_EQ( discovery.mode, ScriptFileSourceMode::ParameterFile );
    ASSERT_EQ( discovery.filePaths.size(), 2u );
    EXPECT_EQ(
        std::filesystem::path( discovery.filePaths[ 0 ] ).filename().string(),
        "control.txt" );
    EXPECT_EQ(
        std::filesystem::path( discovery.filePaths[ 1 ] ).filename().string(),
        "grid.txt" );
}

TEST( ControlManifestReader, DiscoversFilesWhenControlIsAbsent )
{
    TemporaryScriptDirectory testCase;
    testCase.Write( "simu.txt", "string simutask = \"Solve\";\n" );
    testCase.Write( "cfd.txt", "int maxSteps = 20;\n" );
    testCase.Write( "solver.txt", "NsSolver\n" );

    const ScriptFileDiscovery discovery =
        DiscoverCaseScriptFiles( testCase.CasePath() );

    EXPECT_EQ( discovery.mode, ScriptFileSourceMode::Directory );
    ASSERT_EQ( discovery.filePaths.size(), 2u );
    ASSERT_EQ( discovery.diagnostics.size(), 2u );
    EXPECT_NE( discovery.diagnostics[ 1 ].find( "solver.txt" ), std::string::npos );
    EXPECT_EQ(
        std::filesystem::path( discovery.filePaths[ 0 ] ).filename().string(),
        "cfd.txt" );
    EXPECT_EQ(
        std::filesystem::path( discovery.filePaths[ 1 ] ).filename().string(),
        "simu.txt" );
}
