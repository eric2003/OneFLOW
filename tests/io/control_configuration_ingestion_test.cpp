/*---------------------------------------------------------------------------*\\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\\*---------------------------------------------------------------------------*/

#include "ParaFile.h"
#include "DataBase.h"
#include "LogFile.h"
#include "Prj.h"
#include "SimuContextBootstrap.h"
#include "SimuContextCaseState.h"
#include "SimuContext.h"
#include "SimuContextIO.h"
#include "Parallel.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ONEFLOW;

namespace {

class ScopedCaseLogBinding {
public:
    explicit ScopedCaseLogBinding( const std::string & caseDir )
    {
        logFile.SetCaseDir( caseDir );
    }

    ~ScopedCaseLogBinding()
    {
        logFile.ClearCaseDir();
    }

    ScopedCaseLogBinding( const ScopedCaseLogBinding & ) = delete;
    ScopedCaseLogBinding & operator=( const ScopedCaseLogBinding & ) = delete;
};

class TemporaryControlCase {
public:
    TemporaryControlCase()
    {
        static std::atomic<unsigned int> nextId{ 0 };
        const auto id = nextId.fetch_add( 1 );
        const auto timestamp =
            std::chrono::steady_clock::now().time_since_epoch().count();
        path_ = std::filesystem::temp_directory_path()
            / ( "oneflow_control_ingestion_" + std::to_string( timestamp )
                + "_" + std::to_string( id ) );
        std::filesystem::create_directories( path_ / "script" );
    }

    ~TemporaryControlCase()
    {
        RequireGlobalDataBase().RequireDataPara().Clear();
        std::error_code error;
        std::filesystem::remove_all( path_, error );
    }

    std::filesystem::path ScriptPath( const std::string & name ) const
    {
        return path_ / "script" / name;
    }

    void WriteScript( const std::string & name, const std::string & contents ) const
    {
        std::ofstream output( ScriptPath( name ) );
        if ( ! output )
        {
            throw std::runtime_error( "unable to create control script" );
        }
        output << contents;
        if ( ! output )
        {
            throw std::runtime_error( "unable to write control script" );
        }
    }

    void WriteManifest( const std::string & contents ) const
    {
        std::ofstream output( ScriptPath( "control.txt" ) );
        if ( ! output )
        {
            throw std::runtime_error( "unable to create control manifest" );
        }
        output << contents;
        if ( ! output )
        {
            throw std::runtime_error( "unable to write control manifest" );
        }
    }

    std::string CaseDirectory() const { return path_.string(); }

private:
    std::filesystem::path path_;
};

void ClearRuntimeParameters()
{
    RequireGlobalDataBase().RequireDataPara().Clear();
}

class ScopedCaseRuntimeReset {
public:
    ScopedCaseRuntimeReset() = default;

    ~ScopedCaseRuntimeReset() noexcept
    {
        try
        {
            SimuContextDetail::ClearCaseConfigurationState();
        }
        catch ( ... )
        {
        }
    }

    ScopedCaseRuntimeReset( const ScopedCaseRuntimeReset & ) = delete;
    ScopedCaseRuntimeReset & operator=( const ScopedCaseRuntimeReset & ) = delete;
};

void EnsureParallelEnvironment()
{
    static const bool initialized = []()
    {
        SetUpParallelEnvironment();
        return true;
    }();
    (void)initialized;
}


TEST( ControlConfigurationIngestion, MissingScriptPreservesEarlierFileCommits )
{
    TemporaryControlCase testCase;
    testCase.WriteScript( "available.txt", "int ingestion_before_missing = 17;\n" );

    const std::vector< std::string > scriptFiles{
        testCase.ScriptPath( "available.txt" ).string(),
        testCase.ScriptPath( "missing.txt" ).string()
    };

    ClearRuntimeParameters();

    EXPECT_THROW( ReadMultiScriptFiles( scriptFiles ), std::runtime_error );

    const DataBase & database = RequireGlobalDataBase();
    EXPECT_EQ(
        GetDataValue<int>( "ingestion_before_missing", &database ), 17 );
}

TEST( ControlConfigurationIngestion, ParseFailureDoesNotCommitPartialCurrentFile )
{
    TemporaryControlCase testCase;
    testCase.WriteScript( "available.txt", "int ingestion_before_parse_error = 23;\n" );
    testCase.WriteScript(
        "invalid.txt",
        "int ingestion_uncommitted_prefix = 29;\n"
        "int ingestion_invalid_array[missing_array_size] = 1;\n" );

    const std::vector< std::string > scriptFiles{
        testCase.ScriptPath( "available.txt" ).string(),
        testCase.ScriptPath( "invalid.txt" ).string()
    };

    ClearRuntimeParameters();

    EXPECT_THROW( ReadMultiScriptFiles( scriptFiles ), std::runtime_error );

    const DataBase & database = RequireGlobalDataBase();
    EXPECT_EQ(
        GetDataValue<int>( "ingestion_before_parse_error", &database ), 23 );
    EXPECT_THROW(
        GetDataValue<int>( "ingestion_uncommitted_prefix", &database ),
        std::runtime_error );
}


TEST( ControlConfigurationIngestion, CommitFailureCanLeaveEarlierEntriesFromCurrentFile )
{
    TemporaryControlCase testCase;
    testCase.WriteScript(
        "base.txt",
        "int ingestion_conflict_value = 31;\n" );
    testCase.WriteScript(
        "conflict.txt",
        "int ingestion_committed_before_conflict = 37;\n"
        "Real ingestion_conflict_value = 0.5;\n" );

    const std::vector< std::string > scriptFiles{
        testCase.ScriptPath( "base.txt" ).string(),
        testCase.ScriptPath( "conflict.txt" ).string()
    };

    ClearRuntimeParameters();

    EXPECT_THROW( ReadMultiScriptFiles( scriptFiles ), std::runtime_error );

    const DataBase & database = RequireGlobalDataBase();
    EXPECT_EQ(
        GetDataValue<int>( "ingestion_conflict_value", &database ), 31 );
    EXPECT_EQ(
        GetDataValue<int>(
            "ingestion_committed_before_conflict", &database ), 37 );
}



} // namespace

TEST( ControlConfigurationIngestion, ReadsManifestScriptsAndCommitsTypedValuesInOrder )
{
    TemporaryControlCase testCase;
    testCase.WriteManifest(
        "\"override.txt\";\n"
        "\n"
        "\"base.txt\";\n" );
    testCase.WriteScript(
        "base.txt",
        "int ingestion_dimension = 2;\n"
        "Real ingestion_time_step = 0.001;\n"
        "string ingestion_task = \"Grid\";\n" );
    testCase.WriteScript(
        "override.txt",
        "int ingestion_dimension = 3;\n"
        "string ingestion_marker = \"second-file-loaded\";\n" );
    testCase.WriteScript(
        "unlisted.txt", "int ingestion_unlisted = 43;\n" );

    ClearRuntimeParameters();

    std::vector< std::string > scriptFiles;
    ReadScriptFileNameList( testCase.CaseDirectory(), scriptFiles );
    ASSERT_EQ( scriptFiles.size(), 3u );

    // Exercise the production case-level reader composition:
    // manifest discovery followed by ordered multi-file parsing and commit.
    ReadPrjScript( testCase.CaseDirectory() );

    DataBase & mutableDatabase = RequireGlobalDataBase();
    EXPECT_EQ(
        GetDataValue<int>( "ingestion_dimension", &mutableDatabase ), 3 );

    const DataBase & database = mutableDatabase;
    EXPECT_EQ( GetDataValue<int>( "ingestion_dimension", &database ), 3 );
    EXPECT_DOUBLE_EQ(
        GetDataValue<Real>( "ingestion_time_step", &database ), 0.001 );
    EXPECT_EQ(
        GetDataValue<std::string>( "ingestion_task", &database ), "Grid" );
    EXPECT_EQ(
        GetDataValue<std::string>( "ingestion_marker", &database ),
        "second-file-loaded" );
    EXPECT_EQ( GetDataValue<int>( "ingestion_unlisted", &database ), 43 );

    const ScriptFileDiscovery discovery =
        DiscoverCaseScriptFiles( testCase.CaseDirectory() );
    EXPECT_TRUE( std::any_of(
        discovery.diagnostics.begin(), discovery.diagnostics.end(),
        []( const std::string & item ) {
            return item.find( "filename order" ) != std::string::npos;
        } ) );
}


TEST( ControlReaderIntegration, CaseEnvironmentSetupFailureCanBeRolledBack )
{
    TemporaryControlCase testCase;
    ScopedCaseRuntimeReset runtimeReset;
    testCase.WriteManifest(
        "\"z_invalid.txt\";\n"
        "\"a_available.txt\";\n" );
    testCase.WriteScript(
        "a_available.txt",
        "string partial_case_bootstrap_value = \"must_be_rolled_back\";\n" );
    testCase.WriteScript(
        "z_invalid.txt",
        "int invalid_array[missing_array_size] = 1;\n" );

    ClearRuntimeParameters();
    Prj::current_dir = std::filesystem::current_path().string();
    EnsureParallelEnvironment();

    SimuContext context( testCase.CaseDirectory(), false );

    EXPECT_THROW( context.SetupCaseEnvironment(), std::runtime_error );
    EXPECT_FALSE( context.IsEnvironmentReady() );

    // Bootstrap reads actual files in filename order. The malformed second
    // file must roll back the first file and the case-relative bindings.
    EXPECT_TRUE( Prj::prjBaseDir.empty() );
    EXPECT_THROW(
        GetDataValue<std::string>(
            "partial_case_bootstrap_value", &RequireGlobalDataBase() ),
        std::runtime_error );

    // The same context can retry after the failed bootstrap has been rolled back.
    testCase.WriteScript(
        "z_invalid.txt",
        "string recovered_case_bootstrap_value = \"recovered\";\n" );
    EXPECT_NO_THROW( context.SetupCaseEnvironment() );
    EXPECT_TRUE( context.IsEnvironmentReady() );
    EXPECT_EQ(
        GetDataValue<std::string>(
            "recovered_case_bootstrap_value", &RequireGlobalDataBase() ),
        "recovered" );
}

TEST( ControlReaderIntegration, CaseBootstrapAppliesSettingsAndRebindsNextCase )
{
    TemporaryControlCase firstCase;
    TemporaryControlCase secondCase;
    // This guard is destroyed before the temporary case directories. It also
    // runs when a fatal GoogleTest assertion exits the test early.
    ScopedCaseRuntimeReset runtimeReset;

    const std::filesystem::path root =
        std::filesystem::path( firstCase.CaseDirectory() ).parent_path();
    const std::filesystem::path projectRoot = root / "project";
    const std::filesystem::path inputRoot = root / "input";
    const std::filesystem::path outputRoot = root / "output";
    const std::filesystem::path meshInputRoot = root / "mesh-input";
    const std::filesystem::path flowInputRoot = root / "flow-input";
    const std::filesystem::path meshOutputRoot = root / "mesh-output";
    const std::filesystem::path flowOutputRoot = root / "flow-output";

    firstCase.WriteManifest( "\"bootstrap.txt\";\n" );
    firstCase.WriteScript(
        "bootstrap.txt",
        "string projectRoot = \"" + projectRoot.generic_string() + "\";\n"
        "string inputRoot = \"" + inputRoot.generic_string() + "\";\n"
        "string outputRoot = \"" + outputRoot.generic_string() + "\";\n"
        "string meshInputRoot = \"" + meshInputRoot.generic_string() + "\";\n"
        "string flowInputRoot = \"" + flowInputRoot.generic_string() + "\";\n"
        "string meshOutputRoot = \"" + meshOutputRoot.generic_string() + "\";\n"
        "string flowOutputRoot = \"" + flowOutputRoot.generic_string() + "\";\n"
        "string ioPathInfoFileName = \"bootstrap_io_paths.txt\";\n" );
    secondCase.WriteManifest( "" );

    ClearRuntimeParameters();
    Prj::current_dir = std::filesystem::current_path().string();
    EnsureParallelEnvironment();

    SimuContextDetail::InitializeCaseConfiguration( firstCase.CaseDirectory() );

    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::projectBaseDir, projectRoot.string() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetInputBaseDir(), inputRoot.string() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetOutputBaseDir(), outputRoot.string() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetMeshInputBaseDir(), meshInputRoot.string() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetFlowInputBaseDir(), flowInputRoot.string() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetMeshOutputBaseDir(), meshOutputRoot.string() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetFlowOutputBaseDir(), flowOutputRoot.string() ) );
    EXPECT_TRUE( std::filesystem::exists( meshOutputRoot / "bootstrap_io_paths.txt" ) );
    EXPECT_TRUE( std::filesystem::exists( flowOutputRoot / "bootstrap_io_paths.txt" ) );

    // Use the same production cleanup operation as SimuContext::TeardownCase.
    // The process runtime remains active while case-owned state is released.
    SimuContextDetail::ClearCaseConfigurationState();

    EXPECT_THROW(
        GetDataValue<std::string>(
            "projectRoot", &RequireGlobalDataBase() ),
        std::runtime_error );
    EXPECT_TRUE( Prj::prjBaseDir.empty() );

    SimuContextDetail::InitializeCaseConfiguration( secondCase.CaseDirectory() );

    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::prjBaseDir, secondCase.CaseDirectory() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetInputBaseDir(), secondCase.CaseDirectory() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetOutputBaseDir(), secondCase.CaseDirectory() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetMeshInputBaseDir(), secondCase.CaseDirectory() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetFlowInputBaseDir(), secondCase.CaseDirectory() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetMeshOutputBaseDir(), secondCase.CaseDirectory() ) );
    EXPECT_TRUE( Prj::AreSameDirectoryPath( Prj::GetFlowOutputBaseDir(), secondCase.CaseDirectory() ) );
    EXPECT_TRUE( std::filesystem::exists(
        std::filesystem::path( secondCase.CaseDirectory() ) / "results" / "io_paths.txt" ) );

    SimuContextDetail::ClearCaseConfigurationState();
}

TEST( ControlReaderIntegration, ContextSetupPublishesReadinessAfterBootstrap )
{
    TemporaryControlCase testCase;
    ScopedCaseRuntimeReset runtimeReset;

    testCase.WriteManifest( "\"bootstrap.txt\";\n" );
    testCase.WriteScript(
        "bootstrap.txt",
        "int context_bootstrap_dimension = 2;\n" );

    ClearRuntimeParameters();
    EnsureParallelEnvironment();

    SimuContext context( testCase.CaseDirectory(), false );
    EXPECT_FALSE( context.IsEnvironmentReady() );

    EXPECT_NO_THROW( context.SetupCaseEnvironment() );

    EXPECT_TRUE( context.IsEnvironmentReady() );
    EXPECT_EQ(
        GetDataValue<int>(
            "context_bootstrap_dimension", &RequireGlobalDataBase() ),
        2 );
}

TEST( ControlReaderIntegration, FailedContextSetupDoesNotPublishReadiness )
{
    TemporaryControlCase testCase;
    ScopedCaseRuntimeReset runtimeReset;

    testCase.WriteManifest( "\"missing_bootstrap.txt\";\n" );

    ClearRuntimeParameters();
    EnsureParallelEnvironment();

    SimuContext context( testCase.CaseDirectory(), false );
    EXPECT_FALSE( context.IsEnvironmentReady() );

    EXPECT_THROW( context.SetupCaseEnvironment(), std::runtime_error );

    EXPECT_FALSE( context.IsEnvironmentReady() );
}

TEST( ControlReaderIntegration, ReadsCaseScriptsAndWritesDatabaseDump )
{
    TemporaryControlCase testCase;
    testCase.WriteManifest(
        "\"reader_base.txt\";\n"
        "\"reader_override.txt\";\n" );
    testCase.WriteScript(
        "reader_base.txt",
        "int ingestion_reader_dimension = 2;\n"
        "Real ingestion_reader_time_step = 0.002;\n" );
    testCase.WriteScript(
        "reader_override.txt",
        "int ingestion_reader_dimension = 3;\n"
        "string ingestion_reader_marker = \"reader-entry-loaded\";\n" );

    ClearRuntimeParameters();
    EnsureParallelEnvironment();

    {
        ScopedCaseLogBinding logBinding( testCase.CaseDirectory() );
        ReadControlInfo( testCase.CaseDirectory() );
    }

    const DataBase & database = RequireGlobalDataBase();
    EXPECT_EQ(
        GetDataValue<int>( "ingestion_reader_dimension", &database ), 3 );
    EXPECT_DOUBLE_EQ(
        GetDataValue<Real>( "ingestion_reader_time_step", &database ), 0.002 );
    EXPECT_EQ(
        GetDataValue<std::string>(
            "ingestion_reader_marker", &database ),
        "reader-entry-loaded" );

    const std::filesystem::path databaseDump =
        std::filesystem::path( testCase.CaseDirectory() )
        / "log" / "database.log";
    ASSERT_TRUE( std::filesystem::exists( databaseDump ) );

    std::ifstream input( databaseDump );
    ASSERT_TRUE( input.good() );
    const std::string dumpContents{
        std::istreambuf_iterator<char>( input ),
        std::istreambuf_iterator<char>() };
    EXPECT_NE(
        dumpContents.find( "ingestion_reader_dimension" ),
        std::string::npos );
    EXPECT_NE(
        dumpContents.find( "ingestion_reader_marker" ),
        std::string::npos );

    EXPECT_TRUE(
        std::filesystem::exists(
            std::filesystem::path( testCase.CaseDirectory() )
            / "log" / "log0.log" ) );
}


TEST( ControlConfigurationIngestion, ReadsParameterStyleControlAsARealFile )
{
    TemporaryControlCase testCase;
    testCase.WriteScript(
        "control.txt",
        "string simutask = \"Grid\";\n"
        "int dimension = 2;\n" );
    testCase.WriteScript( "grid.txt", "int ingestion_grid_scale = 4;\n" );

    ClearRuntimeParameters();
    ReadPrjScript( testCase.CaseDirectory() );

    const DataBase & database = RequireGlobalDataBase();
    EXPECT_EQ( GetDataValue<std::string>( "simutask", &database ), "Grid" );
    EXPECT_EQ( GetDataValue<int>( "dimension", &database ), 2 );
    EXPECT_EQ( GetDataValue<int>( "ingestion_grid_scale", &database ), 4 );
}

TEST( ControlConfigurationIngestion, DiscoversAndReadsFilesWithoutControl )
{
    TemporaryControlCase testCase;
    testCase.WriteScript(
        "simu.txt", "string simutask = \"Solve\";\n" );
    testCase.WriteScript(
        "cfd.txt", "int ingestion_max_steps = 25;\n" );

    ClearRuntimeParameters();
    ReadPrjScript( testCase.CaseDirectory() );

    const DataBase & database = RequireGlobalDataBase();
    EXPECT_EQ( GetDataValue<std::string>( "simutask", &database ), "Solve" );
    EXPECT_EQ( GetDataValue<int>( "ingestion_max_steps", &database ), 25 );
}
