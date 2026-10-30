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
#include <gtest/gtest.h>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ONEFLOW;

namespace {

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


TEST( ControlConfigurationIngestion, MissingScriptPreservesEarlierFileCommits )
{
    TemporaryControlCase testCase;
    testCase.WriteScript( "available.txt", "int ingestion_before_missing = 17;\\n" );

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
    testCase.WriteScript( "available.txt", "int ingestion_before_parse_error = 23;\\n" );
    testCase.WriteScript(
        "invalid.txt",
        "int ingestion_uncommitted_prefix = 29;\\n"
        "int ingestion_invalid_array[missing_array_size] = 1;\\n" );

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


} // namespace

TEST( ControlConfigurationIngestion, ReadsManifestScriptsAndCommitsTypedValuesInOrder )
{
    TemporaryControlCase testCase;
    testCase.WriteManifest(
        "\"base.txt\";\n"
        "\n"
        "\"override.txt\";\n" );
    testCase.WriteScript(
        "base.txt",
        "int ingestion_dimension = 2;\n"
        "Real ingestion_time_step = 0.001;\n"
        "string ingestion_task = \"Grid\";\n" );
    testCase.WriteScript(
        "override.txt",
        "int ingestion_dimension = 3;\n"
        "string ingestion_marker = \"second-file-loaded\";\n" );

    ClearRuntimeParameters();

    std::vector< std::string > scriptFiles;
    ReadScriptFileNameList( testCase.CaseDirectory(), scriptFiles );
    ASSERT_EQ( scriptFiles.size(), 2u );

    ReadMultiScriptFiles( scriptFiles );

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
}
