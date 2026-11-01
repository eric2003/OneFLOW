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
----------------------------------------------------------------------------*/

#include "DataBase.h"
#include "Prj.h"
#include "SimuImp.h"

#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

using namespace ONEFLOW;

namespace
{

class TemporaryCase
{
public:
    TemporaryCase()
    {
        static std::atomic<unsigned int> nextId{ 0 };
        const auto id = nextId.fetch_add( 1 );
        const auto timestamp =
            std::chrono::steady_clock::now().time_since_epoch().count();
        path_ = std::filesystem::temp_directory_path()
            / ( "oneflow_simu_imp_lifecycle_" + std::to_string( timestamp )
                + "_" + std::to_string( id ) );
        std::filesystem::create_directories( path_ / "script" );
    }

    ~TemporaryCase()
    {
        std::error_code error;
        std::filesystem::remove_all( path_, error );
    }

    void WriteScript( const std::string & name, const std::string & contents )
    {
        std::ofstream output( path_ / "script" / name );
        if ( ! output )
        {
            throw std::runtime_error( "unable to create case script" );
        }
        output << contents;
        if ( ! output )
        {
            throw std::runtime_error( "unable to write case script" );
        }
    }

    void WriteManifest( const std::string & contents )
    {
        std::ofstream output( path_ / "script" / "control.txt" );
        if ( ! output )
        {
            throw std::runtime_error( "unable to create case manifest" );
        }
        output << contents;
        if ( ! output )
        {
            throw std::runtime_error( "unable to write case manifest" );
        }
    }

    std::string Directory() const
    {
        return path_.string();
    }

private:
    std::filesystem::path path_;
};

class ProcessRuntimeGuard
{
public:
    explicit ProcessRuntimeGuard( SimuImp & simulation )
        : simulation_( simulation )
    {
    }

    ~ProcessRuntimeGuard() noexcept
    {
        if ( simulation_.Context().IsProcessEnvironmentActive() )
        {
            try
            {
                simulation_.FinalizeEnvironment();
            }
            catch ( ... )
            {
            }
        }
    }

    ProcessRuntimeGuard( const ProcessRuntimeGuard & ) = delete;
    ProcessRuntimeGuard & operator=( const ProcessRuntimeGuard & ) = delete;

private:
    SimuImp & simulation_;
};

TEST( SimuImpCaseLifecycleIntegration,
      FailedProductionSetupRunsRealCaseTeardownAndKeepsProcessRuntimeActive )
{
    TemporaryCase testCase;
    testCase.WriteManifest(
        "\"available.txt\";\n"
        "\"missing.txt\";\n" );
    testCase.WriteScript(
        "available.txt",
        "string partial_simulation_value = \"must_be_cleared\";\n" );

    SimuImp simulation( testCase.Directory(), false );
    SimuContext & context = simulation.Context();

    // Seed context-owned state so the real TeardownCase reset is observable.
    context.SetTask( TaskEnum::SOLVE_FIELD );
    StringField expandedNames;
    expandedNames.push_back( "UNsSolver" );
    context.SetExpandedSolverNames( expandedNames );

    context.SetupProcessEnvironment();
    ProcessRuntimeGuard runtimeGuard( simulation );

    ASSERT_TRUE( context.IsProcessEnvironmentActive() );
    EXPECT_THROW( simulation.RunCase(), std::runtime_error );

    // These assertions distinguish the production TeardownCase() call from
    // the narrower configuration rollback performed by SetupCaseEnvironment().
    EXPECT_FALSE( context.IsEnvironmentReady() );
    EXPECT_FALSE( context.IsTaskResolved() );
    EXPECT_FALSE( context.HasExpandedSolverNames() );
    EXPECT_EQ( context.TaskName(), "Solve" );
    EXPECT_TRUE( Prj::prjBaseDir.empty() );
    EXPECT_THROW(
        GetDataValue<std::string>(
            "partial_simulation_value", &RequireGlobalDataBase() ),
        std::runtime_error );

    // RunCase owns per-case teardown, not process finalization.
    EXPECT_TRUE( context.IsProcessEnvironmentActive() );
    EXPECT_NO_THROW( simulation.FinalizeEnvironment() );
    EXPECT_FALSE( context.IsProcessEnvironmentActive() );
}

} // namespace
