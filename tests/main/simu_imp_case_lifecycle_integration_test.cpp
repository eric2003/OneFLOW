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

#include "AccelBackend.h"
#include "AccelRuntime.h"
#include "EulerDomain.h"
#include "SimuCaseLifecycle.h"
#include "DataBase.h"
#include "Prj.h"
#include "SimuImp.h"
#include "SimuTask.h"

#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

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

std::vector<std::string> successfulCaseMarkers;
std::vector<int> successfulCaseDimensions;

class CaseMarkerTask final : public ISimuTask
{
public:
    void Execute( SimuContext & /*context*/ ) override
    {
        DataBase & database = RequireGlobalDataBase();
        successfulCaseMarkers.push_back(
            GetDataValue<std::string>( "case_marker", &database ) );
        successfulCaseDimensions.push_back(
            GetDataValue<int>( "dimension", &database ) );
    }
};


class RuntimeBoundTestState final : public EulerDomainState
{
public:
    explicit RuntimeBoundTestState( bool & destroyedWhileRuntimeActive )
        : destroyedWhileRuntimeActive_( destroyedWhileRuntimeActive )
    {
    }

    ~RuntimeBoundTestState() override
    {
        destroyedWhileRuntimeActive_ =
            AccelRuntime::Instance().IsInitialized();
    }

private:
    bool & destroyedWhileRuntimeActive_;
};

bool failProcessBackendFinalize = false;
int processBackendFinalizeCalls = 0;

class RetryableProcessHipBackend final : public AccelBackend
{
public:
    AccelBackendKind Kind() const override { return AccelBackendKind::HIP; }
    const char * Name() const override { return "RETRYABLE_PROCESS_HIP"; }
    bool IsAccelerator() const override { return true; }
    bool SupportsMultiDevice() const override { return true; }

    void Initialize( AccelContext & context ) override
    {
        context.deviceCount = 1;
        context.selectedDevice = 0;
    }

    void Finalize() override
    {
        ++processBackendFinalizeCalls;
        if ( failProcessBackendFinalize )
        {
            throw std::runtime_error(
                "simulated process backend finalization failure" );
        }
    }

    void Synchronize() override {}
    void * Allocate( std::size_t ) override { return nullptr; }
    void Deallocate( void * ) override {}
    void Copy( void *, const void *, std::size_t, AccelCopyKind ) override {}
};

void SetBackendEnvironment( const char * name, const char * value )
{
#ifdef _WIN32
    _putenv_s( name, value );
#else
    setenv( name, value, 1 );
#endif
}

class BackendEnvironmentGuard
{
public:
    explicit BackendEnvironmentGuard( const char * name ) : name_( name )
    {
        const char * value = std::getenv( name );
        hadValue_ = value != nullptr;
        if ( hadValue_ ) oldValue_ = value;
    }

    ~BackendEnvironmentGuard()
    {
        if ( hadValue_ )
        {
            SetBackendEnvironment( name_.c_str(), oldValue_.c_str() );
        }
        else
        {
#ifdef _WIN32
            _putenv_s( name_.c_str(), "" );
#else
            unsetenv( name_.c_str() );
#endif
        }
    }

private:
    std::string name_;
    std::string oldValue_;
    bool hadValue_ = false;
};

class ProcessBackendFailureGuard
{
public:
    ~ProcessBackendFailureGuard()
    {
        failProcessBackendFinalize = false;
    }
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
      ProcessFinalizationReleasesContextBackendStatesBeforeRuntime )
{
    // SimuImp resolves its case path during construction. Initialize the
    // production project environment and pass an absolute temporary path so
    // this test does not depend on an implicit process working directory.
    Prj::Init();
    TemporaryCase testCase;
    SimuImp simulation( testCase.Directory(), false );
    SimuContext & context = simulation.Context();

    context.SetupProcessEnvironment();
    ProcessRuntimeGuard runtimeGuard( simulation );

    ASSERT_TRUE( context.IsProcessEnvironmentActive() );
    ASSERT_TRUE( AccelRuntime::Instance().IsInitialized() );

    bool destroyedWhileRuntimeActive = false;
    const EulerDomainStateKey key{
        0, 0, 0, AccelBackendKind::CPU, -1
    };
    context.AccelStates().Insert(
        key,
        std::make_unique< RuntimeBoundTestState >(
            destroyedWhileRuntimeActive ) );
    ASSERT_EQ( context.AccelStates().Size(), 1u );

    EXPECT_NO_THROW( simulation.FinalizeEnvironment() );

    EXPECT_TRUE( destroyedWhileRuntimeActive );
    EXPECT_EQ( context.AccelStates().Size(), 0u );
    EXPECT_FALSE( AccelRuntime::Instance().IsInitialized() );
    EXPECT_FALSE( context.IsProcessEnvironmentActive() );
}

TEST( SimuImpCaseLifecycleIntegration,
      FailedAcceleratorFinalizationKeepsProcessRuntimeActiveForRetry )
{
    BackendEnvironmentGuard backendSetting( "ONEFLOW_ACCEL_BACKEND" );
    ProcessBackendFailureGuard failureGuard;
    SetBackendEnvironment( "ONEFLOW_ACCEL_BACKEND", "HIP" );
    AccelBackendRegistry::Instance().Register(
        AccelBackendKind::HIP,
        []() { return std::make_unique<RetryableProcessHipBackend>(); } );

    processBackendFinalizeCalls = 0;
    failProcessBackendFinalize = false;

    Prj::Init();
    TemporaryCase testCase;
    SimuImp simulation( testCase.Directory(), false );
    SimuContext & context = simulation.Context();
    context.SetupProcessEnvironment();
    ProcessRuntimeGuard runtimeGuard( simulation );

    ASSERT_TRUE( context.IsProcessEnvironmentActive() );
    ASSERT_TRUE( AccelRuntime::Instance().IsInitialized() );
    ASSERT_STREQ(
        AccelRuntime::Instance().Backend().Name(),
        "RETRYABLE_PROCESS_HIP" );

    failProcessBackendFinalize = true;
    EXPECT_THROW(
        simulation.FinalizeEnvironment(),
        std::runtime_error );

    EXPECT_TRUE( context.IsProcessEnvironmentActive() );
    EXPECT_TRUE( AccelRuntime::Instance().IsInitialized() );
    EXPECT_EQ( processBackendFinalizeCalls, 1 );

    failProcessBackendFinalize = false;
    EXPECT_NO_THROW( simulation.FinalizeEnvironment() );

    EXPECT_FALSE( context.IsProcessEnvironmentActive() );
    EXPECT_FALSE( AccelRuntime::Instance().IsInitialized() );
    EXPECT_EQ( processBackendFinalizeCalls, 2 );
}

TEST( SimuImpCaseLifecycleIntegration,
      ProductionRuntimeContainsCpuBackendRegistration )
{
    EXPECT_TRUE(
        AccelBackendRegistry::Instance().Contains( AccelBackendKind::CPU ) );
}

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


TEST( SimuImpCaseLifecycleIntegration,
      SuccessfulCasesRemainIsolatedWithinOneProcess )
{
    TemporaryCase firstCase;
    TemporaryCase secondCase;

    firstCase.WriteManifest( "\"case.txt\";\n" );
    firstCase.WriteScript(
        "case.txt",
        "string simutask = \"ToyModel\";\n"
        "string case_marker = \"first_case\";\n"
        "int dimension = 2;\n" );

    secondCase.WriteManifest( "\"case.txt\";\n" );
    secondCase.WriteScript(
        "case.txt",
        "string simutask = \"ToyModel\";\n"
        "string case_marker = \"second_case\";\n"
        "int dimension = 3;\n" );

    successfulCaseMarkers.clear();
    successfulCaseDimensions.clear();
    TaskRegistry::Instance().Register(
        "ToyModel",
        []() { return std::make_unique<CaseMarkerTask>(); } );

    Prj::Init();
    SimuContext processContext( std::vector<std::string>{} );
    const std::vector<std::string> caseDirectories{
        firstCase.Directory(), secondCase.Directory()
    };
    std::vector<bool> caseStateCleared;
    std::vector<bool> processRuntimeStayedActive;

    EXPECT_NO_THROW(
        SimuCaseLifecycleDetail::RunCasesWithProcessLifecycle(
            [&]() { processContext.SetupProcessEnvironment(); },
            caseDirectories,
            [&]( const std::string & caseDirectory )
            {
                SimuImp caseSimulation( caseDirectory, false );
                caseSimulation.RunCase();

                const SimuContext & caseContext = caseSimulation.Context();
                bool markerWasCleared = false;
                try
                {
                    (void)GetDataValue<std::string>(
                        "case_marker", &RequireGlobalDataBase() );
                }
                catch ( const std::runtime_error & )
                {
                    markerWasCleared = true;
                }

                caseStateCleared.push_back(
                    ! caseContext.IsEnvironmentReady()
                    && ! caseContext.IsTaskResolved()
                    && ! caseContext.HasExpandedSolverNames()
                    && Prj::prjBaseDir.empty()
                    && markerWasCleared );
                processRuntimeStayedActive.push_back(
                    processContext.IsProcessEnvironmentActive() );
            },
            [&]() { processContext.FinalizeEnvironment(); } ) );

    EXPECT_EQ(
        successfulCaseMarkers,
        ( std::vector<std::string>{ "first_case", "second_case" } ) );
    EXPECT_EQ(
        successfulCaseDimensions,
        ( std::vector<int>{ 2, 3 } ) );
    EXPECT_TRUE( std::filesystem::exists(
        std::filesystem::path( firstCase.Directory() ) / "log" / "database.log" ) );
    EXPECT_TRUE( std::filesystem::exists(
        std::filesystem::path( secondCase.Directory() ) / "log" / "database.log" ) );
    EXPECT_TRUE( std::filesystem::exists(
        std::filesystem::path( firstCase.Directory() ) / "results" / "io_paths.txt" ) );
    EXPECT_TRUE( std::filesystem::exists(
        std::filesystem::path( secondCase.Directory() ) / "results" / "io_paths.txt" ) );
    ASSERT_EQ( caseStateCleared.size(), 2u );
    EXPECT_TRUE( caseStateCleared[ 0 ] );
    EXPECT_TRUE( caseStateCleared[ 1 ] );
    EXPECT_EQ(
        processRuntimeStayedActive,
        ( std::vector<bool>{ true, true } ) );
    EXPECT_FALSE( processContext.IsProcessEnvironmentActive() );
}

} // namespace
