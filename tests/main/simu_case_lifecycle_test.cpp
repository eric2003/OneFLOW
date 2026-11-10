#include <gtest/gtest.h>

#include "SimuCaseLifecycle.h"
#include "SimuProcessLifecycle.h"

#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace ONEFLOW;

TEST( SimuCaseLifecycleTest, RunsCleanupAfterSuccessfulExecution )
{
    std::vector< std::string > events;

    SimuCaseLifecycleDetail::ExecuteWithCleanup(
        [&]() { events.emplace_back( "execute" ); },
        [&]() { events.emplace_back( "cleanup" ); } );

    ASSERT_EQ( events.size(), 2 );
    EXPECT_EQ( events[ 0 ], "execute" );
    EXPECT_EQ( events[ 1 ], "cleanup" );
}

TEST( SimuCaseLifecycleTest, RunsCleanupAndPreservesExecutionFailure )
{
    int cleanups = 0;

    try
    {
        SimuCaseLifecycleDetail::ExecuteWithCleanup(
            []() { throw std::runtime_error( "setup failed" ); },
            [&]() { ++cleanups; } );
        FAIL() << "Expected the execution exception";
    }
    catch ( const std::runtime_error & error )
    {
        EXPECT_EQ( std::string( error.what() ), "setup failed" );
    }

    EXPECT_EQ( cleanups, 1 );
}

TEST( SimuCaseLifecycleTest, CleanupFailureDoesNotReplaceExecutionFailure )
{
    int cleanups = 0;

    try
    {
        SimuCaseLifecycleDetail::ExecuteWithCleanup(
            []() { throw std::runtime_error( "original failure" ); },
            [&]() {
                ++cleanups;
                throw std::logic_error( "cleanup failure" );
            } );
        FAIL() << "Expected the execution exception";
    }
    catch ( const std::runtime_error & error )
    {
        EXPECT_EQ( std::string( error.what() ), "original failure" );
    }

    EXPECT_EQ( cleanups, 1 );
}

TEST( SimuCaseLifecycleTest, SuccessfulExecutionPropagatesCleanupFailure )
{
    EXPECT_THROW(
        SimuCaseLifecycleDetail::ExecuteWithCleanup(
            []() {},
            []() { throw std::logic_error( "cleanup failure" ); } ),
        std::logic_error );
}
TEST( SimuCaseLifecycleTest, RunsEveryCleanupStepAndPreservesFirstFailure )
{
    std::vector< std::string > events;

    try
    {
        SimuCaseLifecycleDetail::ExecuteAllCleanupSteps(
            [&]()
            {
                events.emplace_back( "first" );
                throw std::runtime_error( "first cleanup failure" );
            },
            [&]()
            {
                events.emplace_back( "second" );
                throw std::logic_error( "second cleanup failure" );
            },
            [&]() { events.emplace_back( "third" ); } );
        FAIL() << "Expected the first cleanup exception";
    }
    catch ( const std::runtime_error & error )
    {
        EXPECT_EQ( std::string( error.what() ), "first cleanup failure" );
    }

    const std::vector< std::string > expected = {
        "first", "second", "third"
    };
    EXPECT_EQ( events, expected );
}

TEST( SimuCaseLifecycleTest, ExecutesAllCleanupStepsWhenNoStepFails )
{
    std::vector< int > events;

    SimuCaseLifecycleDetail::ExecuteAllCleanupSteps(
        [&]() { events.push_back( 1 ); },
        [&]() { events.push_back( 2 ); },
        [&]() { events.push_back( 3 ); } );

    const std::vector< int > expected = { 1, 2, 3 };
    EXPECT_EQ( events, expected );
}

TEST( SimuCaseLifecycleTest, SetupRollbackDoesNotRunAfterSuccess )
{
    std::vector< std::string > events;

    SimuCaseLifecycleDetail::ExecuteWithRollback(
        [&]() { events.emplace_back( "setup" ); },
        [&]() { events.emplace_back( "rollback" ); } );

    const std::vector< std::string > expected = { "setup" };
    EXPECT_EQ( events, expected );
}

TEST( SimuCaseLifecycleTest, SetupFailureRunsRollbackAndPreservesFailure )
{
    std::vector< std::string > events;

    try
    {
        SimuCaseLifecycleDetail::ExecuteWithRollback(
            [&]()
            {
                events.emplace_back( "setup" );
                throw std::runtime_error( "setup failed" );
            },
            [&]()
            {
                events.emplace_back( "rollback" );
                throw std::logic_error( "rollback failed" );
            } );
        FAIL() << "Expected the setup exception";
    }
    catch ( const std::runtime_error & error )
    {
        EXPECT_EQ( std::string( error.what() ), "setup failed" );
    }

    const std::vector< std::string > expected = { "setup", "rollback" };
    EXPECT_EQ( events, expected );
}

TEST( SimuCaseLifecycleTest, FinalizesSingleCaseProcessAfterExecutionFailure )
{
    std::vector< std::string > events;

    try
    {
        SimuCaseLifecycleDetail::ExecuteWithProcessLifecycle(
            [&]() { events.emplace_back( "setup" ); },
            [&]()
            {
                events.emplace_back( "case" );
                throw std::runtime_error( "case failed" );
            },
            [&]() { events.emplace_back( "finalize" ); } );
        FAIL() << "Expected the case execution exception";
    }
    catch ( const std::runtime_error & error )
    {
        EXPECT_EQ( std::string( error.what() ), "case failed" );
    }

    const std::vector< std::string > expected = {
        "setup", "case", "finalize"
    };
    EXPECT_EQ( events, expected );
}
TEST( SimuCaseLifecycleTest, RunsMultipleCasesInCommandLineOrder )
{
    const std::vector< std::string > caseDirs = {
        "workspace/grid/mesh-layout/3d/minimal_cube",
        "workspace/grid/mesh_workflow",
        "workspace/grid/mesh-layout/2d/minimal_rectangle",
        "workspace/rae2822grid",
        "workspace/rae2822_roe_sa",
        "workspace/platescript"
    };
    std::vector< std::string > events;

    SimuCaseLifecycleDetail::RunCasesWithProcessLifecycle(
        [&]() { events.emplace_back( "setup" ); },
        caseDirs,
        [&]( const std::string & caseDir ) { events.push_back( caseDir ); },
        [&]() { events.emplace_back( "finalize" ); } );

    std::vector< std::string > expected = { "setup" };
    expected.insert( expected.end(), caseDirs.begin(), caseDirs.end() );
    expected.emplace_back( "finalize" );
    EXPECT_EQ( events, expected );
}

TEST( SimuCaseLifecycleTest, RunsRepeatedCasePathsMoreThanOnce )
{
    const std::string caseDir = "workspace/rae2822grid";
    const std::vector< std::string > caseDirs = { caseDir, caseDir };
    std::vector< std::string > events;

    SimuCaseLifecycleDetail::RunCasesWithProcessLifecycle(
        [&]() { events.emplace_back( "setup" ); },
        caseDirs,
        [&]( const std::string & caseDir ) { events.push_back( caseDir ); },
        [&]() { events.emplace_back( "finalize" ); } );

    const std::vector< std::string > expected = {
        "setup", caseDir, caseDir, "finalize"
    };
    EXPECT_EQ( events, expected );
}

TEST( SimuCaseLifecycleTest, FinalizesProcessAfterCaseFailure )
{
    const std::vector< std::string > caseDirs = { "first", "failing", "last" };
    std::vector< std::string > events;

    try
    {
        SimuCaseLifecycleDetail::RunCasesWithProcessLifecycle(
            [&]() { events.emplace_back( "setup" ); },
            caseDirs,
            [&]( const std::string & caseDir )
            {
                events.push_back( caseDir );
                if ( caseDir == "failing" )
                {
                    throw std::runtime_error( "case failed" );
                }
            },
            [&]() { events.emplace_back( "finalize" ); } );
        FAIL() << "Expected the case execution exception";
    }
    catch ( const std::runtime_error & error )
    {
        EXPECT_EQ( std::string( error.what() ), "case failed" );
    }

    const std::vector< std::string > expected = {
        "setup", "first", "failing", "finalize"
    };
    EXPECT_EQ( events, expected );
}

TEST( SimuCaseLifecycleTest, FinalizesProcessAfterSetupFailure )
{
    std::vector< std::string > events;

    EXPECT_THROW(
        SimuCaseLifecycleDetail::RunCasesWithProcessLifecycle(
            [&]()
            {
                events.emplace_back( "setup" );
                throw std::runtime_error( "process setup failed" );
            },
            std::vector< std::string >{ "first" },
            [&]( const std::string & caseDir ) { events.push_back( caseDir ); },
            [&]() { events.emplace_back( "finalize" ); } ),
        std::runtime_error );

    const std::vector< std::string > expected = { "setup", "finalize" };
    EXPECT_EQ( events, expected );
}

TEST( SimuProcessLifecycleTest, RejectsRepeatedInitialization )
{
    SimuProcessLifecycleDetail::ProcessRuntimeState state;

    state.Start();
    EXPECT_TRUE( state.IsActive() );
    EXPECT_THROW( state.Start(), std::logic_error );
}

TEST( SimuProcessLifecycleTest, FailedCleanupDoesNotCommitFinalizedState )
{
    SimuProcessLifecycleDetail::ProcessRuntimeState state;
    state.Start();

    int laterCleanupCount = 0;
    EXPECT_THROW(
        SimuCaseLifecycleDetail::ExecuteAllCleanupSteps(
            []() { throw std::runtime_error( "accelerator release failed" ); },
            [&]() { ++laterCleanupCount; } ),
        std::runtime_error );

    EXPECT_EQ( laterCleanupCount, 1 );
    EXPECT_TRUE( state.IsActive() );
    EXPECT_FALSE( state.IsFinalized() );
}

TEST( SimuProcessLifecycleTest, SuccessfulCleanupCanCommitFinalizedState )
{
    SimuProcessLifecycleDetail::ProcessRuntimeState state;
    state.Start();

    EXPECT_NO_THROW(
        SimuCaseLifecycleDetail::ExecuteAllCleanupSteps(
            []() {},
            []() {} ) );
    state.MarkFinalized();

    EXPECT_FALSE( state.IsActive() );
    EXPECT_TRUE( state.IsFinalized() );
}

TEST( SimuProcessLifecycleTest, FinalizationIsIdempotent )
{
    SimuProcessLifecycleDetail::ProcessRuntimeState state;

    EXPECT_FALSE( state.IsActive() );
    EXPECT_FALSE( state.IsFinalized() );

    state.Start();
    state.MarkFinalized();
    state.MarkFinalized();

    EXPECT_FALSE( state.IsActive() );
    EXPECT_TRUE( state.IsFinalized() );
    EXPECT_THROW( state.Start(), std::logic_error );
}


TEST( SimuProcessLifecycleTest, MoveConstructionTransfersActiveRuntimeState )
{
    SimuProcessLifecycleDetail::ProcessRuntimeState source;
    source.Start();

    SimuProcessLifecycleDetail::ProcessRuntimeState destination(
        std::move( source ) );

    EXPECT_TRUE( destination.IsActive() );
    EXPECT_FALSE( source.IsActive() );
    EXPECT_FALSE( source.IsFinalized() );
    EXPECT_THROW( source.Start(), std::logic_error );

    destination.MarkFinalized();
    EXPECT_TRUE( destination.IsFinalized() );
}

TEST( SimuProcessLifecycleTest, RuntimeStateCannotBeCopiedOrMoveAssigned )
{
    EXPECT_FALSE( std::is_copy_constructible<
        SimuProcessLifecycleDetail::ProcessRuntimeState >::value );
    EXPECT_FALSE( std::is_copy_assignable<
        SimuProcessLifecycleDetail::ProcessRuntimeState >::value );
    EXPECT_FALSE( std::is_move_assignable<
        SimuProcessLifecycleDetail::ProcessRuntimeState >::value );
}
