#include <gtest/gtest.h>

#include "SimuCaseLifecycle.h"

#include <stdexcept>
#include <string>

using namespace ONEFLOW;

TEST( SimuCaseLifecycleTest, RunsCleanupAfterSuccessfulExecution )
{
    int executions = 0;
    int cleanups = 0;

    SimuCaseLifecycleDetail::ExecuteWithCleanup(
        [&]() { ++executions; },
        [&]() { ++cleanups; } );

    EXPECT_EQ( executions, 1 );
    EXPECT_EQ( cleanups, 1 );
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
