#include <gtest/gtest.h>
#include "SimuContext.h"
#include "SimuTask.h"
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ONEFLOW;

namespace {

int g_noop_executions = 0;

class NoOpTask : public ISimuTask
{
public:
    void Execute( SimuContext& /*ctx*/ ) override { ++g_noop_executions; }
};

// Task that *uses* the context - proves Execute(SimuContext&) plumbing.
class ContextAwareTask : public ISimuTask
{
public:
    void Execute( SimuContext& ctx ) override
    {
        saw_rank_ = ctx.Rank();
        saw_size_ = ctx.Size();
        saw_task_name_ = ctx.TaskName();
        saw_args_count_ = static_cast<int>( ctx.Args().size() );
        ctx.SetTaskByName( "ToyModel" );
        executed_ = true;
    }

    bool executed_ = false;
    int saw_rank_ = -1;
    int saw_size_ = -1;
    std::string saw_task_name_;
    int saw_args_count_ = -1;
};

} // namespace

TEST( SimuContextTest, StoresArgs )
{
    std::vector<std::string> args = { "OneFLOW", "grid.file", "control.file" };
    SimuContext ctx( args );

    ASSERT_EQ( ctx.Args().size(), 3u );
    EXPECT_EQ( ctx.Args()[ 0 ], "OneFLOW" );
    EXPECT_EQ( ctx.Args()[ 2 ], "control.file" );
}

TEST( SimuContextTest, DefaultParallelIsSerial )
{
    SimuContext ctx( std::vector<std::string>{ "OneFLOW" } );
    EXPECT_EQ( ctx.Rank(), 0 );
    EXPECT_EQ( ctx.Size(), 1 );
    EXPECT_FALSE( ctx.IsEnvironmentReady() );
    EXPECT_FALSE( ctx.IsTaskResolved() );
}

TEST( SimuContextTest, SetParallelInfoAcceptsValidRange )
{
    SimuContext ctx( std::vector<std::string>{} );
    ctx.SetParallelInfo( 2, 4 );
    EXPECT_EQ( ctx.Rank(), 2 );
    EXPECT_EQ( ctx.Size(), 4 );
}

TEST( SimuContextTest, SetParallelInfoRejectsBadSize )
{
    SimuContext ctx( std::vector<std::string>{} );
    EXPECT_THROW( ctx.SetParallelInfo( 0, 0 ), std::invalid_argument );
    EXPECT_THROW( ctx.SetParallelInfo( 0, -1 ), std::invalid_argument );
}

TEST( SimuContextTest, SetParallelInfoRejectsBadRank )
{
    SimuContext ctx( std::vector<std::string>{} );
    EXPECT_THROW( ctx.SetParallelInfo( -1, 4 ), std::invalid_argument );
    EXPECT_THROW( ctx.SetParallelInfo( 4, 4 ), std::invalid_argument );
}

TEST( SimuContextTest, SetTaskUpdatesNameAndFlag )
{
    SimuContext ctx( std::vector<std::string>{} );
    ctx.SetTask( TaskEnum::TOY_MODEL );
    EXPECT_TRUE( ctx.IsTaskResolved() );
    EXPECT_EQ( ctx.Task(), TaskEnum::TOY_MODEL );
    EXPECT_EQ( ctx.TaskName(), "ToyModel" );
}

TEST( SimuContextTest, SetTaskByNameRoundTrip )
{
    SimuContext ctx( std::vector<std::string>{} );
    ctx.SetTaskByName( "Grid" );
    EXPECT_EQ( ctx.Task(), TaskEnum::CREATE_GRID );
    EXPECT_EQ( ctx.TaskName(), "Grid" );
    EXPECT_EQ( ctx.TaskName(), TaskEnumToString( TaskEnum::CREATE_GRID ) );
}

TEST( SimuContextTest, SetTaskByNameRejectsUnknown )
{
    SimuContext ctx( std::vector<std::string>{} );
    EXPECT_THROW( ctx.SetTaskByName( "NotARealTask" ), std::invalid_argument );
    EXPECT_FALSE( ctx.IsTaskResolved() );
}

TEST( SimuContextTest, MarkEnvironmentReady )
{
    SimuContext ctx( std::vector<std::string>{} );
    EXPECT_FALSE( ctx.IsEnvironmentReady() );
    ctx.MarkEnvironmentReady( true );
    EXPECT_TRUE( ctx.IsEnvironmentReady() );
    ctx.MarkEnvironmentReady( false );
    EXPECT_FALSE( ctx.IsEnvironmentReady() );
}

TEST( SimuContextTest, OwnsAndClearsAcceleratorStates )
{
    class TestState final : public EulerDomainState
    {
    };

    SimuContext ctx( std::vector<std::string>{} );
    const EulerDomainStateKey key{ 2, 4, 0, AccelBackendKind::CPU };
    ctx.AccelStates().Insert( key, std::make_unique<TestState>() );

    EXPECT_EQ( ctx.AccelStates().Size(), 1u );
    EXPECT_TRUE( ctx.AccelStates().Contains( key ) );

    ctx.ClearAccelStates();
    EXPECT_EQ( ctx.AccelStates().Size(), 0u );
}

TEST( SimuContextTest, InjectedTaskNameWorksWithRegistry )
{
    g_noop_executions = 0;
    auto& reg = TaskRegistry::Instance();
    reg.Register( "ToyModel", []() {
        return std::make_unique<NoOpTask>();
    } );

    SimuContext ctx( std::vector<std::string>{} );
    ctx.SetTaskByName( "ToyModel" );

    auto task = TaskRegistry::Instance().Create( ctx.TaskName() );
    ASSERT_NE( task, nullptr );
    task->Execute( ctx );
    EXPECT_EQ( g_noop_executions, 1 );
}

// ---- phase 2.1: Execute receives SimuContext& ----

TEST( SimuContextTest, ContextAwareTaskSeesInjectedState )
{
    auto& reg = TaskRegistry::Instance();
    reg.Register( "UnitTest_ContextAware", []() {
        return std::make_unique<ContextAwareTask>();
    } );

    std::vector<std::string> args = { "OneFLOW", "a.grd", "b.ctr" };
    SimuContext ctx( args );
    ctx.SetParallelInfo( 3, 8 );
    ctx.SetTaskByName( "Theory" );

    auto task = reg.Create( "UnitTest_ContextAware" );
    ASSERT_NE( task, nullptr );

    // Same call shape as SimuImp::RunSimu: task->Execute( *ctx_ )
    task->Execute( ctx );

    // Downcast only to read test observations (production code uses interface only).
    auto* aware = dynamic_cast<ContextAwareTask*>( task.get() );
    ASSERT_NE( aware, nullptr );
    EXPECT_TRUE( aware->executed_ );
    EXPECT_EQ( aware->saw_rank_, 3 );
    EXPECT_EQ( aware->saw_size_, 8 );
    EXPECT_EQ( aware->saw_task_name_, "Theory" );
    EXPECT_EQ( aware->saw_args_count_, 3 );
    EXPECT_EQ( ctx.TaskName(), "ToyModel" );
}

TEST( SimuContextSolverNames, DefaultHasNoInjectedNames )
{
    std::vector<std::string> args;
    ONEFLOW::SimuContext ctx( args );
    EXPECT_FALSE( ctx.HasExpandedSolverNames() );
}

TEST( SimuContextSolverNames, SetAndClearExpandedNames )
{
    std::vector<std::string> args;
    ONEFLOW::SimuContext ctx( args );

    ONEFLOW::StringField names;
    names.push_back( "UNsSolver" );
    ctx.SetExpandedSolverNames( names );

    ASSERT_TRUE( ctx.HasExpandedSolverNames() );
    EXPECT_EQ( ctx.ExpandedSolverNames().size(), 1u );
    EXPECT_EQ( ctx.ExpandedSolverNames()[ 0 ], "UNsSolver" );

    ctx.ClearExpandedSolverNames();
    EXPECT_FALSE( ctx.HasExpandedSolverNames() );
}


TEST( SimuContextSolverNames, EnsureFillsOnlyWhenEmpty )
{
    std::vector<std::string> args;
    ONEFLOW::SimuContext ctx( args );

    ONEFLOW::StringField first;
    first.push_back( "UNsSolver" );
    ctx.EnsureExpandedSolverNames( first );

    ASSERT_TRUE( ctx.HasExpandedSolverNames() );
    EXPECT_EQ( ctx.ExpandedSolverNames()[ 0 ], "UNsSolver" );

    // Second ensure must not overwrite (production preload vs test injection).
    ONEFLOW::StringField second;
    second.push_back( "UTurbSolver" );
    ctx.EnsureExpandedSolverNames( second );

    ASSERT_EQ( ctx.ExpandedSolverNames().size(), 1u );
    EXPECT_EQ( ctx.ExpandedSolverNames()[ 0 ], "UNsSolver" );
}
