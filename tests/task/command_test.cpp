//command_test.cpp
#include <gtest/gtest.h>

#include "Command.h"
#include "Task.h"
#include "TaskState.h"

#include <memory>
#include <type_traits>

namespace
{

    int g_runCount = 0;
    int g_destroyCount = 0;
    std::vector< int > g_executionLog;

    class CountingTask : public ONEFLOW::Task
    {
    public:
        CountingTask() = default;

        ~CountingTask() override
        {
            ++ g_destroyCount;
        }

        void Run() override
        {
            ++ g_runCount;
        }
    };

    class QueueTestTask : public ONEFLOW::Task
    {
    public:
        explicit QueueTestTask( int id )
            : id_( id )
        {
        }

        ~QueueTestTask() override
        {
            ++ g_destroyCount;
        }

        void Run() override
        {
            g_executionLog.push_back( id_ );
        }

    private:
        int id_;
    };

    class QueueTestCommand : public ONEFLOW::Command
    {
    public:
        explicit QueueTestCommand( int id )
            : id_( id )
        {
        }

        void Execute() override
        {
            g_executionLog.push_back( id_ );
        }

    private:
        int id_;
    };


}

class CommandTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        g_runCount = 0;
        g_destroyCount = 0;
        ONEFLOW::TaskState::task = nullptr;
    }

    void TearDown() override
    {
        ONEFLOW::TaskState::task = nullptr;
    }
};

TEST_F( CommandTest, DefaultCommandHasNoTasks )
{
    ONEFLOW::SimpleCmd command;

    ASSERT_NE( command.GetTaskList(), nullptr );
    EXPECT_TRUE( command.GetTaskList()->empty() );
}

TEST_F( CommandTest, AddRawTaskStoresTask )
{
    ONEFLOW::SimpleCmd command;

    ONEFLOW::Task * task = new ONEFLOW::Task();
    task->taskId = 123;
    task->taskName = "test_task";

    command.AddTask( task );

    ASSERT_EQ( command.GetTaskList()->size(), 1 );

    ONEFLOW::Task * storedTask =
        ( * command.GetTaskList() )[ 0 ];

    ASSERT_NE( storedTask, nullptr );
    EXPECT_EQ( storedTask, task );
    EXPECT_EQ( storedTask->taskId, 123 );
    EXPECT_EQ( storedTask->taskName, "test_task" );
}

TEST_F( CommandTest, AddUniqueTaskTransfersOwnership )
{
    ONEFLOW::SimpleCmd command;

    auto task = std::make_unique< ONEFLOW::Task >();

    ONEFLOW::Task * rawTask = task.get();

    command.AddTask( std::move( task ) );

    EXPECT_EQ( task, nullptr );

    ASSERT_EQ( command.GetTaskList()->size(), 1 );
    EXPECT_EQ( ( * command.GetTaskList() )[ 0 ], rawTask );
}

TEST_F( CommandTest, AddNullRawTaskDoesNothing )
{
    ONEFLOW::SimpleCmd command;

    command.AddTask(
        static_cast< ONEFLOW::Task * >( nullptr )
    );

    EXPECT_TRUE( command.GetTaskList()->empty() );
}

TEST_F( CommandTest, AddNullUniqueTaskDoesNothing )
{
    ONEFLOW::SimpleCmd command;

    std::unique_ptr< ONEFLOW::Task > task;

    command.AddTask( std::move( task ) );

    EXPECT_TRUE( command.GetTaskList()->empty() );
}

TEST_F( CommandTest, SimpleCmdExecutesStoredTask )
{
    ONEFLOW::SimpleCmd command;

    command.AddTask(
        std::make_unique< CountingTask >()
    );

    EXPECT_EQ( g_runCount, 0 );

    command.Execute();

    EXPECT_EQ( g_runCount, 1 );
}

TEST_F( CommandTest, SimpleCmdExecutesTasksInOrder )
{
    class OrderedTask : public ONEFLOW::Task
    {
    public:
        OrderedTask( int id, ONEFLOW::HXVector< int > * order )
            : id_( id ),
            order_( order )
        {
        }

        void Run() override
        {
            order_->push_back( id_ );
        }

    private:
        int id_;
        ONEFLOW::HXVector< int > * order_;
    };

    ONEFLOW::HXVector< int > order;

    ONEFLOW::SimpleCmd command;

    command.AddTask(
        std::make_unique< OrderedTask >( 1, & order )
    );

    command.AddTask(
        std::make_unique< OrderedTask >( 2, & order )
    );

    command.AddTask(
        std::make_unique< OrderedTask >( 3, & order )
    );

    command.Execute();

    ASSERT_EQ( order.size(), 3 );

    EXPECT_EQ( order[ 0 ], 1 );
    EXPECT_EQ( order[ 1 ], 2 );
    EXPECT_EQ( order[ 2 ], 3 );
}

TEST_F( CommandTest, ExecuteSetsCurrentTask )
{
    ONEFLOW::SimpleCmd command;

    auto task = std::make_unique< CountingTask >();
    ONEFLOW::Task * rawTask = task.get();

    command.AddTask( std::move( task ) );

    command.Execute();

    EXPECT_EQ( ONEFLOW::TaskState::task, rawTask );
}

TEST_F( CommandTest, CommandDestroysOwnedTask )
{
    EXPECT_EQ( g_destroyCount, 0 );

    {
        ONEFLOW::SimpleCmd command;

        command.AddTask(
            std::make_unique< CountingTask >()
        );

        EXPECT_EQ( g_destroyCount, 0 );
    }

    EXPECT_EQ( g_destroyCount, 1 );
}

TEST_F( CommandTest, RawTaskOwnershipIsTransferredToCommand )
{
    EXPECT_EQ( g_destroyCount, 0 );

    {
        ONEFLOW::SimpleCmd command;

        command.AddTask(
            new CountingTask()
        );

        EXPECT_EQ( g_destroyCount, 0 );
    }

    EXPECT_EQ( g_destroyCount, 1 );
}

static_assert(
    ! std::is_copy_constructible< ONEFLOW::SimpleCmd >::value,
    "SimpleCmd must not be copy constructible"
    );

static_assert(
    ! std::is_copy_assignable< ONEFLOW::SimpleCmd >::value,
    "SimpleCmd must not be copy assignable"
    );

TEST( CommandQueueTest, InitCreatesEmptyQueue )
{
    ONEFLOW::CMD::Free();
    ONEFLOW::CMD::Init();

    ASSERT_NE( ONEFLOW::CMD::cmdList, nullptr );
    EXPECT_TRUE( ONEFLOW::CMD::cmdList->empty() );

    ONEFLOW::CMD::Free();
}

TEST( CommandQueueTest, AddCmdTransfersOwnership )
{
    g_destroyCount = 0;

    ONEFLOW::CMD::Free();

    {
        auto command = std::make_unique< ONEFLOW::SimpleCmd >();

        command->AddTask(
            std::make_unique< QueueTestTask >( 1 )
        );

        ONEFLOW::CMD::AddCmd( std::move( command ) );

        ASSERT_NE( ONEFLOW::CMD::cmdList, nullptr );
        EXPECT_EQ( ONEFLOW::CMD::cmdList->size(), 1u );

        /*
        * The Task is owned by Command, and Command is owned by CMD.
        */
        EXPECT_EQ( g_destroyCount, 0 );
    }

    /*
    * The original unique_ptr no longer owns the Command.
    */
    EXPECT_EQ( g_destroyCount, 0 );

    ONEFLOW::CMD::Clear();

    /*
    * Clear() destroys Command and therefore its Task.
    */
    EXPECT_EQ( g_destroyCount, 1 );

    ONEFLOW::CMD::Free();
}

TEST( CommandQueueTest, CommandsExecuteInInsertionOrder )
{
    g_executionLog.clear();

    ONEFLOW::CMD::Free();

    auto command1 = std::make_unique< QueueTestCommand >( 1 );
    auto command2 = std::make_unique< QueueTestCommand >( 2 );
    auto command3 = std::make_unique< QueueTestCommand >( 3 );

    ONEFLOW::CMD::AddCmd( std::move( command1 ) );
    ONEFLOW::CMD::AddCmd( std::move( command2 ) );
    ONEFLOW::CMD::AddCmd( std::move( command3 ) );

    ASSERT_EQ( ONEFLOW::CMD::cmdList->size(), 3u );

    ONEFLOW::CMD::ExecuteCmd();

    ASSERT_EQ( g_executionLog.size(), 3u );

    EXPECT_EQ( g_executionLog[ 0 ], 1 );
    EXPECT_EQ( g_executionLog[ 1 ], 2 );
    EXPECT_EQ( g_executionLog[ 2 ], 3 );

    EXPECT_TRUE( ONEFLOW::CMD::cmdList->empty() );

    ONEFLOW::CMD::Free();
}

TEST( CommandQueueTest, EachDispatchCycleCanHaveDifferentCommandList )
{
    g_executionLog.clear();

    /*
    * First iteration:
    *
    * A -> B -> C
    */
    ONEFLOW::CMD::Free();

    ONEFLOW::CMD::AddCmd(
        std::make_unique< QueueTestCommand >( 1 )
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< QueueTestCommand >( 2 )
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< QueueTestCommand >( 3 )
    );

    ONEFLOW::CMD::ExecuteCmd();

    ASSERT_EQ( g_executionLog.size(), 3u );
    EXPECT_EQ( g_executionLog[ 0 ], 1 );
    EXPECT_EQ( g_executionLog[ 1 ], 2 );
    EXPECT_EQ( g_executionLog[ 2 ], 3 );

    /*
    * Second iteration:
    *
    * A -> D
    *
    * The command list is rebuilt dynamically.
    */
    g_executionLog.clear();

    ONEFLOW::CMD::AddCmd(
        std::make_unique< QueueTestCommand >( 1 )
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< QueueTestCommand >( 4 )
    );

    ONEFLOW::CMD::ExecuteCmd();

    ASSERT_EQ( g_executionLog.size(), 2u );
    EXPECT_EQ( g_executionLog[ 0 ], 1 );
    EXPECT_EQ( g_executionLog[ 1 ], 4 );

    /*
    * Third iteration:
    *
    * Only E.
    */
    g_executionLog.clear();

    ONEFLOW::CMD::AddCmd(
        std::make_unique< QueueTestCommand >( 5 )
    );

    ONEFLOW::CMD::ExecuteCmd();

    ASSERT_EQ( g_executionLog.size(), 1u );
    EXPECT_EQ( g_executionLog[ 0 ], 5 );

    ONEFLOW::CMD::Free();
}

TEST( CommandQueueTest, ClearDestroysQueuedCommands )
{
    g_destroyCount = 0;

    ONEFLOW::CMD::Free();

    {
        auto command = std::make_unique< ONEFLOW::SimpleCmd >();

        command->AddTask(
            std::make_unique< QueueTestTask >( 10 )
        );

        ONEFLOW::CMD::AddCmd( std::move( command ) );

        EXPECT_EQ( g_destroyCount, 0 );
        ASSERT_EQ( ONEFLOW::CMD::cmdList->size(), 1u );

        ONEFLOW::CMD::Clear();

        EXPECT_EQ( g_destroyCount, 1 );
        EXPECT_TRUE( ONEFLOW::CMD::cmdList->empty() );
    }

    /*
    * Nothing should be destroyed here because ownership was already
    * released by CMD::Clear().
    */
    EXPECT_EQ( g_destroyCount, 1 );

    ONEFLOW::CMD::Free();
}

TEST( CommandQueueTest, FreeDestroysQueuedCommands )
{
    g_destroyCount = 0;

    ONEFLOW::CMD::Free();

    ONEFLOW::CMD::AddCmd(
        std::make_unique< ONEFLOW::SimpleCmd >()
    );

    EXPECT_EQ( ONEFLOW::CMD::cmdList->size(), 1u );

    ONEFLOW::CMD::Free();

    EXPECT_EQ( ONEFLOW::CMD::cmdList, nullptr );
    EXPECT_EQ( g_destroyCount, 0 );
}

TEST( CommandQueueTest, FreeReleasesCommandsAndTasks )
{
    g_destroyCount = 0;

    ONEFLOW::CMD::Free();

    auto command = std::make_unique< ONEFLOW::SimpleCmd >();

    command->AddTask(
        std::make_unique< QueueTestTask >( 20 )
    );

    ONEFLOW::CMD::AddCmd( std::move( command ) );

    EXPECT_EQ( g_destroyCount, 0 );

    ONEFLOW::CMD::Free();

    EXPECT_EQ( g_destroyCount, 1 );
    EXPECT_EQ( ONEFLOW::CMD::cmdList, nullptr );
}

TEST( CommandQueueTest, RunCmdDoesNotTakeOwnership )
{
    g_executionLog.clear();
    g_destroyCount = 0;

    ONEFLOW::CMD::Free();

    {
        auto command = std::make_unique< ONEFLOW::SimpleCmd >();

        command->AddTask(
            std::make_unique< QueueTestTask >( 30 )
        );

        ONEFLOW::CMD::RunCmd( command.get() );

        ASSERT_EQ( g_executionLog.size(), 1u );
        EXPECT_EQ( g_executionLog[ 0 ], 30 );

        /*
        * RunCmd() only executes the Command.
        * It does not take ownership.
        */
        EXPECT_EQ( g_destroyCount, 0 );
    }

    EXPECT_EQ( g_destroyCount, 1 );

    ONEFLOW::CMD::Free();
}