#include <gtest/gtest.h>

#include "Command.h"
#include "Task.h"
#include "TaskState.h"

#include <memory>
#include <type_traits>
#include <vector>

namespace
{
    int g_taskRunCount = 0;
    int g_taskDestroyCount = 0;

    class CountingTask : public ONEFLOW::Task
    {
    public:
        ~CountingTask() override
        {
            ++ g_taskDestroyCount;
        }

        void Run() override
        {
            ++ g_taskRunCount;
        }
    };

    class OrderedTask : public ONEFLOW::Task
    {
    public:
        explicit OrderedTask( int id )
            : id_( id )
        {
        }

        void Run() override
        {
            order_.push_back( id_ );
        }

        static void ClearOrder()
        {
            order_.clear();
        }

        static const std::vector< int > & GetOrder()
        {
            return order_;
        }

    private:
        int id_;

        static std::vector< int > order_;
    };

    std::vector< int > OrderedTask::order_;

    class TaskStateCheckingTask : public ONEFLOW::Task
    {
    public:
        void Run() override
        {
            observedTask_ = ONEFLOW::TaskState::task;
        }

        static ONEFLOW::Task * GetObservedTask()
        {
            return observedTask_;
        }

        static void ClearObservedTask()
        {
            observedTask_ = nullptr;
        }

    private:
        static ONEFLOW::Task * observedTask_;
    };

    ONEFLOW::Task * TaskStateCheckingTask::observedTask_ = nullptr;
}

class CommandTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        g_taskRunCount = 0;
        g_taskDestroyCount = 0;

        OrderedTask::ClearOrder();
        TaskStateCheckingTask::ClearObservedTask();

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

    EXPECT_EQ( command.GetTaskCount(), 0 );
}

TEST_F( CommandTest, AddRawTaskStoresTask )
{
    auto * task = new CountingTask();

    ONEFLOW::SimpleCmd command;
    command.AddTask( task );

    ASSERT_EQ( command.GetTaskCount(), 1 );
    EXPECT_EQ( command.GetTask( 0 ), task );
}

TEST_F( CommandTest, AddUniqueTaskTransfersOwnership )
{
    auto task = std::make_unique< CountingTask >();

    ONEFLOW::SimpleCmd command;
    command.AddTask( std::move( task ) );

    EXPECT_EQ( task, nullptr );

    ASSERT_EQ( command.GetTaskCount(), 1 );
    EXPECT_NE( command.GetTask( 0 ), nullptr );
}

TEST_F( CommandTest, AddNullRawTaskDoesNothing )
{
    ONEFLOW::SimpleCmd command;

    command.AddTask(
        static_cast< ONEFLOW::Task * >( nullptr )
    );

    EXPECT_EQ( command.GetTaskCount(), 0 );
}

TEST_F( CommandTest, AddNullUniqueTaskDoesNothing )
{
    ONEFLOW::SimpleCmd command;

    std::unique_ptr< ONEFLOW::Task > task;

    command.AddTask( std::move( task ) );

    EXPECT_EQ( command.GetTaskCount(), 0 );
}

TEST_F( CommandTest, SimpleCmdExecutesStoredTask )
{
    ONEFLOW::SimpleCmd command;

    command.AddTask(
        std::make_unique< CountingTask >()
    );

    command.Execute();

    EXPECT_EQ( g_taskRunCount, 1 );
}

TEST_F( CommandTest, SimpleCmdExecutesTasksInOrder )
{
    ONEFLOW::SimpleCmd command;

    command.AddTask(
        std::make_unique< OrderedTask >( 1 )
    );

    command.AddTask(
        std::make_unique< OrderedTask >( 2 )
    );

    command.AddTask(
        std::make_unique< OrderedTask >( 3 )
    );

    command.Execute();

    const auto & order = OrderedTask::GetOrder();

    ASSERT_EQ( order.size(), 3 );

    EXPECT_EQ( order[ 0 ], 1 );
    EXPECT_EQ( order[ 1 ], 2 );
    EXPECT_EQ( order[ 2 ], 3 );
}

TEST_F( CommandTest, ExecuteSetsCurrentTask )
{
    auto task = std::make_unique< TaskStateCheckingTask >();

    ONEFLOW::Task * rawTask = task.get();

    ONEFLOW::SimpleCmd command;
    command.AddTask( std::move( task ) );

    command.Execute();

    EXPECT_EQ(
        TaskStateCheckingTask::GetObservedTask(),
        rawTask
    );
}

TEST_F( CommandTest, CommandDestroysOwnedTask )
{
    {
        ONEFLOW::SimpleCmd command;

        command.AddTask(
            std::make_unique< CountingTask >()
        );

        EXPECT_EQ( g_taskDestroyCount, 0 );
    }

    EXPECT_EQ( g_taskDestroyCount, 1 );
}

TEST_F( CommandTest, RawTaskOwnershipIsTransferredToCommand )
{
    CountingTask * task = new CountingTask();

    {
        ONEFLOW::SimpleCmd command;

        command.AddTask( task );

        EXPECT_EQ( g_taskDestroyCount, 0 );
    }

    EXPECT_EQ( g_taskDestroyCount, 1 );
}

TEST( CommandTypeTest, SimpleCmdIsNotCopyable )
{
    EXPECT_FALSE(
        std::is_copy_constructible< ONEFLOW::SimpleCmd >::value
    );

    EXPECT_FALSE(
        std::is_copy_assignable< ONEFLOW::SimpleCmd >::value
    );
}

TEST_F( CommandTest, ConstAccessorsProvideReadOnlyView )
{
    ONEFLOW::SimpleCmd command;

    command.AddTask(
        std::make_unique< CountingTask >()
    );

    const ONEFLOW::Command * constCommand = & command;

    ASSERT_EQ( constCommand->GetTaskCount(), 1 );
    EXPECT_NE( constCommand->GetTask( 0 ), nullptr );
}

TEST_F( CommandTest, TaskAccessorsMatchOwnedTasks )
{
    ONEFLOW::SimpleCmd command;

    auto task1 = std::make_unique< CountingTask >();
    auto task2 = std::make_unique< CountingTask >();

    ONEFLOW::Task * rawTask1 = task1.get();
    ONEFLOW::Task * rawTask2 = task2.get();

    command.AddTask( std::move( task1 ) );
    command.AddTask( std::move( task2 ) );

    ASSERT_EQ( command.GetTaskCount(), 2 );
    EXPECT_EQ( command.GetTask( 0 ), rawTask1 );
    EXPECT_EQ( command.GetTask( 1 ), rawTask2 );
}