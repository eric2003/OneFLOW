#include <gtest/gtest.h>

#include "Command.h"
#include "Task.h"
#include "TaskState.h"

#include <memory>
#include <vector>

namespace
{
    int g_cmdRunCount = 0;
    int g_cmdTaskDestroyCount = 0;
    int g_cmdDestroyCount = 0;

    class CmdCountingTask : public ONEFLOW::Task
    {
    public:
        ~CmdCountingTask() override
        {
            ++ g_cmdTaskDestroyCount;
        }

        void Run() override
        {
            ++ g_cmdRunCount;
        }
    };

    class CmdOrderedTask : public ONEFLOW::Task
    {
    public:
        explicit CmdOrderedTask( int id )
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

    std::vector< int > CmdOrderedTask::order_;

    class CountingCommand : public ONEFLOW::Command
    {
    public:
        ~CountingCommand() override
        {
            ++ g_cmdDestroyCount;
        }

        void Execute() override
        {
        }
    };

    class ExecutingCommand : public ONEFLOW::Command
    {
    public:
        explicit ExecutingCommand( int id )
            : id_( id )
        {
        }

        void Execute() override
        {
            executionOrder_.push_back( id_ );
        }

        static void ClearExecutionOrder()
        {
            executionOrder_.clear();
        }

        static const std::vector< int > & GetExecutionOrder()
        {
            return executionOrder_;
        }

    private:
        int id_;

        static std::vector< int > executionOrder_;
    };

    std::vector< int > ExecutingCommand::executionOrder_;
}

class CMDTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        ONEFLOW::CMD::Free();

        ONEFLOW::TaskState::task = nullptr;

        g_cmdRunCount = 0;
        g_cmdTaskDestroyCount = 0;
        g_cmdDestroyCount = 0;

        CmdOrderedTask::ClearOrder();
        ExecutingCommand::ClearExecutionOrder();
    }

    void TearDown() override
    {
        ONEFLOW::CMD::Free();

        ONEFLOW::TaskState::task = nullptr;
    }
};

TEST_F( CMDTest, InitCreatesEmptyCommandList )
{
    ONEFLOW::CMD::Init();

    const auto * commandList =
        ONEFLOW::CMD::GetCmdList();

    ASSERT_NE( commandList, nullptr );
    EXPECT_TRUE( commandList->empty() );
}

TEST_F( CMDTest, AddCmdTransfersOwnership )
{
    CountingCommand * command =
        new CountingCommand();

    ONEFLOW::CMD::AddCmd( command );

    const auto * commandList =
        ONEFLOW::CMD::GetCmdList();

    ASSERT_NE( commandList, nullptr );
    EXPECT_EQ( commandList->size(), 1 );


    EXPECT_EQ(
        ( * commandList )[ 0 ],
        command
    );

    EXPECT_EQ( g_cmdDestroyCount, 0 );

    ONEFLOW::CMD::Clear();

    EXPECT_EQ( g_cmdDestroyCount, 1 );
}

TEST_F( CMDTest, AddUniqueCmdTransfersOwnership )
{
    auto command =
        std::make_unique< CountingCommand >();

    ONEFLOW::CMD::AddCmd( std::move( command ) );

    EXPECT_EQ( command, nullptr );

    const auto * commandList =
        ONEFLOW::CMD::GetCmdList();


    ASSERT_NE( commandList, nullptr );
    ASSERT_EQ( commandList->size(), 1 );

    EXPECT_EQ( g_cmdDestroyCount, 0 );

    ONEFLOW::CMD::Clear();

    EXPECT_EQ( g_cmdDestroyCount, 1 );
}

TEST_F( CMDTest, AddNullRawCommandDoesNothing )
{
    const auto * commandList =
        ONEFLOW::CMD::GetCmdList();

    EXPECT_EQ( commandList, nullptr );

    ONEFLOW::CMD::AddCmd(
        static_cast< ONEFLOW::Command * >( nullptr )
    );

    // Adding a null command should not initialize the queue.
    EXPECT_EQ( ONEFLOW::CMD::GetCmdList(), nullptr );
}

TEST_F( CMDTest, AddNullUniqueCommandDoesNothing )
{
    std::unique_ptr< ONEFLOW::Command > command;

    EXPECT_EQ( ONEFLOW::CMD::GetCmdList(), nullptr );

    ONEFLOW::CMD::AddCmd( std::move( command ) );

    // Adding a null command should not initialize the queue.
    EXPECT_EQ( ONEFLOW::CMD::GetCmdList(), nullptr );
}

TEST_F( CMDTest, CommandsExecuteInInsertionOrder )
{
    ONEFLOW::CMD::AddCmd(
        std::make_unique< ExecutingCommand >( 1 )
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< ExecutingCommand >( 2 )
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< ExecutingCommand >( 3 )
    );

    ONEFLOW::CMD::ExecuteCmd();

    const auto & order =
        ExecutingCommand::GetExecutionOrder();

    ASSERT_EQ( order.size(), 3 );

    EXPECT_EQ( order[ 0 ], 1 );
    EXPECT_EQ( order[ 1 ], 2 );
    EXPECT_EQ( order[ 2 ], 3 );
}

TEST_F( CMDTest, ExecuteCmdRunsTask )
{
    auto command =
        std::make_unique< ONEFLOW::SimpleCmd >();

    command->AddTask(
        std::make_unique< CmdCountingTask >()
    );

    ONEFLOW::CMD::AddCmd( std::move( command ) );

    ONEFLOW::CMD::ExecuteCmd();

    EXPECT_EQ( g_cmdRunCount, 1 );
}

TEST_F( CMDTest, ExecuteCmdDestroysOwnedTask )
{
    auto command =
        std::make_unique< ONEFLOW::SimpleCmd >();

    command->AddTask(
        std::make_unique< CmdCountingTask >()
    );

    ONEFLOW::CMD::AddCmd( std::move( command ) );

    EXPECT_EQ( g_cmdTaskDestroyCount, 0 );

    ONEFLOW::CMD::ExecuteCmd();

    EXPECT_EQ( g_cmdTaskDestroyCount, 1 );
}

TEST_F( CMDTest, ExecuteCmdClearsCommandList )
{
    ONEFLOW::CMD::AddCmd(
        std::make_unique< ONEFLOW::SimpleCmd >()
    );

    const auto * commandList =
        ONEFLOW::CMD::GetCmdList();


    ASSERT_NE( commandList, nullptr );
    ASSERT_EQ( commandList->size(), 1 );

    ONEFLOW::CMD::ExecuteCmd();

    ASSERT_NE( commandList, nullptr );
    EXPECT_TRUE( commandList->empty() );
}

TEST_F( CMDTest, ExecuteCmdClearsTaskState )
{
    auto command =
        std::make_unique< ONEFLOW::SimpleCmd >();

    command->AddTask(
        std::make_unique< CmdCountingTask >()
    );

    ONEFLOW::CMD::AddCmd( std::move( command ) );

    ONEFLOW::CMD::ExecuteCmd();

    EXPECT_EQ(
        ONEFLOW::TaskState::task,
        nullptr
    );
}

TEST_F( CMDTest, ClearDestroysQueuedCommands )
{
    ONEFLOW::CMD::AddCmd(
        std::make_unique< CountingCommand >()
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< CountingCommand >()
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< CountingCommand >()
    );

    EXPECT_EQ( g_cmdDestroyCount, 0 );

    ONEFLOW::CMD::Clear();

    EXPECT_EQ( g_cmdDestroyCount, 3 );

    const auto * commandList =
        ONEFLOW::CMD::GetCmdList();

    ASSERT_NE( commandList, nullptr );
    EXPECT_TRUE( commandList->empty() );
}

TEST_F( CMDTest, FreeDestroysQueuedCommands )
{
    ONEFLOW::CMD::AddCmd(
        std::make_unique< CountingCommand >()
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< CountingCommand >()
    );

    EXPECT_EQ( g_cmdDestroyCount, 0 );

    ONEFLOW::CMD::Free();

    EXPECT_EQ( g_cmdDestroyCount, 2 );
    EXPECT_EQ( ONEFLOW::CMD::GetCmdList(), nullptr);
}

TEST_F( CMDTest, FreeReleasesCommandsAndTasks )
{
    auto command =
        std::make_unique< ONEFLOW::SimpleCmd >();

    command->AddTask(
        std::make_unique< CmdCountingTask >()
    );

    ONEFLOW::CMD::AddCmd( std::move( command ) );

    EXPECT_EQ( g_cmdTaskDestroyCount, 0 );

    ONEFLOW::CMD::Free();

    EXPECT_EQ( g_cmdTaskDestroyCount, 1 );
    EXPECT_EQ( ONEFLOW::CMD::GetCmdList(), nullptr);
}

TEST_F( CMDTest, RunCmdDoesNotTakeOwnership )
{
    auto command =
        std::make_unique< CountingCommand >();

    ONEFLOW::Command * rawCommand = command.get();

    ONEFLOW::CMD::RunCmd( rawCommand );

    EXPECT_EQ( g_cmdDestroyCount, 0 );

    command.reset();

    EXPECT_EQ( g_cmdDestroyCount, 1 );
}

TEST_F( CMDTest, FreeCanBeCalledWhenNotInitialized )
{
    EXPECT_NO_THROW(
        ONEFLOW::CMD::Free()
    );

    EXPECT_EQ(
        ONEFLOW::CMD::GetCmdList(),
        nullptr
    );
}

TEST_F( CMDTest, ExecuteCmdBeforeInitDoesNotCrash )
{
    EXPECT_NO_THROW(
        ONEFLOW::CMD::ExecuteCmd()
    );
}

TEST_F( CMDTest, EachDispatchCycleCanHaveDifferentCommandList )
{
    // Dispatch cycle 1: A -> B -> C.
    ONEFLOW::CMD::AddCmd(
        std::make_unique< ExecutingCommand >( 1 )
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< ExecutingCommand >( 2 )
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< ExecutingCommand >( 3 )
    );

    ONEFLOW::CMD::ExecuteCmd();

    {
        const auto & order =
            ExecutingCommand::GetExecutionOrder();

        ASSERT_EQ( order.size(), 3 );

        EXPECT_EQ( order[ 0 ], 1 );
        EXPECT_EQ( order[ 1 ], 2 );
        EXPECT_EQ( order[ 2 ], 3 );
    }

    ExecutingCommand::ClearExecutionOrder();

    // Dispatch cycle 2: A -> D.
    ONEFLOW::CMD::AddCmd(
        std::make_unique< ExecutingCommand >( 1 )
    );

    ONEFLOW::CMD::AddCmd(
        std::make_unique< ExecutingCommand >( 4 )
    );

    ONEFLOW::CMD::ExecuteCmd();

    {
        const auto & order =
            ExecutingCommand::GetExecutionOrder();

        ASSERT_EQ( order.size(), 2 );

        EXPECT_EQ( order[ 0 ], 1 );
        EXPECT_EQ( order[ 1 ], 4 );
    }

    ExecutingCommand::ClearExecutionOrder();

    // Dispatch cycle 3: E only.
    ONEFLOW::CMD::AddCmd(
        std::make_unique< ExecutingCommand >( 5 )
    );

    ONEFLOW::CMD::ExecuteCmd();

    {
        const auto & order =
            ExecutingCommand::GetExecutionOrder();

        ASSERT_EQ( order.size(), 1 );

        EXPECT_EQ( order[ 0 ], 5 );
    }
}

TEST_F( CMDTest, ClearCanBeCalledMultipleTimes )
{
    ONEFLOW::CMD::AddCmd(
        std::make_unique< CountingCommand >()
    );

    ONEFLOW::CMD::Clear();

    EXPECT_EQ( g_cmdDestroyCount, 1 );

    EXPECT_NO_THROW(
        ONEFLOW::CMD::Clear()
    );

    EXPECT_EQ( g_cmdDestroyCount, 1 );
}

TEST_F( CMDTest, GetCmdListProvidesReadOnlyView )
{
    ONEFLOW::CMD::AddCmd(
        std::make_unique< ONEFLOW::SimpleCmd >()
    );

    const ONEFLOW::HXVector< ONEFLOW::Command * > * commandList =
        ONEFLOW::CMD::GetCmdList();

    ASSERT_NE( commandList, nullptr );
    ASSERT_EQ( commandList->size(), 1 );

    EXPECT_NE(
        ( * commandList )[ 0 ],
        nullptr
    );
}