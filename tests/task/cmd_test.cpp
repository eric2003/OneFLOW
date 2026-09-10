//cmd_test.cpp
#include <gtest/gtest.h>

#include "Command.h"
#include "Task.h"
#include "TaskState.h"

#include <memory>

namespace
{

    int g_cmdRunCount = 0;
    int g_cmdTaskDestroyCount = 0;

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

    ASSERT_NE( ONEFLOW::CMD::cmdList, nullptr );
    EXPECT_TRUE( ONEFLOW::CMD::cmdList->empty() );
}

TEST_F( CMDTest, AddCmdStoresCommand )
{
    ONEFLOW::SimpleCmd * command =
        new ONEFLOW::SimpleCmd();

    ONEFLOW::CMD::AddCmd( command );

    ASSERT_NE( ONEFLOW::CMD::cmdList, nullptr );
    ASSERT_EQ( ONEFLOW::CMD::cmdList->size(), 1 );

    EXPECT_EQ(
        ( * ONEFLOW::CMD::cmdList )[ 0 ],
        command
    );

    // ExecuteCmd also deletes the command.
    ONEFLOW::CMD::ExecuteCmd();
}

TEST_F( CMDTest, ExecuteCmdRunsTask )
{
    ONEFLOW::SimpleCmd * command =
        new ONEFLOW::SimpleCmd();

    command->AddTask(
        std::make_unique< CmdCountingTask >()
    );

    ONEFLOW::CMD::AddCmd( command );

    ONEFLOW::CMD::ExecuteCmd();

    EXPECT_EQ( g_cmdRunCount, 1 );
}

TEST_F( CMDTest, ExecuteCmdDestroysOwnedTask )
{
    ONEFLOW::SimpleCmd * command =
        new ONEFLOW::SimpleCmd();

    command->AddTask(
        std::make_unique< CmdCountingTask >()
    );

    ONEFLOW::CMD::AddCmd( command );

    EXPECT_EQ( g_cmdTaskDestroyCount, 0 );

    ONEFLOW::CMD::ExecuteCmd();

    EXPECT_EQ( g_cmdTaskDestroyCount, 1 );
}

TEST_F( CMDTest, ExecuteCmdClearsCommandList )
{
    ONEFLOW::SimpleCmd * command =
        new ONEFLOW::SimpleCmd();

    ONEFLOW::CMD::AddCmd( command );

    ONEFLOW::CMD::ExecuteCmd();

    ASSERT_NE( ONEFLOW::CMD::cmdList, nullptr );
    EXPECT_TRUE( ONEFLOW::CMD::cmdList->empty() );
}

TEST_F( CMDTest, ExecuteCmdClearsDanglingTaskState )
{
    ONEFLOW::SimpleCmd * command =
        new ONEFLOW::SimpleCmd();

    command->AddTask(
        std::make_unique< CmdCountingTask >()
    );

    ONEFLOW::CMD::AddCmd( command );

    ONEFLOW::CMD::ExecuteCmd();

    EXPECT_EQ(
        ONEFLOW::TaskState::task,
        nullptr
    );
}

TEST_F( CMDTest, FreeCanBeCalledWhenNotInitialized )
{
    EXPECT_NO_THROW(
        ONEFLOW::CMD::Free()
    );
}

TEST_F( CMDTest, ExecuteCmdBeforeInitDoesNotCrash )
{
    EXPECT_NO_THROW(
        ONEFLOW::CMD::ExecuteCmd()
    );
}

