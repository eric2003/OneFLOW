// TaskTest.cpp

#include <gtest/gtest.h>
#include "Task.h"

namespace
{

    void DummyAction()
    {
    }

}

TEST( TaskTest, DefaultStateIsSafe )
{
    ONEFLOW::Task task;

    EXPECT_EQ( task.taskId, -1 );
    EXPECT_TRUE( task.taskName.empty() );

    EXPECT_EQ( task.action, nullptr );
    EXPECT_EQ( task.sendAction, nullptr );
    EXPECT_EQ( task.recvAction, nullptr );

    EXPECT_NE( task.dataBook, nullptr );
    EXPECT_NE( task.fileInfo, nullptr );
}

TEST( TaskTest, OwnsDataBookAndFileInfo )
{
    ONEFLOW::Task task;

    EXPECT_NE( task.dataBook, nullptr );
    EXPECT_NE( task.fileInfo, nullptr );
}

TEST( TaskTest, CallbackCanBeAssigned )
{
    ONEFLOW::Task task;

    task.action = &DummyAction;

    EXPECT_NE( task.action, nullptr );
}