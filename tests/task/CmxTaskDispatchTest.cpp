// CmxTaskDispatchTest.cpp
#include <gtest/gtest.h>
#include <vector>
#include <string>
#include "CmxTask.h"
#include "Command.h"
#include "Task.h"
#include "TaskState.h"
#include "SolverState.h"
#include "Message.h"
#include "Register.h"
#include "HXClone.h"

// This test covers the "queued execution" path inside GenerateCmdList:
// when RegisterFactory has no MESG_FUNC-registered class for the given
// message name, AddCmdToList() is used instead, which creates a Task via
// the TASK_FUNC registry slot, binds a file via FILE_FUNC, wraps it in a
// SimpleCmd, and defers execution to CMD::ExecuteCmd().
//
// We deliberately bypass SolverMap by calling SolverState::SetSolverType()
// directly, since SolverMap's source is not yet available. A second test
// suite will cover MultiSolverMultiGridTask's outer loop once SolverMap
// is provided.

namespace
{
    std::vector<std::string> g_callLog;

    // Stub Task-registration class (TASK_FUNC slot): produces a Task whose
    // action() records into g_callLog when Task::Run() invokes it.
    void StubTaskAction()
    {
        g_callLog.push_back( "StubTaskAction" );
    }

    // Minimal HXClone used to occupy the TASK_FUNC registry slot. Its
    // Solve() creates a SimpleTask-like Task and wires up the action
    // function pointer, mirroring what a real *_TaskImp registration
    // class would do (see TaskImp.cpp's DEFINE_DATA_CLASS pattern).
    class StubTaskCreator : public ONEFLOW::HXClone
    {
    public:
        ONEFLOW::HXClone * Clone() const override
        {
            return new StubTaskCreator( *this );
        }
        void Solve() override
        {
            ONEFLOW::Task * task = new ONEFLOW::Task();
            task->action = & StubTaskAction;
            ONEFLOW::TaskState::task = task;
        }
    };
}

class CmxTaskDispatchTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        g_callLog.clear();
        ONEFLOW::MessageMap::Init();
        ONEFLOW::RegisterFactory::Init();
        ONEFLOW::CMD::Init();
        ONEFLOW::SolverState::SetSolverType( /*solverType=*/1 );
    }

    void TearDown() override
    {
        ONEFLOW::MessageMap::Free();
        ONEFLOW::RegisterFactory::FreeMRegister();
        ONEFLOW::CMD::Clear();
        ONEFLOW::CMD::Free();
        g_callLog.clear();
    }
};

