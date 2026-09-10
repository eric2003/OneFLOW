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

TEST_F( CmxTaskDispatchTest, DISABLED_MessageWithNoMesgFuncClassIsQueuedAndRunsViaCmd )
{
    const int SOLVER_TYPE = 1;
    const int TASK_FUNC_SLOT = 3; // matches CmxTask.h's TASK_FUNC constant

    // Register the message name so MessageMap::GetMsgId/GetMsgName work.
    ONEFLOW::MessageMap::Register( "TestMessage" );

    // Wire up the TASK_FUNC registry slot for this solverType, so
    // CreateTask() finds a class and produces our stub Task instead of
    // falling back to `new SimpleTask()`.
    ONEFLOW::RegisterFactory::AddMRegister( SOLVER_TYPE );
    ONEFLOW::MRegister * mRegister = ONEFLOW::RegisterFactory::GetMRegister( SOLVER_TYPE );
    ASSERT_NE( mRegister, nullptr );

    ONEFLOW::StringField fileNames( 5 ); // one slot per msgType (COMM/RECV/MESG/TASK/FILE)
    mRegister->SetSolverFileNames( fileNames );
    // NOTE: RegisterAll()/AllocateData() would normally build these from
    // file names via TextFileParser. Since we want to inject a stub class
    // directly without real files, we register into the TASK_FUNC slot's
    // HXRegister by hand instead of calling RegisterAll().
    ONEFLOW::HXClone::Register( "StubTaskCreator", new StubTaskCreator() );

    // NOTE: this next step exposes a real gap - MRegister::GetRegister(index)
    // returns nullptr until AllocateData() has run, and there is currently
    // no public API to inject a class into a specific HXRegister slot
    // without going through a real file. This test is left as a sketch
    // to surface that gap rather than force a workaround.
    ONEFLOW::HXRegister * taskRegister = mRegister->GetRegister( TASK_FUNC_SLOT );
    ASSERT_NE( taskRegister, nullptr ) << "Requires MRegister to expose a "
        "way to populate a specific HXRegister slot without real files - "
        "see accompanying note.";
}