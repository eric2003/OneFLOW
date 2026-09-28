// TaskRegisterTest.cpp
#include <gtest/gtest.h>
#include <vector>
#include <string>
#include "TaskRegister.h"

// TaskRegister is entirely static (no instantiable "Imp" class like
// ActionMapImp), so tests must reset global state via Free()/Register()
// between cases. This is only safe now that Free() nulls its pointers
// after delete (see the fix in TaskRegister.cpp) -- otherwise each
// test's cleanup would corrupt the next test's state.
//
// CAUTION: if this test binary also links TaskImp.cpp (which contains
// `REGISTER_TASK(RegisterComTask)`, registered via static initialization
// before main() runs), calling Free() here will permanently discard that
// real registration for the remainder of the process. If other code in
// the same binary depends on it, consider isolating these tests into
// their own executable.

namespace
{
    std::vector<std::string> g_callLog;

    void StubTaskA() { g_callLog.push_back( "A" ); }
    void StubTaskB() { g_callLog.push_back( "B" ); }
}

class TaskRegisterTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Start every test from a known-empty registry.
        ONEFLOW::TaskRegister::Free();
        g_callLog.clear();
    }

    void TearDown() override
    {
        ONEFLOW::TaskRegister::Free();
    }
};

TEST_F( TaskRegisterTest, RunOnEmptyRegistryDoesNotCrash )
{
    // No Register() call has happened -> taskList is null.
    // Run() previously dereferenced this null pointer.
    EXPECT_NO_THROW( ONEFLOW::TaskRegister::Run() );
    EXPECT_TRUE( g_callLog.empty() );
}

TEST_F( TaskRegisterTest, RunInvokesRegisteredFunctionsInRegistrationOrder )
{
    ONEFLOW::TaskRegister::Register( &StubTaskA, "StubTaskA" );
    ONEFLOW::TaskRegister::Register( &StubTaskB, "StubTaskB" );

    ONEFLOW::TaskRegister::Run();

    ASSERT_EQ( g_callLog.size(), 2u );
    EXPECT_EQ( g_callLog[0], "A" );
    EXPECT_EQ( g_callLog[1], "B" );
}

TEST_F( TaskRegisterTest, RegisterStoresTaskNamesInOrder )
{
    ONEFLOW::TaskRegister::Register( &StubTaskA, "StubTaskA" );
    ONEFLOW::TaskRegister::Register( &StubTaskB, "StubTaskB" );

    ASSERT_NE( ONEFLOW::TaskRegister::taskNameList, nullptr );
    ASSERT_EQ( ONEFLOW::TaskRegister::taskNameList->size(), 2u );
    EXPECT_EQ( ( * ONEFLOW::TaskRegister::taskNameList )[0], "StubTaskA" );
    EXPECT_EQ( ( * ONEFLOW::TaskRegister::taskNameList )[1], "StubTaskB" );
}

TEST_F( TaskRegisterTest, FreeThenRegisterAgainWorksCorrectly )
{
    // Regression test for the dangling-pointer bug: after Free(), the
    // pointers must be null so that Register() allocates fresh lists
    // instead of pushing into already-freed memory.
    ONEFLOW::TaskRegister::Register( &StubTaskA, "StubTaskA" );
    ONEFLOW::TaskRegister::Free();

    ONEFLOW::TaskRegister::Register( &StubTaskB, "StubTaskB" );
    ONEFLOW::TaskRegister::Run();

    ASSERT_EQ( g_callLog.size(), 1u );
    EXPECT_EQ( g_callLog[0], "B" );
}

TEST_F( TaskRegisterTest, FreeIsIdempotentWhenCalledTwice )
{
    // Regression test for the double-free hazard: calling Free() twice
    // in a row (e.g. once here, once more via Tmp_Free_TaskRegister's
    // destructor at process exit) must not crash.
    ONEFLOW::TaskRegister::Register( &StubTaskA, "StubTaskA" );
    ONEFLOW::TaskRegister::Free();

    EXPECT_NO_THROW( ONEFLOW::TaskRegister::Free() );
}