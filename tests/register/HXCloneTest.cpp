// HXCloneTest.cpp
#include <gtest/gtest.h>
#include "HXClone.h"

namespace
{
    // Minimal concrete HXClone for testing SafeClone/Register without
    // pulling in any real solver/task class.
    class StubClone : public ONEFLOW::HXClone
    {
    public:
        ONEFLOW::HXClone * Clone() const override
        {
            return new StubClone( *this );
        }
    };
}

class HXCloneTest : public ::testing::Test
{
protected:
    void TearDown() override
    {
        // HXClone::classMap has no public Free()/Clear() - see note below.
        // For now, tests must use unique type names to avoid cross-test
        // pollution, since we cannot safely reset classMap here.
    }
};

TEST_F( HXCloneTest, RegisterThenSafeCloneReturnsANewInstance )
{
    ONEFLOW::HXClone::Register( "HXCloneTest_TypeA", new StubClone() );

    ONEFLOW::HXClone * cloned = ONEFLOW::HXClone::SafeClone( "HXCloneTest_TypeA" );

    ASSERT_NE( cloned, nullptr );
    delete cloned; // SafeClone returns a new heap instance; caller owns it
}

TEST_F( HXCloneTest, RegisterIsIdempotentAndDeletesTheDuplicateArgument )
{
    ONEFLOW::HXClone * first = new StubClone();
    ONEFLOW::HXClone::Register( "HXCloneTest_TypeB", first );

    // Registering the same type name again passes ownership of a new
    // instance in, which Register() deletes internally (see original
    // behavior: `delete clone; return iter->second;`). We must not
    // touch `second` after this call except through the registry.
    ONEFLOW::HXClone * second = new StubClone();
    ONEFLOW::HXClone * returned = ONEFLOW::HXClone::Register( "HXCloneTest_TypeB", second );

    EXPECT_EQ( returned, first ); // the original registration wins
}

// NOTE: no test for "unregistered type" (classMap null or type not
// found) because Fatal()'s actual control-flow behavior is unknown to
// us - if it calls exit()/abort(), a test exercising that path would
// kill the whole test binary. Please confirm Fatal's implementation
// before adding coverage for that branch.

TEST_F( HXCloneTest, SafeCloneOnUnregisteredTypeThrows )
{
    EXPECT_THROW(
        ONEFLOW::HXClone::SafeClone( "HXCloneTest_NeverRegistered" ),
        std::runtime_error
    );
}