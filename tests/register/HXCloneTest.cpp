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
        // The registry is process-wide, so tests use unique type names.
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

    // Registering the same type again transfers ownership of the duplicate
    // instance, which Register() discards while preserving the first entry.
    ONEFLOW::HXClone * second = new StubClone();
    ONEFLOW::HXClone * returned = ONEFLOW::HXClone::Register( "HXCloneTest_TypeB", second );

    EXPECT_EQ( returned, first ); // the original registration wins
}

TEST_F( HXCloneTest, SafeCloneOnUnregisteredTypeThrows )
{
    EXPECT_THROW(
        ONEFLOW::HXClone::SafeClone( "HXCloneTest_NeverRegistered" ),
        std::runtime_error
    );
}
