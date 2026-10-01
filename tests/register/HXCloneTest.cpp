// HXCloneTest.cpp
#include <gtest/gtest.h>
#include <memory>
#include "HXClone.h"

namespace
{
    // Minimal concrete HXClone for testing SafeClone/Register without
    // pulling in any real solver/task class.
    class StubClone : public ONEFLOW::HXClone
    {
    public:
        std::unique_ptr< ONEFLOW::HXClone > Clone() const override
        {
            return std::make_unique< StubClone >( *this );
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
    ONEFLOW::HXClone::Register( "HXCloneTest_TypeA", std::make_unique< StubClone >() );

    std::unique_ptr< ONEFLOW::HXClone > cloned =
        ONEFLOW::HXClone::SafeCloneUnique( "HXCloneTest_TypeA" );

    ASSERT_NE( cloned, nullptr );
}

TEST_F( HXCloneTest, RegisterIsIdempotentAndDeletesTheDuplicateArgument )
{
    ONEFLOW::HXClone * first = ONEFLOW::HXClone::Register(
        "HXCloneTest_TypeB", std::make_unique< StubClone >() );

    // Registering the same type again transfers ownership of the duplicate
    // instance, which Register() discards while preserving the first entry.
    ONEFLOW::HXClone * returned = ONEFLOW::HXClone::Register(
        "HXCloneTest_TypeB", std::make_unique< StubClone >() );

    EXPECT_EQ( returned, first ); // the original registration wins
}

TEST_F( HXCloneTest, SafeCloneOnUnregisteredTypeThrows )
{
    EXPECT_THROW(
        ONEFLOW::HXClone::SafeCloneUnique( "HXCloneTest_NeverRegistered" ),
        std::runtime_error
    );
}
