// CategoryTest.cpp
#include <gtest/gtest.h>
#include "Category.h"
// Category::GetCategory currently dereferences map::end() (or a null map
// pointer) when the queried solverType was never registered / Init() was
// never called. These tests characterize the buggy behavior first, then
// pin down the fixed behavior after the patch below.

class CategoryTest : public ::testing::Test
{
protected:
    void TearDown() override
    {
        // Always leave global state clean for the next test, regardless
        // of whether this test called Init() or not.
        ONEFLOW::Category::Free();
    }
};

TEST_F( CategoryTest, GetCategoryOnUnregisteredSolverTypeReturnsSentinel )
{
    ONEFLOW::Category::Init();
    ONEFLOW::Category::AddCategory( /*solverType=*/1, /*category=*/100 );

    // solverType 5 was never registered via AddCategory
    int result = ONEFLOW::Category::GetCategory( 5 );

    EXPECT_EQ( result, -1 );
}

TEST_F( CategoryTest, GetCategoryOnRegisteredSolverTypeReturnsStoredValue )
{
    ONEFLOW::Category::Init();
    ONEFLOW::Category::AddCategory( 2, 200 );

    EXPECT_EQ( ONEFLOW::Category::GetCategory( 2 ), 200 );
}

TEST_F( CategoryTest, GetCategoryBeforeInitDoesNotCrash )
{
    // Category::Init() was never called in this test -> Category::data is null.
    int result = ONEFLOW::Category::GetCategory( 1 );

    EXPECT_EQ( result, -1 );
}

TEST_F( CategoryTest, AddCategoryDoesNotOverwriteExistingEntry )
{
    // AddCategory's existing behavior: only inserts if the key is absent.
    // Pinning this down so the fix below doesn't accidentally change it.
    ONEFLOW::Category::Init();
    ONEFLOW::Category::AddCategory( 3, 300 );
    ONEFLOW::Category::AddCategory( 3, 999 ); // should be ignored

    EXPECT_EQ( ONEFLOW::Category::GetCategory( 3 ), 300 );
}
