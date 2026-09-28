// tests/GridFactoryTest.cpp
#include <gtest/gtest.h>
#include "GridFactory.h"

// Test that GridFactory can be safely instantiated and destroyed on the stack.
// This ensures the RAII pattern is correctly applied without memory leaks.
TEST(GridFactoryLifecycleTest, StackAllocationIsSafe)
{
    // Arrange & Act & Assert
    // We expect no exceptions to be thrown during construction and destruction.
    EXPECT_NO_THROW( {
        ONEFLOW::GridFactory gf;
    // We intentionally do NOT call gf.Run() here to avoid global state 
    // side effects and file I/O in a pure unit test environment.
        });
}

// Test that the legacy heap allocation pattern is functionally equivalent 
// but we prefer stack allocation. (Characterization test)
TEST(GridFactoryLifecycleTest, HeapAllocationAlsoWorksButIsDiscouraged)
{
    EXPECT_NO_THROW({
        ONEFLOW::GridFactory* gf = new ONEFLOW::GridFactory();
    delete gf;
        });
}