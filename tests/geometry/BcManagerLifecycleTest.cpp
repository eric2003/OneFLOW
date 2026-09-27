// tests/BcManagerLifecycleTest.cpp
#include <gtest/gtest.h>
#include "BcRecord.h"

// Verify that BcManager correctly initializes its internal BcRecord unique_ptrs
// and that basic operations (like Update) work without memory leaks.
TEST(BcManagerLifecycleTest, InitializationAndCopyAssignment)
{
    ONEFLOW::BcManager manager;
    
    // 1. Verify internal pointers are correctly initialized (not null)
    ASSERT_NE(manager.bcRecord, nullptr);
    ASSERT_NE(manager.bcRecordNew, nullptr);
    
    // 2. Simulate adding data to the "new" record
    manager.bcRecordNew->bcType.push_back(2); // BC::SOLID_SURFACE
    manager.bcRecordNew->bcNameId.push_back(1);
    
    // 3. Verify the Update() method (which uses operator* for value copy) works correctly
    manager.Update();
    
    // 4. Verify the data was successfully copied to the active record
    EXPECT_EQ(manager.bcRecord->bcType.size(), 1);
    EXPECT_EQ(manager.bcRecord->bcType[0], 2);
    EXPECT_EQ(manager.bcRecord->bcNameId[0], 1);
    
    // 5. When 'manager' goes out of scope, both unique_ptrs are safely destroyed.
}

// Verify that BcManager can be safely instantiated and destroyed on the stack.
TEST(BcManagerLifecycleTest, StackAllocationIsSafe)
{
    EXPECT_NO_THROW({
        ONEFLOW::BcManager manager;
        manager.PreProcess();
    });
}