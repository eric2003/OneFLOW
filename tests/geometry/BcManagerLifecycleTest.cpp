// tests/BcManagerLifecycleTest.cpp
#include <gtest/gtest.h>
#include "BcRecord.h"
#include "FaceTopo.h"

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

TEST(BcManagerLifecycleTest, BoundaryFaceCountApiMatchesCompatibilityWrapper)
{
    // Intent: GetBoundaryFaceCount is preferred; GetNBFace is a compatibility
    // alias that must report the same value after Init and after direct edits.
    ONEFLOW::BcRecord record;
    EXPECT_EQ(record.GetBoundaryFaceCount(), 0);
    EXPECT_EQ(record.GetNBFace(), record.GetBoundaryFaceCount());

    record.Init(4);
    EXPECT_EQ(record.GetBoundaryFaceCount(), 4);
    EXPECT_EQ(record.GetNBFace(), record.GetBoundaryFaceCount());

    record.bcType.push_back(2);
    record.bcNameId.push_back(0);
    EXPECT_EQ(record.GetBoundaryFaceCount(), 5);
    EXPECT_EQ(record.GetNBFace(), record.GetBoundaryFaceCount());
}

TEST(BcManagerLifecycleTest, FaceTopoBoundaryCountDelegatesToBcRecord)
{
    // Intent: FaceTopo boundary-face count is owned by the active BcRecord.
    ONEFLOW::FaceTopo faceTopo;
    faceTopo.SetBoundaryFaceCount(3);
    EXPECT_EQ(faceTopo.GetBoundaryFaceCount(), static_cast<ONEFLOW::HXSize_t>(3));
    EXPECT_EQ(faceTopo.GetNBFaces(), faceTopo.GetBoundaryFaceCount());
    EXPECT_EQ(
        faceTopo.GetBcRecord().GetBoundaryFaceCount(),
        static_cast<int>(faceTopo.GetBoundaryFaceCount()));
}
