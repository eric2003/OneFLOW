// point_locator_test.cpp
#include <gtest/gtest.h>
#include <cmath>

// Include necessary headers from your project
#include "PointLocator.h"
#include "HXDefine.h"

class PointLocatorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Safe initialization for HXVector (inherits from std::vector)
        pmin_ = ONEFLOW::RealField{0.0, 0.0, 0.0};
        pmax_ = ONEFLOW::RealField{10.0, 10.0, 10.0};
        tolerance_ = 1.0e-4; // Explicit tolerance for testing
    }

    ONEFLOW::RealField pmin_;
    ONEFLOW::RealField pmax_;
    ONEFLOW::Real tolerance_;
};

//// Test the critical scenario: Multiple points within tolerance bounding box
//// This previously caused a FATAL ERROR and Stop()
//TEST_F(PointLocatorTest, MultipleNodesInToleranceBox) {
//    ONEFLOW::PointLocator locator;
//
//    // Initialize by passing RealField by reference (perfect match)
//    locator.Initialize(pmin_, pmax_, tolerance_);
//
//    // Add a base point (ID 0)
//    locator.AddPoint(5.0, 5.0, 5.0); 
//
//    // Add another point VERY close to the base point, but still within tolerance
//    // Distance is 0.5e-5, which is < 1.0e-4 tolerance
//    locator.AddPoint(5.0 + 0.5e-5, 5.0, 5.0); // ID 1
//
//    // Add a third point also within tolerance of the base point
//    locator.AddPoint(5.0, 5.0 + 0.8e-5, 5.0); // ID 2
//
//    // Now, query for a point exactly at the base location
//    // It should safely return the ID of the closest point (ID 0) without crashing
//    int foundId = locator.FindPoint(5.0, 5.0, 5.0);
//
//    EXPECT_EQ(foundId, 0);
//    EXPECT_EQ(locator.GetNPoint(), 3); // Ensure all 3 were added as unique points
//}

// Test the critical scenario: Multiple nodes in tolerance bounding box
// This previously caused a FATAL ERROR and Stop()
TEST_F(PointLocatorTest, MultipleNodesInToleranceBox) {
    ONEFLOW::PointLocator locator;

    // Use a larger tolerance for clearer geometric demonstration
    double test_tol = 1.0; 
    locator.Initialize(pmin_, pmax_, test_tol);

    // 1. Add points that are >= tolerance apart, ensuring they are all uniquely added
    locator.AddPoint(0.0, 0.0, 0.0); // ID 0
    locator.AddPoint(1.0, 0.0, 0.0); // ID 1 (Distance to ID 0 is 1.0, not < 1.0, so it's added)
    locator.AddPoint(0.0, 1.0, 0.0); // ID 2 (Distance to ID 0 is 1.0, so it's added)

    // Verify all 3 points were successfully added (Deduplication worked correctly)
    EXPECT_EQ(locator.GetNPoint(), 3);

    // 2. Query a point whose tolerance bounding box covers ALL THREE points
    // Query point: (0.5, 0.5, 0.0)
    // Tolerance box: [-0.5, 1.5] for X and Y. 
    // All three points (0,0,0), (1,0,0), and (0,1,0) fall inside this box.
    int foundId = locator.FindPoint(0.5, 0.5, 0.0);

    // 3. The refactored code should NOT crash. 
    // It should calculate exact distances and return the ID of the closest point.
    // Distance to ID 0: sqrt(0.5^2 + 0.5^2) = sqrt(0.5) ¡Ö 0.707 < 1.0 (Valid)
    // Distance to ID 1: sqrt(0.5^2 + 0.5^2) = sqrt(0.5) ¡Ö 0.707 < 1.0 (Valid)
    // Distance to ID 2: sqrt(0.5^2 + 0.5^2) = sqrt(0.5) ¡Ö 0.707 < 1.0 (Valid)

    // Since distances are equal, it will return the first one it evaluates as the minimum (likely ID 0).
    // The critical success criterion is that it returns a VALID index (0, 1, or 2) instead of crashing.
    EXPECT_TRUE(foundId == 0 || foundId == 1 || foundId == 2);
}

// Test that points outside tolerance are correctly identified as new
TEST_F(PointLocatorTest, PointsOutsideToleranceAreUnique) {
    ONEFLOW::PointLocator locator;
    locator.Initialize(pmin_, pmax_, tolerance_);

    locator.AddPoint(1.0, 1.0, 1.0); // ID 0

    // This point is further away than the tolerance
    int newId = locator.AddPoint(1.0 + tolerance_ * 2.0, 1.0, 1.0);

    EXPECT_EQ(newId, 1);
    EXPECT_EQ(locator.GetNPoint(), 2);
}

// Test finding a non-existent point
TEST_F(PointLocatorTest, FindNonExistentPoint) {
    ONEFLOW::PointLocator locator;
    locator.Initialize(pmin_, pmax_, tolerance_);

    locator.AddPoint(2.0, 2.0, 2.0);

    // Query a point far outside the tolerance
    int foundId = locator.FindPoint(9.0, 9.0, 9.0);

    // Note: If INVALID_INDEX is not found, replace it with -1 or include the header where it's defined (e.g., HXDefine.h)
    EXPECT_EQ(foundId, ONEFLOW::INVALID_INDEX);
}