// adt_tree_test.cpp
#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <cmath>
#include <cstdlib> // For rand() and RAND_MAX

// Include OneFLOW namespaces and headers
// Ensure CMake is configured with the correct include paths
#include "AdtTree.h" 
#include "HXDefine.h"

// ============================================================================
// 1. Test Fixture Definition
// Used to share initialization logic across multiple test cases
// ============================================================================
class AdtTreeTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize the bounding box of the tree: [0, 10] x [0, 10] x [0, 10]
        pmin_ = {0.0, 0.0, 0.0};
        pmax_ = {10.0, 10.0, 10.0};
    }

    std::vector<double> pmin_;
    std::vector<double> pmax_;
};

// ============================================================================
// 2. Specific Test Cases
// ============================================================================

// Test basic insertion and node counting
TEST_F(AdtTreeTest, BasicInsertionAndCount) {
    ONEFLOW::HXAdtTree<int, double> tree(3, pmin_.data(), pmax_.data());

    double pts[5][3] = {
        {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}, {5.0, 5.0, 5.0},
        {8.0, 8.0, 8.0}, {9.0, 9.0, 9.0}
    };

    for (int i = 0; i < 5; ++i) {
        // Note: The refactored AddNode takes ownership of the raw pointer
        tree.AddNode(new ONEFLOW::HXAdtNode<int, double>(3, pts[i], i));
    }

    // GTest assertion: prints expected vs. actual values on failure without aborting subsequent tests
    EXPECT_EQ(tree.nCount(), 5); 
}

// Test the accuracy of region queries
TEST_F(AdtTreeTest, RegionQueryAccuracy) {
    ONEFLOW::HXAdtTree<int, double> tree(3, pmin_.data(), pmax_.data());

    double pts[4][3] = {
        {1.0, 1.0, 1.0},  // ID 0: Inside the query box
        {2.0, 8.0, 2.0},  // ID 1: Outside the query box (Y is too high)
        {3.0, 3.0, 3.0},  // ID 2: Inside the query box
        {9.0, 9.0, 9.0}   // ID 3: Outside the query box
    };

    for (int i = 0; i < 4; ++i) {
        tree.AddNode(new ONEFLOW::HXAdtNode<int, double>(3, pts[i], i));
    }

    // Query box: [0.5, 4.0] for all dimensions
    double qmin[] = {0.5, 0.5, 0.5};
    double qmax[] = {4.0, 4.0, 4.0};

    ONEFLOW::HXAdtTree<int, double>::AdtNodeList results;
    tree.FindNodesInRegion(qmin, qmax, results);

    // Should find exactly 2 points
    EXPECT_EQ(results.size(), 2u); 

    // Verify that the correct point IDs are returned
    bool found0 = false, found2 = false;
    for (auto* node : results) {
        if (node->GetData() == 0) found0 = true;
        if (node->GetData() == 2) found2 = true;
    }
    EXPECT_TRUE(found0) << "Failed to find point ID 0";
    EXPECT_TRUE(found2) << "Failed to find point ID 2";
}

// Test floating-point boundary conditions (Critical for CFD applications)
TEST_F(AdtTreeTest, FloatingPointBoundaryConditions) {
    ONEFLOW::HXAdtTree<int, double> tree(3, pmin_.data(), pmax_.data());

    double pt[] = {5.0, 5.0, 5.0};
    tree.AddNode(new ONEFLOW::HXAdtNode<int, double>(3, pt, 99));

    // Query box exactly touching the point, with a tiny floating-point offset
    double qmin[] = {5.0 - 1e-12, 5.0 - 1e-12, 5.0 - 1e-12};
    double qmax[] = {6.0, 6.0, 6.0};

    ONEFLOW::HXAdtTree<int, double>::AdtNodeList results;
    tree.FindNodesInRegion(qmin, qmax, results);

    EXPECT_EQ(results.size(), 1u);
    EXPECT_EQ(results[0]->GetData(), 99);
}

// Test empty tree query to prevent null pointer dereference
TEST_F(AdtTreeTest, EmptyTreeQuery) {
    ONEFLOW::HXAdtTree<int, double> tree(3, pmin_.data(), pmax_.data());

    double qmin[] = {1.0, 1.0, 1.0};
    double qmax[] = {2.0, 2.0, 2.0};

    ONEFLOW::HXAdtTree<int, double>::AdtNodeList results;

    // Should not crash and should return an empty list
    EXPECT_NO_THROW(tree.FindNodesInRegion(qmin, qmax, results));
    EXPECT_TRUE(results.empty());
}

// Test large-scale data and memory safety (Verifies unique_ptr centralized management)
TEST_F(AdtTreeTest, StressTestAndMemorySafety) {
    ONEFLOW::HXAdtTree<int, double> tree(3, pmin_.data(), pmax_.data());

    const int NUM_POINTS = 10000;
    for (int i = 0; i < NUM_POINTS; ++i) {
        double pt[] = {
            static_cast<double>(rand()) / RAND_MAX * 10.0,
            static_cast<double>(rand()) / RAND_MAX * 10.0,
            static_cast<double>(rand()) / RAND_MAX * 10.0
        };
        tree.AddNode(new ONEFLOW::HXAdtNode<int, double>(3, pt, i));
    }

    EXPECT_EQ(tree.nCount(), NUM_POINTS);

    // When this TEST_F ends, 'tree' is automatically destroyed.
    // If the refactored memory management (unique_ptr) has flaws, 
    // this will trigger AddressSanitizer errors or a stack overflow crash.
}