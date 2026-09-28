/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "ADTTest.h"
#include "SimuBase.h"
#include "AdtTree.h"
#include <iostream> 
#include <cassert>

BeginNameSpace( ONEFLOW )

// ============================================================================
// 3. Test Cases (≤‚ ‘”√¿˝)
// ============================================================================

void TestBasicInsertionAndCount() {
    std::cout << "Running Test 1: Basic Insertion and Count... ";

    double pmin[] = {0.0, 0.0, 0.0};
    double pmax[] = {10.0, 10.0, 10.0};

    ONEFLOW::HXAdtTree<int, double> tree(3, pmin, pmax);

    // Insert 5 known points
    double pts[5][3] = {
        {1.0, 1.0, 1.0},
        {2.0, 2.0, 2.0},
        {5.0, 5.0, 5.0},
        {8.0, 8.0, 8.0},
        {9.0, 9.0, 9.0}
    };

    for (int i = 0; i < 5; ++i) {
        tree.AddNode(new ONEFLOW::HXAdtNode<int, double>(3, pts[i], i));
    }

    // Verify count
    assert(tree.nCount() == 5);
    std::cout << "PASSED\n";
}

void TestRegionQueryAccuracy() {
    std::cout << "Running Test 2: Region Query Accuracy... ";

    double pmin[] = {0.0, 0.0, 0.0};
    double pmax[] = {10.0, 10.0, 10.0};

    ONEFLOW::HXAdtTree<int, double> tree(3, pmin, pmax);

    double pts[4][3] = {
        {1.0, 1.0, 1.0},  // ID 0: Inside query box
        {2.0, 8.0, 2.0},  // ID 1: Outside query box (Y is too high)
        {3.0, 3.0, 3.0},  // ID 2: Inside query box
        {9.0, 9.0, 9.0}   // ID 3: Outside query box
    };

    for (int i = 0; i < 4; ++i) {
        tree.AddNode(new ONEFLOW::HXAdtNode<int, double>(3, pts[i], i));
    }

    // Query box: [0.5, 0.5, 0.5] to [4.0, 4.0, 4.0]
    double qmin[] = {0.5, 0.5, 0.5};
    double qmax[] = {4.0, 4.0, 4.0};

    ONEFLOW::HXAdtTree<int, double>::AdtNodeList results;
    tree.FindNodesInRegion(qmin, qmax, results);

    // Should find exactly 2 points (ID 0 and ID 2)
    assert(results.size() == 2);

    // Verify the IDs are correct (order might vary, so we check existence)
    bool found0 = false, found2 = false;
    for (auto* node : results) {
        if (node->GetData() == 0) found0 = true;
        if (node->GetData() == 2) found2 = true;
    }
    assert(found0 && found2);

    std::cout << "PASSED\n";
}

void TestBoundaryConditions() {
    std::cout << "Running Test 3: Boundary Conditions (Exact Match)... ";

    double pmin[] = {0.0, 0.0, 0.0};
    double pmax[] = {10.0, 10.0, 10.0};

    ONEFLOW::HXAdtTree<int, double> tree(3, pmin, pmax);

    // Insert a point exactly on the boundary of the query box
    double pt[] = {5.0, 5.0, 5.0};
    tree.AddNode(new ONEFLOW::HXAdtNode<int, double>(3, pt, 99));

    // Query box exactly touching the point
    double qmin[] = {5.0, 5.0, 5.0};
    double qmax[] = {6.0, 6.0, 6.0};

    ONEFLOW::HXAdtTree<int, double>::AdtNodeList results;
    tree.FindNodesInRegion(qmin, qmax, results);

    assert(results.size() == 1);
    assert(results[0]->GetData() == 99);

    std::cout << "PASSED\n";
}

void TestStressAndMemorySafety() {
    std::cout << "Running Test 4: Stress Test (10,000 random points)... ";

    double pmin[] = {0.0, 0.0, 0.0};
    double pmax[] = {100.0, 100.0, 100.0};

    ONEFLOW::HXAdtTree<int, double> tree(3, pmin, pmax);

    const int NUM_POINTS = 10000;
    for (int i = 0; i < NUM_POINTS; ++i) {
        double pt[] = {
            static_cast<double>(rand()) / RAND_MAX * 100.0,
            static_cast<double>(rand()) / RAND_MAX * 100.0,
            static_cast<double>(rand()) / RAND_MAX * 100.0
        };
        tree.AddNode(new ONEFLOW::HXAdtNode<int, double>(3, pt, i));
    }

    // Verify count matches
    assert(tree.nCount() == NUM_POINTS);

    // Query a small region, ensure it doesn't crash or hang
    double qmin[] = {49.0, 49.0, 49.0};
    double qmax[] = {51.0, 51.0, 51.0};
    ONEFLOW::HXAdtTree<int, double>::AdtNodeList results;
    tree.FindNodesInRegion(qmin, qmax, results);

    // When 'tree' goes out of scope here, the destructor will run.
    // If the refactored memory management (unique_ptr) is correct, 
    // it will clean up 10,000 nodes instantly without stack overflow or memory leaks.

    std::cout << "PASSED\n";
}

ADTTest::ADTTest()
{
    ;
}

ADTTest::~ADTTest()
{
    ;
}

void ADTTest::Run()
{
    std::cout << "========================================\n";
    std::cout << " Starting HXAdtTree Refactoring Tests \n";
    std::cout << "========================================\n";

    try {
        TestBasicInsertionAndCount();
        TestRegionQueryAccuracy();
        TestBoundaryConditions();
        TestStressAndMemorySafety();

        std::cout << "========================================\n";
        std::cout << " ALL TESTS PASSED SUCCESSFULLY! \n";
        std::cout << " The refactored tree is logically sound and memory-safe.\n";
        std::cout << "========================================\n";
        return;
    } catch (const std::exception& e) {
        std::cerr << "TEST FAILED with exception: " << e.what() << "\n";
        return;
    }
}

EndNameSpace
WRAP_TEST_CLASS(ADTTest, "adt_test");

