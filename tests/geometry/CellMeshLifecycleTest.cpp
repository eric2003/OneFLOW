// tests/CellMeshLifecycleTest.cpp
#include <gtest/gtest.h>
#include "CellMesh.h"

// Verify that CellMesh correctly manages its internal value-type cellTopo field
// without memory leaks or crashes, mimicking the behavior of Alloc().
TEST(CellMeshLifecycleTest, ValueSemanticsForCellTopo)
{
    ONEFLOW::CellMesh cellMesh;
    
    // 1. Verify cellTopo is accessible and empty initially
    EXPECT_EQ(cellMesh.GetCellTopo().GetNumberOfCells(), 0);
    
    // 2. Simulate internal resizing (mimicking UnsGrid::ReadGrid)
    cellMesh.GetCellTopo().Alloc(100);
    EXPECT_EQ(cellMesh.GetCellTopo().GetNumberOfCells(), 100);
    
    // 3. cellMesh goes out of scope here. 
    // The CellTopo value type is safely cleaned up by RAII.
}

// Verify that CellMesh can be safely instantiated and destroyed on the stack.
TEST(CellMeshLifecycleTest, StackAllocationIsSafe)
{
    EXPECT_NO_THROW({
        ONEFLOW::CellMesh cellMesh;
        cellMesh.GetCellTopo().Alloc(10);
    });
}