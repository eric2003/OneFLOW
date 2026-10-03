// tests/GridElemLifecycleTest.cpp
#include <gtest/gtest.h>
#include "ElemFeature.h"
#include "PointManager.h"
#include "FaceSolver.h"

// Verify that the internal components of GridElem can be safely 
// instantiated and destroyed on the stack as value types.
TEST(GridElemLifecycleTest, ValueSemanticsForInternalComponents)
{
    EXPECT_NO_THROW({
        // 1. ElemFeature (eTypes is a value type from Step 4.2)
        ONEFLOW::ElemFeature feature;
        feature.eTypes.push_back(12); // HEXA_8
        EXPECT_EQ(feature.eTypes.size(), 1);
        
        // 2. MeshPointManager
        ONEFLOW::MeshPointManager manager;
        manager.localToGlobal.push_back(0);
        EXPECT_EQ(manager.localToGlobal.size(), 1);
        
        // 3. FaceSolver (faceTopo is internally managed)
        ONEFLOW::FaceSolver solver;
        EXPECT_NE(solver.faceTopo, nullptr);

        // ElemFeature no longer stores a FaceSolver observer.
        // ScanElements receives its dependency explicitly when needed.

        // When these objects go out of scope, their destructors run safely.
        // This proves that GridElem's internal components can use value semantics.
    });
}