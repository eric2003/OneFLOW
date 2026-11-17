// tests/ElemFeatureLifecycleTest.cpp
#include <gtest/gtest.h>
#include "ElemFeature.h"

// Verify that ElemFeature correctly manages its internal value-type eTypes field
// without memory leaks or crashes, mimicking the behavior of SetElementTypeAndNode.
TEST(ElemFeatureLifecycleTest, ValueSemanticsForETypes)
{
    ONEFLOW::ElemFeature feature;
    
    // 1. Verify eTypes is accessible and empty initially
    EXPECT_EQ(feature.eTypes.size(), 0);
    
    // 2. Simulate internal resizing (mimicking CgnsSection::SetElementTypeAndNode)
    feature.eTypes.push_back(12); // HEXA_8
    feature.eTypes.push_back(12);
    feature.eTypes.push_back(5);  // PYRA_5
    
    EXPECT_EQ(feature.eTypes.size(), 3);
    EXPECT_EQ(feature.eTypes[0], 12);
    EXPECT_EQ(feature.eTypes[2], 5);
    
    // 3. feature goes out of scope here. 
    // The IntField value type is safely cleaned up by RAII.
}

// Verify that ElemFeature can be safely instantiated and destroyed on the stack.
TEST(ElemFeatureLifecycleTest, StackAllocationIsSafe)
{
    EXPECT_NO_THROW({
        ONEFLOW::ElemFeature feature;
        feature.eTypes.resize(1000);
    });
}