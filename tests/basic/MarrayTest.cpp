#include <gtest/gtest.h>

#include "Marray.h"

namespace
{
using ONEFLOW::HXSize_t;
using ONEFLOW::Marray;

TEST(MarrayTest, ConstructorPreservesBothSizeTypeExtents)
{
    const HXSize_t equationCount = 2;
    const HXSize_t valueCount = 3;

    Marray< int > values( equationCount, valueCount );

    EXPECT_EQ(values.GetNEqu(), equationCount);
    ASSERT_EQ(values[ 0 ].size(), valueCount);
    ASSERT_EQ(values[ 1 ].size(), valueCount);
}
}
