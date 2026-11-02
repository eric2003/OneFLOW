#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "HXVector.h"

namespace
{
using ONEFLOW::HXVector;

TEST(HXVectorTest, InitializerListAssignmentReplacesExistingContents)
{
    HXVector< int > values = { 9, 8, 7, 6 };

    values = { 0, 1, 0 };

    ASSERT_EQ(values.size(), 3u);
    EXPECT_EQ(values[ 0 ], 0);
    EXPECT_EQ(values[ 1 ], 1);
    EXPECT_EQ(values[ 2 ], 0);
}

TEST(HXVectorTest, EmptyInitializerListAssignmentClearsContents)
{
    HXVector< int > values = { 1, 2, 3 };

    values = {};

    EXPECT_TRUE(values.empty());
}

TEST(HXVectorTest, ScalarAssignmentStillFillsExistingElements)
{
    HXVector< int > values = { 1, 2, 3 };

    values = 7;

    ASSERT_EQ(values.size(), 3u);
    EXPECT_EQ(values[ 0 ], 7);
    EXPECT_EQ(values[ 1 ], 7);
    EXPECT_EQ(values[ 2 ], 7);
}

TEST(HXVectorTest, StandardVectorCopyAssignmentRemainsSupported)
{
    HXVector< int > values = { 1, 2, 3 };
    const std::vector< int > source = { 4, 5 };

    values = source;

    ASSERT_EQ(values.size(), 2u);
    EXPECT_EQ(values[ 0 ], 4);
    EXPECT_EQ(values[ 1 ], 5);
}

TEST(HXVectorTest, MoveOnlyElementTypeRemainsSupported)
{
    HXVector< std::unique_ptr< int > > values;
    values.push_back(std::make_unique< int >(42));

    ASSERT_EQ(values.size(), 1u);
    ASSERT_NE(values[ 0 ], nullptr);
    EXPECT_EQ(*values[ 0 ], 42);
}
}
