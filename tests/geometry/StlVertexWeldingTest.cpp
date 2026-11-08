#include "PointLocator.h"

#include <gtest/gtest.h>

namespace {

void InitializeLocator(ONEFLOW::PointLocator& locator, double tolerance)
{
    ONEFLOW::RealField minimum(3), maximum(3);
    for (int axis = 0; axis < 3; ++axis)
    {
        minimum[axis] = -10.0;
        maximum[axis] = 10.0;
    }
    locator.Initialize(minimum, maximum, tolerance);
}

TEST(StlVertexWeldingTest, MatchesPointsAcrossRoundingCellBoundary)
{
    ONEFLOW::PointLocator locator;
    InitializeLocator(locator, 0.1);
    const int first = locator.AddPoint(0.049, 0.0, 0.0);
    const int second = locator.AddPoint(0.051, 0.0, 0.0);
    EXPECT_EQ(first, second);
    EXPECT_EQ(locator.GetNPoint(), 1);
}

TEST(StlVertexWeldingTest, DoesNotMergeDistantPointsInSameRoundedCell)
{
    ONEFLOW::PointLocator locator;
    InitializeLocator(locator, 1.0);
    const int first = locator.AddPoint(0.49, 0.49, 0.0);
    const int second = locator.AddPoint(-0.49, -0.49, 0.0);
    EXPECT_NE(first, second);
    EXPECT_EQ(locator.GetNPoint(), 2);
}

TEST(StlVertexWeldingTest, UsesRepresentativeRatherThanTransitiveClosure)
{
    ONEFLOW::PointLocator locator;
    InitializeLocator(locator, 1.0);
    const int a = locator.AddPoint(0.0, 0.0, 0.0);
    const int b = locator.AddPoint(0.75, 0.0, 0.0);
    const int c = locator.AddPoint(1.5, 0.0, 0.0);
    EXPECT_EQ(a, b);
    EXPECT_NE(a, c);
    EXPECT_EQ(locator.GetNPoint(), 2);
}

} // namespace
