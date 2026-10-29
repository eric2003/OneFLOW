#include <gtest/gtest.h>

#include "ULimiterPolicy.h"

using namespace ONEFLOW;

TEST( LimiterPolicy, BarthReturnsZeroForOpposingFaceDirection )
{
    EXPECT_DOUBLE_EQ( BarthFunction( 0.5, -1.0, 1.0, -0.1 ), 0.0 );
}

TEST( LimiterPolicy, BarthConstrainsPositiveIncrement )
{
    EXPECT_DOUBLE_EQ( BarthFunction( 2.0, -1.0, 1.0, 1.0 ), 0.5 );
}

TEST( LimiterPolicy, BarthConstrainsNegativeIncrement )
{
    EXPECT_DOUBLE_EQ( BarthFunction( -2.0, -1.0, 1.0, 1.0 ), 0.5 );
}

TEST( LimiterPolicy, BarthPreservesIncrementInsideBounds )
{
    EXPECT_DOUBLE_EQ( BarthFunction( 0.25, -1.0, 1.0, 1.0 ), 1.0 );
}

TEST( LimiterPolicy, VenkatReturnsOneForNearZeroIncrement )
{
    EXPECT_DOUBLE_EQ( VencatFunction( 0.0, -1.0, 1.0, 1.0e-8 ), 1.0 );
}

TEST( LimiterPolicy, VenkatUsesUpperBoundForPositiveIncrement )
{
    const Real eps = 0.1;
    const Real expected = ( 1.0 + eps + 4.0 ) / ( 1.0 + 8.0 + 2.0 + eps );
    EXPECT_NEAR( VencatFunction( 2.0, -1.0, 1.0, eps ), expected, 1.0e-12 );
}

TEST( LimiterPolicy, VenkatUsesLowerBoundForNegativeIncrement )
{
    const Real eps = 0.1;
    const Real expected = ( 1.0 + eps + 4.0 ) / ( 1.0 + 8.0 + 2.0 + eps );
    EXPECT_NEAR( VencatFunction( -2.0, -1.0, 1.0, eps ), expected, 1.0e-12 );
}
