#include <gtest/gtest.h>

#include "ULimiterPolicy.h"
#include "HXMath.h"

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
    const Real expected = ( 9.0 + eps + 12.0 ) / ( 9.0 + 8.0 + 6.0 + eps );
    EXPECT_NEAR( VencatFunction( -2.0, -3.0, 1.0, eps ), expected, 1.0e-12 );
}

TEST( LimiterPolicy, EricVariantUsesItsOwnCrossTermCoefficient )
{
    const Real eps = 0.1;
    const Real expected = ( 1.0 + eps + 2.0 ) / ( 1.0 + 8.0 + 2.0 + eps );
    EXPECT_NEAR( VencatEric( 1.0, 2.0, eps ), expected, 1.0e-12 );
}

TEST( LimiterPolicy, EricVariantSelectsLowerBoundForNegativeIncrement )
{
    const Real eps = 0.1;
    const Real expected = ( 9.0 + eps + 6.0 ) / ( 9.0 + 8.0 + 6.0 + eps );
    EXPECT_NEAR( VencatFunctionEric( -2.0, -3.0, 1.0, eps ), expected, 1.0e-12 );
}

TEST( LimiterPolicy, VenkatRegularizationUsesUnitFloorForNarrowGlobalRange )
{
    const Real eps = VencatRegularizationEpsilon( 0.1, -0.2, 0.2 );
    EXPECT_DOUBLE_EQ( eps, 0.01 + SMALL );
}

TEST( LimiterPolicy, VenkatRegularizationScalesWithGlobalOffsetRange )
{
    const Real eps = VencatRegularizationEpsilon( 0.1, -2.0, 3.0 );
    EXPECT_DOUBLE_EQ( eps, 0.25 + SMALL );
}

TEST( LimiterPolicy, ReconstructFaceIncrementIsGradientDotDisplacement )
{
    EXPECT_DOUBLE_EQ(
        ReconstructFaceIncrement( 1.0, 2.0, 3.0, 0.5, 0.0, 0.0 ), 0.5 );
    EXPECT_DOUBLE_EQ(
        ReconstructFaceIncrement( 1.0, 2.0, 3.0, 0.0, 1.0, 0.0 ), 2.0 );
}

TEST( LimiterPolicy, FaceSideDirectionDotUsesOrientationSign )
{
    // Unit normal +x, displacement +x of length 2 -> +1 with orientation +1.
    EXPECT_NEAR(
        FaceSideDirectionDot( 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0 ), 1.0, 1.0e-12 );
    EXPECT_NEAR(
        FaceSideDirectionDot( 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, -1.0 ), -1.0, 1.0e-12 );
}

TEST( LimiterPolicy, TightenLimiterKeepsStricterBound )
{
    EXPECT_DOUBLE_EQ( TightenLimiter( 1.0, 0.3 ), 0.3 );
    EXPECT_DOUBLE_EQ( TightenLimiter( 0.2, 0.5 ), 0.2 );
}

TEST( LimiterPolicy, ApplyLimitedFaceReconstructionAddsScaledIncrement )
{
    EXPECT_DOUBLE_EQ(
        ApplyLimitedFaceReconstruction( 10.0, 0.5, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0 ),
        11.0 );
}

TEST( LimiterPolicy, BlendFaceLimiterRespectsBothSidesMode )
{
    EXPECT_DOUBLE_EQ( BlendFaceLimiter( 0.8, 0.3, false ), 0.8 );
    EXPECT_DOUBLE_EQ( BlendFaceLimiter( 0.8, 0.3, true ), 0.3 );
}
