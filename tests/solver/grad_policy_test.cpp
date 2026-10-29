#include <gtest/gtest.h>

#include "UGradPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( GradPolicy, ArithmeticFaceAverageIsHalfSum )
{
    EXPECT_DOUBLE_EQ( ArithmeticFaceAverage( 2.0, 4.0 ), 3.0 );
}

TEST( GradPolicy, InverseDistanceWeightsSumToOne )
{
    Real cl = 0.0, cr = 0.0;
    // left center farther (2) than right (1) -> left weight smaller
    InverseDistanceFaceWeights( 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, cl, cr );
    EXPECT_NEAR( cl + cr, 1.0, 1.0e-12 );
    EXPECT_NEAR( cl, 1.0 / 3.0, 1.0e-12 );
    EXPECT_NEAR( cr, 2.0 / 3.0, 1.0e-12 );
}

TEST( GradPolicy, InverseDistanceSymmetricWhenEqualDistance )
{
    Real cl = 0.0, cr = 0.0;
    InverseDistanceFaceWeights( 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, cl, cr );
    EXPECT_NEAR( cl, 0.5, 1.0e-12 );
    EXPECT_NEAR( cr, 0.5, 1.0e-12 );
}

TEST( GradPolicy, WeightedFaceAverageMatchesWeights )
{
    EXPECT_DOUBLE_EQ( WeightedFaceAverage( 10.0, 20.0, 0.25, 0.75 ), 17.5 );
}

TEST( GradPolicy, AreaWeightedFaceNormalScalesUnitNormal )
{
    Real snx, sny, snz;
    AreaWeightedFaceNormal( 0.0, 1.0, 0.0, 3.0, snx, sny, snz );
    EXPECT_DOUBLE_EQ( snx, 0.0 );
    EXPECT_DOUBLE_EQ( sny, 3.0 );
    EXPECT_DOUBLE_EQ( snz, 0.0 );
}

TEST( GradPolicy, ScaleGradientByVolumeDividesComponents )
{
    Real gx = 2.0, gy = 4.0, gz = 6.0;
    ScaleGradientByVolume( 2.0, gx, gy, gz );
    EXPECT_DOUBLE_EQ( gx, 1.0 );
    EXPECT_DOUBLE_EQ( gy, 2.0 );
    EXPECT_DOUBLE_EQ( gz, 3.0 );
}

TEST( GradPolicy, CopyGradientCopiesAllComponents )
{
    Real dx = 0, dy = 0, dz = 0;
    CopyGradient( 1.5, -2.5, 3.5, dx, dy, dz );
    EXPECT_DOUBLE_EQ( dx, 1.5 );
    EXPECT_DOUBLE_EQ( dy, -2.5 );
    EXPECT_DOUBLE_EQ( dz, 3.5 );
}

TEST( GradPolicy, AccumulateFaceGradientAppliesSign )
{
    Real gx = 0, gy = 0, gz = 0;
    AccumulateFaceGradient( 2.0, 0.0, 0.0, 3.0, 1.0, gx, gy, gz );
    EXPECT_DOUBLE_EQ( gx, 6.0 );
    AccumulateFaceGradient( 2.0, 0.0, 0.0, 3.0, -1.0, gx, gy, gz );
    EXPECT_DOUBLE_EQ( gx, 0.0 );
}

TEST( GradPolicy, FaceSkewAngleDegreesIsZeroForParallel )
{
    // projectedLength equal to |d| -> asin(1) = 90 degrees
    EXPECT_NEAR( FaceSkewAngleDegrees( 1.0, 0.0, 0.0, 1.0 ), 90.0, 1.0e-10 );
    // zero projection -> 0 degrees
    EXPECT_NEAR( FaceSkewAngleDegrees( 1.0, 0.0, 0.0, 0.0 ), 0.0, 1.0e-10 );
}

TEST( GradPolicy, FaceSkewAngleDegreesClampsArgument )
{
    // |proj| > |d| would exceed asin domain without clamp
    EXPECT_NEAR( FaceSkewAngleDegrees( 1.0, 0.0, 0.0, 2.0 ), 90.0, 1.0e-10 );
    EXPECT_NEAR( FaceSkewAngleDegrees( 1.0, 0.0, 0.0, -2.0 ), -90.0, 1.0e-10 );
}

TEST( GradPolicy, UnitizeDisplacementMakesUnitVector )
{
    Real dx = 3.0, dy = 0.0, dz = 4.0, inv = 0.0;
    UnitizeDisplacement( dx, dy, dz, inv );
    EXPECT_NEAR( inv, 0.2, 1.0e-12 );
    EXPECT_NEAR( dx, 0.6, 1.0e-12 );
    EXPECT_NEAR( dy, 0.0, 1.0e-12 );
    EXPECT_NEAR( dz, 0.8, 1.0e-12 );
}

TEST( GradPolicy, CorrectFaceGradientMatchesEndpointDifference )
{
    // Unit direction +x, invLength=1; qLeft=0, qRight=2 -> target=2
    Real gx = 0.0, gy = 0.0, gz = 0.0;
    CorrectFaceGradient( 0.0, 2.0, gx, gy, gz, 1.0, 0.0, 0.0, 1.0 );
    EXPECT_NEAR( gx, 2.0, 1.0e-12 );
    EXPECT_NEAR( gy, 0.0, 1.0e-12 );
}
