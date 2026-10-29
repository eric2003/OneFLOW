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
