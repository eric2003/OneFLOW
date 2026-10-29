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
