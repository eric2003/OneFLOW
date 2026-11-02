#include <gtest/gtest.h>

#include "TurbInvFluxPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbInvFluxPolicy, SplitRelativeNormalVelocityPositive )
{
    Real vp, vm;
    SplitRelativeNormalVelocity( 4.0, vp, vm );
    EXPECT_DOUBLE_EQ( vp, 4.0 );
    EXPECT_DOUBLE_EQ( vm, 0.0 );
}

TEST( TurbInvFluxPolicy, SplitRelativeNormalVelocityNegative )
{
    Real vp, vm;
    SplitRelativeNormalVelocity( -3.0, vp, vm );
    EXPECT_DOUBLE_EQ( vp, 0.0 );
    EXPECT_DOUBLE_EQ( vm, -3.0 );
}

TEST( TurbInvFluxPolicy, TurbConvectiveDensityWeightSaIsUnity )
{
    EXPECT_DOUBLE_EQ( TurbConvectiveDensityWeight( 1.0, 99.0 ), 1.0 );
}

TEST( TurbInvFluxPolicy, TurbConvectiveDensityWeightKwIsDensity )
{
    EXPECT_DOUBLE_EQ( TurbConvectiveDensityWeight( 0.0, 1.25 ), 1.25 );
}

TEST( TurbInvFluxPolicy, UpwindScalarConvectiveFluxLeftOnly )
{
    EXPECT_DOUBLE_EQ(
        UpwindScalarConvectiveFlux( 2.0, 1.0, 3.0, 0.0, 1.0, 9.0 ), 6.0 );
}
