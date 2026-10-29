#include <gtest/gtest.h>

#include "NsInvFluxPolicy.h"

using namespace ONEFLOW;

TEST( InvFluxPolicy, RelativeFaceNormalVelocitySubtractsGridSpeed )
{
    EXPECT_DOUBLE_EQ(
        RelativeFaceNormalVelocity( 1.0, 0.0, 0.0, 0.5, 2.0, 0.0, 0.0 ), 1.5 );
}

TEST( InvFluxPolicy, EulerConvectiveFluxMatchesMassAndPressureTerms )
{
    Real fr, fx, fy, fz, fe;
    EulerConvectiveFlux(
        1.0, 2.0, 0.0, 0.0, 3.0, 10.0,
        2.0, 1.0, 0.0, 0.0, 0.0,
        fr, fx, fy, fz, fe );
    EXPECT_DOUBLE_EQ( fr, 2.0 );
    EXPECT_DOUBLE_EQ( fx, 2.0 * 2.0 + 3.0 );
    EXPECT_DOUBLE_EQ( fy, 0.0 );
    EXPECT_DOUBLE_EQ( fz, 0.0 );
    EXPECT_DOUBLE_EQ( fe, 2.0 * 10.0 );
}

TEST( InvFluxPolicy, EulerConvectiveFluxZeroNormalVelocityIsPressureOnly )
{
    Real fr, fx, fy, fz, fe;
    EulerConvectiveFlux(
        1.2, 0.0, 0.0, 0.0, 5.0, 8.0,
        0.0, 0.0, 1.0, 0.0, 0.0,
        fr, fx, fy, fz, fe );
    EXPECT_DOUBLE_EQ( fr, 0.0 );
    EXPECT_DOUBLE_EQ( fx, 0.0 );
    EXPECT_DOUBLE_EQ( fy, 5.0 );
    EXPECT_DOUBLE_EQ( fz, 0.0 );
    EXPECT_DOUBLE_EQ( fe, 0.0 );
}

TEST( InvFluxPolicy, RoeAverageWeightsMatchDensityRatio )
{
    Real ratio = 0.0, coef = 0.0;
    RoeAverageWeights( 1.0, 4.0, ratio, coef );
    EXPECT_DOUBLE_EQ( ratio, 2.0 );
    EXPECT_DOUBLE_EQ( coef, 1.0 / 3.0 );
}

TEST( InvFluxPolicy, RoeAverageComponentIsWeightedMean )
{
    EXPECT_DOUBLE_EQ( RoeAverageComponent( 1.0, 4.0, 2.0, 1.0 / 3.0 ), 3.0 );
}

TEST( InvFluxPolicy, RoeAverageStateUsesGeometricDensity )
{
    Real rm, um, vm, wm, pm, hm, gama;
    RoeAverageState(
        1.0, 0.0, 0.0, 0.0, 1.0, 3.5, 1.4,
        4.0, 3.0, 0.0, 0.0, 4.0, 7.0, 1.4,
        rm, um, vm, wm, pm, hm, gama );
    EXPECT_DOUBLE_EQ( rm, 2.0 );
    EXPECT_DOUBLE_EQ( um, 2.0 );
    EXPECT_DOUBLE_EQ( gama, 1.4 );
}
