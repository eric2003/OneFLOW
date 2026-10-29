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
