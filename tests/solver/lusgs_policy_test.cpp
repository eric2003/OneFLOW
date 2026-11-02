#include <gtest/gtest.h>

#include "NsLusgsPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( LusgsPolicy, SplitFaceEigenvaluesSymmetricAtZeroSign )
{
    Real l1, l2, l3;
    LusgsSplitFaceEigenvalues( 1.0, 2.0, 0.0, 9.0, 9.0, 9.0, l1, l2, l3 );
    EXPECT_DOUBLE_EQ( l1, 0.5 );
    EXPECT_DOUBLE_EQ( l2, 1.5 );
    EXPECT_DOUBLE_EQ( l3, -0.5 );
}

TEST( LusgsPolicy, LinearizationWeightsMatchFormula )
{
    Real x1, x2;
    LusgsLinearizationWeights( 1.0, 3.0, 0.0, 4.0, 2.0, x1, x2 );
    EXPECT_DOUBLE_EQ( x1, ( 1.0 + 1.0 - 3.0 - 0.0 ) / ( 4.0 + 4.0 ) );
    EXPECT_DOUBLE_EQ( x2, ( 3.0 - 0.0 ) / ( 2.0 + 2.0 ) );
}

TEST( LusgsPolicy, FluxIncrementScalesWithArea )
{
    Real fr, fx, fy, fz, fe;
    LusgsFluxIncrementFiveEq(
        2.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 3.0,
        1.0, 0.0, 0.0, 0.0, 0.0,
        fr, fx, fy, fz, fe );
    EXPECT_DOUBLE_EQ( fr, 6.0 );
}

// Intermediate: lmd = 1/2 (lambda_phys + sign * faceAbs)
TEST( LusgsPolicy, SplitFaceEigenvaluesWithPositiveSign )
{
    Real l1, l2, l3;
    // vn=1, c=2 => phys 1, 3, -1; sign=+1; faceAbs all 0.4
    LusgsSplitFaceEigenvalues( 1.0, 2.0, 1.0, 0.4, 0.4, 0.4, l1, l2, l3 );
    EXPECT_NEAR( l1, 0.5 * ( 1.0 + 0.4 ), 1.0e-12 );
    EXPECT_NEAR( l2, 0.5 * ( 3.0 + 0.4 ), 1.0e-12 );
    EXPECT_NEAR( l3, 0.5 * ( -1.0 + 0.4 ), 1.0e-12 );
}

TEST( LusgsPolicy, SplitFaceEigenvaluesWithNegativeSign )
{
    Real l1, l2, l3;
    LusgsSplitFaceEigenvalues( 1.0, 2.0, -1.0, 0.4, 0.4, 0.4, l1, l2, l3 );
    EXPECT_NEAR( l1, 0.5 * ( 1.0 - 0.4 ), 1.0e-12 );
    EXPECT_NEAR( l2, 0.5 * ( 3.0 - 0.4 ), 1.0e-12 );
    EXPECT_NEAR( l3, 0.5 * ( -1.0 - 0.4 ), 1.0e-12 );
}

TEST( LusgsPolicy, LinearizationWeightsMoreSamples )
{
    Real x1, x2;
    // lmd1=2,lmd2=4,lmd3=0, c2=9, c=3
    // x1 = (2+2-4-0)/(9+9) = 0/18 = 0
    // x2 = (4-0)/(3+3) = 4/6
    LusgsLinearizationWeights( 2.0, 4.0, 0.0, 9.0, 3.0, x1, x2 );
    EXPECT_NEAR( x1, 0.0, 1.0e-12 );
    EXPECT_NEAR( x2, 4.0 / 6.0, 1.0e-12 );
}
