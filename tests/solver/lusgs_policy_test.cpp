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
