#include <gtest/gtest.h>

#include "UResidualPolicy.h"

using namespace ONEFLOW;

TEST( ResidualDeposit, BoundaryAffectsOnlyLeft )
{
    Real left = 10.0;
    Real right = 20.0;
    DepositOrientedFaceFlux( 3.0, /*isBoundaryFace=*/true, left, right );
    EXPECT_DOUBLE_EQ( left, 7.0 );
    EXPECT_DOUBLE_EQ( right, 20.0 );
}

TEST( ResidualDeposit, InteriorOppositeSigns )
{
    Real left = 0.0;
    Real right = 0.0;
    DepositOrientedFaceFlux( 5.0, /*isBoundaryFace=*/false, left, right );
    EXPECT_DOUBLE_EQ( left, -5.0 );
    EXPECT_DOUBLE_EQ( right, 5.0 );
}

TEST( ResidualDeposit, ZeroFluxIsNoOp )
{
    Real left = 1.5;
    Real right = -2.5;
    DepositOrientedFaceFlux( 0.0, false, left, right );
    EXPECT_DOUBLE_EQ( left, 1.5 );
    EXPECT_DOUBLE_EQ( right, -2.5 );
}
