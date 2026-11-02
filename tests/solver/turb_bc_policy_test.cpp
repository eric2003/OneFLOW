#include <gtest/gtest.h>
#include <cmath>

#include "TurbBcPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbBcPolicy, SaWallImageNegates )
{
    EXPECT_DOUBLE_EQ( TurbWallSaImage( 3.0 ), -3.0 );
}

TEST( TurbBcPolicy, WallOmegaMenter )
{
    // 60 * 2 / (1 * 0.075 * 4 * 1) = 120 / 0.3 = 400
    EXPECT_DOUBLE_EQ(
        TurbWallOmegaMenter( 2.0, 1.0, 0.075, 4.0, 1.0 ), 400.0 );
}

TEST( TurbBcPolicy, KineticEnergyImage )
{
    EXPECT_DOUBLE_EQ( TurbWallKineticEnergyImage( 1.5 ), -1.5 );
}

TEST( TurbBcPolicy, FarfieldUsesInteriorWhenOutflow )
{
    EXPECT_TRUE( TurbFarfieldUsesInterior( 1.0 ) );
    EXPECT_TRUE( TurbFarfieldUsesInterior( 0.0 ) );
    EXPECT_FALSE( TurbFarfieldUsesInterior( -0.1 ) );
}

// Intermediate: omega_w = 60 mu / (rho beta1 d^2 Re) on a d grid.
// Ref constants mu=rho=Re=1, beta1=0.075 => omega = 800 / d^2.
TEST( TurbBcPolicy, WallOmegaMenterCurveMatchesClosedForm )
{
    const Real mu = 1.0;
    const Real rho = 1.0;
    const Real beta1 = 0.075;
    const Real re = 1.0;
    const Real ds[] = { 0.05, 0.1, 0.2, 0.5, 1.0, 2.0 };
    for ( Real d : ds )
    {
        const Real d2 = d * d;
        const Real expected = 60.0 * mu / ( rho * beta1 * d2 * re );
        EXPECT_NEAR(
            TurbWallOmegaMenter( mu, rho, beta1, d2, re ),
            expected,
            1.0e-12 );
    }
}
