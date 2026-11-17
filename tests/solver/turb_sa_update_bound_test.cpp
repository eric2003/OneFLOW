#include <gtest/gtest.h>

#include "TurbSaPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( SaUpdateBound, UpperBoundFormula )
{
    // 10 * 2 / (1 * 0.5) = 40
    EXPECT_DOUBLE_EQ(
        SaWorkingVariableUpperBound( 10.0, 2.0, 1.0, 0.5 ), 40.0 );
}

TEST( SaUpdateBound, NonNegativeKeepsPositive )
{
    EXPECT_DOUBLE_EQ( SaNonNegativeWorkingVariable( 1.5, 0.1 ), 1.5 );
}

TEST( SaUpdateBound, NonNegativeReplacesNegative )
{
    EXPECT_DOUBLE_EQ( SaNonNegativeWorkingVariable( -2.0, 0.3 ), 0.3 );
}

// Intermediate: SaWorkingVariableUpperBound on an fv1 grid.
// Ref: rmax=10, mu=2, rho=1; bound = 20 / fv1.
TEST( SaUpdateBound, SaWorkingVariableUpperBoundCurveMatchesClosedForm )
{
    const Real rmax = 10.0;
    const Real mu = 2.0;
    const Real density = 1.0;
    const Real fv1s[] = { 0.25, 0.5, 1.0, 2.0 };
    for ( Real fv1 : fv1s )
    {
        EXPECT_NEAR(
            SaWorkingVariableUpperBound( rmax, mu, density, fv1 ),
            rmax * mu / ( density * fv1 ),
            1.0e-12 );
    }
}

TEST( SaUpdateBound, SaWorkingVariableUpperBoundInverselyScalesWithFv1 )
{
    const Real b1 = SaWorkingVariableUpperBound( 10.0, 2.0, 1.0, 0.5 );
    const Real b2 = SaWorkingVariableUpperBound( 10.0, 2.0, 1.0, 1.0 );
    EXPECT_NEAR( b1, 2.0 * b2, 1.0e-12 );
}
