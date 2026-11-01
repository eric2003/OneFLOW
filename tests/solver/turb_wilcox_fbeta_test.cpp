#include <gtest/gtest.h>
#include <cmath>

#include "TurbSstPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( WilcoxFbeta, FkIsOneWhenXkNonPositive )
{
    EXPECT_DOUBLE_EQ( WilcoxFkFromXk( 0.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( WilcoxFkFromXk( -1.0 ), 1.0 );
}

TEST( WilcoxFbeta, FkFormulaWhenXkPositive )
{
    const Real xk = 0.1;
    const Real xk2 = xk * xk;
    const Real expected = ( 1.0 + 680.0 * xk2 ) / ( 1.0 + 400.0 * xk2 );
    EXPECT_DOUBLE_EQ( WilcoxFkFromXk( xk ), expected );
}

TEST( WilcoxFbeta, Fw1998AtZeroIsOne )
{
    EXPECT_DOUBLE_EQ( WilcoxFw1998( 0.0 ), 1.0 );
}

TEST( WilcoxFbeta, Fw2006AtZeroIsOne )
{
    EXPECT_DOUBLE_EQ( WilcoxFw2006( 0.0 ), 1.0 );
}

TEST( WilcoxFbeta, Fw1998ApproachesSevenEighths )
{
    // large xw -> 70/80 = 0.875
    EXPECT_NEAR( WilcoxFw1998( 1.0e6 ), 0.875, 1.0e-6 );
}

// Intermediate: fk(xk) on a sample grid vs closed form.
// Proves WilcoxFkFromXk; does not prove full omega-equation residual.
TEST( WilcoxFbeta, FkXkCurveMatchesClosedForm )
{
    const Real samples[] = { -0.5, -0.1, 0.0, 0.05, 0.1, 0.2, 0.5, 1.0 };
    for ( Real xk : samples )
    {
        Real expected = 1.0;
        if ( xk > 0.0 )
        {
            const Real xk2 = xk * xk;
            expected = ( 1.0 + 680.0 * xk2 ) / ( 1.0 + 400.0 * xk2 );
        }
        EXPECT_NEAR( WilcoxFkFromXk( xk ), expected, 1.0e-12 );
    }
}
