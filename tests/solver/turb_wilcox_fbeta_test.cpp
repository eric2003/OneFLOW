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

// Intermediate: fk(xk) vs closed form on a sample grid.
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

// Intermediate: fw1998(xw)=(1+70 xw)/(1+80 xw)
TEST( WilcoxFbeta, Fw1998CurveMatchesClosedForm )
{
    const Real samples[] = { 0.0, 0.1, 0.5, 1.0, 2.0, 10.0 };
    for ( Real xw : samples )
    {
        const Real expected = ( 1.0 + 70.0 * xw ) / ( 1.0 + 80.0 * xw );
        EXPECT_NEAR( WilcoxFw1998( xw ), expected, 1.0e-12 );
    }
}

// Intermediate: fw2006(xw)=(1+85 xw)/(1+100 xw)
TEST( WilcoxFbeta, Fw2006CurveMatchesClosedForm )
{
    const Real samples[] = { 0.0, 0.1, 0.5, 1.0, 2.0, 10.0 };
    for ( Real xw : samples )
    {
        const Real expected = ( 1.0 + 85.0 * xw ) / ( 1.0 + 100.0 * xw );
        EXPECT_NEAR( WilcoxFw2006( xw ), expected, 1.0e-12 );
    }
}

// Intermediate: vortex-stretching invariant |xw| from strain/vorticity tensors.
// Closed-form matches WilcoxVortexStretchingInvariant; not a full residual.
TEST( WilcoxFbeta, VortexStretchingZeroWhenQuiet )
{
    EXPECT_DOUBLE_EQ(
        WilcoxVortexStretchingInvariant(
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
        0.0 );
}

TEST( WilcoxFbeta, VortexStretchingMatchesClosedFormSample )
{
    // w12=1, others 0; s11=2, s22=3, rest 0
    // xw = - w12^2 (s11+s22) = -1*(2+3) = -5, abs = 5
    EXPECT_DOUBLE_EQ(
        WilcoxVortexStretchingInvariant(
            1.0, 0.0, 0.0,
            2.0, 3.0, 0.0, 0.0, 0.0, 0.0 ),
        5.0 );
}

TEST( WilcoxFbeta, VortexStretchingCrossTerms )
{
    // w12=0,w13=1,w23=1; s12=2, rest strain 0 except s13=0,s23=0
    // xw = - (w13^2 (s11+s33) + w23^2 (s22+s33)) + 2*(-w13*w23*s12 + ...)
    //    = - (0 + 0) + 2*(-1*1*2) = -4, abs = 4
    EXPECT_DOUBLE_EQ(
        WilcoxVortexStretchingInvariant(
            0.0, 1.0, 1.0,
            0.0, 0.0, 0.0, 2.0, 0.0, 0.0 ),
        4.0 );
}

// Intermediate: WilcoxXkFromCrossTerm on an omega grid.
// Ref: cross=8; xk = 8 / omega^3.
TEST( TurbWilcoxFbeta, WilcoxXkFromCrossTermCurveMatchesClosedForm )
{
    const Real cross = 8.0;
    const Real omegas[] = { 0.5, 1.0, 2.0, 4.0 };
    for ( Real omega : omegas )
    {
        EXPECT_NEAR(
            WilcoxXkFromCrossTerm( cross, omega ),
            cross / ( omega * omega * omega ),
            1.0e-12 );
    }
}

TEST( TurbWilcoxFbeta, WilcoxXkFeedsFk )
{
    const Real xk = WilcoxXkFromCrossTerm( 8.0, 2.0 ); // 1.0
    EXPECT_NEAR( xk, 1.0, 1.0e-12 );
    const Real xk2 = xk * xk;
    const Real expected = ( 1.0 + 680.0 * xk2 ) / ( 1.0 + 400.0 * xk2 );
    EXPECT_NEAR( WilcoxFkFromXk( xk ), expected, 1.0e-12 );
}
