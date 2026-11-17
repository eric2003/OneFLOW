#include <gtest/gtest.h>

#include "TurbSstPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSstResidual, SourceKIsProdMinusDiss )
{
    EXPECT_DOUBLE_EQ( SstSourceK( 5.0, 2.0 ), 3.0 );
}

TEST( TurbSstResidual, SourceWIncludesCrossDiff )
{
    EXPECT_DOUBLE_EQ( SstSourceW( 4.0, 1.0, 0.5 ), 3.5 );
}

TEST( TurbSstResidual, DiagonalK )
{
    EXPECT_DOUBLE_EQ( SstDiagonalK( 1.0, 0.09, 2.0 ), -0.36 );
}

TEST( TurbSstResidual, NegativeSourcePart )
{
    EXPECT_DOUBLE_EQ( SstNegativeSourcePart( 3.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( SstNegativeSourcePart( -4.0 ), -4.0 );
}

TEST( TurbSstResidual, DiagonalWWithoutCrossDiff )
{
    // -2 * fBeta * beta * omega - 0/(rho*omega)
    EXPECT_DOUBLE_EQ( SstDiagonalW( 1.0, 0.075, 2.0, 0.0, 1.0 ), -0.3 );
}

TEST( TurbSstResidual, DiagonalWWithCrossDiff )
{
    // -2*1*0.075*2 - 0.4/(1*2) = -0.3 - 0.2 = -0.5
    EXPECT_DOUBLE_EQ( SstDiagonalW( 1.0, 0.075, 2.0, 0.4, 1.0 ), -0.5 );
}

// Intermediate: SstSourceK on a production grid.
// Ref: eps=2; Sk = Pk - 2.
TEST( TurbSstResidual, SstSourceKCurveMatchesClosedForm )
{
    const Real eps = 2.0;
    const Real prods[] = { 0.0, 1.0, 2.0, 3.0, 5.0, 10.0 };
    for ( Real p : prods )
    {
        EXPECT_NEAR( SstSourceK( p, eps ), p - eps, 1.0e-12 );
    }
}

TEST( TurbSstResidual, SstSourceKAssembledFromDissipationAndLimiter )
{
    const Real eps = SstDisskStandard( 1.0, 0.09, 1.0, 2.0, 3.0 );
    const Real pLim = SstLimitProdk( 20.0, 10.0, eps );
    EXPECT_NEAR( SstSourceK( pLim, eps ), pLim - eps, 1.0e-12 );
}

// Intermediate: SstSourceW on a cross-diffusion grid.
// Ref: Pw=4, Dw=1; Sw = 3 + CD.
TEST( TurbSstResidual, SstSourceWCurveMatchesClosedForm )
{
    const Real prodw = 4.0;
    const Real dissw = 1.0;
    const Real cds[] = { -1.0, -0.5, 0.0, 0.5, 1.0, 2.0 };
    for ( Real cd : cds )
    {
        EXPECT_NEAR(
            SstSourceW( prodw, dissw, cd ),
            prodw - dissw + cd,
            1.0e-12 );
    }
}

TEST( TurbSstResidual, SstSourceWAssembledFromDissw )
{
    const Real Dw = SstDissw( 1.0, 0.075, 1.0, 2.0 );
    // Dw = 0.3; Pw=1, CD=0 => Sw = 0.7
    EXPECT_NEAR( SstSourceW( 1.0, Dw, 0.0 ), 1.0 - Dw, 1.0e-12 );
}

// Intermediate: SstDiagonalK on an omega grid.
// Ref: fBetas=1, betas=0.09; dk = -0.18 * omega.
TEST( TurbSstResidual, SstDiagonalKCurveMatchesClosedForm )
{
    const Real fBetas = 1.0;
    const Real betas = 0.09;
    const Real omegas[] = { 0.5, 1.0, 2.0, 3.0, 5.0, 10.0 };
    for ( Real omega : omegas )
    {
        EXPECT_NEAR(
            SstDiagonalK( fBetas, betas, omega ),
            -2.0 * fBetas * betas * omega,
            1.0e-12 );
    }
}

TEST( TurbSstResidual, SstDiagonalKAlwaysNonPositive )
{
    EXPECT_LE( SstDiagonalK( 1.0, 0.09, 1.0 ), 0.0 );
    EXPECT_LE( SstDiagonalK( 1.0, 0.09, 10.0 ), 0.0 );
}

// Intermediate: SstDiagonalW on an omega grid (no cross-diffusion).
// Ref: fBeta=1, beta=0.075, rho=1, |CD|=0; dw = -0.15 * omega.
TEST( TurbSstResidual, SstDiagonalWCurveMatchesClosedForm )
{
    const Real fBeta = 1.0;
    const Real beta = 0.075;
    const Real density = 1.0;
    const Real omegas[] = { 0.5, 1.0, 2.0, 3.0, 5.0, 10.0 };
    for ( Real omega : omegas )
    {
        EXPECT_NEAR(
            SstDiagonalW( fBeta, beta, omega, 0.0, density ),
            -2.0 * fBeta * beta * omega,
            1.0e-12 );
    }
}

TEST( TurbSstResidual, SstDiagonalWAddsCrossDiffTerm )
{
    // -2*1*0.075*2 - 0.4/(1*2) = -0.3 - 0.2 = -0.5
    EXPECT_NEAR(
        SstDiagonalW( 1.0, 0.075, 2.0, 0.4, 1.0 ),
        -0.5,
        1.0e-12 );
}

// Intermediate: SstNegativeSourcePart on a source grid.
// n = 0.5 * (S - |S|) = min(S, 0).
TEST( TurbSstResidual, SstNegativeSourcePartCurveMatchesClosedForm )
{
    const Real sources[] = { -4.0, -1.0, 0.0, 1.0, 3.0, 10.0 };
    for ( Real s : sources )
    {
        const Real expected = ( s < 0.0 ) ? s : 0.0;
        EXPECT_NEAR( SstNegativeSourcePart( s ), expected, 1.0e-12 );
    }
}

TEST( TurbSstResidual, SstNegativeSourcePartFromSourceK )
{
    const Real sk = SstSourceK( 1.0, 3.0 ); // -2
    EXPECT_NEAR( SstNegativeSourcePart( sk ), -2.0, 1.0e-12 );
    EXPECT_NEAR( SstNegativeSourcePart( SstSourceK( 5.0, 1.0 ) ), 0.0, 1.0e-12 );
}
