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
