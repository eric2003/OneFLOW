#include <gtest/gtest.h>
#include <cmath>

#include "TurbSstPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSstPolicy, ProdkStrainOnly )
{
    EXPECT_DOUBLE_EQ( SstProdkStrainOnly( 2.0, 3.0, 0.5 ), 3.0 );
}

TEST( TurbSstPolicy, ProdkBoussinesqZeroDiv )
{
    EXPECT_DOUBLE_EQ(
        SstProdkBoussinesq( 2.0, 4.0, 0.0, 1.0, 1.0, 1.0 ), 8.0 );
}

TEST( TurbSstPolicy, DisskStandard )
{
    EXPECT_DOUBLE_EQ(
        SstDisskStandard( 1.0, 0.09, 1.0, 2.0, 3.0 ), 0.09 * 2.0 * 3.0 );
}

TEST( TurbSstPolicy, LimitProdkCaps )
{
    EXPECT_DOUBLE_EQ( SstLimitProdk( 10.0, 2.0, 3.0 ), 6.0 );
}

TEST( TurbSstPolicy, CrossDiffScalar )
{
    EXPECT_DOUBLE_EQ(
        SstCrossDiffScalar( 1.0, 0.5, 2.0, 1.0 ), 2.0 );
}

TEST( TurbSstPolicy, CdkwWithFloor )
{
    EXPECT_DOUBLE_EQ( SstCdkwWithFloor( 1.0e-30, 1.0e-20 ), 1.0e-20 );
}

TEST( TurbSstPolicy, BlendingArg1TakesMinOfMax )
{
    EXPECT_DOUBLE_EQ( SstBlendingArg1( 1.0, 3.0, 2.0 ), 2.0 );
    EXPECT_DOUBLE_EQ( SstBlendingArg1( 1.0, 0.5, 2.0 ), 1.0 );
}

TEST( TurbSstPolicy, F1FromArg1AtZeroIsZero )
{
    EXPECT_DOUBLE_EQ( SstF1FromArg1( 0.0 ), 0.0 );
}

TEST( TurbSstPolicy, F1FromArg1LargeApproachesOne )
{
    EXPECT_NEAR( SstF1FromArg1( 2.0 ), 1.0, 1.0e-6 );
}

TEST( TurbSstPolicy, GradientDotProduct )
{
    EXPECT_DOUBLE_EQ(
        SstGradientDotProduct( 1.0, 0.0, 0.0, 2.0, 0.0, 0.0 ), 2.0 );
}

TEST( TurbSstPolicy, BlendIsConvexCombination )
{
    // phi = F1*phi1 + (1-F1)*phi2
    EXPECT_DOUBLE_EQ( SstBlend( 0.25, 4.0, 8.0 ), 7.0 );
    EXPECT_DOUBLE_EQ( SstBlend( 1.0, 4.0, 8.0 ), 4.0 );
    EXPECT_DOUBLE_EQ( SstBlend( 0.0, 4.0, 8.0 ), 8.0 );
}

// Intermediate verification: F1(arg1) vs closed form tanh(arg1^4).
TEST( TurbSstPolicy, F1Arg1CurveMatchesClosedForm )
{
    const Real args[] = { 0.0, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0 };
    for ( Real a : args )
    {
        const Real expected = std::tanh( std::pow( a, 4 ) );
        EXPECT_NEAR( SstF1FromArg1( a ), expected, 1.0e-12 );
    }
}

// Intermediate: SstDisskStandard on an omega grid.
// Ref: fBetas=1, betas=0.09, rho=1, k=2; eps = 0.18 * omega.
TEST( TurbSstPolicy, SstDisskStandardCurveMatchesClosedForm )
{
    const Real fBetas = 1.0;
    const Real betas = 0.09;
    const Real density = 1.0;
    const Real k = 2.0;
    const Real omegas[] = { 0.5, 1.0, 2.0, 3.0, 5.0, 10.0 };
    for ( Real omega : omegas )
    {
        EXPECT_NEAR(
            SstDisskStandard( fBetas, betas, density, k, omega ),
            fBetas * betas * density * k * omega,
            1.0e-12 );
    }
}

TEST( TurbSstPolicy, SstDisskStandardLinearInKineticEnergy )
{
    EXPECT_NEAR(
        SstDisskStandard( 1.0, 0.09, 1.0, 4.0, 2.0 ),
        2.0 * SstDisskStandard( 1.0, 0.09, 1.0, 2.0, 2.0 ),
        1.0e-12 );
}

// Intermediate: SstLimitProdk on a production grid.
// Ref: pklim=2, eps=3 => cap=6.
TEST( TurbSstPolicy, SstLimitProdkCurveMatchesClosedForm )
{
    const Real pklim = 2.0;
    const Real eps = 3.0;
    const Real cap = pklim * eps;
    const Real prods[] = { 0.0, 1.0, 3.0, 6.0, 8.0, 12.0 };
    for ( Real p : prods )
    {
        EXPECT_NEAR(
            SstLimitProdk( p, pklim, eps ),
            ( p < cap ? p : cap ),
            1.0e-12 );
    }
}

TEST( TurbSstPolicy, SstLimitProdkUsesDissipationFromSstDisskStandard )
{
    const Real eps = SstDisskStandard( 1.0, 0.09, 1.0, 2.0, 3.0 );
    // eps = 0.54; pklim=10 => cap=5.4; production 20 limited to 5.4
    EXPECT_NEAR( SstLimitProdk( 20.0, 10.0, eps ), 10.0 * eps, 1.0e-12 );
}
