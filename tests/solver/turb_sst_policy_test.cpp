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

// Intermediate: SstDissw on an omega grid.
// Ref: fBeta=1, beta=0.075, rho=1; Dw = 0.075 * omega^2.
TEST( TurbSstPolicy, SstDisswCurveMatchesClosedForm )
{
    const Real fBeta = 1.0;
    const Real beta = 0.075;
    const Real density = 1.0;
    const Real omegas[] = { 0.5, 1.0, 2.0, 3.0, 5.0, 10.0 };
    for ( Real omega : omegas )
    {
        EXPECT_NEAR(
            SstDissw( fBeta, beta, density, omega ),
            fBeta * beta * density * omega * omega,
            1.0e-12 );
    }
}

TEST( TurbSstPolicy, SstDisswQuadraticInOmega )
{
    const Real d1 = SstDissw( 1.0, 0.075, 1.0, 2.0 );
    const Real d2 = SstDissw( 1.0, 0.075, 1.0, 4.0 );
    EXPECT_NEAR( d2, 4.0 * d1, 1.0e-12 );
}

// Intermediate: SstCrossDiffMenter on an omega grid.
// Ref: F1=0, rho=1, sigw2=0.856, cross=1; CD = 1.712 / omega.
TEST( TurbSstPolicy, SstCrossDiffMenterCurveMatchesClosedForm )
{
    const Real f1 = 0.0;
    const Real density = 1.0;
    const Real sigw2 = 0.856;
    const Real cross = 1.0;
    const Real omegas[] = { 0.5, 1.0, 2.0, 4.0, 8.0 };
    for ( Real omega : omegas )
    {
        EXPECT_NEAR(
            SstCrossDiffMenter( f1, density, sigw2, cross, omega ),
            2.0 * ( 1.0 - f1 ) * density * sigw2 * cross / omega,
            1.0e-12 );
    }
}

TEST( TurbSstPolicy, SstCrossDiffMenterVanishesWhenF1IsOne )
{
    EXPECT_NEAR(
        SstCrossDiffMenter( 1.0, 1.0, 0.856, 10.0, 2.0 ),
        0.0,
        1.0e-12 );
}

// Intermediate: SstProdwMenter on a Pk grid.
// Ref: alphaw=0.44, rho=1, mut=0.1, Re=1; Pw = 4.4 * Pk.
TEST( TurbSstPolicy, SstProdwMenterCurveMatchesClosedForm )
{
    const Real alphaw = 0.44;
    const Real density = 1.0;
    const Real mut = 0.1;
    const Real reynolds = 1.0;
    const Real pks[] = { 0.0, 0.5, 1.0, 2.0, 5.0 };
    for ( Real pk : pks )
    {
        EXPECT_NEAR(
            SstProdwMenter( alphaw, density, mut, pk, reynolds ),
            alphaw * density / mut * pk * reynolds,
            1.0e-12 );
    }
}

TEST( TurbSstPolicy, SstProdwMenterLinearInProdk )
{
    const Real p1 = SstProdwMenter( 0.44, 1.0, 0.1, 1.0, 1.0 );
    const Real p2 = SstProdwMenter( 0.44, 1.0, 0.1, 2.0, 1.0 );
    EXPECT_NEAR( p2, 2.0 * p1, 1.0e-12 );
}

// Intermediate: SstProdkBoussinesq on a strain-rate-squared grid (div=0).
// Ref: mut=2, invRe=1, rho=1, k=1, div=0; Pk = 2 * S2.
TEST( TurbSstPolicy, SstProdkBoussinesqCurveMatchesClosedForm )
{
    const Real mut = 2.0;
    const Real invRe = 1.0;
    const Real s2s[] = { 0.0, 0.5, 1.0, 2.0, 4.0 };
    for ( Real s2 : s2s )
    {
        EXPECT_NEAR(
            SstProdkBoussinesq( mut, s2, 0.0, invRe, 1.0, 1.0 ),
            mut * s2 * invRe,
            1.0e-12 );
    }
}

TEST( TurbSstPolicy, SstProdkBoussinesqSubtractsDivergenceTerm )
{
    // mut=0, div=3, rho=1, k=2 => Pk = -2/3 * 1 * 2 * 3 = -4
    EXPECT_NEAR(
        SstProdkBoussinesq( 0.0, 0.0, 3.0, 1.0, 1.0, 2.0 ),
        -4.0,
        1.0e-12 );
}

// Intermediate: SstProdkStrainOnly on a strain-rate-squared grid.
// Ref: mut=2, invRe=0.5; Pk = S2.
TEST( TurbSstPolicy, SstProdkStrainOnlyCurveMatchesClosedForm )
{
    const Real mut = 2.0;
    const Real invRe = 0.5;
    const Real s2s[] = { 0.0, 0.5, 1.0, 2.0, 4.0, 8.0 };
    for ( Real s2 : s2s )
    {
        EXPECT_NEAR(
            SstProdkStrainOnly( mut, s2, invRe ),
            mut * s2 * invRe,
            1.0e-12 );
    }
}

TEST( TurbSstPolicy, SstProdkStrainOnlyLinearInEddyViscosity )
{
    EXPECT_NEAR(
        SstProdkStrainOnly( 4.0, 3.0, 0.5 ),
        2.0 * SstProdkStrainOnly( 2.0, 3.0, 0.5 ),
        1.0e-12 );
}

// Intermediate: SstProdkVorticity on a vorticity-squared grid (div=0).
// Ref: mut=2, invRe=1, rho=1, k=1, div=0; Pk = 2 * W2.
TEST( TurbSstPolicy, SstProdkVorticityCurveMatchesClosedForm )
{
    const Real mut = 2.0;
    const Real invRe = 1.0;
    const Real w2s[] = { 0.0, 0.5, 1.0, 2.0, 4.0 };
    for ( Real w2 : w2s )
    {
        EXPECT_NEAR(
            SstProdkVorticity( mut, w2, invRe, 1.0, 1.0, 0.0 ),
            mut * w2 * invRe,
            1.0e-12 );
    }
}

TEST( TurbSstPolicy, SstProdkVorticitySubtractsDivergenceTerm )
{
    // mut=0, div=3, rho=1, k=2 => Pk = -2/3 * 1 * 2 * 3 = -4
    EXPECT_NEAR(
        SstProdkVorticity( 0.0, 0.0, 1.0, 1.0, 2.0, 3.0 ),
        -4.0,
        1.0e-12 );
}
