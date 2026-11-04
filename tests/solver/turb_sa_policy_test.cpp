#include <gtest/gtest.h>
#include <cmath>

#include "TurbSaPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSaPolicy, SaChiIsProduct )
{
    EXPECT_NEAR( SaChi( 2.0, 0.5 ), 1.0 + SMALL, 1.0e-12 );
}

TEST( TurbSaPolicy, SaFv1AtLargeChiApproachesOne )
{
    const Real chi = 100.0;
    const Real cv13 = 7.1 * 7.1 * 7.1;
    EXPECT_NEAR( SaFv1( chi, cv13 ), 1.0, 1.0e-3 );
}

TEST( TurbSaPolicy, SaFv2AtZeroChiIsOne )
{
    EXPECT_DOUBLE_EQ( SaFv2( 0.0, 0.0 ), 1.0 );
}

TEST( TurbSaPolicy, SaModifiedStrainSimpleBranch )
{
    EXPECT_DOUBLE_EQ( SaModifiedStrainRate( 2.0, 1.0, 0.7, 0.9 ), 3.0 );
}

TEST( TurbSaPolicy, SaRRatioIsClamped )
{
    EXPECT_DOUBLE_EQ( SaRRatio( 100.0, 1.0 ), 10.0 );
}

TEST( TurbSaPolicy, SaGAndFwAtRZero )
{
    EXPECT_DOUBLE_EQ( SaG( 0.0, 0.3 ), 0.0 );
}

TEST( TurbSaPolicy, SaFt2DisabledIsZero )
{
    EXPECT_DOUBLE_EQ( SaFt2( 1.0, 1.2, 0.5, false ), 0.0 );
}

TEST( TurbSaPolicy, SaEddyViscosityIsRhoNuFv1 )
{
    EXPECT_DOUBLE_EQ( SaEddyViscosity( 1.2, 0.5, 0.8 ), 0.48 );
}

// Source-term components are checked against independently evaluated
// closed-form reference values from the equations documented in sa_source.rst.
TEST( TurbSaPolicy, SaProductionMatchesClosedForm )
{
    EXPECT_NEAR( SaProduction( 0.1355, 0.2, 3.0, 0.4 ), 0.13008, 1.0e-14 );
}

TEST( TurbSaPolicy, SaDestructionMatchesClosedForm )
{
    EXPECT_NEAR( SaDestruction( 3.0, 0.5, 0.1355, 0.2, 0.4 ), 0.58916, 1.0e-14 );
}

TEST( TurbSaPolicy, SaDiffusionMatchesClosedForm )
{
    EXPECT_NEAR( SaDiffusion( 1.3, 0.25, 0.5 ), 0.1625, 1.0e-14 );
}

TEST( TurbSaPolicy, SaLimitedDiffusionSelectsTheSmallerContribution )
{
    EXPECT_DOUBLE_EQ( SaLimitedDiffusion( 0.4, 0.5, 0.3 ), 0.15 );
    EXPECT_DOUBLE_EQ( SaLimitedDiffusion( 0.1, 0.5, 0.3 ), 0.1 );
}

TEST( TurbSaPolicy, SaSourceResidualCombinesProductionDestructionAndDiffusion )
{
    EXPECT_NEAR( SaSourceResidual( 0.13008, 0.58916, 0.1625 ), -0.29658, 1.0e-14 );
}

TEST( TurbSaPolicy, SaNegativePartMatchesPiecewiseDefinition )
{
    EXPECT_DOUBLE_EQ( SaNegativePart( 2.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( SaNegativePart( 0.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( SaNegativePart( -2.0 ), 2.0 );
}

// Intermediate verification: fv1(chi) vs closed-form theory (not a full CFD case).
// cv1 = 7.1 is the standard SA constant (NASA TMR / Spalart-Allmaras).
TEST( TurbSaPolicy, Fv1ChiCurveMatchesClosedForm )
{
    const Real cv1 = 7.1;
    const Real cv13 = cv1 * cv1 * cv1;
    const Real chis[] = { 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0 };
    for ( Real chi : chis )
    {
        const Real chi3 = chi * chi * chi;
        const Real expected = chi3 / ( chi3 + cv13 );
        EXPECT_NEAR( SaFv1( chi, cv13 ), expected, 1.0e-12 );
    }
}

// Intermediate verification: fw(r) vs published closed form (cw2=0.3, cw3=2).
// Proves pure SaG/SaFw match the SA wall-destruction function algebra.
// Does not prove full transport residual or wall-resolved CFD accuracy.
TEST( TurbSaPolicy, FwRCurveMatchesClosedForm )
{
    const Real cw2 = 0.3;
    const Real cw3 = 2.0;
    const Real cw36 = std::pow( cw3, 6 );
    const Real rs[] = { 0.0, 0.25, 0.5, 1.0, 2.0, 10.0 };
    for ( Real r : rs )
    {
        const Real g = r + cw2 * ( std::pow( r, 6 ) - r );
        const Real expectedFw =
            g * std::pow( ( one + cw36 ) / ( std::pow( g, 6 ) + cw36 ), one / six );
        EXPECT_NEAR( SaG( r, cw2 ), g, 1.0e-12 );
        EXPECT_NEAR( SaFw( g, cw36 ), expectedFw, 1.0e-12 );
    }
}

// Intermediate: fv2(chi) with standard cv1=7.1 via fv1(chi).
// Proves SaFv2 algebra; does not prove full modified-strain residual.
TEST( TurbSaPolicy, Fv2ChiCurveMatchesClosedForm )
{
    const Real cv1 = 7.1;
    const Real cv13 = cv1 * cv1 * cv1;
    const Real chis[] = { 0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0 };
    for ( Real chi : chis )
    {
        const Real fv1 = SaFv1( chi, cv13 );
        const Real expected = one - chi / ( one + chi * fv1 );
        EXPECT_NEAR( SaFv2( chi, fv1 ), expected, 1.0e-12 );
    }
}

// Modified strain: regular branch and limited branch (c2=0.7, c3=0.9).
TEST( TurbSaPolicy, ModifiedStrainRegularBranch )
{
    const Real omega = 2.0;
    const Real sBar = -0.5; // > -c2*omega = -1.4
    EXPECT_DOUBLE_EQ( SaModifiedStrainRate( omega, sBar, 0.7, 0.9 ), 1.5 );
}

TEST( TurbSaPolicy, ModifiedStrainLimitedBranch )
{
    const Real omega = 1.0;
    const Real sBar = -1.0; // < -c2*omega = -0.7
    const Real c2 = 0.7;
    const Real c3 = 0.9;
    const Real term1 = c2 * c2 * omega + c3 * sBar;
    const Real term2 = ( c3 - two * c2 ) * omega - sBar;
    const Real expected = omega + omega * term1 / term2;
    EXPECT_NEAR( SaModifiedStrainRate( omega, sBar, c2, c3 ), expected, 1.0e-12 );
}

TEST( TurbSaPolicy, CompleteCellClosureReportsIntermediateValues )
{
    SaEddyViscosityInput input;
    input.density = 1.2;
    input.workingVariable = 0.5;
    input.laminarViscosity = 0.3;
    input.cv1Cubed = 7.1 * 7.1 * 7.1;
    input.maxViscosityRatio = 100.0;

    const SaEddyViscosityResult result = SaComputeEddyViscosity( input );
    const Real expectedChi = 0.5 * ( 1.2 / ( 0.3 + SMALL ) );
    const Real expectedFv1 =
        expectedChi * expectedChi * expectedChi /
        ( expectedChi * expectedChi * expectedChi + input.cv1Cubed );
    const Real expectedUncapped = 1.2 * 0.5 * expectedFv1;

    EXPECT_NEAR( result.chi, expectedChi, 1.0e-14 );
    EXPECT_NEAR( result.fv1, expectedFv1, 1.0e-14 );
    EXPECT_NEAR( result.uncappedEddyViscosity, expectedUncapped, 1.0e-14 );
    EXPECT_NEAR( result.eddyViscosity, expectedUncapped, 1.0e-14 );
}

TEST( TurbSaPolicy, CompleteCellClosureAppliesSharedViscosityCap )
{
    SaEddyViscosityInput input;
    input.density = 1.0;
    input.workingVariable = 100.0;
    input.laminarViscosity = 0.5;
    input.cv1Cubed = 7.1 * 7.1 * 7.1;
    input.maxViscosityRatio = 2.0;

    const SaEddyViscosityResult result = SaComputeEddyViscosity( input );
    EXPECT_GT( result.uncappedEddyViscosity, 1.0 );
    EXPECT_DOUBLE_EQ( result.eddyViscosity, 1.0 );
}

// Intermediate: SaProduction on a nu_tilde grid.
// Ref: cb1=0.1355, ft2=0.2, S=3 => P = 0.3252 * nu_tilde.
TEST( TurbSaPolicy, SaProductionCurveMatchesClosedForm )
{
    const Real cb1 = 0.1355;
    const Real ft2 = 0.2;
    const Real S = 3.0;
    const Real nus[] = { 0.0, 0.1, 0.4, 1.0, 2.0, 5.0 };
    for ( Real nu : nus )
    {
        EXPECT_NEAR(
            SaProduction( cb1, ft2, S, nu ),
            cb1 * ( 1.0 - ft2 ) * S * nu,
            1.0e-12 );
    }
}

TEST( TurbSaPolicy, SaProductionVanishesWhenFt2IsOne )
{
    EXPECT_NEAR( SaProduction( 0.1355, 1.0, 3.0, 10.0 ), 0.0, 1.0e-12 );
}

// Intermediate: SaDestruction on a nu_tilde_rs grid.
// Ref: cw1k=3, fw=0.5, cb1=0.1355, ft2=0.2.
TEST( TurbSaPolicy, SaDestructionCurveMatchesClosedForm )
{
    const Real cw1k = 3.0;
    const Real fw = 0.5;
    const Real cb1 = 0.1355;
    const Real ft2 = 0.2;
    const Real nus[] = { 0.0, 0.1, 0.4, 1.0, 2.0, 5.0 };
    for ( Real nu : nus )
    {
        EXPECT_NEAR(
            SaDestruction( cw1k, fw, cb1, ft2, nu ),
            ( cw1k * fw - cb1 * ft2 ) * nu,
            1.0e-12 );
    }
}

TEST( TurbSaPolicy, SaDestructionCanBeNegativeWhenFt2Dominates )
{
    // cw1k*fw = 0.1, cb1*ft2 = 0.1355 => negative prefactor
    EXPECT_LT( SaDestruction( 0.2, 0.5, 0.1355, 1.0, 1.0 ), 0.0 );
}

// Intermediate: SaSourceResidual on a diffusion grid.
// Ref: P=0.13008, D=0.58916; S = P - D + Diff.
TEST( TurbSaPolicy, SaSourceResidualCurveMatchesClosedForm )
{
    const Real P = 0.13008;
    const Real D = 0.58916;
    const Real diffs[] = { 0.0, 0.1, 0.1625, 0.5, 1.0 };
    for ( Real diff : diffs )
    {
        EXPECT_NEAR(
            SaSourceResidual( P, D, diff ),
            P - D + diff,
            1.0e-12 );
    }
}

TEST( TurbSaPolicy, SaSourceResidualAssembledFromPieces )
{
    const Real P = SaProduction( 0.1355, 0.2, 3.0, 0.4 );
    const Real D = SaDestruction( 3.0, 0.5, 0.1355, 0.2, 0.4 );
    const Real Diff = SaLimitedDiffusion(
        SaDiffusion( 1.3, 0.25, 0.5 ), 0.5, P );
    EXPECT_NEAR(
        SaSourceResidual( P, D, Diff ),
        P - D + Diff,
        1.0e-12 );
}

// Intermediate: SaDiffusion on a |grad nu|^2 grid.
// Ref: cb2OverSigma=1.3, invRe=0.5; Diff = 0.65 * g2.
TEST( TurbSaPolicy, SaDiffusionCurveMatchesClosedForm )
{
    const Real cb2OverSigma = 1.3;
    const Real invRe = 0.5;
    const Real g2s[] = { 0.0, 0.25, 1.0, 2.0, 4.0 };
    for ( Real g2 : g2s )
    {
        EXPECT_NEAR(
            SaDiffusion( cb2OverSigma, g2, invRe ),
            cb2OverSigma * g2 * invRe,
            1.0e-12 );
    }
}

TEST( TurbSaPolicy, SaDiffusionLinearInInvReynolds )
{
    EXPECT_NEAR(
        SaDiffusion( 1.3, 0.25, 1.0 ),
        2.0 * SaDiffusion( 1.3, 0.25, 0.5 ),
        1.0e-12 );
}

// Intermediate: SaLimitedDiffusion on a diffusion grid.
// Ref: rprod=0.5, P=0.4 => cap=0.2.
TEST( TurbSaPolicy, SaLimitedDiffusionCurveMatchesClosedForm )
{
    const Real rprod = 0.5;
    const Real production = 0.4;
    const Real cap = rprod * production;
    const Real diffs[] = { 0.0, 0.05, 0.1, 0.2, 0.4, 1.0 };
    for ( Real d : diffs )
    {
        EXPECT_NEAR(
            SaLimitedDiffusion( d, rprod, production ),
            ( d < cap ? d : cap ),
            1.0e-12 );
    }
}

TEST( TurbSaPolicy, SaLimitedDiffusionUsesSaDiffusionInput )
{
    const Real raw = SaDiffusion( 1.3, 0.25, 0.5 ); // 0.1625
    // rprod*P = 0.5*0.1 = 0.05 < 0.1625 => limited to 0.05
    EXPECT_NEAR( SaLimitedDiffusion( raw, 0.5, 0.1 ), 0.05, 1.0e-12 );
}

// Intermediate: SaEddyViscosity on a nu_tilde grid.
// Ref: rho=1.2, fv1=0.8; mut = 0.96 * nu.
TEST( TurbSaPolicy, SaEddyViscosityCurveMatchesClosedForm )
{
    const Real density = 1.2;
    const Real fv1 = 0.8;
    const Real nus[] = { 0.0, 0.25, 0.5, 1.0, 2.0, 5.0 };
    for ( Real nu : nus )
    {
        EXPECT_NEAR(
            SaEddyViscosity( density, nu, fv1 ),
            density * nu * fv1,
            1.0e-12 );
    }
}

TEST( TurbSaPolicy, SaEddyViscosityComposedFromFv1 )
{
    const Real chi = SaChi( 2.0, 0.5 );
    const Real fv1 = SaFv1( chi, 7.1 * 7.1 * 7.1 );
    EXPECT_NEAR(
        SaEddyViscosity( 1.0, 2.0, fv1 ),
        1.0 * 2.0 * fv1,
        1.0e-12 );
}
