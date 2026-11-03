#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>

#include "NsInvFluxPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( InvFluxPolicy, RelativeFaceNormalVelocitySubtractsGridSpeed )
{
    EXPECT_DOUBLE_EQ(
        RelativeFaceNormalVelocity( 1.0, 0.0, 0.0, 0.5, 2.0, 0.0, 0.0 ), 1.5 );
}

TEST( InvFluxPolicy, EulerConvectiveFluxMatchesMassAndPressureTerms )
{
    Real fr, fx, fy, fz, fe;
    EulerConvectiveFlux(
        1.0, 2.0, 0.0, 0.0, 3.0, 10.0,
        2.0, 1.0, 0.0, 0.0, 0.0,
        fr, fx, fy, fz, fe );
    EXPECT_DOUBLE_EQ( fr, 2.0 );
    EXPECT_DOUBLE_EQ( fx, 2.0 * 2.0 + 3.0 );
    EXPECT_DOUBLE_EQ( fy, 0.0 );
    EXPECT_DOUBLE_EQ( fz, 0.0 );
    EXPECT_DOUBLE_EQ( fe, 2.0 * 10.0 );
}

TEST( InvFluxPolicy, EulerConvectiveFluxZeroNormalVelocityIsPressureOnly )
{
    Real fr, fx, fy, fz, fe;
    EulerConvectiveFlux(
        1.2, 0.0, 0.0, 0.0, 5.0, 8.0,
        0.0, 0.0, 1.0, 0.0, 0.0,
        fr, fx, fy, fz, fe );
    EXPECT_DOUBLE_EQ( fr, 0.0 );
    EXPECT_DOUBLE_EQ( fx, 0.0 );
    EXPECT_DOUBLE_EQ( fy, 5.0 );
    EXPECT_DOUBLE_EQ( fz, 0.0 );
    EXPECT_DOUBLE_EQ( fe, 0.0 );
}

TEST( InvFluxPolicy, RoeAverageWeightsMatchDensityRatio )
{
    Real ratio = 0.0, coef = 0.0;
    RoeAverageWeights( 1.0, 4.0, ratio, coef );
    EXPECT_DOUBLE_EQ( ratio, 2.0 );
    EXPECT_DOUBLE_EQ( coef, 1.0 / 3.0 );
}

TEST( InvFluxPolicy, RoeAverageComponentIsWeightedMean )
{
    EXPECT_DOUBLE_EQ( RoeAverageComponent( 1.0, 4.0, 2.0, 1.0 / 3.0 ), 3.0 );
}

TEST( InvFluxPolicy, RoeAverageStateUsesGeometricDensity )
{
    Real rm, um, vm, wm, pm, hm, gama;
    RoeAverageState(
        1.0, 0.0, 0.0, 0.0, 1.0, 3.5, 1.4,
        4.0, 3.0, 0.0, 0.0, 4.0, 7.0, 1.4,
        rm, um, vm, wm, pm, hm, gama );
    EXPECT_DOUBLE_EQ( rm, 2.0 );
    EXPECT_DOUBLE_EQ( um, 2.0 );
    EXPECT_DOUBLE_EQ( gama, 1.4 );
}

TEST( InvFluxPolicy, TotalEnthalpyAddsKinetic )
{
    EXPECT_DOUBLE_EQ( TotalEnthalpy( 3.0, 8.0 ), 7.0 );
}

TEST( InvFluxPolicy, RoePressureUsesTotalEnthalpyClosure )
{
    // Legacy UNsLusgs face state: rho=2, H=6.5, |u|^2=4, gamma=1.4.
    // This gives 18/7; the weighted static pressure for the same pair is 3.
    EXPECT_NEAR(
        RoePressureFromTotalEnthalpy( 2.0, 6.5, 4.0, 1.4 ),
        18.0 / 7.0,
        1.0e-14 );
}

TEST( InvFluxPolicy, LusgsFacePressureModeSelectsLegacyEnthalpyClosure )
{
    EXPECT_NEAR(
        LusgsFacePressureFromMode( 0, 3.0, 2.0, 6.5, 4.0, 1.4 ),
        18.0 / 7.0,
        1.0e-14 );
}

TEST( InvFluxPolicy, LusgsFacePressureModeSelectsRoeWeightedPressure )
{
    EXPECT_DOUBLE_EQ(
        LusgsFacePressureFromMode( 1, 3.0, 2.0, 6.5, 4.0, 1.4 ),
        3.0 );
}

TEST( InvFluxPolicy, LusgsFacePressureModeRejectsUnknownMode )
{
    EXPECT_THROW(
        LusgsFacePressureFromMode( 2, 3.0, 2.0, 6.5, 4.0, 1.4 ),
        std::invalid_argument );
}

TEST( InvFluxPolicy, IdealGasSoundSpeedMatchesIsentropicRelation )
{
    EXPECT_NEAR( IdealGasSoundSpeed( 1.4, 3.5, 0.0 ), std::sqrt( 1.4 ), 1.0e-12 );
}

TEST( InvFluxPolicy, RoeAbsoluteEigenvaluesTakeAbsoluteValues )
{
    Real a1, a2, a3;
    RoeAbsoluteEigenvalues( -2.0, 1.0, a1, a2, a3 );
    EXPECT_DOUBLE_EQ( a1, 2.0 );
    EXPECT_DOUBLE_EQ( a2, 1.0 );
    EXPECT_DOUBLE_EQ( a3, 3.0 );
}

TEST( InvFluxPolicy, HartenEntropyFixRaisesFloorWhenEnabled )
{
    const Real absVn = 0.1;
    const Real absVnPlusC = 2.0;
    const Real absVnMinusC = 1.5;
    const Real centropy1 = 0.5;
    const Real centropy2 = 0.25;

    Real f1, f2, f3;
    ApplyHartenEntropyFix(
        absVn, absVnPlusC, absVnMinusC, centropy1, centropy2, true, f1, f2, f3 );

    const Real maxEigenvalue = MAX( absVnPlusC, absVnMinusC );
    const Real m1 = maxEigenvalue * centropy1;
    const Real m2 = maxEigenvalue * centropy2;
    EXPECT_DOUBLE_EQ( f1, MAX( m1, absVn ) );
    EXPECT_DOUBLE_EQ( f2, MAX( m2, absVnPlusC ) );
    EXPECT_DOUBLE_EQ( f3, MAX( m2, absVnMinusC ) );
}

TEST( InvFluxPolicy, HartenEntropyFixIsNoOpWhenDisabled )
{
    Real f1, f2, f3;
    ApplyHartenEntropyFix( 0.1, 2.0, 1.5, 0.5, 0.25, false, f1, f2, f3 );
    EXPECT_DOUBLE_EQ( f1, 0.1 );
    EXPECT_DOUBLE_EQ( f2, 2.0 );
    EXPECT_DOUBLE_EQ( f3, 1.5 );
}

TEST( InvFluxPolicy, FMSplit1IsUpwindMach )
{
    EXPECT_DOUBLE_EQ( FMSplit1( 2.0, 1.0 ), 2.0 );
    EXPECT_DOUBLE_EQ( FMSplit1( -2.0, 1.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( FMSplit1( -2.0, -1.0 ), -2.0 );
    EXPECT_DOUBLE_EQ( FMSplit1( 2.0, -1.0 ), 0.0 );
}

TEST( InvFluxPolicy, FMSplit2IsQuadraticNearZero )
{
    EXPECT_DOUBLE_EQ( FMSplit2( 0.0, 1.0 ), 0.25 );
    EXPECT_DOUBLE_EQ( FMSplit2( 0.0, -1.0 ), -0.25 );
}

TEST( InvFluxPolicy, FMSplit4FallsBackOutsideSonic )
{
    EXPECT_DOUBLE_EQ( FMSplit4( 1.5, 0.125, 1.0 ), FMSplit1( 1.5, 1.0 ) );
    EXPECT_DOUBLE_EQ( FMSplit4( -1.5, 0.125, -1.0 ), FMSplit1( -1.5, -1.0 ) );
}

TEST( InvFluxPolicy, FMSplit4SubsonicMatchesPolynomialForm )
{
    const Real mach = 0.5;
    const Real beta = 0.125;
    const Real ipn = 1.0;
    const Real expected =
        FMSplit2( mach, ipn ) * ( one - ipn * sixteen * beta * FMSplit2( mach, - ipn ) );
    EXPECT_DOUBLE_EQ( FMSplit4( mach, beta, ipn ), expected );
}

TEST( InvFluxPolicy, FPSplit5SupersonicIsStep )
{
    EXPECT_DOUBLE_EQ( FPSplit5( 1.5, 0.1875, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( FPSplit5( -1.5, 0.1875, 1.0 ), 0.0 );
}

TEST( InvFluxPolicy, FPSplit5SubsonicMatchesPolynomialForm )
{
    const Real mach = 0.5;
    const Real alpha = 0.1875;
    const Real ipn = 1.0;
    const Real expected =
        FMSplit2( mach, ipn ) * ( ( two * ipn - mach ) - ipn * sixteen * alpha * mach * FMSplit2( mach, - ipn ) );
    EXPECT_DOUBLE_EQ( FPSplit5( mach, alpha, ipn ), expected );
}

TEST( InvFluxPolicy, IdealGasStaticEnthalpyMatchesGammaLaw )
{
    // gamma=1.4, p=1, rho=1 -> h = 3.5
    EXPECT_DOUBLE_EQ( IdealGasStaticEnthalpy( 1.0, 1.0, 1.4 ), 3.5 );
}

TEST( InvFluxPolicy, IdealGasSpecificInternalEnergyAddsKinetic )
{
    // e_int = 1/(0.4) = 2.5, + 0.5*4 = 2 -> 4.5
    EXPECT_DOUBLE_EQ(
        IdealGasSpecificInternalEnergy( 1.0, 1.0, 1.4, 4.0 ), 4.5 );
}

TEST( InvFluxPolicy, IdealGasSoundSpeedFromPressure )
{
    EXPECT_NEAR(
        IdealGasSoundSpeedFromPressure( 1.0, 1.0, 1.4 ),
        std::sqrt( 1.4 ),
        1.0e-12 );
}

TEST( InvFluxPolicy, IdealGasPrimToConservativeMatchesManual )
{
    Real qr, qu, qv, qw, qe;
    IdealGasPrimToConservative( 1.0, 2.0, 0.0, 0.0, 1.0, 1.4, qr, qu, qv, qw, qe );
    EXPECT_DOUBLE_EQ( qr, 1.0 );
    EXPECT_DOUBLE_EQ( qu, 2.0 );
    EXPECT_DOUBLE_EQ( qv, 0.0 );
    EXPECT_DOUBLE_EQ( qw, 0.0 );
    // e = 1/0.4 + 0.5*4 = 2.5+2 = 4.5, rho*e = 4.5
    EXPECT_DOUBLE_EQ( qe, 4.5 );
}

TEST( InvFluxPolicy, FaceMachNumberIsVnOverC )
{
    EXPECT_DOUBLE_EQ( FaceMachNumber( 3.0, 1.5 ), 2.0 );
}

TEST( InvFluxPolicy, IdealGasConservativeToPrimInvertsPrimToConservative )
{
    Real qr, qu, qv, qw, qe;
    IdealGasPrimToConservative( 2.0, 1.0, 0.0, 0.0, 4.0, 1.4, qr, qu, qv, qw, qe );
    Real rho, u, v, w, p;
    IdealGasConservativeToPrim( qr, qu, qv, qw, qe, 1.4, rho, u, v, w, p );
    EXPECT_NEAR( rho, 2.0, 1.0e-12 );
    EXPECT_NEAR( u, 1.0, 1.0e-12 );
    EXPECT_NEAR( p, 4.0, 1.0e-12 );
}

TEST( InvFluxPolicy, RoeDissipationWaveCoefficientsMatchFormula )
{
    Real xi1, xi2;
    RoeDissipationWaveCoefficients( 1.0, 3.0, 0.5, 4.0, 2.0, xi1, xi2 );
    // xi1 = (2*1 - 3 - 0.5) / 8 = (-1.5)/8
    EXPECT_DOUBLE_EQ( xi1, ( 2.0 * 1.0 - 3.0 - 0.5 ) / ( 2.0 * 4.0 ) );
    EXPECT_DOUBLE_EQ( xi2, ( 3.0 - 0.5 ) / ( 2.0 * 2.0 ) );
}

TEST( InvFluxPolicy, RoeContinuityNormalJumpMatchesDotProduct )
{
    EXPECT_DOUBLE_EQ(
        RoeContinuityNormalJump( 2.0, 1.0, 0.0, 0.0, 3.0, 4.0, 0.0, 0.0 ),
        2.0 * 3.0 - 4.0 );
}

TEST( InvFluxPolicy, RoeTotalEnthalpyJumpMatchesFormula )
{
    // u=1,v=0,w=0,gama=1.4; dq only energy +1 -> dh = -ae * (0 - 1) = ae = 0.4
    EXPECT_NEAR(
        RoeTotalEnthalpyJump( 1.0, 0.0, 0.0, 1.4, 0.0, 0.0, 0.0, 0.0, 1.0 ),
        0.4,
        1.0e-12 );
}

TEST( InvFluxPolicy, RoeDissipationCouplingTermsMatchFormula )
{
    Real t1, t2;
    RoeDissipationCouplingTerms( 1.0, 2.0, 4.0, 0.5, 0.25, t1, t2 );
    EXPECT_DOUBLE_EQ( t1, 1.0 * 0.5 + 2.0 * 0.25 );
    EXPECT_DOUBLE_EQ( t2, 4.0 * 2.0 * 0.5 + 1.0 * 0.25 );
}

TEST( InvFluxPolicy, ApplyRoeDissipationFiveEqSubtractsContinuityTerm )
{
    Real fr = 10.0, fx = 0.0, fy = 0.0, fz = 0.0, fe = 0.0;
    ApplyRoeDissipationFiveEq(
        1.0, 2.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        3.0, 0.0, 0.0, 0.0, 0.0,
        fr, fx, fy, fz, fe );
    // fr -= (1*3 - 2) = 1 -> 9
    EXPECT_DOUBLE_EQ( fr, 9.0 );
}

TEST( InvFluxPolicy, IdealGasBilateralSoundSpeedsAveragesSides )
{
    Real cl, cr, cm;
    IdealGasBilateralSoundSpeeds( 1.0, 1.0, 1.4, 1.0, 1.0, 1.4, cl, cr, cm );
    EXPECT_NEAR( cl, std::sqrt( 1.4 ), 1.0e-12 );
    EXPECT_NEAR( cr, std::sqrt( 1.4 ), 1.0e-12 );
    EXPECT_NEAR( cm, std::sqrt( 1.4 ), 1.0e-12 );
}


// Intermediate verification: AUSM FMSplit4 / FPSplit5 vs closed form on a Mach grid.
// Proves pure inline splits in NsInvFluxPolicy.h; does not prove full AUSM face flux.
TEST( InvFluxPolicy, FMSplit4MachCurveMatchesClosedForm )
{
    const Real beta = 0.125;
    const Real ipn = 1.0;
    const Real machs[] = { -2.0, -1.5, -1.0, -0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0 };
    for ( Real mach : machs )
    {
        Real expected = 0.0;
        if ( ABS( mach ) >= one )
        {
            expected = FMSplit1( mach, ipn );
        }
        else
        {
            expected = FMSplit2( mach, ipn )
                * ( one - ipn * sixteen * beta * FMSplit2( mach, - ipn ) );
        }
        EXPECT_NEAR( FMSplit4( mach, beta, ipn ), expected, 1.0e-12 );
    }
}

TEST( InvFluxPolicy, FPSplit5MachCurveMatchesClosedForm )
{
    const Real alpha = 0.1875;
    const Real ipn = 1.0;
    const Real machs[] = { -2.0, -1.5, -1.0, -0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0 };
    for ( Real mach : machs )
    {
        Real expected = 0.0;
        const Real macha = ABS( mach );
        if ( macha >= one )
        {
            expected = half * ( one + ipn * mach / ( macha + SMALL ) );
        }
        else
        {
            expected = FMSplit2( mach, ipn )
                * ( ( two * ipn - mach )
                    - ipn * sixteen * alpha * mach * FMSplit2( mach, - ipn ) );
        }
        EXPECT_NEAR( FPSplit5( mach, alpha, ipn ), expected, 1.0e-12 );
    }
}

TEST( InvFluxPolicy, FMSplit4NegativeIpnIsMirror )
{
    // ipn = -1 branch at M=0.5 should equal -FMSplit4(+0.5) with +ipn for mass-like split symmetry
    const Real beta = 0.125;
    EXPECT_NEAR(
        FMSplit4( 0.5, beta, -1.0 ),
        - FMSplit4( -0.5, beta, 1.0 ),
        1.0e-12 );
}

// Intermediate: Harten floors |lambda| on an |vn| grid (enabled).
// Ref: |vn+c|=2, |vn-c|=1.5, centropy1=0.5, centropy2=0.25 => maxEig=2, m1=1, m2=0.5.
TEST( InvFluxPolicy, HartenEntropyFixCurveMatchesClosedForm )
{
    const Real absVnPlusC = 2.0;
    const Real absVnMinusC = 1.5;
    const Real c1 = 0.5;
    const Real c2 = 0.25;
    const Real absVns[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0 };
    for ( Real absVn : absVns )
    {
        Real f1, f2, f3;
        ApplyHartenEntropyFix(
            absVn, absVnPlusC, absVnMinusC, c1, c2, true, f1, f2, f3 );
        const Real maxEig = MAX( absVnPlusC, absVnMinusC );
        const Real m1 = maxEig * c1;
        const Real m2 = maxEig * c2;
        EXPECT_NEAR( f1, MAX( m1, absVn ), 1.0e-12 );
        EXPECT_NEAR( f2, MAX( m2, absVnPlusC ), 1.0e-12 );
        EXPECT_NEAR( f3, MAX( m2, absVnMinusC ), 1.0e-12 );
    }
}

TEST( InvFluxPolicy, HartenEntropyFixDisabledLeavesInputs )
{
    const Real absVns[] = { 0.0, 0.1, 1.0, 5.0 };
    for ( Real absVn : absVns )
    {
        Real f1, f2, f3;
        ApplyHartenEntropyFix( absVn, 2.0, 1.5, 0.5, 0.25, false, f1, f2, f3 );
        EXPECT_DOUBLE_EQ( f1, absVn );
        EXPECT_DOUBLE_EQ( f2, 2.0 );
        EXPECT_DOUBLE_EQ( f3, 1.5 );
    }
}

// Intermediate: RoeAbsoluteEigenvalues on a vn grid with fixed sound speed.
TEST( InvFluxPolicy, RoeAbsoluteEigenvaluesCurveMatchesClosedForm )
{
    const Real c = 1.0;
    const Real vns[] = { -3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0 };
    for ( Real vn : vns )
    {
        Real a1, a2, a3;
        RoeAbsoluteEigenvalues( vn, c, a1, a2, a3 );
        EXPECT_NEAR( a1, ABS( vn ), 1.0e-12 );
        EXPECT_NEAR( a2, ABS( vn + c ), 1.0e-12 );
        EXPECT_NEAR( a3, ABS( vn - c ), 1.0e-12 );
    }
}

TEST( InvFluxPolicy, RoeAbsoluteEigenvaluesSonicPoints )
{
    // At vn = +/- c the acoustic branch hits zero before absolute value.
    Real a1, a2, a3;
    RoeAbsoluteEigenvalues( 1.0, 1.0, a1, a2, a3 );
    EXPECT_DOUBLE_EQ( a1, 1.0 );
    EXPECT_DOUBLE_EQ( a2, 2.0 );
    EXPECT_DOUBLE_EQ( a3, 0.0 );
    RoeAbsoluteEigenvalues( -1.0, 1.0, a1, a2, a3 );
    EXPECT_DOUBLE_EQ( a1, 1.0 );
    EXPECT_DOUBLE_EQ( a2, 0.0 );
    EXPECT_DOUBLE_EQ( a3, 2.0 );
}

// Intermediate: RoeDissipationWaveCoefficients on an |vn| grid.
// Ref: a2=3, a3=0.5, c=2 (c^2=4).
TEST( InvFluxPolicy, RoeDissipationWaveCoefficientsCurveMatchesClosedForm )
{
    const Real a2 = 3.0;
    const Real a3 = 0.5;
    const Real c = 2.0;
    const Real c2 = 4.0;
    const Real a1s[] = { 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0 };
    for ( Real a1 : a1s )
    {
        Real xi1, xi2;
        RoeDissipationWaveCoefficients( a1, a2, a3, c2, c, xi1, xi2 );
        EXPECT_NEAR( xi1, ( 2.0 * a1 - a2 - a3 ) / ( 2.0 * c2 ), 1.0e-12 );
        EXPECT_NEAR( xi2, ( a2 - a3 ) / ( 2.0 * c ), 1.0e-12 );
    }
}

TEST( InvFluxPolicy, RoeDissipationWaveCoefficientXi2IndependentOfAbsVn )
{
    Real x1a, x2a, x1b, x2b;
    RoeDissipationWaveCoefficients( 0.0, 3.0, 0.5, 4.0, 2.0, x1a, x2a );
    RoeDissipationWaveCoefficients( 4.0, 3.0, 0.5, 4.0, 2.0, x1b, x2b );
    EXPECT_DOUBLE_EQ( x2a, x2b );
    EXPECT_NE( x1a, x1b );
}

// Intermediate: RoeDissipationCouplingTerms on a Delta-H grid.
// Ref: dc=2, c^2=4, xi1=0.5, xi2=0.25.
TEST( InvFluxPolicy, RoeDissipationCouplingTermsCurveMatchesClosedForm )
{
    const Real dc = 2.0;
    const Real c2 = 4.0;
    const Real xi1 = 0.5;
    const Real xi2 = 0.25;
    const Real dhs[] = { -4.0, -1.0, 0.0, 0.5, 1.0, 2.0, 4.0 };
    for ( Real dh : dhs )
    {
        Real t1, t2;
        RoeDissipationCouplingTerms( dh, dc, c2, xi1, xi2, t1, t2 );
        EXPECT_NEAR( t1, dh * xi1 + dc * xi2, 1.0e-12 );
        EXPECT_NEAR( t2, c2 * dc * xi1 + dh * xi2, 1.0e-12 );
    }
}

TEST( InvFluxPolicy, RoeDissipationCouplingTerm2HasAcousticWeight )
{
    // term2 acoustic piece is c2*dc*xi1 independent of dh.
    Real t1a, t2a, t1b, t2b;
    RoeDissipationCouplingTerms( 0.0, 2.0, 4.0, 0.5, 0.25, t1a, t2a );
    RoeDissipationCouplingTerms( 3.0, 2.0, 4.0, 0.5, 0.25, t1b, t2b );
    EXPECT_NEAR( t2a - 0.0 * 0.25, 4.0 * 2.0 * 0.5, 1.0e-12 );
    EXPECT_NEAR( t2b - 3.0 * 0.25, 4.0 * 2.0 * 0.5, 1.0e-12 );
}

// Intermediate: RoeAverageWeights on a density-ratio grid.
TEST( InvFluxPolicy, RoeAverageWeightsCurveMatchesClosedForm )
{
    const Real rhoL = 1.0;
    const Real rhoRs[] = { 0.25, 0.5, 1.0, 2.0, 4.0, 9.0, 16.0 };
    for ( Real rhoR : rhoRs )
    {
        Real ratio = 0.0, coef = 0.0;
        RoeAverageWeights( rhoL, rhoR, ratio, coef );
        EXPECT_NEAR( ratio, std::sqrt( rhoR / rhoL ), 1.0e-12 );
        EXPECT_NEAR( coef, 1.0 / ( 1.0 + ratio ), 1.0e-12 );
    }
}

TEST( InvFluxPolicy, RoeAverageWeightsSymmetricAtUnitRatio )
{
    Real ratio = 0.0, coef = 0.0;
    RoeAverageWeights( 2.0, 2.0, ratio, coef );
    EXPECT_DOUBLE_EQ( ratio, 1.0 );
    EXPECT_DOUBLE_EQ( coef, 0.5 );
}

// Intermediate: IdealGasStaticEnthalpy on a p/rho grid (gamma=1.4).
TEST( InvFluxPolicy, IdealGasStaticEnthalpyCurveMatchesClosedForm )
{
    const Real gama = 1.4;
    const Real coef = gama / ( gama - 1.0 ); // 3.5
    const Real density = 1.0;
    const Real pressures[] = { 0.2, 0.5, 1.0, 1.5, 2.0, 4.0 };
    for ( Real p : pressures )
    {
        EXPECT_NEAR(
            IdealGasStaticEnthalpy( density, p, gama ),
            coef * ( p / density ),
            1.0e-12 );
    }
}

TEST( InvFluxPolicy, IdealGasStaticEnthalpyScalesWithPressureDensity )
{
    // Doubling p and rho keeps h constant.
    EXPECT_DOUBLE_EQ(
        IdealGasStaticEnthalpy( 1.0, 2.0, 1.4 ),
        IdealGasStaticEnthalpy( 2.0, 4.0, 1.4 ) );
}

// Intermediate: IdealGasSoundSpeed on a total-enthalpy grid (u^2=0).
TEST( InvFluxPolicy, IdealGasSoundSpeedCurveMatchesClosedForm )
{
    const Real gama = 1.4;
    const Real speedSquared = 0.0;
    const Real Hs[] = { 0.5, 1.0, 2.0, 3.5, 4.0, 5.0 };
    for ( Real H : Hs )
    {
        const Real expected = sqrt( ( gama - 1.0 ) * H );
        EXPECT_NEAR( IdealGasSoundSpeed( gama, H, speedSquared ), expected, 1.0e-12 );
    }
}

TEST( InvFluxPolicy, IdealGasSoundSpeedSubtractsKineticEnergy )
{
    // H=3.5, |u|^2=1 -> static part 3.0? Wait: H - 0.5*|u|^2 = 3.5 - 0.5 = 3.0
    // c = sqrt(0.4 * 3.0) = sqrt(1.2)
    EXPECT_NEAR(
        IdealGasSoundSpeed( 1.4, 3.5, 1.0 ),
        sqrt( 0.4 * 3.0 ),
        1.0e-12 );
}

// Intermediate: RoePressureFromTotalEnthalpy on an H grid (rho=1, u^2=0).
TEST( InvFluxPolicy, RoePressureFromTotalEnthalpyCurveMatchesClosedForm )
{
    const Real density = 1.0;
    const Real speedSquared = 0.0;
    const Real gama = 1.4;
    const Real fac = ( gama - 1.0 ) / gama;
    const Real Hs[] = { 0.5, 1.0, 2.0, 3.5, 4.0, 5.0 };
    for ( Real H : Hs )
    {
        EXPECT_NEAR(
            RoePressureFromTotalEnthalpy( density, H, speedSquared, gama ),
            density * H * fac,
            1.0e-12 );
    }
}

TEST( InvFluxPolicy, RoePressureFromTotalEnthalpyRemovesKinetic )
{
    EXPECT_NEAR(
        RoePressureFromTotalEnthalpy( 2.0, 6.5, 4.0, 1.4 ),
        2.0 * 4.5 * ( 0.4 / 1.4 ),
        1.0e-12 );
}

// Intermediate: IdealGasSpecificInternalEnergy on a kinetic grid.
// Ref: rho=1, p=1, gamma=1.4 -> e_int = 2.5.
TEST( InvFluxPolicy, IdealGasSpecificInternalEnergyCurveMatchesClosedForm )
{
    const Real density = 1.0;
    const Real pressure = 1.0;
    const Real gama = 1.4;
    const Real eInt = ( pressure / density ) / ( gama - 1.0 );
    const Real u2s[] = { 0.0, 1.0, 2.0, 4.0, 8.0, 16.0 };
    for ( Real u2 : u2s )
    {
        EXPECT_NEAR(
            IdealGasSpecificInternalEnergy( density, pressure, gama, u2 ),
            eInt + 0.5 * u2,
            1.0e-12 );
    }
}

TEST( InvFluxPolicy, IdealGasSpecificInternalEnergyScalesWithPressure )
{
    // Double p at fixed rho and u=0 doubles e_int.
    EXPECT_NEAR(
        IdealGasSpecificInternalEnergy( 1.0, 2.0, 1.4, 0.0 ),
        2.0 * IdealGasSpecificInternalEnergy( 1.0, 1.0, 1.4, 0.0 ),
        1.0e-12 );
}

// Intermediate: IdealGasSoundSpeedFromPressure on a pressure grid.
// Ref: rho=1, gamma=1.4.
TEST( InvFluxPolicy, IdealGasSoundSpeedFromPressureCurveMatchesClosedForm )
{
    const Real density = 1.0;
    const Real gama = 1.4;
    const Real pressures[] = { 0.25, 0.5, 1.0, 1.4, 2.0, 4.0 };
    for ( Real p : pressures )
    {
        EXPECT_NEAR(
            IdealGasSoundSpeedFromPressure( density, p, gama ),
            sqrt( gama * p / density ),
            1.0e-12 );
    }
}

TEST( InvFluxPolicy, IdealGasSoundSpeedFromPressureScalesWithSqrtPressure )
{
    const Real c1 = IdealGasSoundSpeedFromPressure( 1.0, 1.0, 1.4 );
    const Real c4 = IdealGasSoundSpeedFromPressure( 1.0, 4.0, 1.4 );
    EXPECT_NEAR( c4, 2.0 * c1, 1.0e-12 );
}

// Intermediate: FaceMachNumber on a vn grid (sound speed fixed).
TEST( InvFluxPolicy, FaceMachNumberCurveMatchesClosedForm )
{
    const Real c = 1.5;
    const Real vns[] = { -3.0, -1.5, 0.0, 0.75, 1.5, 3.0 };
    for ( Real vn : vns )
    {
        EXPECT_NEAR( FaceMachNumber( vn, c ), vn / c, 1.0e-12 );
    }
}

TEST( InvFluxPolicy, FaceMachNumberInverselyScalesWithSoundSpeed )
{
    EXPECT_NEAR(
        FaceMachNumber( 3.0, 1.5 ),
        2.0 * FaceMachNumber( 3.0, 3.0 ),
        1.0e-12 );
}

// Intermediate: RelativeFaceNormalVelocity on a grid-speed grid.
// Ref: n=(1,0,0), u=(2,0,0); vn = 2 - vfn.
TEST( InvFluxPolicy, RelativeFaceNormalVelocityCurveMatchesClosedForm )
{
    const Real vfns[] = { -1.0, -0.5, 0.0, 0.5, 1.0, 2.0 };
    for ( Real vfn : vfns )
    {
        EXPECT_NEAR(
            RelativeFaceNormalVelocity( 1.0, 0.0, 0.0, vfn, 2.0, 0.0, 0.0 ),
            2.0 - vfn,
            1.0e-12 );
    }
}

TEST( InvFluxPolicy, RelativeFaceNormalVelocityUsesFullNormal )
{
    // n=(0,1,0), u=(0,3,0), vfn=1 => vn=2
    EXPECT_NEAR(
        RelativeFaceNormalVelocity( 0.0, 1.0, 0.0, 1.0, 0.0, 3.0, 0.0 ),
        2.0,
        1.0e-12 );
}
