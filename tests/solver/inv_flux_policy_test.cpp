#include <gtest/gtest.h>
#include <cmath>

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
