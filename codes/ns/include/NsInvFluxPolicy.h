#pragma once

#include "HXDefine.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm euler-convective-flux
Real RelativeFaceNormalVelocity(
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real u, Real v, Real w );

void EulerConvectiveFlux(
    Real rho, Real u, Real v, Real w, Real p, Real h,
    Real vnRel,
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real & fluxRho,
    Real & fluxMx,
    Real & fluxMy,
    Real & fluxMz,
    Real & fluxE );

void RoeAverageWeights(
    Real densityLeft, Real densityRight,
    Real & ratio, Real & coef );

Real RoeAverageComponent(
    Real left, Real right, Real ratio, Real coef );

void RoeAverageState(
    Real densityLeft, Real uLeft, Real vLeft, Real wLeft,
    Real pLeft, Real hLeft, Real gamaLeft,
    Real densityRight, Real uRight, Real vRight, Real wRight,
    Real pRight, Real hRight, Real gamaRight,
    Real & densityMid, Real & uMid, Real & vMid, Real & wMid,
    Real & pMid, Real & hMid, Real & gamaMid );

Real TotalEnthalpy( Real internalEnthalpy, Real speedSquared );

// @doc-algorithm roe-average
// @doc-stage implicit-face-pressure
// @doc-pure
// @doc-formula p = rho * (H - |u|^2/2) * (gamma - 1) / gamma
// Reconstruct the LUSGS face pressure from the Roe-averaged total enthalpy.
Real RoePressureFromTotalEnthalpy(
    Real density, Real totalEnthalpy, Real speedSquared, Real gama );

Real IdealGasSoundSpeed(
    Real gama, Real totalEnthalpy, Real speedSquared );

void RoeAbsoluteEigenvalues(
    Real vnRel, Real soundSpeed,
    Real & absVn, Real & absVnPlusC, Real & absVnMinusC );

void ApplyHartenEntropyFix(
    Real absVn, Real absVnPlusC, Real absVnMinusC,
    Real centropy1, Real centropy2, bool enabled,
    Real & fixedVn, Real & fixedVnPlusC, Real & fixedVnMinusC );

// AUSM splits: inline so there is exactly one definition (no LNK2005).
// @doc-algorithm ausm-split
inline Real FMSplit1( const Real & mach, const Real & ipn )
{
    return half * ( mach + ipn * ABS( mach ) );
}

inline Real FMSplit2( const Real & mach, const Real & ipn )
{
    return ipn * fourth * SQR( mach + ipn );
}

inline Real FMSplit4( const Real & mach, const Real & beta, const Real & ipn )
{
    if ( ABS( mach ) >= one )
    {
        return FMSplit1( mach, ipn );
    }
    return FMSplit2( mach, ipn ) * ( one - ipn * sixteen * beta * FMSplit2( mach, - ipn ) );
}

inline Real FPSplit5( const Real & mach, const Real & alpha, const Real & ipn )
{
    const Real macha = ABS( mach );
    if ( macha >= one )
    {
        return half * ( one + ipn * mach / ( macha + SMALL ) );
    }
    return FMSplit2( mach, ipn ) * ( ( two * ipn - mach ) - ipn * sixteen * alpha * mach * FMSplit2( mach, - ipn ) );
}


// @doc-algorithm ideal-gas-thermo
// @doc-stage static-enthalpy
// @doc-pure
// @doc-formula h = (gamma/(gamma-1)) * (p/rho)
Real IdealGasStaticEnthalpy( Real density, Real pressure, Real gama );

// @doc-algorithm ideal-gas-thermo
// @doc-stage specific-internal-energy
// @doc-pure
// @doc-formula e = p/(rho*(gamma-1)) + |v|^2/2
Real IdealGasSpecificInternalEnergy(
    Real density, Real pressure, Real gama, Real speedSquared );

// @doc-algorithm ideal-gas-thermo
// @doc-stage sound-speed-from-pressure
// @doc-pure
// @doc-formula c = sqrt(gamma * p / rho)
Real IdealGasSoundSpeedFromPressure( Real density, Real pressure, Real gama );


// @doc-algorithm ideal-gas-thermo
// @doc-stage prim-to-conservative
// @doc-pure
// Five-equation ideal-gas map (species handled by caller).
void IdealGasPrimToConservative(
    Real density, Real u, Real v, Real w, Real pressure, Real gama,
    Real & qRho, Real & qRhoU, Real & qRhoV, Real & qRhoW, Real & qRhoE );

// @doc-algorithm ideal-gas-thermo
// @doc-stage face-mach
// @doc-pure
// @doc-formula M = vn_rel / c
Real FaceMachNumber( Real vnRel, Real soundSpeed );


// @doc-algorithm ideal-gas-thermo
// @doc-stage conservative-to-prim
// @doc-pure
void IdealGasConservativeToPrim(
    Real qRho, Real qRhoU, Real qRhoV, Real qRhoW, Real qRhoE, Real gama,
    Real & density, Real & u, Real & v, Real & w, Real & pressure );

// @doc-algorithm roe-dissipation
// @doc-stage wave-coefficients
// @doc-pure
// @doc-formula xi1 = (2*a1-a2-a3)/(2*c^2), xi2 = (a2-a3)/(2*c)
void RoeDissipationWaveCoefficients(
    Real absVnFixed, Real absVnPlusCFixed, Real absVnMinusCFixed,
    Real soundSpeedSquared, Real soundSpeed,
    Real & xi1, Real & xi2 );

// @doc-algorithm roe-dissipation
// @doc-stage continuity-normal-jump
// @doc-pure
// @doc-formula dc = vn * dq_rho - n * dq_{rho v}
Real RoeContinuityNormalJump(
    Real vnFlow,
    Real xfn, Real yfn, Real zfn,
    Real dqRho, Real dqRhoU, Real dqRhoV, Real dqRhoW );


// @doc-algorithm roe-dissipation
// @doc-stage total-enthalpy-jump
// @doc-pure
Real RoeTotalEnthalpyJump(
    Real u, Real v, Real w, Real gama,
    Real dqRho, Real dqRhoU, Real dqRhoV, Real dqRhoW, Real dqRhoE );

// @doc-algorithm roe-dissipation
// @doc-stage coupling-terms
// @doc-pure
// @doc-formula term1 = dh*xi1 + dc*xi2, term2 = c2*dc*xi1 + dh*xi2
void RoeDissipationCouplingTerms(
    Real dh, Real dc, Real soundSpeedSquared, Real xi1, Real xi2,
    Real & term1, Real & term2 );

// @doc-algorithm roe-dissipation
// @doc-stage apply-five-eq
// @doc-pure
// Subtracts Roe dissipation from physical flux (in place).
void ApplyRoeDissipationFiveEq(
    Real meig1, Real term1, Real term2,
    Real uMid, Real vMid, Real wMid, Real hMid, Real vnFlow,
    Real xfn, Real yfn, Real zfn,
    Real dqRho, Real dqRhoU, Real dqRhoV, Real dqRhoW, Real dqRhoE,
    Real & fluxRho, Real & fluxMx, Real & fluxMy, Real & fluxMz, Real & fluxE );


void IdealGasBilateralSoundSpeeds(
    Real densityLeft, Real pressureLeft, Real gamaLeft,
    Real densityRight, Real pressureRight, Real gamaRight,
    Real & soundLeft, Real & soundRight, Real & soundAverage );

EndNameSpace
