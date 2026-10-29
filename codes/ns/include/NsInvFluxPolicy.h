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

EndNameSpace
