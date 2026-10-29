#include "NsInvFluxPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

Real RelativeFaceNormalVelocity(
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real u, Real v, Real w )
{
    return xfn * u + yfn * v + zfn * w - vfn;
}

void EulerConvectiveFlux(
    Real rho, Real u, Real v, Real w, Real p, Real h,
    Real vnRel,
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real & fluxRho,
    Real & fluxMx,
    Real & fluxMy,
    Real & fluxMz,
    Real & fluxE )
{
    const Real massFlux = rho * vnRel;
    fluxRho = massFlux;
    fluxMx  = massFlux * u + xfn * p;
    fluxMy  = massFlux * v + yfn * p;
    fluxMz  = massFlux * w + zfn * p;
    fluxE   = massFlux * h + vfn * p;
}

void RoeAverageWeights(
    Real densityLeft, Real densityRight,
    Real & ratio, Real & coef )
{
    ratio = sqrt( densityRight / densityLeft );
    coef = 1.0 / ( 1.0 + ratio );
}

Real RoeAverageComponent(
    Real left, Real right, Real ratio, Real coef )
{
    return ( left + right * ratio ) * coef;
}

void RoeAverageState(
    Real densityLeft, Real uLeft, Real vLeft, Real wLeft,
    Real pLeft, Real hLeft, Real gamaLeft,
    Real densityRight, Real uRight, Real vRight, Real wRight,
    Real pRight, Real hRight, Real gamaRight,
    Real & densityMid, Real & uMid, Real & vMid, Real & wMid,
    Real & pMid, Real & hMid, Real & gamaMid )
{
    Real ratio = 0.0;
    Real coef = 0.0;
    RoeAverageWeights( densityLeft, densityRight, ratio, coef );
    densityMid = sqrt( densityLeft * densityRight );
    uMid = RoeAverageComponent( uLeft, uRight, ratio, coef );
    vMid = RoeAverageComponent( vLeft, vRight, ratio, coef );
    wMid = RoeAverageComponent( wLeft, wRight, ratio, coef );
    pMid = RoeAverageComponent( pLeft, pRight, ratio, coef );
    hMid = RoeAverageComponent( hLeft, hRight, ratio, coef );
    gamaMid = RoeAverageComponent( gamaLeft, gamaRight, ratio, coef );
}

Real TotalEnthalpy( Real internalEnthalpy, Real speedSquared )
{
    return internalEnthalpy + half * speedSquared;
}

Real IdealGasSoundSpeed(
    Real gama, Real totalEnthalpy, Real speedSquared )
{
    const Real c2 = ( gama - one ) * ( totalEnthalpy - half * speedSquared );
    return sqrt( ABS( c2 ) );
}

void RoeAbsoluteEigenvalues(
    Real vnRel, Real soundSpeed,
    Real & absVn, Real & absVnPlusC, Real & absVnMinusC )
{
    absVn = ABS( vnRel );
    absVnPlusC = ABS( vnRel + soundSpeed );
    absVnMinusC = ABS( vnRel - soundSpeed );
}

void ApplyHartenEntropyFix(
    Real absVn, Real absVnPlusC, Real absVnMinusC,
    Real centropy1, Real centropy2, bool enabled,
    Real & fixedVn, Real & fixedVnPlusC, Real & fixedVnMinusC )
{
    fixedVn = absVn;
    fixedVnPlusC = absVnPlusC;
    fixedVnMinusC = absVnMinusC;
    if ( ! enabled )
    {
        return;
    }
    const Real maxEigenvalue = MAX( absVnPlusC, absVnMinusC );
    const Real m1 = maxEigenvalue * centropy1;
    const Real m2 = maxEigenvalue * centropy2;
    fixedVn = MAX( m1, absVn );
    fixedVnPlusC = MAX( m2, absVnPlusC );
    fixedVnMinusC = MAX( m2, absVnMinusC );
}

EndNameSpace
