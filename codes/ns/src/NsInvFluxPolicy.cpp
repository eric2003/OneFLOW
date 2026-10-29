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


Real IdealGasStaticEnthalpy( Real density, Real pressure, Real gama )
{
    return ( gama / ( gama - one ) ) * ( pressure / density );
}

Real IdealGasSpecificInternalEnergy(
    Real density, Real pressure, Real gama, Real speedSquared )
{
    return ( pressure / density ) / ( gama - one ) + half * speedSquared;
}

Real IdealGasSoundSpeedFromPressure( Real density, Real pressure, Real gama )
{
    return sqrt( ABS( gama * pressure / density ) );
}


void IdealGasPrimToConservative(
    Real density, Real u, Real v, Real w, Real pressure, Real gama,
    Real & qRho, Real & qRhoU, Real & qRhoV, Real & qRhoW, Real & qRhoE )
{
    const Real speedSquared = u * u + v * v + w * w;
    const Real em = IdealGasSpecificInternalEnergy(
        density, pressure, gama, speedSquared );
    qRho = density;
    qRhoU = density * u;
    qRhoV = density * v;
    qRhoW = density * w;
    qRhoE = density * em;
}

Real FaceMachNumber( Real vnRel, Real soundSpeed )
{
    return vnRel / soundSpeed;
}


void IdealGasConservativeToPrim(
    Real qRho, Real qRhoU, Real qRhoV, Real qRhoW, Real qRhoE, Real gama,
    Real & density, Real & u, Real & v, Real & w, Real & pressure )
{
    density = ABS( qRho );
    const Real rd = 1.0 / density;
    u = qRhoU * rd;
    v = qRhoV * rd;
    w = qRhoW * rd;
    const Real speedSquared = u * u + v * v + w * w;
    const Real reint = qRhoE - half * density * speedSquared;
    pressure = ( gama - one ) * reint;
}

void RoeDissipationWaveCoefficients(
    Real absVnFixed, Real absVnPlusCFixed, Real absVnMinusCFixed,
    Real soundSpeedSquared, Real soundSpeed,
    Real & xi1, Real & xi2 )
{
    xi1 = ( two * absVnFixed - absVnPlusCFixed - absVnMinusCFixed )
        / ( two * soundSpeedSquared );
    xi2 = ( absVnPlusCFixed - absVnMinusCFixed ) / ( two * soundSpeed );
}

Real RoeContinuityNormalJump(
    Real vnFlow,
    Real xfn, Real yfn, Real zfn,
    Real dqRho, Real dqRhoU, Real dqRhoV, Real dqRhoW )
{
    return vnFlow * dqRho - xfn * dqRhoU - yfn * dqRhoV - zfn * dqRhoW;
}


Real RoeTotalEnthalpyJump(
    Real u, Real v, Real w, Real gama,
    Real dqRho, Real dqRhoU, Real dqRhoV, Real dqRhoW, Real dqRhoE )
{
    const Real ae = gama - one;
    const Real speedSquared = u * u + v * v + w * w;
    const Real af = half * ae * speedSquared;
    Real dh = - ae * ( u * dqRhoU + v * dqRhoV + w * dqRhoW - dqRhoE );
    dh += af * dqRho;
    return dh;
}

void RoeDissipationCouplingTerms(
    Real dh, Real dc, Real soundSpeedSquared, Real xi1, Real xi2,
    Real & term1, Real & term2 )
{
    term1 = dh * xi1 + dc * xi2;
    term2 = soundSpeedSquared * dc * xi1 + dh * xi2;
}

void ApplyRoeDissipationFiveEq(
    Real meig1, Real term1, Real term2,
    Real uMid, Real vMid, Real wMid, Real hMid, Real vnFlow,
    Real xfn, Real yfn, Real zfn,
    Real dqRho, Real dqRhoU, Real dqRhoV, Real dqRhoW, Real dqRhoE,
    Real & fluxRho, Real & fluxMx, Real & fluxMy, Real & fluxMz, Real & fluxE )
{
    fluxRho -= ( meig1 * dqRho - term1 );
    fluxMx  -= ( meig1 * dqRhoU - uMid * term1 + xfn * term2 );
    fluxMy  -= ( meig1 * dqRhoV - vMid * term1 + yfn * term2 );
    fluxMz  -= ( meig1 * dqRhoW - wMid * term1 + zfn * term2 );
    fluxE   -= ( meig1 * dqRhoE - hMid * term1 + vnFlow * term2 );
}


void IdealGasBilateralSoundSpeeds(
    Real densityLeft, Real pressureLeft, Real gamaLeft,
    Real densityRight, Real pressureRight, Real gamaRight,
    Real & soundLeft, Real & soundRight, Real & soundAverage )
{
    soundLeft = IdealGasSoundSpeedFromPressure( densityLeft, pressureLeft, gamaLeft );
    soundRight = IdealGasSoundSpeedFromPressure( densityRight, pressureRight, gamaRight );
    soundAverage = half * ( soundLeft + soundRight );
}

EndNameSpace
