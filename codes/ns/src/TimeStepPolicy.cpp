#include "TimeStepPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

Real InviscidFaceSpectralRadius(
    Real vnRelMid, Real soundSpeedMid, Real faceArea )
{
    return half * faceArea * ( ABS( vnRelMid ) + soundSpeedMid );
}

Real ViscousFaceSpectralRadiusModel1(
    Real densityMid,
    Real laminarViscosity, Real turbulentViscosity,
    Real invPrandtlLaminar, Real invPrandtlTurbulent,
    Real gama, Real reynolds, Real faceArea )
{
    const Real c1 = 4.0 / 3.0 * ( laminarViscosity + turbulentViscosity );
    const Real c2 = gama * ( laminarViscosity * invPrandtlLaminar
                           + turbulentViscosity * invPrandtlTurbulent );
    const Real c3 = two * MAX( c1, c2 ) / ( reynolds * densityMid );
    return SQR( faceArea ) * c3;
}

Real ViscousFaceSpectralRadiusModel2(
    Real densityMid, Real totalViscosity,
    Real cellDistance, Real reynolds, Real faceArea )
{
    const Real c1 = two * totalViscosity
        / ( densityMid * cellDistance * reynolds + SMALL );
    return half * c1 * faceArea;
}

Real ViscousFaceSpectralRadiusModel3(
    Real densityMid, Real reynolds, Real faceArea )
{
    return SQR( faceArea ) / ( reynolds * densityMid + SMALL );
}

Real CellTimeStepFromSpectrum(
    Real cfl, Real cellVolume, Real spectralRadius )
{
    return cfl * cellVolume / spectralRadius;
}

Real CombineInvAndVisTimeStep(
    Real invTimeStep, Real visTimeStep )
{
    return invTimeStep * visTimeStep / ( invTimeStep + visTimeStep );
}

EndNameSpace
