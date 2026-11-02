#include "NsVisFluxPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

void StokesStressTensor(
    Real viscosity,
    Real dudx, Real dudy, Real dudz,
    Real dvdx, Real dvdy, Real dvdz,
    Real dwdx, Real dwdy, Real dwdz,
    Real & txx, Real & tyy, Real & tzz,
    Real & txy, Real & txz, Real & tyz )
{
    const Real divv2p3 = two3rd * ( dudx + dvdy + dwdz );
    txx = viscosity * ( two * dudx - divv2p3 );
    tyy = viscosity * ( two * dvdy - divv2p3 );
    tzz = viscosity * ( two * dwdz - divv2p3 );
    txy = viscosity * ( dudy + dvdx );
    txz = viscosity * ( dudz + dwdx );
    tyz = viscosity * ( dvdz + dwdy );
}

void ProjectStressMomentumFlux(
    Real xfn, Real yfn, Real zfn,
    Real txx, Real tyy, Real tzz,
    Real txy, Real txz, Real tyz,
    Real & fluxMx, Real & fluxMy, Real & fluxMz )
{
    fluxMx = xfn * txx + yfn * txy + zfn * txz;
    fluxMy = xfn * txy + yfn * tyy + zfn * tyz;
    fluxMz = xfn * txz + yfn * tyz + zfn * tzz;
}

Real ViscousEnergyFlux(
    Real u, Real v, Real w,
    Real fluxMx, Real fluxMy, Real fluxMz,
    Real heatFluxNormal )
{
    return u * fluxMx + v * fluxMy + w * fluxMz + heatFluxNormal;
}

Real EffectiveThermalConductivity(
    Real laminarViscosity, Real turbulentViscosity,
    Real invPrandtlLaminar, Real invPrandtlTurbulent,
    Real specificHeat )
{
    return ( laminarViscosity * invPrandtlLaminar
           + turbulentViscosity * invPrandtlTurbulent ) * specificHeat;
}

Real ConductiveHeatFluxNormal(
    Real thermalConductivity, Real temperatureNormalGrad )
{
    return thermalConductivity * temperatureNormalGrad;
}

Real FaceNormalDot(
    Real xfn, Real yfn, Real zfn,
    Real fx, Real fy, Real fz )
{
    return xfn * fx + yfn * fy + zfn * fz;
}


Real ViscousFaceFluxCoefficient( Real invReynolds, Real faceArea )
{
    return - invReynolds * faceArea;
}

void DepositFaceScalarToCells(
    Real faceValue, Real & cellLeft, Real & cellRight )
{
    cellLeft += faceValue;
    cellRight += faceValue;
}

EndNameSpace
