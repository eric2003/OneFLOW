/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

#include "HXTypeBasic.h"

BeginNameSpace( ONEFLOW )

// Viscous flux arithmetic, shared by every execution space.
//
// These are the same formulas as codes/ns/src/NsVisFluxPolicy.cpp, written as
// header-inline so a device kernel and the host chain evaluate identical
// expressions instead of two transcriptions that can drift.  The host policy
// file stays as the reference for the general (stateful) host path; this header
// is what the device port compiles against, and a contract test compares the
// two.
//
// Nothing here allocates, reads a global, or touches a container: every input
// is a scalar.  That is what makes the same code valid inside a HIP kernel,
// a CUDA kernel, or a plain host loop.
namespace ViscousFluxArithmetic
{

// Stress tensor for a Newtonian fluid, with the -2/3 divergence correction.
// Mirrors NsVisFluxPolicy::StokesStressTensor exactly, including the sign of
// the diagonal term.
inline void StokesStress(
    Real viscosity,
    Real dudx, Real dudy, Real dudz,
    Real dvdx, Real dvdy, Real dvdz,
    Real dwdx, Real dwdy, Real dwdz,
    Real & txx, Real & tyy, Real & tzz,
    Real & txy, Real & txz, Real & tyz )
{
    const Real twoThirds = Real( 2 ) / Real( 3 );
    const Real divergence = dudx + dvdy + dwdz;
    const Real diagonal = twoThirds * divergence;
    txx = viscosity * ( Real( 2 ) * dudx - diagonal );
    tyy = viscosity * ( Real( 2 ) * dvdy - diagonal );
    tzz = viscosity * ( Real( 2 ) * dwdz - diagonal );
    txy = viscosity * ( dudy + dvdx );
    txz = viscosity * ( dudz + dwdx );
    tyz = viscosity * ( dvdz + dwdy );
}

// Project the stress tensor onto the face normal: flux_momentum = n . tau.
inline void ProjectStressMomentum(
    Real xNormal, Real yNormal, Real zNormal,
    Real txx, Real tyy, Real tzz,
    Real txy, Real txz, Real tyz,
    Real & fluxX, Real & fluxY, Real & fluxZ )
{
    fluxX = xNormal * txx + yNormal * txy + zNormal * txz;
    fluxY = xNormal * txy + yNormal * tyy + zNormal * tyz;
    fluxZ = xNormal * txz + yNormal * tyz + zNormal * tzz;
}

// Energy flux carried by the viscous momentum flux plus conduction.
inline Real ViscousEnergy(
    Real u, Real v, Real w,
    Real fluxX, Real fluxY, Real fluxZ,
    Real heatFluxNormal )
{
    return u * fluxX + v * fluxY + w * fluxZ + heatFluxNormal;
}

// Conductivity from the two viscosities and their inverse Prandtl numbers.
inline Real ThermalConductivity(
    Real laminarViscosity, Real turbulentViscosity,
    Real invPrandtlLaminar, Real invPrandtlTurbulent,
    Real specificHeat )
{
    return ( laminarViscosity * invPrandtlLaminar
           + turbulentViscosity * invPrandtlTurbulent ) * specificHeat;
}

// Normal conduction flux, q = k * dT/dn.
inline Real ConductiveHeatFlux( Real conductivity, Real temperatureNormalGrad )
{
    return conductivity * temperatureNormalGrad;
}

inline Real NormalDot( Real xNormal, Real yNormal, Real zNormal,
                       Real fx, Real fy, Real fz )
{
    return xNormal * fx + yNormal * fy + zNormal * fz;
}

// The flux is assembled in a stress-positive sign and scaled once per face.
// Mirrors NsVisFluxPolicy::ViscousFaceFluxCoefficient, which carries the minus.
inline Real FaceFluxCoefficient( Real inverseReynolds, Real faceArea )
{
    return - inverseReynolds * faceArea;
}

// Nondimensional Sutherland law, copied term for term from
// model/src/Sutherland.cpp::SutherlandViscosityNondim.
//
// The host writes sqrt(t^3) rather than pow(t, 1.5).  The two agree
// mathematically but not in the last bits, so the device path copies this exact
// form to stay inside the contract rather than being merely close.  Note there
// is no separate reference viscosity: this is the nondimensional form and the
// reference is absorbed into the Reynolds number, so inventing one would be a
// silent scale error.
inline Real SutherlandViscosityNondim(
    Real temperatureNondim, Real sutherlandTsNondim )
{
    const Real t3 = temperatureNondim * temperatureNondim * temperatureNondim;
    return sqrt( t3 )
        * ( Real( 1 ) + sutherlandTsNondim )
        / ( temperatureNondim + sutherlandTsNondim );
}

}

EndNameSpace
