/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

#include "HXTypeBasic.h"

#include <cmath>

// Host/device callable marker.
//
// These functions are the one implementation that both the host chain and the
// HIP kernels compile.  Without the attribute the device compiler rejects the
// call from a __device__ body ("reference to __host__ function ... in
// __device__ function"), and duplicating the formulas into device-only copies
// is exactly the transcription drift this header exists to prevent.
//
// The guard is written on __HIPCC__ rather than __CUDACC__ so the same header
// also works in a CUDA translation unit; on a plain host compiler the macro is
// empty and the functions stay ordinary inline functions.
#if defined( __HIPCC__ ) || defined( __CUDACC__ )
#define ONEFLOW_VISCOUS_HD __host__ __device__
#else
#define ONEFLOW_VISCOUS_HD
#endif

BeginNameSpace( ONEFLOW )

// Viscous mean-flow flux arithmetic, shared by every execution space.
//
// This is the same decomposition the production host chain performs, written as
// header-inline scalar functions so the HIP kernel and the host reference
// compile the same expressions instead of two transcriptions that can drift.
// The host sources stay the reference:
//
//   codes/uns/src/UNsVisFlux.cpp:262      UNsVisFlux::PrepareFaceValue
//   codes/usolver/src/VisGrad.cpp:86      VisGradGeom::CalcGradCoef
//   codes/usolver/src/VisGrad.cpp:144     VisGrad::AverFaceValue
//   codes/usolver/src/VisGrad.cpp:152     VisGrad::CorrectFaceGrad
//   codes/usolver/src/VisGrad.cpp:160     VisGrad::CalcNormalGrad
//   codes/usolver/src/UGradPolicy.cpp:86  CorrectFaceGradient
//   codes/ns/src/NsVisFluxPolicy.cpp:6    StokesStressTensor
//   codes/ns/src/NsVisFluxPolicy.cpp:23   ProjectStressMomentumFlux
//   codes/ns/src/NsVisFluxPolicy.cpp:34   ViscousEnergyFlux
//   codes/ns/src/NsVisFluxPolicy.cpp:42   EffectiveThermalConductivity
//   codes/ns/src/NsVisFluxPolicy.cpp:51   ConductiveHeatFluxNormal
//   codes/ns/src/NsVisFluxPolicy.cpp:65   ViscousFaceFluxCoefficient
//   codes/uns/src/UNsVisFlux.cpp:439      UNsVisFlux::UpdateFaceVisFlux
//   codes/uns/src/UNsVisFlux.cpp:162      UNsVisFlux::SaveHeatFlux
//   codes/model/src/Sutherland.cpp:65     SutherlandViscosityNondim
//
// Every function takes scalars or fixed-size arrays and reads no global, so a
// host translation unit and a device translation unit can both call it.  The
// host policy headers (NsVisFluxPolicy.h, UGradPolicy.h) cannot serve that role:
// they are ordinary host declarations with no device annotation, and depending
// on the toolchain either fail to compile from a __device__ body or compile as
// separate host functions the kernel cannot reach.
//
// Sutherland uses the nondimensional sqrt( t^3 ) form, not pow( t, 1.5 ).  The
// two differ in the last bits because pow takes a different libm path, and the
// contract test compares at 1e-14 relative.
// ---------------------------------------------------------------------------

// Mirror of SutherlandViscosityNondim (codes/model/src/Sutherland.cpp:65).
ONEFLOW_VISCOUS_HD inline Real SutherlandViscosityNondimDevice(
    Real temperatureNondim, Real sutherlandTsNondim )
{
    const Real t3 = temperatureNondim * temperatureNondim * temperatureNondim;
    return sqrt( t3 ) * ( 1.0 + sutherlandTsNondim )
        / ( temperatureNondim + sutherlandTsNondim );
}

// Mirror of UnitizeDisplacement (codes/usolver/src/UGradPolicy.cpp:76).
// Returns 1 / |d|; d is replaced by the unit vector.
ONEFLOW_VISCOUS_HD inline Real ViscousUnitizeDisplacement(
    Real & dx, Real & dy, Real & dz )
{
    const Real invLength = 1.0 / sqrt( dx * dx + dy * dy + dz * dz );
    dx *= invLength;
    dy *= invLength;
    dz *= invLength;
    return invLength;
}

// Mirror of CorrectFaceGradient (codes/usolver/src/UGradPolicy.cpp:86), which
// VisGrad::CorrectFaceGrad calls through its CorrectGrad wrapper.
//
// Project the averaged gradient onto the cell-centre unit direction, replace
// that component with the finite difference across the face, and add the
// correction back along the same direction.
ONEFLOW_VISCOUS_HD inline void CorrectViscousFaceGradient(
    Real valueLeft, Real valueRight,
    Real & dqdx, Real & dqdy, Real & dqdz,
    Real unitDx, Real unitDy, Real unitDz, Real invLength )
{
    const Real projected = dqdx * unitDx + dqdy * unitDy + dqdz * unitDz;
    const Real target = ( valueRight - valueLeft ) * invLength;
    const Real df = target - projected;
    dqdx += df * unitDx;
    dqdy += df * unitDy;
    dqdz += df * unitDz;
}

// Stokes stress tensor (codes/ns/src/NsVisFluxPolicy.cpp:6).
//
//   tau_ij = mu * ( du_i/dx_j + du_j/dx_i - 2/3 delta_ij div v )
//
// The -2/3 div v term is inside the parentheses on the host, so it is inside
// here too; it is not folded into a separately computed pressure correction.
ONEFLOW_VISCOUS_HD inline void ViscousStokesStress(
    Real viscosity,
    Real dudx, Real dudy, Real dudz,
    Real dvdx, Real dvdy, Real dvdz,
    Real dwdx, Real dwdy, Real dwdz,
    Real & txx, Real & tyy, Real & tzz,
    Real & txy, Real & txz, Real & tyz )
{
    const Real divv2p3 = ( 2.0 / 3.0 ) * ( dudx + dvdy + dwdz );
    txx = viscosity * ( 2.0 * dudx - divv2p3 );
    tyy = viscosity * ( 2.0 * dvdy - divv2p3 );
    tzz = viscosity * ( 2.0 * dwdz - divv2p3 );
    txy = viscosity * ( dudy + dvdx );
    txz = viscosity * ( dudz + dwdx );
    tyz = viscosity * ( dvdz + dwdy );
}

// Project the stress tensor onto the face normal (NsVisFluxPolicy.cpp:23).
ONEFLOW_VISCOUS_HD inline void ViscousProjectStressMomentumFlux(
    Real xfn, Real yfn, Real zfn,
    Real txx, Real tyy, Real tzz,
    Real txy, Real txz, Real tyz,
    Real & fluxMx, Real & fluxMy, Real & fluxMz )
{
    fluxMx = xfn * txx + yfn * txy + zfn * txz;
    fluxMy = xfn * txy + yfn * tyy + zfn * tyz;
    fluxMz = xfn * txz + yfn * tyz + zfn * tzz;
}

// Work done by the viscous momentum flux plus the normal heat flux
// (NsVisFluxPolicy.cpp:34).
ONEFLOW_VISCOUS_HD inline Real ViscousEnergyFluxDevice(
    Real u, Real v, Real w,
    Real fluxMx, Real fluxMy, Real fluxMz,
    Real heatFluxNormal )
{
    return u * fluxMx + v * fluxMy + w * fluxMz + heatFluxNormal;
}

// k_eff = ( mu_l / Pr_l + mu_t / Pr_t ) * cp (NsVisFluxPolicy.cpp:42).
// The Prandtl arguments are the host's *inverse* numbers (nscom.oprl/oprt).
ONEFLOW_VISCOUS_HD inline Real ViscousEffectiveThermalConductivity(
    Real laminarViscosity, Real turbulentViscosity,
    Real inversePrandtlLaminar, Real inversePrandtlTurbulent,
    Real specificHeat )
{
    return ( laminarViscosity * inversePrandtlLaminar
           + turbulentViscosity * inversePrandtlTurbulent ) * specificHeat;
}

// ConductiveHeatFluxNormal (NsVisFluxPolicy.cpp:51).
ONEFLOW_VISCOUS_HD inline Real ViscousConductiveHeatFluxNormal(
    Real thermalConductivity, Real temperatureNormalGrad )
{
    return thermalConductivity * temperatureNormalGrad;
}

// ViscousFaceFluxCoefficient (NsVisFluxPolicy.cpp:65): -1/Re * A.
ONEFLOW_VISCOUS_HD inline Real ViscousFaceFluxCoefficientDevice(
    Real inverseReynolds, Real faceArea )
{
    return - inverseReynolds * faceArea;
}

// Wall heat flux as UNsVisFlux::SaveHeatFlux (codes/uns/src/UNsVisFlux.cpp:162)
// writes it: minus 1/Re times the face-normal conduction flux.
ONEFLOW_VISCOUS_HD inline Real ViscousWallHeatFlux(
    Real inverseReynolds, Real normalHeatFlux )
{
    return - inverseReynolds * normalHeatFlux;
}

EndNameSpace
