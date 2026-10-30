#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm stokes-stress
// @doc-stage tensor
// @doc-pure
void StokesStressTensor(
    Real viscosity,
    Real dudx, Real dudy, Real dudz,
    Real dvdx, Real dvdy, Real dvdz,
    Real dwdx, Real dwdy, Real dwdz,
    Real & txx, Real & tyy, Real & tzz,
    Real & txy, Real & txz, Real & tyz );

// @doc-algorithm stokes-stress
// @doc-stage momentum-flux
// @doc-pure
void ProjectStressMomentumFlux(
    Real xfn, Real yfn, Real zfn,
    Real txx, Real tyy, Real tzz,
    Real txy, Real txz, Real tyz,
    Real & fluxMx, Real & fluxMy, Real & fluxMz );

// @doc-algorithm stokes-stress
// @doc-stage energy-flux
// @doc-pure
Real ViscousEnergyFlux(
    Real u, Real v, Real w,
    Real fluxMx, Real fluxMy, Real fluxMz,
    Real heatFluxNormal );

// @doc-algorithm heat-flux
// @doc-stage conductivity
// @doc-pure
// @doc-formula k/cp = mu_l/Pr_l + mu_t/Pr_t  (returns k when multiplied by cp)
Real EffectiveThermalConductivity(
    Real laminarViscosity, Real turbulentViscosity,
    Real invPrandtlLaminar, Real invPrandtlTurbulent,
    Real specificHeat );

// @doc-algorithm heat-flux
// @doc-stage normal-conductive
// @doc-pure
// @doc-formula q_n = k * dT/dn
Real ConductiveHeatFluxNormal(
    Real thermalConductivity, Real temperatureNormalGrad );

// @doc-algorithm heat-flux
// @doc-stage face-normal-dot
// @doc-pure
Real FaceNormalDot(
    Real xfn, Real yfn, Real zfn,
    Real fx, Real fy, Real fz );


// @doc-algorithm viscous-face-scale
// @doc-stage nondim-coefficient
// @doc-pure
// @doc-formula coeff = - (1/Re) * A
Real ViscousFaceFluxCoefficient( Real invReynolds, Real faceArea );

// @doc-algorithm spectral-radius
// @doc-stage bilateral-deposit
// @doc-pure
void DepositFaceScalarToCells(
    Real faceValue, Real & cellLeft, Real & cellRight );

EndNameSpace
