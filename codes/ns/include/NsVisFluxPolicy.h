#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm stokes-stress
// @doc-stage tensor
// @doc-pure
// @doc-formula tau_ij = mu (dui/dxj + duj/dxi) - (2/3) mu (div v) delta_ij
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
// @doc-formula F_mom = tau · n
void ProjectStressMomentumFlux(
    Real xfn, Real yfn, Real zfn,
    Real txx, Real tyy, Real tzz,
    Real txy, Real txz, Real tyz,
    Real & fluxMx, Real & fluxMy, Real & fluxMz );

// @doc-algorithm stokes-stress
// @doc-stage energy-flux
// @doc-pure
// @doc-formula F_E = v · F_mom + q_n
Real ViscousEnergyFlux(
    Real u, Real v, Real w,
    Real fluxMx, Real fluxMy, Real fluxMz,
    Real heatFluxNormal );

EndNameSpace
