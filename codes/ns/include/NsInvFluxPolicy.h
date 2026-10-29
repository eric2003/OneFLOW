#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm euler-convective-flux
// @doc-stage relative-normal-velocity
// @doc-pure
// @doc-formula vn_rel = n·v - v_fn
Real RelativeFaceNormalVelocity(
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real u, Real v, Real w );

// @doc-algorithm euler-convective-flux
// @doc-stage physical-flux
// @doc-pure
// @doc-formula F = m*[1,u,v,w,h]^T + p*[0,nx,ny,nz,v_fn]^T, m = rho * vn_rel
void EulerConvectiveFlux(
    Real rho, Real u, Real v, Real w, Real p, Real h,
    Real vnRel,
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real & fluxRho,
    Real & fluxMx,
    Real & fluxMy,
    Real & fluxMz,
    Real & fluxE );


// @doc-algorithm roe-average
// @doc-stage weights
// @doc-pure
// @doc-formula r = sqrt(rho_R/rho_L), c = 1/(1+r)
void RoeAverageWeights(
    Real densityLeft, Real densityRight,
    Real & ratio, Real & coef );

// @doc-algorithm roe-average
// @doc-stage component
// @doc-pure
// @doc-formula mid = (left + right * r) * c
Real RoeAverageComponent(
    Real left, Real right, Real ratio, Real coef );

// @doc-algorithm roe-average
// @doc-stage state
// @doc-pure
// Density uses geometric mean; other scalars use RoeAverageComponent.
void RoeAverageState(
    Real densityLeft, Real uLeft, Real vLeft, Real wLeft,
    Real pLeft, Real hLeft, Real gamaLeft,
    Real densityRight, Real uRight, Real vRight, Real wRight,
    Real pRight, Real hRight, Real gamaRight,
    Real & densityMid, Real & uMid, Real & vMid, Real & wMid,
    Real & pMid, Real & hMid, Real & gamaMid );

EndNameSpace
