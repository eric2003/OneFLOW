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

EndNameSpace
