#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm turb-spectrum
// @doc-stage inviscid-radius
// @doc-pure
// @doc-formula lambda_inv = (1/2) * A * |vn|
Real TurbInviscidFaceSpectralRadius( Real faceArea, Real absVnMid );

// @doc-algorithm turb-spectrum
// @doc-stage inviscid-matrix
// @doc-pure
// @doc-formula M_L += (A/2)(-vn_L - |vn|), M_R += (A/2)(vn_R - |vn|)
void TurbInviscidFaceMatrixContrib(
    Real faceArea, Real vnLeft, Real vnRight, Real absVnMid,
    Real & matrixLeftAdd, Real & matrixRightAdd );

// @doc-algorithm turb-spectrum
// @doc-stage sa-effective-nu
// @doc-pure
// @doc-formula nu_et = mu_l/rho + nu_tilde
Real SaEffectiveKinematicViscosity(
    Real laminarViscosity, Real density, Real nuTilde );

// @doc-algorithm turb-spectrum
// @doc-stage diffusive-radius
// @doc-pure
// @doc-formula lambda_diff = (1/Re) * (1/sigma) * nu_eff * A^2 / V
Real TurbDiffusiveFaceSpectralRadius(
    Real invReynolds, Real invSigma, Real nuEff,
    Real faceAreaSquared, Real invCellVolume );

EndNameSpace
