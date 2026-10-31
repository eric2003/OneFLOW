#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm turb-vis-flux
// @doc-stage face-average
// @doc-pure
Real ArithmeticMean( Real a, Real b );

// @doc-algorithm turb-vis-flux
// @doc-stage sa-com-vis
// @doc-pure
// @doc-formula mu_eff = (nu_l + nu_t) / sigma
Real SaFaceEffectiveDiffusivity(
    Real laminarKinematicViscosity,
    Real workingVariableAverage,
    Real invSigma );

// @doc-algorithm turb-vis-flux
// @doc-stage sst-com-vis
// @doc-pure
// @doc-formula Gamma = mu_l + sigma * mu_t
Real SstFaceEffectiveDiffusivity(
    Real laminarViscosity, Real turbulentViscosity, Real sigma );

// @doc-algorithm turb-vis-flux
// @doc-stage scalar-flux
// @doc-pure
// @doc-formula F = - (1/Re) * Gamma * (dq/dn) * A
Real TurbScalarDiffusiveFlux(
    Real invReynolds, Real effectiveDiffusivity,
    Real normalGradient, Real faceArea );

EndNameSpace
