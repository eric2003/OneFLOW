#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm turb-inv-flux
// @doc-stage upwind-split
// @doc-pure
// @doc-formula vn+ = 1/2 (vn+|vn|), vn- = 1/2 (vn-|vn|)
void SplitRelativeNormalVelocity(
    Real vnRel, Real & vnPlus, Real & vnMinus );

// @doc-algorithm turb-inv-flux
// @doc-stage density-weight
// @doc-pure
// @doc-formula w = alpha + (1-alpha)*rho  (alpha=1 -> unity for SA)
Real TurbConvectiveDensityWeight( Real rhoCoefFlag, Real density );

// @doc-algorithm turb-inv-flux
// @doc-stage scalar-flux
// @doc-pure
// @doc-formula F = (vn+ * w_L) phi_L + (vn- * w_R) phi_R
Real UpwindScalarConvectiveFlux(
    Real vnPlus, Real weightLeft, Real phiLeft,
    Real vnMinus, Real weightRight, Real phiRight );

EndNameSpace
