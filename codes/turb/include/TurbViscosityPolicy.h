#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm turbulence-viscosity-cap
// @doc-stage ratio-cap
// @doc-pure
// @doc-formula mu_t = min(mu_t, r_max * mu_l)
// Shared eddy-viscosity cap used by SA, SST, and other turbulence closures.
Real TurbCappedEddyViscosity(
    Real eddyViscosity, Real laminarViscosity, Real maxViscosityRatio );

// Reject negative or NaN density values before evaluating a turbulence closure.
bool TurbDensityIsInvalid( Real density );

EndNameSpace
