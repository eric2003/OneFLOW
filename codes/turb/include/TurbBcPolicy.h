#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm turb-bc
// @doc-stage wall-sa
// @doc-pure
// @doc-formula nu_ghost = - nu_interior  (Dirichlet 0 via linear image)
Real TurbWallSaImage( Real interiorWorkingVariable );

// @doc-algorithm turb-bc
// @doc-stage wall-omega
// @doc-pure
// @doc-formula omega_w = 60 mu_l / (rho beta1 d^2 Re)
Real TurbWallOmegaMenter(
    Real laminarViscosity, Real density, Real beta1,
    Real wallDistanceSquared, Real reynolds );

// @doc-algorithm turb-bc
// @doc-stage wall-k-image
// @doc-pure
// @doc-formula k_ghost = - k_interior
Real TurbWallKineticEnergyImage( Real interiorKineticEnergy );


// @doc-algorithm turb-bc
// @doc-stage farfield-direction
// @doc-pure
// Outward normal velocity relative to the face: use interior state when
// fluid leaves the domain (vn >= 0); otherwise use freestream turbulence.
bool TurbFarfieldUsesInterior( Real outwardNormalVelocity );

EndNameSpace
