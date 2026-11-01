#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// Shared eddy-viscosity cap used by SA, SST, and other turbulence closures.
Real TurbCappedEddyViscosity(
    Real eddyViscosity, Real laminarViscosity, Real maxViscosityRatio );

EndNameSpace
