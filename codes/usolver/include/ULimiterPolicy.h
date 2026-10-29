#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// dq is the predicted face increment; minv/maxv are allowed offset bounds.
// dot is the directional factor used by the Barth face-side policy.
Real BarthFunction( Real dq, Real minv, Real maxv, Real dot );

// eps is the already-scaled regularization term used by the rational policy.
Real VencatFunction( Real dq, Real minv, Real maxv, Real eps );
Real VencatFunctionEric( Real dq, Real minv, Real maxv, Real eps );

// Build the squared global-range regularization term used by VencatFunction.
Real VencatRegularizationEpsilon( Real coefficient, Real minOffsetRange, Real maxOffsetRange );

// Direct rational forms retained as separate numerical variants.
Real Vencat( Real x, Real y, Real eps );
Real VencatEric( Real x, Real y, Real eps );

EndNameSpace
