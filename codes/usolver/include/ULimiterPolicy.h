#pragma once

#include "HXDefine.h"
#include "HXMath.h"

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


// @doc-algorithm slope-limiter-face
// @doc-stage reconstruct
// @doc-pure
// @doc-formula dq_f = g · (x_f - x_c)
// Unlimited linear reconstruction from cell gradient to face center.
Real ReconstructFaceIncrement(
    Real dqdx, Real dqdy, Real dqdz,
    Real dx, Real dy, Real dz );

// @doc-algorithm slope-limiter-face
// @doc-stage barth-direction
// @doc-pure
// @doc-formula s = orientation * (n · d) / (|d| + SMALL)
// orientation is +1 on the left side, -1 on the right side.
Real FaceSideDirectionDot(
    Real xfn, Real yfn, Real zfn,
    Real dx, Real dy, Real dz,
    Real orientationSign );

// @doc-algorithm slope-limiter-face
// @doc-stage tighten
// @doc-pure
// Keep the stricter (smaller) positive limiter bound.
inline Real TightenLimiter( Real current, Real faceValue )
{
    return MIN( current, faceValue );
}


// @doc-algorithm limited-face-reconstruction
// @doc-stage apply
// @doc-pure
// @doc-formula q_f = q0 + phi * (g · d)
Real ApplyLimitedFaceReconstruction(
    Real q0,
    Real phi,
    Real dqdx, Real dqdy, Real dqdz,
    Real dx, Real dy, Real dz );

// @doc-algorithm limited-face-reconstruction
// @doc-stage blend-phi
// @doc-pure
// Face reconstruction uses the cell limiter, or min(left,right) when both-sides mode is on.
Real BlendFaceLimiter( Real phiPrimary, Real phiNeighbor, bool useMinOfBothSides );

EndNameSpace
