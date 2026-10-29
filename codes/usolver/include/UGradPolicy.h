#pragma once

#include "HXDefine.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm green-gauss-gradient
// @doc-stage face-average
// @doc-pure
// @doc-formula q_f = (1/2) * (q_L + q_R)
Real ArithmeticFaceAverage( Real qLeft, Real qRight );

// @doc-algorithm green-gauss-gradient
// @doc-stage inverse-distance-weights
// @doc-pure
// @doc-formula c_L = d_R / (d_L + d_R + eps), c_R = d_L / (d_L + d_R + eps)
// d_* is distance from cell center to face center.
void InverseDistanceFaceWeights(
    Real dxl, Real dyl, Real dzl,
    Real dxr, Real dyr, Real dzr,
    Real & weightLeft, Real & weightRight );

// @doc-algorithm green-gauss-gradient
// @doc-stage weighted-face-value
// @doc-pure
Real WeightedFaceAverage( Real qLeft, Real qRight, Real weightLeft, Real weightRight );

// @doc-algorithm green-gauss-gradient
// @doc-stage area-weighted-normal
// @doc-pure
// @doc-formula S = n * A
void AreaWeightedFaceNormal(
    Real xfn, Real yfn, Real zfn, Real area,
    Real & snx, Real & sny, Real & snz );

EndNameSpace
