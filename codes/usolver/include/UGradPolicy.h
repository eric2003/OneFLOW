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


// @doc-algorithm green-gauss-gradient
// @doc-stage scale-by-volume
// @doc-pure
// @doc-formula grad = grad / V
void ScaleGradientByVolume(
    Real volume,
    Real & dqdx, Real & dqdy, Real & dqdz );

// @doc-algorithm green-gauss-gradient
// @doc-stage boundary-ghost-copy
// @doc-pure
// Boundary ghost cells inherit the owner (left) cell gradient.
void CopyGradient(
    Real srcDqdx, Real srcDqdy, Real srcDqdz,
    Real & dstDqdx, Real & dstDqdy, Real & dstDqdz );

// @doc-algorithm green-gauss-gradient
// @doc-stage face-flux-accumulate
// @doc-pure
// @doc-formula grad += sign * S * q_f
void AccumulateFaceGradient(
    Real snx, Real sny, Real snz, Real faceValue, Real sign,
    Real & dqdx, Real & dqdy, Real & dqdz );


// @doc-algorithm vis-face-geometry
// @doc-stage skew-angle
// @doc-pure
// @doc-formula angle_deg = asin(clamp(proj / |d|, -1, 1)) * 180 / pi
Real FaceSkewAngleDegrees(
    Real dx, Real dy, Real dz, Real projectedLength );

// @doc-algorithm vis-face-geometry
// @doc-stage unitize
// @doc-pure
// @doc-formula u = d / |d|, invLength = 1 / |d|
void UnitizeDisplacement(
    Real & dx, Real & dy, Real & dz, Real & invLength );


// @doc-algorithm vis-face-geometry
// @doc-stage correct-grad
// @doc-pure
// @doc-formula g <- g + ((fr-fl)/|LR| - g·u) * u
void CorrectFaceGradient(
    Real qLeft, Real qRight,
    Real & dqdx, Real & dqdy, Real & dqdz,
    Real unitDx, Real unitDy, Real unitDz, Real invLength );

EndNameSpace
