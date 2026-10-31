#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm turb-lusgs
// @doc-stage matrix-dq
// @doc-pure
// @doc-formula df = M * dq
Real TurbLusgsMatrixTimesDq( Real matrixDiag, Real deltaQ );

// @doc-algorithm turb-lusgs
// @doc-stage solve-delta
// @doc-pure
// @doc-formula dq = rhs / lambda
Real TurbLusgsSolveByRadius( Real rightHandSide, Real spectralRadius );

// @doc-algorithm turb-lusgs
// @doc-stage residual-diff
// @doc-pure
// @doc-formula r = a - b
Real TurbLusgsDifference( Real a, Real b );

// @doc-algorithm turb-lusgs
// @doc-stage sweep-norm
// @doc-pure
// @doc-formula contrib = (dq - dq0)^2
Real TurbLusgsSweepNormContribution( Real dqNew, Real dqOld );

EndNameSpace
