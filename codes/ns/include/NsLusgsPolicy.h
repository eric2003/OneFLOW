#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm lusgs-sweep
// @doc-stage split-eigenvalues
// @doc-pure
// @doc-formula lmd_i = 1/2 (lambda_i + sign * |lambda|_face,i)
void LusgsSplitFaceEigenvalues(
    Real vnRel, Real soundSpeed, Real signOfMatrix,
    Real faceAbs1, Real faceAbs2, Real faceAbs3,
    Real & lmd1, Real & lmd2, Real & lmd3 );

// @doc-algorithm lusgs-sweep
// @doc-stage linearization-weights
// @doc-pure
// @doc-formula x1 = (2*lmd1-lmd2-lmd3)/(2*c^2), x2 = (lmd2-lmd3)/(2*c)
void LusgsLinearizationWeights(
    Real lmd1, Real lmd2, Real lmd3,
    Real soundSpeedSquared, Real soundSpeed,
    Real & x1, Real & x2 );

// @doc-algorithm lusgs-sweep
// @doc-stage five-eq-increment
// @doc-pure
// Approximate inviscid flux Jacobian times dq, scaled by face area.
void LusgsFluxIncrementFiveEq(
    Real lmd1, Real term1, Real term2,
    Real u, Real v, Real w, Real totalEnthalpy, Real vnFluid,
    Real xfn, Real yfn, Real zfn, Real faceArea,
    Real dqRho, Real dqRhoU, Real dqRhoV, Real dqRhoW, Real dqRhoE,
    Real & dfRho, Real & dfMx, Real & dfMy, Real & dfMz, Real & dfE );

EndNameSpace
