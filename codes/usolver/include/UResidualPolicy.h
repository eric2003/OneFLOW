#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm residual-assembly
// @doc-stage deposit-face-flux
// @doc-pure
// @doc-formula R_L -= F; if interior: R_R += F
// Face flux is oriented from left cell to right cell. Boundary faces have no
// independent right residual contribution in this assembly.
inline void DepositOrientedFaceFlux(
    Real flux,
    bool isBoundaryFace,
    Real & leftResidual,
    Real & rightResidual )
{
    leftResidual -= flux;
    if ( ! isBoundaryFace )
    {
        rightResidual += flux;
    }
}

EndNameSpace
