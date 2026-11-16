/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

#include "AccelViews.h"

BeginNameSpace( ONEFLOW )

// The mean-flow gradient planes, addressed by the same two subscripts the
// solver's MRField uses: [velocityComponent][direction].
//
// This struct exists so the field-to-view assembly can be exercised without a
// live grid.  The assembly it feeds had a real defect that no existing test
// could catch: the SA velocity-gradient slot order was transposed
// ([dudx,dvdx,dwdx,...] instead of [dudx,dudy,dudz,...]).  Every contract test
// constructs the view by hand, so the kernels were verified while the code that
// actually fills the view in the solver was not.  The failure mode is silent:
// it compiles, the contract tests pass, and the residual is simply wrong.
struct SaMeanFlowGradientFields
{
    // gradient[ component ][ direction ].
    const Real * gradient[ 3 ][ 3 ] = {};
};

// Fill the nine velocityGradient slots in the order the kernels read them.
//
// Host order is dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
// (UTurbSrcFlux::PrepareCellValue1Equ:341-349), and
// HipTurbulenceBackend.hip:115-125 reads exactly that.  Slot 0..2 is the
// x-velocity's three derivatives; slot 3..5 the y-velocity's; slot 6..8 the
// z-velocity's.
inline void FillSaVelocityGradient(
    const SaMeanFlowGradientFields & fields, TurbSaStateView & state )
{
    for ( int component = 0; component < 3; ++ component )
    {
        for ( int direction = 0; direction < 3; ++ direction )
        {
            state.velocityGradient[ component * 3 + direction ] =
                fields.gradient[ component ][ direction ];
        }
    }
}

EndNameSpace
