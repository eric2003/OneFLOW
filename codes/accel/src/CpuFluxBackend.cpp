/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\---------------------------------------------------------------------------*/

#include "CpuFluxBackend.h"

#include <cmath>
#include <stdexcept>

namespace ONEFLOW
{

void CpuFluxBackend::CalcInvFlux(
    const FaceStateView & state,
    FaceFluxView & flux,
    int )
{
    if ( state.nFaces != flux.nFaces || flux.nEquations != 1 )
    {
        throw std::invalid_argument( "CPU flux view dimensions are inconsistent." );
    }

    for ( int face = 0; face < state.nFaces; ++ face )
    {
        const Real normalVelocity = state.xNormal[ face ];
        const Real positive = 0.5 * ( normalVelocity + std::abs( normalVelocity ) );
        const Real negative = 0.5 * ( normalVelocity - std::abs( normalVelocity ) );
        flux.values[ face ] =
            ( state.qLeft[ face ] * positive + state.qRight[ face ] * negative )
            * state.faceArea[ face ];
    }
}

void CpuFluxBackend::AddFaceFlux(
    const FaceFluxView & flux,
    const FaceConnectivityView & connectivity,
    ResidualView & residual )
{
    if ( flux.nFaces != connectivity.nFaces || residual.nEquations != 1 )
    {
        throw std::invalid_argument( "CPU residual view dimensions are inconsistent." );
    }

    for ( int face = 0; face < connectivity.nBoundaryFaces; ++ face )
    {
        residual.values[ connectivity.leftCell[ face ] ] -= flux.values[ face ];
    }
    for ( int face = connectivity.nBoundaryFaces; face < connectivity.nFaces; ++ face )
    {
        residual.values[ connectivity.leftCell[ face ] ] -= flux.values[ face ];
        residual.values[ connectivity.rightCell[ face ] ] += flux.values[ face ];
    }
}

}
