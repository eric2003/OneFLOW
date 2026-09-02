/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\---------------------------------------------------------------------------*/

#include "HipFluxBackend.h"

#include "AccelRuntime.h"
#include "DeviceBuffer.h"

#include <hip/hip_runtime.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace ONEFLOW
{
namespace
{

void CheckHip( hipError_t status, const char * operation )
{
    if ( status == hipSuccess ) return;
    throw std::runtime_error(
        std::string( operation ) + ": " + hipGetErrorString( status ) );
}

__global__ void CalcInvFluxKernel(
    const Real * qLeft,
    const Real * qRight,
    const Real * xNormal,
    const Real * area,
    Real * flux,
    int nFaces )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nFaces ) return;

    const Real normalVelocity = xNormal[ face ];
    const Real magnitude = normalVelocity >= 0.0 ? normalVelocity : -normalVelocity;
    const Real positive = 0.5 * ( normalVelocity + magnitude );
    const Real negative = 0.5 * ( normalVelocity - magnitude );
    flux[ face ] =
        ( qLeft[ face ] * positive + qRight[ face ] * negative ) * area[ face ];
}

__global__ void AddFaceFluxKernel(
    const Real * flux,
    const int * leftCell,
    const int * rightCell,
    int nBoundaryFaces,
    int nFaces,
    int nCells,
    Real * residual )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;

    Real value = 0.0;
    for ( int face = 0; face < nBoundaryFaces; ++ face )
    {
        if ( leftCell[ face ] == cell ) value -= flux[ face ];
    }
    for ( int face = nBoundaryFaces; face < nFaces; ++ face )
    {
        if ( leftCell[ face ] == cell ) value -= flux[ face ];
        if ( rightCell[ face ] == cell ) value += flux[ face ];
    }
    residual[ cell ] += value;
}

}

void HipFluxBackend::CalcInvFlux(
    const FaceStateView & state,
    FaceFluxView & flux,
    int )
{
    if ( state.nFaces != flux.nFaces || flux.nEquations != 1 )
    {
        throw std::invalid_argument( "HIP flux view dimensions are inconsistent." );
    }

    DeviceBuffer< Real > qLeft( state.nFaces );
    DeviceBuffer< Real > qRight( state.nFaces );
    DeviceBuffer< Real > xNormal( state.nFaces );
    DeviceBuffer< Real > area( state.nFaces );
    DeviceBuffer< Real > deviceFlux( state.nFaces );

    qLeft.CopyFromHost( state.qLeft, state.nFaces );
    qRight.CopyFromHost( state.qRight, state.nFaces );
    xNormal.CopyFromHost( state.xNormal, state.nFaces );
    area.CopyFromHost( state.faceArea, state.nFaces );

    const int blockSize = 256;
    const int gridSize = ( state.nFaces + blockSize - 1 ) / blockSize;
    hipLaunchKernelGGL(
        CalcInvFluxKernel,
        dim3( gridSize ),
        dim3( blockSize ),
        0,
        0,
        qLeft.Data(),
        qRight.Data(),
        xNormal.Data(),
        area.Data(),
        deviceFlux.Data(),
        state.nFaces );
    CheckHip( hipGetLastError(), "CalcInvFluxKernel launch" );
    CheckHip( hipDeviceSynchronize(), "CalcInvFluxKernel synchronize" );
    deviceFlux.CopyToHost( flux.values, flux.nFaces );
}

void HipFluxBackend::AddFaceFlux(
    const FaceFluxView & flux,
    const FaceConnectivityView & connectivity,
    ResidualView & residual )
{
    if ( flux.nFaces != connectivity.nFaces || residual.nEquations != 1 )
    {
        throw std::invalid_argument( "HIP residual view dimensions are inconsistent." );
    }

    DeviceBuffer< Real > deviceFlux( flux.nFaces );
    DeviceBuffer< int > leftCell( connectivity.nFaces );
    DeviceBuffer< int > rightCell( connectivity.nFaces );
    DeviceBuffer< Real > deviceResidual( residual.nCells );

    deviceFlux.CopyFromHost( flux.values, flux.nFaces );
    leftCell.CopyFromHost( connectivity.leftCell, connectivity.nFaces );
    rightCell.CopyFromHost( connectivity.rightCell, connectivity.nFaces );
    deviceResidual.CopyFromHost( residual.values, residual.nCells );

    const int blockSize = 256;
    const int gridSize = ( residual.nCells + blockSize - 1 ) / blockSize;
    hipLaunchKernelGGL(
        AddFaceFluxKernel,
        dim3( gridSize ),
        dim3( blockSize ),
        0,
        0,
        deviceFlux.Data(),
        leftCell.Data(),
        rightCell.Data(),
        connectivity.nBoundaryFaces,
        connectivity.nFaces,
        residual.nCells,
        deviceResidual.Data() );
    CheckHip( hipGetLastError(), "AddFaceFluxKernel launch" );
    CheckHip( hipDeviceSynchronize(), "AddFaceFluxKernel synchronize" );
    deviceResidual.CopyToHost( residual.values, residual.nCells );
}

}
