/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

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

// ---------------------------------------------------------------------------
// Scalar convection kernel (nEquations == 1)
// ---------------------------------------------------------------------------

__global__ void CalcInvFluxScalarKernel(
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

__global__ void AddFaceFluxScalarKernel(
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

// ---------------------------------------------------------------------------
// Euler Rusanov kernel (nEquations >= 3)
//
// Layout: equation-major: data[eq * nFaces + face]
// Input:  conserved variables Q = [rho, rho*u, rho*v, rho*w, rho*E]
// Output: numerical flux F = 0.5*(F_L+F_R) - 0.5*|lambda|*(Q_R-Q_L) * area
// ---------------------------------------------------------------------------

struct Primitive
{
    Real density;
    Real u, v, w;
    Real pressure;
    Real enthalpy;  // total enthalpy H = (rho*E + p) / rho
    Real soundSpeed;
};

__device__ void ConservedToPrimitive(
    const Real * q, int face, int nFaces, int nEq,
    Real gamma, Primitive & p )
{
    const Real gamm1 = gamma - 1.0;

    p.density = q[ 0 * nFaces + face ];
    const Real rhou = q[ 1 * nFaces + face ];
    p.u = rhou / p.density;
    Real ke = 0.5 * p.u * p.u;

    if ( nEq >= 4 )
    {
        const Real rhov = q[ 2 * nFaces + face ];
        const Real rhow = q[ 3 * nFaces + face ];
        p.v = rhov / p.density;
        p.w = rhow / p.density;
        ke = 0.5 * ( p.u * p.u + p.v * p.v + p.w * p.w );
    }
    else
    {
        p.v = 0.0;
        p.w = 0.0;
    }

    const Real rhoE = q[ ( nEq - 1 ) * nFaces + face ];
    p.pressure = gamm1 * ( rhoE - ke * p.density );
    p.enthalpy = ( rhoE + p.pressure ) / p.density;
    p.soundSpeed = sqrt( fabs( gamma * p.pressure / p.density ) );
}

__device__ void PhysicalFlux(
    const Primitive & p,
    Real nx, Real ny, Real nz, Real vfn,
    Real * f, int nEq )
{
    const Real normalVelocity = p.u * nx + p.v * ny + p.w * nz;
    const Real massFlux = p.density * normalVelocity;

    f[ 0 ] = massFlux;
    f[ 1 ] = massFlux * p.u + p.pressure * nx;
    f[ 2 ] = massFlux * p.v + p.pressure * ny;
    if ( nEq >= 4 )
    {
        f[ 3 ] = massFlux * p.w + p.pressure * nz;
        f[ 4 ] = massFlux * p.enthalpy + p.pressure * vfn;
    }
    else
    {
        f[ 2 ] = massFlux * p.enthalpy + p.pressure * vfn;
    }
}

__device__ Real MaxWaveSpeed(
    const Primitive & p,
    Real nx, Real ny, Real nz, Real vfn )
{
    const Real normalVelocity = p.u * nx + p.v * ny + p.w * nz;
    return fabs( normalVelocity - vfn ) + p.soundSpeed;
}

__global__ void CalcInvFluxEulerKernel(
    const Real * qLeft,
    const Real * qRight,
    const Real * xNormal,
    const Real * yNormal,
    const Real * zNormal,
    const Real * meshVelocityNormal,
    const Real * area,
    Real * flux,
    int nFaces,
    int nEq,
    Real gamma )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nFaces ) return;

    const Real nx = xNormal ? xNormal[ face ] : 1.0;
    const Real ny = yNormal ? yNormal[ face ] : 0.0;
    const Real nz = zNormal ? zNormal[ face ] : 0.0;
    const Real vfn = meshVelocityNormal ? meshVelocityNormal[ face ] : 0.0;
    const Real faceArea = area ? area[ face ] : 1.0;

    Primitive pL, pR;
    ConservedToPrimitive( qLeft,  face, nFaces, nEq, gamma, pL );
    ConservedToPrimitive( qRight, face, nFaces, nEq, gamma, pR );

    Real fL[ 5 ], fR[ 5 ];
    PhysicalFlux( pL, nx, ny, nz, vfn, fL, nEq );
    PhysicalFlux( pR, nx, ny, nz, vfn, fR, nEq );

    const Real lambdaL = MaxWaveSpeed( pL, nx, ny, nz, vfn );
    const Real lambdaR = MaxWaveSpeed( pR, nx, ny, nz, vfn );
    const Real lambdaMax = lambdaL > lambdaR ? lambdaL : lambdaR;

    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const int idx = eq * nFaces + face;
        flux[ idx ] = faceArea * (
            0.5 * ( fL[ eq ] + fR[ eq ] )
            - 0.5 * lambdaMax * ( qRight[ idx ] - qLeft[ idx ] ) );
    }
}

__global__ void AddFaceFluxEulerKernel(
    const Real * flux,
    const int * leftCell,
    const int * rightCell,
    int nBoundaryFaces,
    int nFaces,
    int nCells,
    int nEq,
    Real * residual )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;

    for ( int eq = 0; eq < nEq; ++ eq )
    {
        Real value = 0.0;
        for ( int face = 0; face < nBoundaryFaces; ++ face )
        {
            if ( leftCell[ face ] == cell )
                value -= flux[ eq * nFaces + face ];
        }
        for ( int face = nBoundaryFaces; face < nFaces; ++ face )
        {
            if ( leftCell[ face ] == cell )
                value -= flux[ eq * nFaces + face ];
            if ( rightCell[ face ] == cell )
                value += flux[ eq * nFaces + face ];
        }
        residual[ eq * nCells + cell ] += value;
    }
}

} // namespace

void HipFluxBackend::CalcInvFlux(
    const FaceStateView & state,
    FaceFluxView & flux,
    int )
{
    if ( state.nFaces != flux.nFaces || state.nEquations != flux.nEquations )
    {
        throw std::invalid_argument( "HIP flux view dimensions are inconsistent." );
    }
    ValidateFaceStateView( state );

    const int nFaces = state.nFaces;
    const int nEq = state.nEquations;
    const int totalValues = nFaces * nEq;

    if ( nEq == 1 )
    {
        DeviceBuffer< Real > qLeft( nFaces );
        DeviceBuffer< Real > qRight( nFaces );
        DeviceBuffer< Real > xNormal( nFaces );
        DeviceBuffer< Real > area( nFaces );
        DeviceBuffer< Real > deviceFlux( nFaces );

        qLeft.CopyFromHost( state.qLeft, nFaces );
        qRight.CopyFromHost( state.qRight, nFaces );
        if ( state.xNormal ) xNormal.CopyFromHost( state.xNormal, nFaces );
        if ( state.faceArea ) area.CopyFromHost( state.faceArea, nFaces );

        const int blockSize = 256;
        const int gridSize = ( nFaces + blockSize - 1 ) / blockSize;
        hipLaunchKernelGGL(
            CalcInvFluxScalarKernel,
            dim3( gridSize ),
            dim3( blockSize ),
            0,
            0,
            qLeft.Data(),
            qRight.Data(),
            state.xNormal ? xNormal.Data() : nullptr,
            state.faceArea ? area.Data() : nullptr,
            deviceFlux.Data(),
            nFaces );
        CheckHip( hipGetLastError(), "CalcInvFluxScalarKernel launch" );
        CheckHip( hipDeviceSynchronize(), "CalcInvFluxScalarKernel synchronize" );
        deviceFlux.CopyToHost( flux.values, nFaces );
    }
    else if ( nEq >= 3 )
    {
        DeviceBuffer< Real > qLeft( totalValues );
        DeviceBuffer< Real > qRight( totalValues );
        DeviceBuffer< Real > xNormal( nFaces );
        DeviceBuffer< Real > yNormal( nFaces );
        DeviceBuffer< Real > zNormal( nFaces );
        DeviceBuffer< Real > meshVelocityNormal( nFaces );
        DeviceBuffer< Real > area( nFaces );
        DeviceBuffer< Real > deviceFlux( totalValues );

        qLeft.CopyFromHost( state.qLeft, totalValues );
        qRight.CopyFromHost( state.qRight, totalValues );
        if ( state.xNormal ) xNormal.CopyFromHost( state.xNormal, nFaces );
        if ( state.yNormal ) yNormal.CopyFromHost( state.yNormal, nFaces );
        if ( state.zNormal ) zNormal.CopyFromHost( state.zNormal, nFaces );
        if ( state.meshVelocityNormal )
            meshVelocityNormal.CopyFromHost( state.meshVelocityNormal, nFaces );
        if ( state.faceArea ) area.CopyFromHost( state.faceArea, nFaces );

        const int blockSize = 256;
        const int gridSize = ( nFaces + blockSize - 1 ) / blockSize;
        hipLaunchKernelGGL(
            CalcInvFluxEulerKernel,
            dim3( gridSize ),
            dim3( blockSize ),
            0,
            0,
            qLeft.Data(),
            qRight.Data(),
            state.xNormal ? xNormal.Data() : nullptr,
            state.yNormal ? yNormal.Data() : nullptr,
            state.zNormal ? zNormal.Data() : nullptr,
            state.meshVelocityNormal ? meshVelocityNormal.Data() : nullptr,
            state.faceArea ? area.Data() : nullptr,
            deviceFlux.Data(),
            nFaces,
            nEq,
            state.gamma );
        CheckHip( hipGetLastError(), "CalcInvFluxEulerKernel launch" );
        CheckHip( hipDeviceSynchronize(), "CalcInvFluxEulerKernel synchronize" );
        deviceFlux.CopyToHost( flux.values, totalValues );
    }
    else
    {
        throw std::invalid_argument(
            "HIP flux backend: unsupported number of equations." );
    }
}

void HipFluxBackend::AddFaceFlux(
    const FaceFluxView & flux,
    const FaceConnectivityView & connectivity,
    ResidualView & residual )
{
    if ( flux.nFaces != connectivity.nFaces
         || flux.nEquations != residual.nEquations )
    {
        throw std::invalid_argument( "HIP residual view dimensions are inconsistent." );
    }
    ValidateFaceConnectivityView( connectivity );

    const int nFaces = flux.nFaces;
    const int nCells = residual.nCells;
    const int nEq = flux.nEquations;
    const int totalFlux = nFaces * nEq;
    const int totalResidual = nCells * nEq;

    if ( nEq == 1 )
    {
        DeviceBuffer< Real > deviceFlux( nFaces );
        DeviceBuffer< int > leftCell( nFaces );
        DeviceBuffer< int > rightCell( nFaces );
        DeviceBuffer< Real > deviceResidual( nCells );

        deviceFlux.CopyFromHost( flux.values, nFaces );
        leftCell.CopyFromHost( connectivity.leftCell, nFaces );
        rightCell.CopyFromHost( connectivity.rightCell, nFaces );
        deviceResidual.CopyFromHost( residual.values, nCells );

        const int blockSize = 256;
        const int gridSize = ( nCells + blockSize - 1 ) / blockSize;
        hipLaunchKernelGGL(
            AddFaceFluxScalarKernel,
            dim3( gridSize ),
            dim3( blockSize ),
            0,
            0,
            deviceFlux.Data(),
            leftCell.Data(),
            rightCell.Data(),
            connectivity.nBoundaryFaces,
            nFaces,
            nCells,
            deviceResidual.Data() );
        CheckHip( hipGetLastError(), "AddFaceFluxScalarKernel launch" );
        CheckHip( hipDeviceSynchronize(), "AddFaceFluxScalarKernel synchronize" );
        deviceResidual.CopyToHost( residual.values, nCells );
    }
    else
    {
        DeviceBuffer< Real > deviceFlux( totalFlux );
        DeviceBuffer< int > leftCell( nFaces );
        DeviceBuffer< int > rightCell( nFaces );
        DeviceBuffer< Real > deviceResidual( totalResidual );

        deviceFlux.CopyFromHost( flux.values, totalFlux );
        leftCell.CopyFromHost( connectivity.leftCell, nFaces );
        rightCell.CopyFromHost( connectivity.rightCell, nFaces );
        deviceResidual.CopyFromHost( residual.values, totalResidual );

        const int blockSize = 256;
        const int gridSize = ( nCells + blockSize - 1 ) / blockSize;
        hipLaunchKernelGGL(
            AddFaceFluxEulerKernel,
            dim3( gridSize ),
            dim3( blockSize ),
            0,
            0,
            deviceFlux.Data(),
            leftCell.Data(),
            rightCell.Data(),
            connectivity.nBoundaryFaces,
            nFaces,
            nCells,
            nEq,
            deviceResidual.Data() );
        CheckHip( hipGetLastError(), "AddFaceFluxEulerKernel launch" );
        CheckHip( hipDeviceSynchronize(), "AddFaceFluxEulerKernel synchronize" );
        deviceResidual.CopyToHost( residual.values, totalResidual );
    }
}

}
