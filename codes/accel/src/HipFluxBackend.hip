/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "HipFluxBackend.h"

#include "AccelRuntime.h"
#include "EulerCpuAdapter.h"
#include "DeviceBuffer.h"
#include "StageProfiler.h"

#include <hip/hip_runtime.h>

#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>

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

std::unique_ptr< HipFluxBackend > sharedHipFluxBackend;

class HipNs3DBackendState final : public Ns3DBackendState
{
public:
    explicit HipNs3DBackendState( int deviceIdValue )
        : deviceId_( deviceIdValue )
    {
    }

    AccelBackendKind Backend() const noexcept override
    {
        return AccelBackendKind::HIP;
    }

    int DeviceId() const noexcept override
    {
        return deviceId_;
    }

    HipGradientStorage & Gradient() noexcept
    {
        return gradient_;
    }

private:
    int deviceId_ = -1;
    HipGradientStorage gradient_;
};


// Build the cell-to-face CSR used by the cell-based production kernels.  The
// face arrays stay the source of truth; this is only a deterministic
// enumeration of each owned cell's incident faces (boundary faces appear on
// the left cell only; internal faces appear on both sides).
void EnsureCellFaceConnectivity(
    HipGradientStorage & storage,
    const int * leftCell,
    const int * rightCell,
    int nCells,
    int nFaces,
    int nBoundaryFaces )
{
    if ( storage.cellFaceValid
         && storage.cachedCellFaceLeft == leftCell
         && storage.cachedCellFaceRight == rightCell
         && storage.cachedCellFaceCells == nCells
         && storage.cachedCellFaceFaces == nFaces
         && storage.cachedCellFaceBoundaryFaces == nBoundaryFaces )
    {
        return;
    }
    if ( nCells <= 0 || nFaces <= 0 || leftCell == nullptr
         || rightCell == nullptr )
    {
        throw std::invalid_argument(
            "cell-face connectivity requires a non-empty face partition" );
    }
    std::vector< int > offsets( static_cast< std::size_t >( nCells ) + 1, 0 );
    for ( int face = 0; face < nFaces; ++ face )
    {
        const int left = leftCell[ face ];
        if ( left >= 0 && left < nCells ) ++ offsets[ left + 1 ];
        const int right = rightCell[ face ];
        if ( right >= 0 && right < nCells ) ++ offsets[ right + 1 ];
    }
    for ( int cell = 0; cell < nCells; ++ cell )
    {
        offsets[ cell + 1 ] += offsets[ cell ];
    }
    std::vector< int > ids( static_cast< std::size_t >( offsets[ nCells ] ) );
    std::vector< int > cursor( offsets.begin(), offsets.end() - 1 );
    for ( int face = 0; face < nFaces; ++ face )
    {
        const int left = leftCell[ face ];
        if ( left >= 0 && left < nCells ) ids[ cursor[ left ] ++ ] = face;
        const int right = rightCell[ face ];
        if ( right >= 0 && right < nCells ) ids[ cursor[ right ] ++ ] = face;
    }
    storage.cellFaceOffsets.CopyFromHost(
        offsets.data(), static_cast< int >( offsets.size() ) );
    storage.cellFaceIds.CopyFromHost(
        ids.data(), static_cast< int >( ids.size() ) );
    storage.cachedCellFaceLeft = leftCell;
    storage.cachedCellFaceRight = rightCell;
    storage.cachedCellFaceCells = nCells;
    storage.cachedCellFaceFaces = nFaces;
    storage.cachedCellFaceBoundaryFaces = nBoundaryFaces;
    storage.cellFaceValid = true;
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
    const unsigned char * boundaryMask,
    int nBoundaryFaces,
    int nFaces,
    int nCells,
    Real * residual )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;

    Real value = 0.0;
    for ( int face = 0; face < nFaces; ++ face )
    {
        const bool boundary = boundaryMask != nullptr
            ? boundaryMask[ face ] != 0
            : face < nBoundaryFaces;
        if ( leftCell[ face ] == cell ) value -= flux[ face ];
        if ( ! boundary && rightCell[ face ] == cell )
        {
            value += flux[ face ];
        }
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
    Real gamma,
    int scheme )
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

    Real fL[ 5 ] = {};
    Real fR[ 5 ] = {};
    Real lambdaMax = 0.0;
    if ( scheme == 1 )
    {
        const Real ratio = sqrt( pR.density / pL.density );
        const Real coef = 1.0 / ( 1.0 + ratio );
        const Real densityMean = sqrt( pL.density * pR.density );
        const Real uMean = ( pL.u + pR.u * ratio ) * coef;
        const Real vMean = ( pL.v + pR.v * ratio ) * coef;
        const Real wMean = ( pL.w + pR.w * ratio ) * coef;
        const Real pressureMean =
            ( pL.pressure + pR.pressure * ratio ) * coef;
        const Real normalVelocityMean =
            nx * uMean + ny * vMean + nz * wMean - vfn;
        const Real soundSpeedMean = sqrt(
            fabs( gamma * pressureMean / densityMean ) );
        lambdaMax = fabs( normalVelocityMean ) + soundSpeedMean;

        const Real normalVelocityLeft =
            nx * pL.u + ny * pL.v + nz * pL.w - vfn;
        const Real normalVelocityRight =
            nx * pR.u + ny * pR.v + nz * pR.w - vfn;
        const Real massFluxLeft = pL.density * normalVelocityLeft;
        const Real massFluxRight = pR.density * normalVelocityRight;
        fL[ 0 ] = massFluxLeft;
        fL[ 1 ] = massFluxLeft * pL.u + nx * pL.pressure;
        fR[ 0 ] = massFluxRight;
        fR[ 1 ] = massFluxRight * pR.u + nx * pR.pressure;
        if ( nEq >= 4 )
        {
            fL[ 2 ] = massFluxLeft * pL.v + ny * pL.pressure;
            fL[ 3 ] = massFluxLeft * pL.w + nz * pL.pressure;
            fL[ 4 ] = massFluxLeft * pL.enthalpy + vfn * pL.pressure;
            fR[ 2 ] = massFluxRight * pR.v + ny * pR.pressure;
            fR[ 3 ] = massFluxRight * pR.w + nz * pR.pressure;
            fR[ 4 ] = massFluxRight * pR.enthalpy + vfn * pR.pressure;
        }
        else
        {
            fL[ 2 ] = massFluxLeft * pL.enthalpy + vfn * pL.pressure;
            fR[ 2 ] = massFluxRight * pR.enthalpy + vfn * pR.pressure;
        }
    }
    else
    {
        PhysicalFlux( pL, nx, ny, nz, vfn, fL, nEq );
        PhysicalFlux( pR, nx, ny, nz, vfn, fR, nEq );
        const Real lambdaL = MaxWaveSpeed( pL, nx, ny, nz, vfn );
        const Real lambdaR = MaxWaveSpeed( pR, nx, ny, nz, vfn );
        lambdaMax = lambdaL > lambdaR ? lambdaL : lambdaR;
    }

    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const int idx = eq * nFaces + face;
        flux[ idx ] = faceArea * (
            0.5 * ( fL[ eq ] + fR[ eq ] )
            - 0.5 * lambdaMax * ( qRight[ idx ] - qLeft[ idx ] ) );
    }
}

__device__ void PrimitiveToConservedFace(
    const Real * primitive, int face, int nFaces, int nEq,
    Real gamma, Real * q )
{
    const Real density = primitive[ 0 * nFaces + face ];
    const Real u = primitive[ 1 * nFaces + face ];
    const Real v = nEq >= 4 ? primitive[ 2 * nFaces + face ] : 0.0;
    const Real w = nEq >= 4 ? primitive[ 3 * nFaces + face ] : 0.0;
    const Real pressure = primitive[ ( nEq == 3 ? 2 : 4 ) * nFaces + face ];
    const Real kineticEnergy = 0.5 * density *
        ( u * u + v * v + w * w );

    q[ 0 ] = density;
    q[ 1 ] = density * u;
    if ( nEq == 3 )
    {
        q[ 2 ] = pressure / ( gamma - 1.0 ) + kineticEnergy;
    }
    else
    {
        q[ 2 ] = density * v;
        q[ 3 ] = density * w;
        q[ 4 ] = pressure / ( gamma - 1.0 ) + kineticEnergy;
    }
}

__device__ void ConservedValuesToPrimitive(
    const Real * q, int nEq, Real gamma, Primitive & p )
{
    const Real gamm1 = gamma - 1.0;
    p.density = q[ 0 ];
    p.u = q[ 1 ] / p.density;
    p.v = 0.0;
    p.w = 0.0;
    Real kinetic = 0.5 * p.u * p.u;
    if ( nEq >= 4 )
    {
        p.v = q[ 2 ] / p.density;
        p.w = q[ 3 ] / p.density;
        kinetic = 0.5 * ( p.u * p.u + p.v * p.v + p.w * p.w );
    }
    const Real rhoE = q[ nEq - 1 ];
    p.pressure = gamm1 * ( rhoE - kinetic * p.density );
    p.enthalpy = ( rhoE + p.pressure ) / p.density;
    p.soundSpeed = sqrt( fabs( gamma * p.pressure / p.density ) );
}

__global__ void CalcAndScatterPrimitiveEulerKernel(
    const Real * primitiveLeft,
    const Real * primitiveRight,
    const Real * xNormal,
    const Real * yNormal,
    const Real * zNormal,
    const Real * meshVelocityNormal,
    const Real * area,
    const int * leftCell,
    const int * rightCell,
    const unsigned char * boundaryMask,
    int nBoundaryFaces,
    int nFaces,
    int nCells,
    int nEq,
    Real gamma,
    int scheme,
    Real * flux,
    Real * residual )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nFaces ) return;

    const Real nx = xNormal ? xNormal[ face ] : 1.0;
    const Real ny = yNormal ? yNormal[ face ] : 0.0;
    const Real nz = zNormal ? zNormal[ face ] : 0.0;
    const Real vfn = meshVelocityNormal ? meshVelocityNormal[ face ] : 0.0;
    const Real faceArea = area ? area[ face ] : 1.0;

    Real qL[ 5 ] = {};
    Real qR[ 5 ] = {};
    PrimitiveToConservedFace(
        primitiveLeft, face, nFaces, nEq, gamma, qL );
    PrimitiveToConservedFace(
        primitiveRight, face, nFaces, nEq, gamma, qR );

    Primitive pL, pR;
    ConservedValuesToPrimitive( qL, nEq, gamma, pL );
    ConservedValuesToPrimitive( qR, nEq, gamma, pR );

    Real fL[ 5 ] = {};
    Real fR[ 5 ] = {};
    Real lambdaMax = 0.0;
    if ( scheme == 1 )
    {
        const Real ratio = sqrt( pR.density / pL.density );
        const Real coef = 1.0 / ( 1.0 + ratio );
        const Real densityMean = sqrt( pL.density * pR.density );
        const Real uMean = ( pL.u + pR.u * ratio ) * coef;
        const Real vMean = ( pL.v + pR.v * ratio ) * coef;
        const Real wMean = ( pL.w + pR.w * ratio ) * coef;
        const Real pressureMean =
            ( pL.pressure + pR.pressure * ratio ) * coef;
        const Real normalVelocityMean =
            nx * uMean + ny * vMean + nz * wMean - vfn;
        const Real soundSpeedMean = sqrt(
            fabs( gamma * pressureMean / densityMean ) );
        lambdaMax = fabs( normalVelocityMean ) + soundSpeedMean;

        const Real normalVelocityLeft =
            nx * pL.u + ny * pL.v + nz * pL.w - vfn;
        const Real normalVelocityRight =
            nx * pR.u + ny * pR.v + nz * pR.w - vfn;
        const Real massFluxLeft = pL.density * normalVelocityLeft;
        const Real massFluxRight = pR.density * normalVelocityRight;
        fL[ 0 ] = massFluxLeft;
        fL[ 1 ] = massFluxLeft * pL.u + nx * pL.pressure;
        fR[ 0 ] = massFluxRight;
        fR[ 1 ] = massFluxRight * pR.u + nx * pR.pressure;
        if ( nEq >= 4 )
        {
            fL[ 2 ] = massFluxLeft * pL.v + ny * pL.pressure;
            fL[ 3 ] = massFluxLeft * pL.w + nz * pL.pressure;
            fL[ 4 ] = massFluxLeft * pL.enthalpy + vfn * pL.pressure;
            fR[ 2 ] = massFluxRight * pR.v + ny * pR.pressure;
            fR[ 3 ] = massFluxRight * pR.w + nz * pR.pressure;
            fR[ 4 ] = massFluxRight * pR.enthalpy + vfn * pR.pressure;
        }
        else
        {
            fL[ 2 ] = massFluxLeft * pL.enthalpy + vfn * pL.pressure;
            fR[ 2 ] = massFluxRight * pR.enthalpy + vfn * pR.pressure;
        }
    }
    else
    {
        PhysicalFlux( pL, nx, ny, nz, vfn, fL, nEq );
        PhysicalFlux( pR, nx, ny, nz, vfn, fR, nEq );
        const Real lambdaL = MaxWaveSpeed( pL, nx, ny, nz, vfn );
        const Real lambdaR = MaxWaveSpeed( pR, nx, ny, nz, vfn );
        lambdaMax = lambdaL > lambdaR ? lambdaL : lambdaR;
    }

    // Write the face flux only.  The cell-based gather kernel consumes it;
    // the fused face-parallel atomic scatter was the main 3-D bottleneck.
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        flux[ eq * nFaces + face ] = faceArea * (
            0.5 * ( fL[ eq ] + fR[ eq ] )
            - 0.5 * lambdaMax * ( qR[ eq ] - qL[ eq ] ) );
    }
}

__global__ void GradientZeroKernel( Real * gradient, int nValues )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index < nValues ) gradient[ index ] = 0.0;
}

__global__ void GradientAccumulateKernel(
    const Real * q, const Real * xFace, const Real * yFace,
    const Real * zFace, const Real * xNormal,
    const Real * yNormal, const Real * zNormal, const Real * faceArea,
    const Real * xCell, const Real * yCell,
    const Real * zCell, const int * leftCell, const int * rightCell,
    int nFaces, int nBoundaryFaces, int nCells, int nEq,
    Real * gradient )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nFaces ) return;
    const int left = leftCell[ face ];
    const int right = rightCell[ face ];
    const Real dxl = xFace[ face ] - xCell[ left ];
    const Real dyl = yFace[ face ] - yCell[ left ];
    const Real dzl = zFace[ face ] - zCell[ left ];
    const Real dxr = xFace[ face ] - xCell[ right ];
    const Real dyr = yFace[ face ] - yCell[ right ];
    const Real dzr = zFace[ face ] - zCell[ right ];
    const Real delt1 = sqrt( dxl * dxl + dyl * dyl + dzl * dzl );
    const Real delt2 = sqrt( dxr * dxr + dyr * dyr + dzr * dzr );
    const Real delta = 1.0 / ( delt1 + delt2 + 1.0e-40 );
    const Real cl = delt2 * delta;
    const Real cr = delt1 * delta;
    const Real ax = xNormal[ face ] * faceArea[ face ];
    const Real ay = yNormal[ face ] * faceArea[ face ];
    const Real az = zNormal[ face ] * faceArea[ face ];
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const Real value = cl * q[ eq * ( nCells + nBoundaryFaces ) + left ]
            + cr * q[ eq * ( nCells + nBoundaryFaces ) + right ];
        const int base = eq * 3 * ( nCells + nBoundaryFaces );
        atomicAdd( &gradient[ base + left ], ax * value );
        atomicAdd( &gradient[ base + ( nCells + nBoundaryFaces ) + left ], ay * value );
        atomicAdd( &gradient[ base + 2 * ( nCells + nBoundaryFaces ) + left ], az * value );
        if ( face >= nBoundaryFaces )
        {
            atomicAdd( &gradient[ base + right ], -ax * value );
            atomicAdd( &gradient[ base + ( nCells + nBoundaryFaces ) + right ], -ay * value );
            atomicAdd( &gradient[ base + 2 * ( nCells + nBoundaryFaces ) + right ], -az * value );
        }
    }
}

__global__ void GradientNormalizeKernel(
    Real * gradient, const Real * volume, int nCells,
    int nBoundaryFaces, int nEq )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;
    const int totalCells = nCells + nBoundaryFaces;
    const Real inverseVolume = 1.0 / volume[ cell ];
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const int base = eq * 3 * totalCells;
        gradient[ base + cell ] *= inverseVolume;
        gradient[ base + totalCells + cell ] *= inverseVolume;
        gradient[ base + 2 * totalCells + cell ] *= inverseVolume;
    }
}

__global__ void GradientCopyBoundaryKernel(
    Real * gradient, const int * leftCell, const int * rightCell,
    int nCells, int nBoundaryFaces, int nEq )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nBoundaryFaces ) return;
    const int totalCells = nCells + nBoundaryFaces;
    const int left = leftCell[ face ];
    const int ghost = rightCell[ face ];
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const int base = eq * 3 * totalCells;
        gradient[ base + ghost ] = gradient[ base + left ];
        gradient[ base + totalCells + ghost ] =
            gradient[ base + totalCells + left ];
        gradient[ base + 2 * totalCells + ghost ] =
            gradient[ base + 2 * totalCells + left ];
    }
}

__device__ void AtomicMinPositiveReal( Real * address, Real value )
{
    auto * bits = reinterpret_cast< unsigned long long * >( address );
    const unsigned long long valueBits =
        static_cast< unsigned long long >( __double_as_longlong( value ) );
    unsigned long long previous = *bits;
    while ( valueBits < previous )
    {
        const unsigned long long assumed = previous;
        previous = atomicCAS( bits, assumed, valueBits );
        if ( previous == assumed ) break;
    }
}

__global__ void GradientCellKernel(
    const Real * q, Real * gradient,
    const Real * xFace, const Real * yFace, const Real * zFace,
    const Real * xNormal, const Real * yNormal, const Real * zNormal,
    const Real * faceArea, const Real * xCell, const Real * yCell,
    const Real * zCell, const Real * volume,
    const int * leftCell, const int * rightCell,
    const int * cellFaceOffsets, const int * cellFaceIds,
    int nCells, int nBoundaryFaces, int nFaces, int nEq )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;
    const int totalCells = nCells + nBoundaryFaces;
    Real accumulated[ 5 * 3 ] = {};
    for ( int slot = cellFaceOffsets[ cell ];
          slot < cellFaceOffsets[ cell + 1 ]; ++ slot )
    {
        const int face = cellFaceIds[ slot ];
        const int left = leftCell[ face ];
        const int right = rightCell[ face ];
        const Real dxl = xFace[ face ] - xCell[ left ];
        const Real dyl = yFace[ face ] - yCell[ left ];
        const Real dzl = zFace[ face ] - zCell[ left ];
        const Real dxr = xFace[ face ] - xCell[ right ];
        const Real dyr = yFace[ face ] - yCell[ right ];
        const Real dzr = zFace[ face ] - zCell[ right ];
        const Real delt1 = sqrt( dxl * dxl + dyl * dyl + dzl * dzl );
        const Real delt2 = sqrt( dxr * dxr + dyr * dyr + dzr * dzr );
        const Real delta = 1.0 / ( delt1 + delt2 + 1.0e-40 );
        const Real cl = delt2 * delta;
        const Real cr = delt1 * delta;
        const Real ax = xNormal[ face ] * faceArea[ face ];
        const Real ay = yNormal[ face ] * faceArea[ face ];
        const Real az = zNormal[ face ] * faceArea[ face ];
        const Real sign = left == cell ? 1.0 : -1.0;
        for ( int eq = 0; eq < nEq; ++ eq )
        {
            const Real value = cl * q[ eq * totalCells + left ]
                + cr * q[ eq * totalCells + right ];
            accumulated[ eq * 3 ] += sign * ax * value;
            accumulated[ eq * 3 + 1 ] += sign * ay * value;
            accumulated[ eq * 3 + 2 ] += sign * az * value;
        }
    }
    const Real inverseVolume = 1.0 / volume[ cell ];
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const int base = eq * 3 * totalCells;
        gradient[ base + cell ] = accumulated[ eq * 3 ] * inverseVolume;
        gradient[ base + totalCells + cell ] =
            accumulated[ eq * 3 + 1 ] * inverseVolume;
        gradient[ base + 2 * totalCells + cell ] =
            accumulated[ eq * 3 + 2 ] * inverseVolume;
    }
}

__global__ void PackGradientFacesKernel(
    const Real * gradient, const int * cellIds, int nFaces, int nEq,
    int totalCells, Real * packed )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total = nFaces * 3 * nEq;
    if ( index >= total ) return;
    const int face = index % nFaces;
    const int plane = index / nFaces;
    const int component = plane / nEq;
    const int equation = plane - component * nEq;
    packed[ index ] = gradient[
        equation * 3 * totalCells + component * totalCells
        + cellIds[ face ] ];
}

__global__ void UnpackGradientFacesKernel(
    Real * gradient, const int * ghostCellIds, int nFaces, int nEq,
    int totalCells, const Real * packed )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total = nFaces * 3 * nEq;
    if ( index >= total ) return;
    const int face = index % nFaces;
    const int plane = index / nFaces;
    const int component = plane / nEq;
    const int equation = plane - component * nEq;
    gradient[
        equation * 3 * totalCells + component * totalCells
        + ghostCellIds[ face ] ] = packed[ index ];
}

// Gather the ghost slice of every equation into one contiguous buffer so the
// host can fetch it with a single blocking copy instead of one per equation.
__global__ void PackGhostCellsKernel(
    const Real * cellState, int nCells, int totalCells, int nEquations,
    int nGhostCells, Real * packed )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total = nEquations * nGhostCells;
    if ( index >= total ) return;
    const int equation = index / nGhostCells;
    const int ghost = index - equation * nGhostCells;
    packed[ index ] = cellState[ equation * totalCells + nCells + ghost ];
}

// Gather an arbitrary set of owned cells into an [equation][cell] buffer so a
// host reader needs one blocking D2H instead of nEquations*cellCount copies.
__global__ void PackStateCellsKernel(
    const Real * cellState, const int * cellIds, int cellCount,
    int nEquations, int totalCells, Real * packed )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total = nEquations * cellCount;
    if ( index >= total ) return;
    const int equation = index / cellCount;
    const int slot = index - equation * cellCount;
    packed[ index ] =
        cellState[ equation * totalCells + cellIds[ slot ] ];
}

__global__ void InviscidSpectrumZeroKernel(
    Real * spectrum, int nCells )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell < nCells ) spectrum[ cell ] = 0.0;
}

__global__ void InviscidSpectrumAccumulateKernel(
    const Real * q,
    const Real * xNormal,
    const Real * yNormal,
    const Real * zNormal,
    const Real * meshVelocityNormal,
    const Real * faceArea,
    const int * leftCell,
    const int * rightCell,
    int nCells,
    int nGhostCells,
    int nFaces,
    Real gamma,
    Real * spectrum,
    int * status )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nFaces ) return;

    const int totalCells = nCells + nGhostCells;
    const int left = leftCell[ face ];
    const int right = rightCell[ face ];
    if ( left < 0 || left >= nCells
         || right < 0 || right >= totalCells )
    {
        atomicExch( status, 1 );
        return;
    }

    const Real rhoLeft = q[ 0 * totalCells + left ];
    const Real uLeft = q[ 1 * totalCells + left ];
    const Real vLeft = q[ 2 * totalCells + left ];
    const Real wLeft = q[ 3 * totalCells + left ];
    const Real pressureLeft = q[ 4 * totalCells + left ];
    const Real rhoRight = q[ 0 * totalCells + right ];
    const Real uRight = q[ 1 * totalCells + right ];
    const Real vRight = q[ 2 * totalCells + right ];
    const Real wRight = q[ 3 * totalCells + right ];
    const Real pressureRight = q[ 4 * totalCells + right ];
    if ( rhoLeft <= 0.0 || rhoRight <= 0.0
         || pressureLeft <= 0.0 || pressureRight <= 0.0 )
    {
        atomicExch( status, 1 );
        return;
    }

    const Real nx = xNormal[ face ];
    const Real ny = yNormal[ face ];
    const Real nz = zNormal[ face ];
    const Real vfn = meshVelocityNormal[ face ];
    const Real vnLeft =
        nx * uLeft + ny * vLeft + nz * wLeft - vfn;
    const Real vnRight =
        nx * uRight + ny * vRight + nz * wRight - vfn;
    const Real densityMean = 0.5 * ( rhoLeft + rhoRight );
    const Real pressureMean = 0.5 * ( pressureLeft + pressureRight );
    const Real soundSpeed =
        sqrt( gamma * pressureMean / densityMean );
    const Real normalVelocity = 0.5 * ( vnLeft + vnRight );
    const Real faceSpectrum = 0.5 * faceArea[ face ]
        * ( fabs( normalVelocity ) + soundSpeed );
    if ( ! isfinite( faceSpectrum ) || faceSpectrum < 0.0 )
    {
        atomicExch( status, 1 );
        return;
    }

    atomicAdd( &spectrum[ left ], faceSpectrum );
    if ( right < nCells )
        atomicAdd( &spectrum[ right ], faceSpectrum );
}

__global__ void LocalInviscidTimeStepKernel(
    const Real * spectrum,
    const Real * volume,
    Real cfl,
    int nCells,
    Real * timeStep,
    Real * minimum,
    int * status )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;

    const Real value = cfl * volume[ cell ] / spectrum[ cell ];
    if ( ! isfinite( value ) || value <= 0.0 )
    {
        atomicExch( status, 1 );
        return;
    }
    timeStep[ cell ] = value;
    AtomicMinPositiveReal( minimum, value );
}

__global__ void CapLocalTimeStepKernel(
    Real * timeStep,
    const Real * minimum,
    Real maxTimeRatio,
    int nCells )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;
    const Real maximum = maxTimeRatio * minimum[ 0 ];
    if ( timeStep[ cell ] > maximum ) timeStep[ cell ] = maximum;
}

__device__ Real BoundaryPower( Real value, Real exponent, Real gamma )
{
    if ( gamma > 1.399999999999 && gamma < 1.400000000001
         && exponent > 2.499999999999 && exponent < 2.500000000001 )
    {
        return value * value * sqrt( value );
    }
    return pow( value, exponent );
}

__global__ void UpdateBoundaryPrimitiveKernel(
    Real * q,
    Real * bcQ,
    const int * leftCell,
    const int * rightCell,
    const int * boundaryType,
    const Real * xNormal,
    const Real * yNormal,
    const Real * zNormal,
    const Real * meshVelocityNormal,
    const Real * inflow,
    int faceOuterNormal,
    int nCells,
    int nBoundaryFaces,
    int nEq,
    Real gamma )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nBoundaryFaces ) return;

    const int left = leftCell[ face ];
    const int ghost = rightCell[ face ];
    const int totalCells = nCells + nBoundaryFaces;
    const int type = boundaryType[ face ];

    // Interfaces, periodic faces, and overset faces are owned by an external
    // exchange path. Never overwrite their ghost values locally.
    if ( type < 0 || type == 99 || type >= 1000 ) return;

    Real rho = q[ 0 * totalCells + left ];
    Real u = q[ 1 * totalCells + left ];
    Real v = q[ 2 * totalCells + left ];
    Real w = q[ 3 * totalCells + left ];
    Real pressure = q[ 4 * totalCells + left ];
    Real ghostValues[ 5 ] = { rho, u, v, w, pressure };

    if ( type == 5 ) // INFLOW
    {
        for ( int eq = 0; eq < nEq; ++ eq )
            ghostValues[ eq ] = inflow[ eq ];
    }
    else if ( type == 1 || type == 6 || type == 7
              || type == 71 || type == 72 || type == 73 ) // extrapolation/outflow/pole
    {
        // The CPU oracle uses the owner-cell primitive state for these models.
    }
    else if ( type == 2 || type == 3 ) // inviscid wall/symmetry
    {
        const Real nx = xNormal[ face ] * faceOuterNormal;
        const Real ny = yNormal[ face ] * faceOuterNormal;
        const Real nz = zNormal[ face ] * faceOuterNormal;
        const Real vn = nx * u + ny * v + nz * w
            - meshVelocityNormal[ face ];
        ghostValues[ 1 ] = u - 2.0 * nx * vn;
        ghostValues[ 2 ] = v - 2.0 * ny * vn;
        ghostValues[ 3 ] = w - 2.0 * nz * vn;
    }
    else if ( type == 4 ) // FARFIELD
    {
        const Real rref = inflow[ 0 ];
        const Real uref = inflow[ 1 ];
        const Real vref = inflow[ 2 ];
        const Real wref = inflow[ 3 ];
        const Real pref = inflow[ 4 ];
        const Real nx = xNormal[ face ] * faceOuterNormal;
        const Real ny = yNormal[ face ] * faceOuterNormal;
        const Real nz = zNormal[ face ] * faceOuterNormal;
        const Real vnin = nx * u + ny * v + nz * w
            - meshVelocityNormal[ face ];
        const Real vnref = nx * uref + ny * vref + nz * wref
            - meshVelocityNormal[ face ];
        const Real cref = sqrt( fabs( gamma * pref / rref ) );
        const Real cin = sqrt( fabs( gamma * pressure / rho ) );
        const Real gamm1 = gamma - 1.0;
        const Real velocity = sqrt( u * u + v * v + w * w );

        if ( velocity > cin )
        {
            if ( vnin < 0.0 )
            {
                for ( int eq = 0; eq < nEq; ++ eq )
                    ghostValues[ eq ] = inflow[ eq ];
            }
        }
        else
        {
            const Real riemp = vnin + 2.0 * cin / gamm1;
            const Real riemm = vnref - 2.0 * cref / gamm1;
            const Real vnb = 0.5 * ( riemp + riemm );
            const Real cb = 0.25 * ( riemp - riemm ) * gamm1;
            Real entropy;
            Real vtx;
            Real vty;
            Real vtz;
            if ( vnb >= 0.0 )
            {
                entropy = pressure / BoundaryPower( rho, gamma, gamma );
                vtx = u - nx * vnin;
                vty = v - ny * vnin;
                vtz = w - nz * vnin;
            }
            else
            {
                entropy = pref / BoundaryPower( rref, gamma, gamma );
                vtx = uref - nx * vnref;
                vty = vref - ny * vnref;
                vtz = wref - nz * vnref;
            }
            const Real rhoBoundary = BoundaryPower(
                cb * cb / ( entropy * gamma ), 1.0 / gamm1, gamma );
            const Real ub = vtx + nx * vnb;
            const Real vb = vty + ny * vnb;
            const Real wb = vtz + nz * vnb;
            ghostValues[ 0 ] = rhoBoundary;
            ghostValues[ 1 ] = ub;
            ghostValues[ 2 ] = vb;
            ghostValues[ 3 ] = wb;
            ghostValues[ 4 ] = cb * cb * rhoBoundary / gamma;
        }
    }
    else
    {
        return;
    }

    for ( int eq = 0; eq < nEq; ++ eq )
        q[ eq * totalCells + ghost ] = ghostValues[ eq ];

    if ( type == 2 && bcQ != nullptr )
    {
        // NsBcSolver::SymmetryBc writes the face state as the average of the
        // owner and reflected ghost values; the limiter uses this for both
        // reconstructed sides on solid surfaces.
        for ( int eq = 0; eq < nEq; ++ eq )
        {
            bcQ[ eq * nBoundaryFaces + face ] =
                0.5 * ( q[ eq * totalCells + left ] + ghostValues[ eq ] );
        }
    }
}


__global__ void ReconstructFaceEulerKernel(
    const Real * q,
    const Real * gradient,
    const Real * xFace,
    const Real * yFace,
    const Real * zFace,
    const Real * xCell,
    const Real * yCell,
    const Real * zCell,
    const int * leftCell,
    const int * rightCell,
    const unsigned char * boundaryMask,
    const ReconstructionBoundaryOperation * boundaryOperation,
    const Real * bcQ,
    int bcQStride,
    int nCells,
    int nBoundaryFaces,
    int nFaces,
    int nEq,
    int physicality,
    int densityComponent,
    int pressureComponent,
    Real * qf1,
    Real * qf2 )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nFaces ) return;

    const int totalCells = nCells + nBoundaryFaces;
    const int left = leftCell[ face ];
    const int right = rightCell[ face ];
    Real leftCandidate[ 5 ] = {};
    Real rightCandidate[ 5 ] = {};
    const Real dxl = xFace[ face ] - xCell[ left ];
    const Real dyl = yFace[ face ] - yCell[ left ];
    const Real dzl = zFace[ face ] - zCell[ left ];
    const Real dxr = xFace[ face ] - xCell[ right ];
    const Real dyr = yFace[ face ] - yCell[ right ];
    const Real dzr = zFace[ face ] - zCell[ right ];

    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const int qBase = eq * totalCells;
        const int gradientBase = eq * 3 * totalCells;
        const Real qLeft = q[ qBase + left ];
        const Real qRight = q[ qBase + right ];
        leftCandidate[ eq ] = qLeft
            + gradient[ gradientBase + left ] * dxl
            + gradient[ gradientBase + totalCells + left ] * dyl
            + gradient[ gradientBase + 2 * totalCells + left ] * dzl;
        rightCandidate[ eq ] = qRight
            + gradient[ gradientBase + right ] * dxr
            + gradient[ gradientBase + totalCells + right ] * dyr
            + gradient[ gradientBase + 2 * totalCells + right ] * dzr;
    }

    bool leftPhysical = true;
    bool rightPhysical = true;
    if ( physicality != 0 )
    {
        leftPhysical = leftCandidate[ densityComponent ] > 0.0
            && leftCandidate[ pressureComponent ] > 0.0;
        rightPhysical = rightCandidate[ densityComponent ] > 0.0
            && rightCandidate[ pressureComponent ] > 0.0;
    }
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const int qBase = eq * totalCells;
        qf1[ eq * nFaces + face ] =
            leftPhysical ? leftCandidate[ eq ] : q[ qBase + left ];
        qf2[ eq * nFaces + face ] =
            rightPhysical ? rightCandidate[ eq ] : q[ qBase + right ];
    }

    const bool boundary = boundaryMask != nullptr
        ? boundaryMask[ face ] != 0
        : face < nBoundaryFaces;
    if ( ! boundary ) return;

    ReconstructionBoundaryOperation operation =
        ReconstructionBoundaryOperation::Average;
    if ( boundaryOperation != nullptr ) operation = boundaryOperation[ face ];
    if ( operation == ReconstructionBoundaryOperation::Preserve ) return;
    if ( operation == ReconstructionBoundaryOperation::SolidOverride )
    {
        for ( int eq = 0; eq < nEq; ++ eq )
        {
            const Real value = bcQ[ eq * bcQStride + face ];
            qf1[ eq * nFaces + face ] = value;
            qf2[ eq * nFaces + face ] = value;
        }
        return;
    }
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const int qBase = eq * totalCells;
        const Real value = 0.5 * ( q[ qBase + left ] + q[ qBase + right ] );
        qf1[ eq * nFaces + face ] = value;
        qf2[ eq * nFaces + face ] = value;
    }
}

__global__ void AddFaceFluxEulerKernel(
    const Real * flux,
    const int * leftCell,
    const int * rightCell,
    const unsigned char * boundaryMask,
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
        for ( int face = 0; face < nFaces; ++ face )
        {
            const bool boundary = boundaryMask != nullptr
                ? boundaryMask[ face ] != 0
                : face < nBoundaryFaces;
            if ( leftCell[ face ] == cell )
                value -= flux[ eq * nFaces + face ];
            if ( ! boundary && rightCell[ face ] == cell )
                value += flux[ eq * nFaces + face ];
        }
        residual[ eq * nCells + cell ] += value;
    }
}

// Scatter each face contribution to its owner cells.  The previous kernel
// assigned one thread per cell and scanned every face, which made the
// residual update O(nCells*nFaces) for an unstructured 3D grid.
__global__ void ScatterFaceFluxEulerKernel(
    const Real * flux,
    const int * leftCell,
    const int * rightCell,
    const unsigned char * boundaryMask,
    int nBoundaryFaces,
    int nFaces,
    int nCells,
    int nEq,
    Real * residual )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nFaces ) return;

    const int left = leftCell[ face ];
    const int right = rightCell[ face ];
    const bool boundary = boundaryMask != nullptr
        ? boundaryMask[ face ] != 0
        : face < nBoundaryFaces;
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const Real faceFlux = flux[ eq * nFaces + face ];
        atomicAdd( &residual[ eq * nCells + left ], -faceFlux );
        if ( ! boundary && right >= 0 && right < nCells )
            atomicAdd( &residual[ eq * nCells + right ], faceFlux );
    }
}

__global__ void GatherFaceFluxEulerKernel(
    const Real * flux,
    const int * leftCell,
    const int * rightCell,
    const unsigned char * boundaryMask,
    const int * cellFaceOffsets,
    const int * cellFaceIds,
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
        for ( int slot = cellFaceOffsets[ cell ];
              slot < cellFaceOffsets[ cell + 1 ]; ++ slot )
        {
            const int face = cellFaceIds[ slot ];
            const Real faceFlux = flux[ eq * nFaces + face ];
            if ( leftCell[ face ] == cell )
            {
                value -= faceFlux;
            }
            else if ( rightCell[ face ] == cell )
            {
                const bool boundary = boundaryMask != nullptr
                    ? boundaryMask[ face ] != 0
                    : face < nBoundaryFaces;
                if ( ! boundary ) value += faceFlux;
            }
        }
        residual[ eq * nCells + cell ] += value;
    }
}

__global__ void ScaleResidualKernel(
    Real * residual,
    const Real * timeStep,
    const Real * cellVolume,
    int nCells,
    int nEquations,
    Real rkCoefficient )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total = nCells * nEquations;
    if ( index >= total ) return;

    const int cell = index % nCells;
    residual[ index ] *=
        ( timeStep[ cell ] / cellVolume[ cell ] ) * rkCoefficient;
}

__global__ void UpdatePrimitiveStateKernel(
    Real * primitive,
    const Real * residual,
    int nCells,
    int totalCells,
    Real gamma,
    int * status )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;

    const Real density = primitive[ 0 * totalCells + cell ];
    const Real u = primitive[ 1 * totalCells + cell ];
    const Real v = primitive[ 2 * totalCells + cell ];
    const Real w = primitive[ 3 * totalCells + cell ];
    const Real pressure = primitive[ 4 * totalCells + cell ];
    const Real kinetic =
        0.5 * density * ( u * u + v * v + w * w );

    const Real updatedDensity =
        density + residual[ 0 * nCells + cell ];
    const Real updatedMomentumX =
        density * u + residual[ 1 * nCells + cell ];
    const Real updatedMomentumY =
        density * v + residual[ 2 * nCells + cell ];
    const Real updatedMomentumZ =
        density * w + residual[ 3 * nCells + cell ];
    const Real updatedEnergy =
        pressure / ( gamma - 1.0 ) + kinetic
        + residual[ 4 * nCells + cell ];

    const Real inverseDensity = 1.0 / updatedDensity;
    const Real updatedU = updatedMomentumX * inverseDensity;
    const Real updatedV = updatedMomentumY * inverseDensity;
    const Real updatedW = updatedMomentumZ * inverseDensity;
    const Real updatedPressure = ( gamma - 1.0 )
        * ( updatedEnergy
            - 0.5 * updatedDensity
                * ( updatedU * updatedU
                    + updatedV * updatedV
                    + updatedW * updatedW ) );

    if ( ! isfinite( updatedDensity ) || ! isfinite( updatedU )
         || ! isfinite( updatedV ) || ! isfinite( updatedW )
         || ! isfinite( updatedPressure )
         || updatedDensity <= 0.0 || updatedPressure <= 0.0 )
    {
        atomicExch( status, 1 );
        return;
    }

    primitive[ 0 * totalCells + cell ] = updatedDensity;
    primitive[ 1 * totalCells + cell ] = updatedU;
    primitive[ 2 * totalCells + cell ] = updatedV;
    primitive[ 3 * totalCells + cell ] = updatedW;
    primitive[ 4 * totalCells + cell ] = updatedPressure;
}

__global__ void PackBoundaryKernel(
    const Real * state,
    const int * cellIds,
    int cellCount,
    int interfaceCount,
    int nEquations,
    Real * packed )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int valueCount = interfaceCount * nEquations;
    if ( index >= valueCount ) return;
    const int equation = index / interfaceCount;
    const int interfaceId = index - equation * interfaceCount;
    packed[ index ] =
        state[ equation * cellCount + cellIds[ interfaceId ] ];
}

__global__ void UnpackBoundaryKernel(
    const Real * packed,
    const int * cellIds,
    int cellCount,
    int interfaceCount,
    int nEquations,
    Real * state )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int valueCount = interfaceCount * nEquations;
    if ( index >= valueCount ) return;
    const int equation = index / interfaceCount;
    const int interfaceId = index - equation * interfaceCount;
    state[ equation * cellCount + cellIds[ interfaceId ] ] =
        packed[ index ];
}

} // namespace

std::size_t HipBoundaryEndpointKeyHash::operator()(
    const HipBoundaryEndpointKey & key ) const noexcept
{
    std::size_t result = std::hash< int >{}( key.localZone );
    result ^= std::hash< int >{}( key.neighborZone )
        + static_cast< std::size_t >( 0x9e3779b9 )
        + ( result << 6 ) + ( result >> 2 );
    result ^= std::hash< int >{}( key.ghostLevel )
        + static_cast< std::size_t >( 0x9e3779b9 )
        + ( result << 6 ) + ( result >> 2 );
    return result;
}

std::unique_ptr< Ns3DBackendState > CreateHipNs3DBackendState(
    int deviceId )
{
    if ( deviceId < 0 )
    {
        throw std::invalid_argument(
            "cannot create HIP Ns3D state with an invalid device" );
    }
    return std::make_unique< HipNs3DBackendState >( deviceId );
}

void HipFluxBackend::BindState(
    const void * cacheKey,
    Ns3DDeviceState & state )
{
    if ( cacheKey == nullptr )
    {
        throw std::invalid_argument(
            "cannot bind HIP Ns3D state without a cache key" );
    }
    if ( state.key.backend != AccelBackendKind::HIP
         || state.key.deviceId
             != AccelRuntime::Instance().Context().selectedDevice )
    {
        throw std::invalid_argument(
            "HIP Ns3D state key does not match the active device" );
    }
    auto * hipState =
        dynamic_cast< HipNs3DBackendState * >( state.BackendState() );
    if ( hipState == nullptr
         || hipState->DeviceId() != state.key.deviceId )
    {
        throw std::invalid_argument(
            "HIP Ns3D state has no matching backend ownership" );
    }
    const auto found = boundStates_.find( cacheKey );
    if ( found != boundStates_.end() && found->second != &state )
    {
        throw std::logic_error(
            "HIP grid cache key is already bound to another Ns3D state" );
    }
    boundStates_[ cacheKey ] = &state;
}

void HipFluxBackend::UnbindState( const void * cacheKey ) noexcept
{
    if ( cacheKey != nullptr ) boundStates_.erase( cacheKey );
}

namespace
{

// Pinned host staging for the production gradient mirror.  hipMemcpy to or
// from pageable memory stages through a driver buffer and measured at roughly
// 200 MB/s on this platform, which made the per-stage gradient mirror the
// single largest remaining loop cost.  The buffer is allocated once per thread
// and grown on demand; it is intentionally not freed during static teardown,
// when the runtime may already be gone.
class PinnedStagingBuffer
{
public:
    Real * Get( std::size_t count )
    {
        if ( count > capacity_ )
        {
            if ( data_ != nullptr )
            {
                CheckHip( hipHostFree( data_ ), "pinned staging free" );
            }
            void * allocated = nullptr;
            CheckHip(
                hipHostMalloc(
                    & allocated, count * sizeof( Real ),
                    hipHostMallocDefault ),
                "pinned staging allocate" );
            data_ = static_cast< Real * >( allocated );
            capacity_ = count;
        }
        return data_;
    }

private:
    Real * data_ = nullptr;
    std::size_t capacity_ = 0;
};

PinnedStagingBuffer & LocalPinnedStaging()
{
    static thread_local PinnedStagingBuffer buffer;
    return buffer;
}

bool ResidentContractCheckRequested()
{
    const char * value =
        std::getenv( "ONEFLOW_UNS_HIP_RESIDENT_CONTRACT" );
    return value != nullptr && value[ 0 ] == '1';
}

int ResidentContractRank()
{
    const char * names[] = {
        "OMPI_COMM_WORLD_RANK", "MV2_COMM_WORLD_RANK", "SLURM_PROCID" };
    for ( const char * name : names )
    {
        const char * value = std::getenv( name );
        if ( value != nullptr && value[ 0 ] != '\0' )
            return std::atoi( value );
    }
    return 0;
}

// Consume the resident-state generation contract before a device kernel reads
// the resident primitive state.  Owned-cell generation changes force a full
// upload; ghost-only host mutations refresh just the ghost region.
ResidentStateSync RefreshResidentPrimitiveState(
    HipGradientStorage & storage,
    Ns3DDeviceState * owner,
    const void * ownerToken,
    std::uint64_t topologyGeneration,
    std::uint64_t fieldGeneration,
    int nCells,
    int nGhostCells,
    int nEquations,
    const Real * const * hostComponents )
{
    const int totalCells = nCells + nGhostCells;
    storage.residentContract.Bind(
        ownerToken, topologyGeneration, nCells, nGhostCells, nEquations );
    bool hostPointersChanged = false;
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        if ( hostComponents[ equation ] == nullptr )
        {
            throw std::invalid_argument(
                "resident primitive state host component is null" );
        }
        if ( storage.cachedStateQ[ equation ] != hostComponents[ equation ] )
            hostPointersChanged = true;
    }
    // Without an owning state there is no observable field generation, so the
    // original "refresh from the host on every call" fallback is preserved.
    ResidentStateSync sync = owner == nullptr
        ? ResidentStateSync::FullField
        : storage.residentContract.RequiredSync( fieldGeneration );
    if ( hostPointersChanged ) sync = ResidentStateSync::FullField;
    if ( sync == ResidentStateSync::FullField
         && owner != nullptr
         && storage.residentContract.Resident()
         && ! storage.residentContract.HostMirrorCurrent() )
    {
        // The device advanced the field and nobody pulled it back.  Uploading
        // the host mirror now would overwrite the newest device state with an
        // old field, so fail loudly instead: the host writer must refresh the
        // mirror first (ResidentStateContract::NoteHostMirrorRefresh).
        throw std::logic_error(
            "host primitive mirror is stale; refresh it from the device "
            "before uploading" );
    }
    if ( sync == ResidentStateSync::FullField )
    {
        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            AccelRuntime::Instance().Backend().Copy(
                storage.cellState.Data()
                    + static_cast< std::size_t >( equation ) * totalCells,
                hostComponents[ equation ],
                static_cast< std::size_t >( totalCells ) * sizeof( Real ),
                AccelCopyKind::HostToDevice );
        }
        storage.residentContract.NoteFullFieldSync( fieldGeneration );
    }
    else if ( sync == ResidentStateSync::GhostCells )
    {
        if ( nGhostCells > 0 )
        {
            for ( int equation = 0; equation < nEquations; ++ equation )
            {
                AccelRuntime::Instance().Backend().Copy(
                    storage.cellState.Data()
                        + static_cast< std::size_t >( equation ) * totalCells
                        + nCells,
                    hostComponents[ equation ] + nCells,
                    static_cast< std::size_t >( nGhostCells ) * sizeof( Real ),
                    AccelCopyKind::HostToDevice );
            }
        }
        storage.residentContract.NoteGhostSync();
    }
    for ( int equation = 0; equation < nEquations; ++ equation )
        storage.cachedStateQ[ equation ] = hostComponents[ equation ];
    storage.cachedEquations = nEquations;
    return sync;
}

void VerifyResidentPrimitiveState(
    HipGradientStorage & storage,
    int nCells,
    int nGhostCells,
    int nEquations,
    const Real * const * hostComponents,
    Real tolerance,
    Real * maxAbsoluteDifference,
    long * mismatchedValues,
    long * ownedMismatches = nullptr,
    long * ghostMismatches = nullptr,
    int * firstEquation = nullptr,
    int * firstCell = nullptr,
    Real * firstDeviceValue = nullptr,
    Real * firstHostValue = nullptr )
{
    const int totalCells = nCells + nGhostCells;
    Real maximum = 0.0;
    long mismatches = 0;
    long owned = 0;
    long ghost = 0;
    int firstEq = -1;
    int firstCellIndex = -1;
    Real firstDevice = 0.0;
    Real firstHost = 0.0;
    std::vector< Real > deviceValues(
        static_cast< std::size_t >( totalCells ) );
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        AccelRuntime::Instance().Backend().Copy(
            deviceValues.data(),
            storage.cellState.Data()
                + static_cast< std::size_t >( equation ) * totalCells,
            static_cast< std::size_t >( totalCells ) * sizeof( Real ),
            AccelCopyKind::DeviceToHost );
        const Real * host = hostComponents[ equation ];
        for ( int cell = 0; cell < totalCells; ++ cell )
        {
            const Real difference =
                std::fabs( deviceValues[ cell ] - host[ cell ] );
            if ( difference > maximum ) maximum = difference;
            if ( difference > tolerance )
            {
                ++mismatches;
                if ( cell < nCells ) ++owned; else ++ghost;
                if ( firstEq < 0 )
                {
                    firstEq = equation;
                    firstCellIndex = cell;
                    firstDevice = deviceValues[ cell ];
                    firstHost = host[ cell ];
                }
            }
        }
    }
    if ( maxAbsoluteDifference != nullptr )
        *maxAbsoluteDifference = maximum;
    if ( mismatchedValues != nullptr ) *mismatchedValues = mismatches;
    if ( ownedMismatches != nullptr ) *ownedMismatches = owned;
    if ( ghostMismatches != nullptr ) *ghostMismatches = ghost;
    if ( firstEquation != nullptr ) *firstEquation = firstEq;
    if ( firstCell != nullptr ) *firstCell = firstCellIndex;
    if ( firstDeviceValue != nullptr ) *firstDeviceValue = firstDevice;
    if ( firstHostValue != nullptr ) *firstHostValue = firstHost;
}

void ReportResidentStateContract(
    std::uint64_t topologyGeneration,
    std::uint64_t fieldGeneration,
    ResidentStateSync sync,
    int nCells,
    int nGhostCells,
    int nEquations,
    long mismatches,
    Real maximumDifference,
    long ownedMismatches = 0,
    long ghostMismatches = 0,
    int firstEquation = -1,
    int firstCell = -1,
    Real firstDeviceValue = 0.0,
    Real firstHostValue = 0.0 )
{
    static int printed = 0;
    const bool passed = mismatches == 0;
    if ( ! passed || printed < 3 )
    {
        const char * syncName =
            sync == ResidentStateSync::FullField ? "full"
            : sync == ResidentStateSync::GhostCells ? "ghost" : "none";
        std::cout << "RESIDENT_STATE_CONTRACT rank=" << ResidentContractRank()
                  << " status=" << ( passed ? "PASS" : "FAIL" )
                  << " sync=" << syncName
                  << " topology_generation=" << topologyGeneration
                  << " field_generation=" << fieldGeneration
                  << " cells=" << nCells
                  << " ghost_cells=" << nGhostCells
                  << " equations=" << nEquations
                  << " mismatched_values=" << mismatches
                  << " owned_mismatches=" << ownedMismatches
                  << " ghost_mismatches=" << ghostMismatches
                  << " first_equation=" << firstEquation
                  << " first_cell=" << firstCell
                  << " first_device=" << firstDeviceValue
                  << " first_host=" << firstHostValue
                  << " max_abs_diff=" << maximumDifference
                  << std::endl;
        ++printed;
    }
}

}

void HipFluxBackend::CalcGradient( const CellGradientView & view )
{
    ValidateCellGradientView( view );
    HipGradientStorage * storage = &fallbackGradientStorage_;
    Ns3DDeviceState * owner = nullptr;
    const auto binding = boundStates_.find( view.cacheKey );
    if ( binding != boundStates_.end() )
    {
        owner = binding->second;
        auto * hipState =
            dynamic_cast< HipNs3DBackendState * >( owner->BackendState() );
        if ( hipState == nullptr )
        {
            throw std::logic_error(
                "bound HIP Ns3D state lost its backend ownership" );
        }
        storage = &hipState->Gradient();
    }
    const int totalCells = view.nCells + view.nGhostCells;
    const int totalGradient = view.nEquations * 3 * totalCells;
    storage->cellState.Resize( view.nEquations * totalCells );
    const void * ownerToken =
        owner != nullptr ? owner->OwnerToken() : view.cacheKey;
    const std::uint64_t fieldGeneration =
        owner != nullptr ? owner->FieldGeneration() : 0;
    const std::uint64_t topologyGeneration =
        owner != nullptr ? owner->TopologyGeneration() : 0;
    const ResidentStateSync sync = RefreshResidentPrimitiveState(
        * storage, owner, ownerToken, topologyGeneration, fieldGeneration,
        view.nCells, view.nGhostCells, view.nEquations, view.q );
    if ( owner != nullptr && ResidentContractCheckRequested() )
    {
        // The contract diagnostic compares the resident field against the host
        // mirror.  With a device-authoritative state the mirror is stale by
        // design, so pull it back first; otherwise the check would flag the
        // intended one-sided ownership as a mismatch.
        if ( ! storage->residentContract.HostMirrorCurrent() )
        {
            for ( int equation = 0; equation < view.nEquations; ++ equation )
            {
                AccelRuntime::Instance().Backend().Copy(
                    const_cast< Real * >( view.q[ equation ] ),
                    storage->cellState.Data()
                        + static_cast< std::size_t >( equation ) * totalCells,
                    static_cast< std::size_t >( totalCells ) * sizeof( Real ),
                    AccelCopyKind::DeviceToHost );
                storage->cachedStateQ[ equation ] = view.q[ equation ];
            }
            storage->residentContract.NoteHostMirrorRefresh(
                fieldGeneration );
        }
        Real maximumDifference = 0.0;
        long mismatches = 0;
        long ownedMismatches = 0;
        long ghostMismatches = 0;
        int firstEquation = -1;
        int firstCell = -1;
        Real firstDevice = 0.0;
        Real firstHost = 0.0;
        VerifyResidentPrimitiveState(
            * storage, view.nCells, view.nGhostCells, view.nEquations,
            view.q, 0.0, & maximumDifference, & mismatches,
            & ownedMismatches, & ghostMismatches, & firstEquation,
            & firstCell, & firstDevice, & firstHost );
        ReportResidentStateContract(
            topologyGeneration, fieldGeneration, sync, view.nCells,
            view.nGhostCells, view.nEquations, mismatches,
            maximumDifference, ownedMismatches, ghostMismatches,
            firstEquation, firstCell, firstDevice, firstHost );
    }
    const bool cacheMatches = storage->geometryCacheValid
        && storage->cachedOwnerToken == ownerToken
        && storage->cachedTopologyGeneration == topologyGeneration
        && storage->cachedXFace == view.xFace
        && storage->cachedYFace == view.yFace
        && storage->cachedZFace == view.zFace
        && storage->cachedXNormal == view.xNormal
        && storage->cachedYNormal == view.yNormal
        && storage->cachedZNormal == view.zNormal
        && storage->cachedFaceArea == view.faceArea
        && storage->cachedXCell == view.xCell
        && storage->cachedYCell == view.yCell
        && storage->cachedZCell == view.zCell
        && storage->cachedVolume == view.cellVolume
        && storage->cachedLeftCell == view.leftCell
        && storage->cachedRightCell == view.rightCell
        && storage->cachedCells == view.nCells
        && storage->cachedGhostCells == view.nGhostCells
        && storage->cachedFaces == view.nFaces
        && storage->cachedBoundaryFaces == view.nBoundaryFaces;
    if ( ! cacheMatches )
    {
        storage->xFaceCenter.CopyFromHost( view.xFace, view.nFaces );
        storage->yFaceCenter.CopyFromHost( view.yFace, view.nFaces );
        storage->zFaceCenter.CopyFromHost( view.zFace, view.nFaces );
        storage->xNormal.CopyFromHost( view.xNormal, view.nFaces );
        storage->yNormal.CopyFromHost( view.yNormal, view.nFaces );
        storage->zNormal.CopyFromHost( view.zNormal, view.nFaces );
        storage->faceArea.CopyFromHost( view.faceArea, view.nFaces );
        storage->xCellCenter.CopyFromHost( view.xCell, totalCells );
        storage->yCellCenter.CopyFromHost( view.yCell, totalCells );
        storage->zCellCenter.CopyFromHost( view.zCell, totalCells );
        storage->cellVolume.CopyFromHost( view.cellVolume, view.nCells );
        storage->leftCell.CopyFromHost( view.leftCell, view.nFaces );
        storage->rightCell.CopyFromHost( view.rightCell, view.nFaces );
        storage->cachedOwnerToken = ownerToken;
        storage->cachedTopologyGeneration = topologyGeneration;
        storage->cachedXFace = view.xFace;
        storage->cachedYFace = view.yFace;
        storage->cachedZFace = view.zFace;
        storage->cachedXNormal = view.xNormal;
        storage->cachedYNormal = view.yNormal;
        storage->cachedZNormal = view.zNormal;
        storage->cachedFaceArea = view.faceArea;
        storage->cachedXCell = view.xCell;
        storage->cachedYCell = view.yCell;
        storage->cachedZCell = view.zCell;
        storage->cachedVolume = view.cellVolume;
        storage->cachedLeftCell = view.leftCell;
        storage->cachedRightCell = view.rightCell;
        storage->cachedCells = view.nCells;
        storage->cachedGhostCells = view.nGhostCells;
        storage->cachedFaces = view.nFaces;
        storage->cachedBoundaryFaces = view.nBoundaryFaces;
        storage->geometryCacheValid = true;
    }
    storage->cellGradient.Resize( totalGradient );
    EnsureCellFaceConnectivity(
        *storage, view.leftCell, view.rightCell,
        view.nCells, view.nFaces, view.nBoundaryFaces );
    const int blockSize = 256;
    { ScopedStageTimer kernelTimer( "gradient_kernel" );
    hipLaunchKernelGGL(
        GradientCellKernel,
        dim3( ( view.nCells + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0, storage->cellState.Data(),
        storage->cellGradient.Data(),
        storage->xFaceCenter.Data(), storage->yFaceCenter.Data(),
        storage->zFaceCenter.Data(), storage->xNormal.Data(),
        storage->yNormal.Data(), storage->zNormal.Data(),
        storage->faceArea.Data(), storage->xCellCenter.Data(),
        storage->yCellCenter.Data(), storage->zCellCenter.Data(),
        storage->cellVolume.Data(),
        storage->leftCell.Data(), storage->rightCell.Data(),
        storage->cellFaceOffsets.Data(), storage->cellFaceIds.Data(),
        view.nCells, view.nBoundaryFaces, view.nFaces,
        view.nEquations );
    if ( view.nBoundaryFaces > 0 )
    {
        hipLaunchKernelGGL(
            GradientCopyBoundaryKernel,
            dim3( ( view.nBoundaryFaces + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0, storage->cellGradient.Data(),
            storage->leftCell.Data(), storage->rightCell.Data(),
            view.nCells, view.nBoundaryFaces, view.nEquations );
    }
    CheckHip( hipGetLastError(), "gradient kernel launch" );
    CheckHip( hipDeviceSynchronize(), "gradient kernel synchronize" );
    }
    if ( ! view.downloadToHost ) return;
    // One contiguous device->host transfer plus a host-side scatter instead of
    // 3*neq small transfers with a synchronisation each.  The device layout is
    // [equation][dimension][cell], so the staging buffer is a straight copy of
    // the resident gradient block.
    {
        ScopedStageTimer timer( "gradient_download" );
        const std::size_t totalValues =
            static_cast< std::size_t >( totalGradient );
        const Real * staging =
            LocalPinnedStaging().Get( totalValues );
        AccelRuntime::Instance().Backend().Copy(
            const_cast< Real * >( staging ), storage->cellGradient.Data(),
            totalValues * sizeof( Real ), AccelCopyKind::DeviceToHost );
        for ( int equation = 0; equation < view.nEquations; ++ equation )
        {
            const std::size_t base =
                static_cast< std::size_t >( equation ) * 3 * totalCells;
            std::copy(
                staging + base, staging + base + totalCells,
                view.dqdx[ equation ] );
            std::copy(
                staging + base + totalCells,
                staging + base + 2 * totalCells,
                view.dqdy[ equation ] );
            std::copy(
                staging + base + 2 * totalCells,
                staging + base + 3 * totalCells,
                view.dqdz[ equation ] );
        }
    }
}

void HipFluxBackend::CalcLocalInviscidTimeStep(
    const CellTimeStepView & view )
{
    ValidateCellTimeStepView( view );
    HipGradientStorage * storage = &fallbackGradientStorage_;
    Ns3DDeviceState * owner = nullptr;
    const auto binding = boundStates_.find( view.cacheKey );
    if ( binding != boundStates_.end() )
    {
        owner = binding->second;
        auto * hipState =
            dynamic_cast< HipNs3DBackendState * >( owner->BackendState() );
        if ( hipState == nullptr )
        {
            throw std::logic_error(
                "bound HIP Ns3D state lost its backend ownership" );
        }
        storage = &hipState->Gradient();
    }

    const int totalCells = view.nCells + view.nGhostCells;
    storage->cellState.Resize( view.nEquations * totalCells );
    const void * ownerToken =
        owner != nullptr ? owner->OwnerToken() : view.cacheKey;
    const std::uint64_t fieldGeneration =
        owner != nullptr ? owner->FieldGeneration() : 0;
    const std::uint64_t topologyGeneration =
        owner != nullptr ? owner->TopologyGeneration() : 0;
    // A resident timestep may only reuse q when the generation/ghost contract
    // still matches; otherwise refresh through the same path as gradient so
    // the device never fabricates a timestep from stale halo values.
    if ( view.preferResidentState )
    {
        RefreshResidentPrimitiveState(
            * storage, owner, ownerToken, topologyGeneration,
            fieldGeneration, view.nCells, view.nGhostCells,
            view.nEquations, view.q );
    }
    else
    {
        for ( int equation = 0; equation < view.nEquations; ++ equation )
        {
            AccelRuntime::Instance().Backend().Copy(
                storage->cellState.Data()
                    + static_cast< std::size_t >( equation ) * totalCells,
                view.q[ equation ],
                static_cast< std::size_t >( totalCells ) * sizeof( Real ),
                AccelCopyKind::HostToDevice );
            storage->cachedStateQ[ equation ] = view.q[ equation ];
        }
        storage->residentContract.Bind(
            ownerToken, topologyGeneration, view.nCells,
            view.nGhostCells, view.nEquations );
        storage->residentContract.NoteFullFieldSync( fieldGeneration );
        storage->cachedEquations = view.nEquations;
    }

    const bool geometryMatches =
        storage->timeStepGeometryCacheValid
        && storage->cachedOwnerToken == ownerToken
        && storage->cachedTopologyGeneration == topologyGeneration
        && storage->cachedTimeStepXNormal == view.xNormal
        && storage->cachedTimeStepYNormal == view.yNormal
        && storage->cachedTimeStepZNormal == view.zNormal
        && storage->cachedTimeStepMeshVelocityNormal
            == view.meshVelocityNormal
        && storage->cachedTimeStepFaceArea == view.faceArea
        && storage->cachedTimeStepVolume == view.cellVolume
        && storage->cachedTimeStepLeftCell == view.leftCell
        && storage->cachedTimeStepRightCell == view.rightCell
        && storage->cachedTimeStepCells == view.nCells
        && storage->cachedTimeStepGhostCells == view.nGhostCells
        && storage->cachedTimeStepFaces == view.nFaces
        && storage->cachedTimeStepBoundaryFaces == view.nBoundaryFaces;
    if ( ! geometryMatches )
    {
        storage->xNormal.CopyFromHost( view.xNormal, view.nFaces );
        storage->yNormal.CopyFromHost( view.yNormal, view.nFaces );
        storage->zNormal.CopyFromHost( view.zNormal, view.nFaces );
        storage->meshVelocityNormal.CopyFromHost(
            view.meshVelocityNormal, view.nFaces );
        storage->faceArea.CopyFromHost( view.faceArea, view.nFaces );
        storage->cellVolume.CopyFromHost(
            view.cellVolume, view.nCells );
        storage->leftCell.CopyFromHost( view.leftCell, view.nFaces );
        storage->rightCell.CopyFromHost( view.rightCell, view.nFaces );
        storage->cachedOwnerToken = ownerToken;
        storage->cachedTopologyGeneration = topologyGeneration;
        storage->cachedTimeStepXNormal = view.xNormal;
        storage->cachedTimeStepYNormal = view.yNormal;
        storage->cachedTimeStepZNormal = view.zNormal;
        storage->cachedTimeStepMeshVelocityNormal =
            view.meshVelocityNormal;
        storage->cachedTimeStepFaceArea = view.faceArea;
        storage->cachedTimeStepVolume = view.cellVolume;
        storage->cachedTimeStepLeftCell = view.leftCell;
        storage->cachedTimeStepRightCell = view.rightCell;
        storage->cachedTimeStepCells = view.nCells;
        storage->cachedTimeStepGhostCells = view.nGhostCells;
        storage->cachedTimeStepFaces = view.nFaces;
        storage->cachedTimeStepBoundaryFaces = view.nBoundaryFaces;
        storage->timeStepGeometryCacheValid = true;
    }

    storage->inviscidSpectrum.Resize( view.nCells );
    storage->timeStep.Resize( view.nCells );
    const Real initialMinimum = std::numeric_limits< Real >::max();
    int status = 0;
    storage->minTimeStep.CopyFromHost( & initialMinimum, 1 );
    storage->timeStepStatus.CopyFromHost( & status, 1 );

    const int blockSize = 256;
    hipLaunchKernelGGL(
        InviscidSpectrumZeroKernel,
        dim3( ( view.nCells + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage->inviscidSpectrum.Data(), view.nCells );
    hipLaunchKernelGGL(
        InviscidSpectrumAccumulateKernel,
        dim3( ( view.nFaces + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage->cellState.Data(), storage->xNormal.Data(),
        storage->yNormal.Data(), storage->zNormal.Data(),
        storage->meshVelocityNormal.Data(), storage->faceArea.Data(),
        storage->leftCell.Data(), storage->rightCell.Data(),
        view.nCells, view.nGhostCells, view.nFaces, view.gamma,
        storage->inviscidSpectrum.Data(),
        storage->timeStepStatus.Data() );
    hipLaunchKernelGGL(
        LocalInviscidTimeStepKernel,
        dim3( ( view.nCells + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage->inviscidSpectrum.Data(), storage->cellVolume.Data(),
        view.cfl, view.nCells, storage->timeStep.Data(),
        storage->minTimeStep.Data(), storage->timeStepStatus.Data() );
    CheckHip( hipGetLastError(), "local timestep kernel launch" );
    storage->timeStepStatus.CopyToHost( & status, 1 );
    if ( status != 0 )
    {
        throw std::runtime_error(
            "HIP local timestep produced an invalid spectral radius" );
    }
    if ( view.maxTimeRatio > 0.0 )
    {
        hipLaunchKernelGGL(
            CapLocalTimeStepKernel,
            dim3( ( view.nCells + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0, storage->timeStep.Data(),
            storage->minTimeStep.Data(), view.maxTimeRatio,
            view.nCells );
        CheckHip( hipGetLastError(), "local timestep cap kernel launch" );
    }
    CheckHip(
        hipDeviceSynchronize(), "local timestep kernel synchronize" );

    storage->cachedTimeStepHost = view.timeStep;
    storage->residentTimeStepReady = true;
    if ( view.downloadToHost )
    {
        storage->timeStep.CopyToHost( view.timeStep, view.nCells );
        if ( view.inviscidSpectrum != nullptr )
        {
            storage->inviscidSpectrum.CopyToHost(
                view.inviscidSpectrum, view.nCells );
        }
    }
}

void HipFluxBackend::UpdateBoundaryState( const BoundaryStateView & view )
{
    ValidateBoundaryStateView( view );
    HipGradientStorage * storage = &fallbackGradientStorage_;
    Ns3DDeviceState * owner = nullptr;
    const auto binding = boundStates_.find( view.cacheKey );
    if ( binding != boundStates_.end() )
    {
        owner = binding->second;
        auto * hipState =
            dynamic_cast< HipNs3DBackendState * >( owner->BackendState() );
        if ( hipState == nullptr )
        {
            throw std::logic_error(
                "bound HIP Ns3D state lost its backend ownership" );
        }
        storage = &hipState->Gradient();
    }

    const int totalCells = view.nCells + view.nGhostCells;
    {
    ScopedStageTimer resizeTimer( "bnd_resize" );
    storage->cellState.Resize(
        static_cast< std::size_t >( view.nEquations ) * totalCells );
    }
    const void * ownerToken =
        owner != nullptr ? owner->OwnerToken() : view.cacheKey;
    const std::uint64_t fieldGeneration =
        owner != nullptr ? owner->FieldGeneration() : 0;
    const std::uint64_t topologyGeneration =
        owner != nullptr ? owner->TopologyGeneration() : 0;
    // The device boundary model rewrites ghost cells itself.  It may reuse the
    // resident q only while the generation/ghost contract says it is current.
    if ( view.preferResidentState )
    {
        ScopedStageTimer refreshTimer( "bnd_refresh" );
        RefreshResidentPrimitiveState(
            * storage, owner, ownerToken, topologyGeneration,
            fieldGeneration, view.nCells, view.nGhostCells,
            view.nEquations, view.q );
    }
    else
    {
        for ( int equation = 0; equation < view.nEquations; ++ equation )
        {
            AccelRuntime::Instance().Backend().Copy(
                storage->cellState.Data()
                    + static_cast< std::size_t >( equation ) * totalCells,
                view.q[ equation ],
                static_cast< std::size_t >( totalCells ) * sizeof( Real ),
                AccelCopyKind::HostToDevice );
            storage->cachedStateQ[ equation ] = view.q[ equation ];
        }
        storage->residentContract.Bind(
            ownerToken, topologyGeneration, view.nCells,
            view.nGhostCells, view.nEquations );
        storage->residentContract.NoteFullFieldSync( fieldGeneration );
        storage->cachedEquations = view.nEquations;
    }

    const bool boundaryGeometryMatches =
        storage->boundaryGeometryCacheValid
        && storage->cachedBoundaryOwnerToken == ownerToken
        && storage->cachedBoundaryGeometryGeneration == topologyGeneration
        && storage->cachedBoundaryXNormal == view.xNormal
        && storage->cachedBoundaryYNormal == view.yNormal
        && storage->cachedBoundaryZNormal == view.zNormal
        && storage->cachedBoundaryMeshVelocityNormal
            == view.meshVelocityNormal
        && storage->cachedBoundaryLeftCell == view.leftCell
        && storage->cachedBoundaryRightCell == view.rightCell
        && storage->cachedBoundaryType == view.boundaryType
        && storage->cachedBoundaryFaces == view.nBoundaryFaces;
    if ( ! boundaryGeometryMatches )
    {
        ScopedStageTimer geometryTimer( "bnd_geometry" );
        storage->boundaryXNormal.CopyFromHost(
            view.xNormal, view.nBoundaryFaces );
        storage->boundaryYNormal.CopyFromHost(
            view.yNormal, view.nBoundaryFaces );
        storage->boundaryZNormal.CopyFromHost(
            view.zNormal, view.nBoundaryFaces );
        storage->boundaryMeshVelocityNormal.CopyFromHost(
            view.meshVelocityNormal, view.nBoundaryFaces );
        storage->boundaryLeftCell.CopyFromHost(
            view.leftCell, view.nBoundaryFaces );
        storage->boundaryRightCell.CopyFromHost(
            view.rightCell, view.nBoundaryFaces );
        storage->boundaryType.CopyFromHost(
            view.boundaryType, view.nBoundaryFaces );
        storage->cachedBoundaryOwnerToken = ownerToken;
        storage->cachedBoundaryGeometryGeneration = topologyGeneration;
        storage->cachedBoundaryXNormal = view.xNormal;
        storage->cachedBoundaryYNormal = view.yNormal;
        storage->cachedBoundaryZNormal = view.zNormal;
        storage->cachedBoundaryMeshVelocityNormal =
            view.meshVelocityNormal;
        storage->cachedBoundaryLeftCell = view.leftCell;
        storage->cachedBoundaryRightCell = view.rightCell;
        storage->cachedBoundaryType = view.boundaryType;
        storage->cachedBoundaryFaces = view.nBoundaryFaces;
        storage->boundaryGeometryCacheValid = true;
        static unsigned long long geometryMisses = 0;
        ++ geometryMisses;
        if ( geometryMisses == 1 )
        {
            std::cout << "HIP_BOUNDARY_GEOMETRY_REUPLOAD first miss"
                      << " valid=" << storage->boundaryGeometryCacheValid
                      << std::endl;
        }
    }
    // Inflow is a solver field and may change between calls; keep this copy
    // unconditional even when the static boundary metadata is cached.
    {
    ScopedStageTimer staticTimer( "bnd_static_upload" );
    storage->inflow.CopyFromHost( view.inflow, view.nEquations );
    storage->boundaryFaceState.Resize(
        view.nEquations * view.nBoundaryFaces );
    }

    const int blockSize = 256;
    {
    ScopedStageTimer kernelTimer( "bnd_kernel" );
    hipLaunchKernelGGL(
        UpdateBoundaryPrimitiveKernel,
        dim3( ( view.nBoundaryFaces + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0, storage->cellState.Data(),
        storage->boundaryFaceState.Data(),
        storage->boundaryLeftCell.Data(),
        storage->boundaryRightCell.Data(),
        storage->boundaryType.Data(),
        storage->boundaryXNormal.Data(),
        storage->boundaryYNormal.Data(),
        storage->boundaryZNormal.Data(),
        storage->boundaryMeshVelocityNormal.Data(),
        storage->inflow.Data(), view.faceOuterNormal, view.nCells,
        view.nBoundaryFaces, view.nEquations, view.gamma );
    CheckHip( hipGetLastError(), "boundary kernel launch" );
    // No synchronise here: the packed ghost D2H below is blocking on the same
    // stream, so it already waits for this kernel.
    }

    {
    ScopedStageTimer ghostCopyTimer( "bnd_ghost_copy" );
    // Both copies below used to run once per equation (two blocking D2H each,
    // thirty per step).  Pack on the device, fetch once, scatter on the host.
    const int nGhost = view.nBoundaryFaces;
    const std::size_t ghostValues =
        static_cast< std::size_t >( view.nEquations ) * nGhost;
    const bool mirrorGhostQ = ! view.deviceAuthoritativeHostMirror;
    const bool mirrorBcQ = view.bcQ[ 0 ] != nullptr;
    if ( ghostValues > 0 && ( mirrorGhostQ || mirrorBcQ ) )
    {
        const int blockSizePack = 256;
        if ( mirrorGhostQ )
        {
            storage->boundaryGhostPack.Resize( ghostValues );
            hipLaunchKernelGGL(
                PackGhostCellsKernel,
                dim3( ( static_cast< int >( ghostValues ) + blockSizePack - 1 )
                      / blockSizePack ),
                dim3( blockSizePack ), 0, 0,
                storage->cellState.Data(), view.nCells, totalCells,
                view.nEquations, nGhost,
                storage->boundaryGhostPack.Data() );
            CheckHip( hipGetLastError(), "ghost pack kernel launch" );
            Real * ghostStaging = LocalPinnedStaging().Get( ghostValues );
            AccelRuntime::Instance().Backend().Copy(
                ghostStaging, storage->boundaryGhostPack.Data(),
                ghostValues * sizeof( Real ), AccelCopyKind::DeviceToHost );
            for ( int equation = 0; equation < view.nEquations; ++ equation )
            {
                std::copy(
                    ghostStaging
                        + static_cast< std::size_t >( equation ) * nGhost,
                    ghostStaging
                        + static_cast< std::size_t >( equation + 1 ) * nGhost,
                    view.q[ equation ] + view.nCells );
            }
        }
        if ( mirrorBcQ )
        {
            Real * faceStaging = LocalPinnedStaging().Get( ghostValues );
            AccelRuntime::Instance().Backend().Copy(
                faceStaging, storage->boundaryFaceState.Data(),
                ghostValues * sizeof( Real ), AccelCopyKind::DeviceToHost );
            for ( int equation = 0; equation < view.nEquations; ++ equation )
            {
                std::copy(
                    faceStaging + static_cast< std::size_t >( equation ) * nGhost,
                    faceStaging + static_cast< std::size_t >( equation + 1 ) * nGhost,
                    view.bcQ[ equation ] );
            }
        }
    }
    }
    // The device boundary kernel wrote the ghost cells.  When they were
    // mirrored, both sides now reflect the same ghost epoch; on the
    // device-authoritative path only the device is current.
    storage->residentContract.NoteDeviceGhostWrite(
        ! view.deviceAuthoritativeHostMirror );
}


void HipFluxBackend::ReconstructFaceValues(
    const CellFaceReconstructionView & view,
    bool downloadToHost )
{
    ValidateCellFaceReconstructionView( view );
    if ( view.limiterMode != ReconstructionLimiterMode::Disabled )
    {
        throw std::invalid_argument(
            "HIP device reconstruction currently requires limiter off" );
    }
    const void * cacheKey = view.cacheKey != nullptr
        ? view.cacheKey : view.ownerToken;
    HipGradientStorage * storage = &fallbackGradientStorage_;
    Ns3DDeviceState * owner = nullptr;
    const auto binding = boundStates_.find( cacheKey );
    if ( binding != boundStates_.end() )
    {
        owner = binding->second;
        auto * hipState =
            dynamic_cast< HipNs3DBackendState * >( owner->BackendState() );
        if ( hipState == nullptr )
        {
            throw std::logic_error(
                "bound HIP Ns3D state lost its backend ownership" );
        }
        storage = &hipState->Gradient();
    }
    const std::uint64_t topologyGeneration = owner != nullptr
        ? owner->TopologyGeneration() : view.topologyGeneration;
    if ( storage->cellState.Size() <
             static_cast< std::size_t >( view.nEquations )
                 * ( view.nCells + view.nGhostCells )
         || storage->cellGradient.Size() <
             static_cast< std::size_t >( view.nEquations ) * 3
                 * ( view.nCells + view.nGhostCells )
         || ! storage->geometryCacheValid
         || storage->cachedFaces != view.nFaces
         || storage->cachedBoundaryFaces != view.nBoundaryFaces )
    {
        throw std::logic_error(
            "HIP reconstruction requires a preceding compatible gradient" );
    }

    const int totalCells = view.nCells + view.nGhostCells;
    const int totalFaces = view.nFaces;
    const int totalValues = view.nEquations * totalFaces;
    if ( view.syncGradientGhostsFromHost && view.nGhostCells > 0 )
    {
        // Cross-zone interface ghosts: the host mirror holds the exchanged
        // neighbor gradient after Grad::SwapBcGrad.  Refresh only the ghost
        // slice of the state-owned device gradient (three components per
        // equation) before reconstructing interface faces.  The source is a
        // pinned staging copy so the upload does not go through the driver's
        // pageable staging buffer.
        ScopedStageTimer syncTimer( "halo_ghost_sync" );
        const int nGhostCells = view.nGhostCells;
        const std::size_t ghostValues =
            static_cast< std::size_t >( view.nEquations ) * 3 * nGhostCells;
        Real * pinned = LocalPinnedStaging().Get( ghostValues );
        std::size_t offset = 0;
        for ( int equation = 0; equation < view.nEquations; ++ equation )
        {
            const Real * hostGradient[ 3 ] = {
                view.dqdx[ equation ], view.dqdy[ equation ],
                view.dqdz[ equation ] };
            for ( int component = 0; component < 3; ++ component )
            {
                std::copy(
                    hostGradient[ component ] + view.nCells,
                    hostGradient[ component ] + view.nCells + nGhostCells,
                    pinned + offset );
                offset += static_cast< std::size_t >( nGhostCells );
            }
        }
        offset = 0;
        for ( int equation = 0; equation < view.nEquations; ++ equation )
        {
            const std::size_t base =
                static_cast< std::size_t >( equation ) * 3 * totalCells;
            for ( int component = 0; component < 3; ++ component )
            {
                AccelRuntime::Instance().Backend().Copy(
                    storage->cellGradient.Data() + base
                        + static_cast< std::size_t >( component ) * totalCells
                        + view.nCells,
                    pinned + offset,
                    static_cast< std::size_t >( nGhostCells ) * sizeof( Real ),
                    AccelCopyKind::HostToDevice );
                offset += static_cast< std::size_t >( nGhostCells );
            }
        }
    }
    storage->qf1.Resize( totalValues );
    storage->qf2.Resize( totalValues );
    const bool operationChanged =
        storage->cachedBoundaryOperationHost != view.boundaryOperation
        || storage->cachedBoundaryOperationCount
               != static_cast< std::size_t >( totalFaces )
        || storage->cachedBoundaryTopologyGeneration
               != topologyGeneration;
    storage->boundaryOperation.Resize( totalFaces );
    if ( operationChanged && view.boundaryOperation != nullptr )
    {
        storage->boundaryOperation.CopyFromHost(
            view.boundaryOperation, totalFaces );
    }
    else if ( operationChanged )
    {
        std::vector< ReconstructionBoundaryOperation > defaultOperation(
            totalFaces, ReconstructionBoundaryOperation::Preserve );
        storage->boundaryOperation.CopyFromHost(
            defaultOperation.data(), totalFaces );
    }
    if ( operationChanged )
    {
        storage->cachedBoundaryOperationHost = view.boundaryOperation;
        storage->cachedBoundaryOperationCount = totalFaces;
        storage->cachedBoundaryTopologyGeneration =
            topologyGeneration;
    }
    if ( view.boundaryMask != nullptr )
    {
        storage->boundaryMask.CopyFromHost(
            view.boundaryMask, totalFaces );
        storage->boundaryMaskValid = true;
    }
    else
    {
        storage->boundaryMask.Reset();
        storage->boundaryMaskValid = false;
    }
    const bool hasSolid = view.hasSolidBoundary;
    if ( hasSolid )
    {
        // The host boundary-state buffer and the device copy share the same
        // equation-major [equation][face] layout, so one contiguous transfer
        // replaces nEqu blocking copies.  This transfer is the whole
        // reconstruction imbalance on the 8-zone case: only wall-adjacent zones
        // have solid faces, and they pay nEqu*totalFaces values per stage
        // (5.6 ms/step on those ranks, 0 on the others).
        ScopedStageTimer bcqTimer( "recon_bcq" );
        const std::size_t bcQValues =
            static_cast< std::size_t >( view.nEquations )
            * view.nBoundaryFaces;
        storage->bcQ.Resize( bcQValues );
        AccelRuntime::Instance().Backend().Copy(
            storage->bcQ.Data(), view.bcQ[ 0 ],
            bcQValues * sizeof( Real ),
            AccelCopyKind::HostToDevice );
    }
    else
    {
        storage->bcQ.Reset();
    }
    storage->cachedEquations = view.nEquations;
    storage->cachedFieldGeneration = view.fieldGeneration;

    const int blockSize = 256;
    ScopedStageTimer reconTimer( "recon_kernel" );
    hipLaunchKernelGGL(
        ReconstructFaceEulerKernel,
        dim3( ( totalFaces + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage->cellState.Data(), storage->cellGradient.Data(),
        storage->xFaceCenter.Data(), storage->yFaceCenter.Data(),
        storage->zFaceCenter.Data(), storage->xCellCenter.Data(),
        storage->yCellCenter.Data(), storage->zCellCenter.Data(),
        storage->leftCell.Data(), storage->rightCell.Data(),
        storage->boundaryMaskValid ? storage->boundaryMask.Data() : nullptr,
        storage->boundaryOperation.Data(),
        hasSolid ? storage->bcQ.Data() : nullptr,
        view.nBoundaryFaces,
        view.nCells, view.nBoundaryFaces, view.nFaces, view.nEquations,
        view.physicality == ReconstructionPhysicalityPolicy::None ? 0 : 1,
        view.densityComponent, view.pressureComponent,
        storage->qf1.Data(), storage->qf2.Data() );
    CheckHip( hipGetLastError(), "reconstruction kernel launch" );
    CheckHip( hipDeviceSynchronize(), "reconstruction kernel synchronize" );

    if ( downloadToHost )
    {
        for ( int equation = 0; equation < view.nEquations; ++ equation )
        {
            AccelRuntime::Instance().Backend().Copy(
                view.qLeft[ equation ], storage->qf1.Data() + equation * totalFaces,
                totalFaces * sizeof( Real ), AccelCopyKind::DeviceToHost );
            AccelRuntime::Instance().Backend().Copy(
                view.qRight[ equation ], storage->qf2.Data() + equation * totalFaces,
                totalFaces * sizeof( Real ), AccelCopyKind::DeviceToHost );
        }
    }
}

void HipFluxBackend::CalcAndAddReconstructedFaceFlux(
    const FaceConnectivityView & connectivity,
    ResidualView & residual,
    int scheme,
    Real gamma,
    FaceFluxView * hostFlux,
    const void * cacheKey,
    bool downloadResidual )
{
    ValidateFaceConnectivityView( connectivity );
    if ( residual.nCells <= 0 || residual.nEquations < 3 || gamma <= 1.0 )
    {
        throw std::invalid_argument(
            "invalid reconstructed HIP residual request" );
    }
    const void * key = cacheKey;
    if ( key == nullptr )
    {
        throw std::invalid_argument(
            "reconstructed HIP flux requires a bound cache key" );
    }
    auto binding = boundStates_.find( key );
    HipGradientStorage * storage = &fallbackGradientStorage_;
    if ( binding != boundStates_.end() )
    {
        auto * hipState = dynamic_cast< HipNs3DBackendState * >(
            binding->second->BackendState() );
        if ( hipState == nullptr )
            throw std::logic_error( "bound HIP state is not HIP-owned" );
        storage = &hipState->Gradient();
    }
    const int nFaces = storage->cachedFaces;
    const int nEq = storage->cachedEquations;
    const int nCells = residual.nCells;
    if ( nFaces <= 0 || nEq < 3 || residual.nEquations != nEq
         || connectivity.nFaces != nFaces
         || storage->qf1.Size() < static_cast< std::size_t >( nFaces * nEq )
         || storage->qf2.Size() < static_cast< std::size_t >( nFaces * nEq )
         || storage->cachedBoundaryFaces != connectivity.nBoundaryFaces )
    {
        throw std::invalid_argument(
            "reconstructed HIP flux dimensions do not match state" );
    }
    if ( storage->cachedLeftCell != connectivity.leftCell
         || storage->cachedRightCell != connectivity.rightCell )
    {
        throw std::invalid_argument(
            "reconstructed HIP flux connectivity does not match state" );
    }
    if ( hostFlux != nullptr
         && ( hostFlux->nFaces != nFaces || hostFlux->nEquations != nEq
              || hostFlux->values == nullptr ) )
    {
        throw std::invalid_argument( "invalid reconstructed HIP host flux" );
    }

    const int totalResidual = nCells * nEq;
    storage->residual.Resize( totalResidual );
    if ( residual.values != nullptr )
    {
        storage->residual.CopyFromHost(
            residual.values, totalResidual );
    }
    else
    {
        for ( int equation = 0; equation < nEq; ++ equation )
        {
            AccelRuntime::Instance().Backend().Copy(
                storage->residual.Data() + equation * nCells,
                residual.components[ equation ], nCells * sizeof( Real ),
                AccelCopyKind::HostToDevice );
        }
    }
    storage->residualCells = nCells;
    storage->residualEquations = nEq;
    storage->residualValid = true;
    storage->residualMetricsValid = false;
    storage->faceFlux.Resize( nFaces * nEq );
    EnsureCellFaceConnectivity(
        *storage, connectivity.leftCell, connectivity.rightCell,
        nCells, nFaces, connectivity.nBoundaryFaces );

    const unsigned char * deviceBoundaryMask = nullptr;
    if ( connectivity.boundaryMask != nullptr )
    {
        boundaryMask_.CopyFromHost(
            connectivity.boundaryMask, nFaces );
        deviceBoundaryMask = boundaryMask_.Data();
    }
    const int blockSize = 256;
    { ScopedStageTimer kernelTimer( "HIP kernel" );
    hipLaunchKernelGGL(
        CalcAndScatterPrimitiveEulerKernel,
        dim3( ( nFaces + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage->qf1.Data(), storage->qf2.Data(),
        storage->xNormal.Data(), storage->yNormal.Data(),
        storage->zNormal.Data(), nullptr, storage->faceArea.Data(),
        storage->leftCell.Data(), storage->rightCell.Data(),
        deviceBoundaryMask, connectivity.nBoundaryFaces, nFaces, nCells,
        nEq, gamma, scheme,
        storage->faceFlux.Data(),
        storage->residual.Data() );
    hipLaunchKernelGGL(
        GatherFaceFluxEulerKernel,
        dim3( ( nCells + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage->faceFlux.Data(),
        storage->leftCell.Data(), storage->rightCell.Data(),
        deviceBoundaryMask,
        storage->cellFaceOffsets.Data(), storage->cellFaceIds.Data(),
        connectivity.nBoundaryFaces, nFaces, nCells, nEq,
        storage->residual.Data() );
    CheckHip( hipGetLastError(),
              "reconstructed primitive flux kernel launch" );
    CheckHip( hipDeviceSynchronize(),
              "reconstructed primitive flux kernel synchronize" );
    }
    if ( downloadResidual && residual.values != nullptr )
    {
        storage->residual.CopyToHost( residual.values, totalResidual );
    }
    else if ( downloadResidual )
    {
        for ( int equation = 0; equation < nEq; ++ equation )
        {
            AccelRuntime::Instance().Backend().Copy(
                residual.components[ equation ],
                storage->residual.Data() + equation * nCells,
                nCells * sizeof( Real ), AccelCopyKind::DeviceToHost );
        }
    }
    if ( hostFlux != nullptr )
    {
        storage->faceFlux.CopyToHost( hostFlux->values, nFaces * nEq );
    }
}

void HipFluxBackend::ScaleCurrentResidual(
    const CellStateUpdateView & view )
{
    ValidateCellStateUpdateView( view );
    const auto binding = boundStates_.find( view.cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP state update requires a bound Ns3D state" );
    }
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error(
            "bound HIP state update is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    if ( ! storage.residualValid
         || storage.residualCells != view.nCells
         || storage.residualEquations != view.nEquations
         || storage.cachedCells != view.nCells
         || storage.cachedGhostCells != view.nGhostCells
         || storage.cachedVolume != view.cellVolume
         || storage.cellVolume.Size()
             < static_cast< std::size_t >( view.nCells ) )
    {
        throw std::logic_error(
            "HIP state update requires a compatible resident residual" );
    }

    const bool residentTimeStepUsable =
        view.preferResidentTimeStep
        && storage.residentTimeStepReady
        && storage.cachedTimeStepHost == view.timeStep
        && storage.cachedTimeStepCells == view.nCells;
    if ( ! residentTimeStepUsable )
    {
        storage.timeStep.CopyFromHost( view.timeStep, view.nCells );
        storage.cachedTimeStepHost = view.timeStep;
        storage.cachedTimeStepCells = view.nCells;
        storage.residentTimeStepReady = true;
    }
    const int total = view.nCells * view.nEquations;
    const int blockSize = 256;
    hipLaunchKernelGGL(
        ScaleResidualKernel,
        dim3( ( total + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage.residual.Data(), storage.timeStep.Data(),
        storage.cellVolume.Data(), view.nCells, view.nEquations,
        view.rkCoefficient );
    CheckHip( hipGetLastError(), "state residual scale kernel launch" );
    CheckHip(
        hipDeviceSynchronize(),
        "state residual scale kernel synchronize" );
    storage.residualMetricsValid = true;
}

void HipFluxBackend::UpdatePrimitiveState(
    const CellStateUpdateView & view,
    bool downloadToHost )
{
    ValidateCellStateUpdateView( view );
    const auto binding = boundStates_.find( view.cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP primitive update requires a bound Ns3D state" );
    }
    Ns3DDeviceState * owner = binding->second;
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        owner->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error(
            "bound HIP primitive update is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const int totalCells = view.nCells + view.nGhostCells;
    if ( ! storage.residualValid
         || storage.residualCells != view.nCells
         || storage.residualEquations != view.nEquations
         || storage.cellState.Size()
             < static_cast< std::size_t >(
                 totalCells * view.nEquations ) )
    {
        throw std::logic_error(
            "HIP primitive update requires resident state and residual" );
    }

    int status = 0;
    storage.stateUpdateStatus.CopyFromHost( & status, 1 );
    const int blockSize = 256;
    { ScopedStageTimer updateTimer( "state_update" );
    hipLaunchKernelGGL(
        UpdatePrimitiveStateKernel,
        dim3( ( view.nCells + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage.cellState.Data(), storage.residual.Data(),
        view.nCells, totalCells, view.gamma,
        storage.stateUpdateStatus.Data() );
    CheckHip( hipGetLastError(), "primitive state update kernel launch" );
    storage.stateUpdateStatus.CopyToHost( & status, 1 );
    if ( status != 0 )
    {
        throw std::runtime_error(
            "HIP primitive state update produced a non-physical state" );
    }

    if ( downloadToHost )
    {
        for ( int equation = 0; equation < view.nEquations; ++ equation )
        {
            AccelRuntime::Instance().Backend().Copy(
                view.primitive[ equation ],
                storage.cellState.Data() + equation * totalCells,
                view.nCells * sizeof( Real ),
                AccelCopyKind::DeviceToHost );
        }
    }
    }
    owner->MarkFieldChanged();
    storage.residualValid = false;
    // The contract may only claim "device matches host" when the owned-cell
    // mirror was actually refreshed and the solver is not about to write the
    // host state itself.  Otherwise force a full refresh on the next read.
    // A device-authoritative update keeps the field on the device and marks the
    // host mirror stale; host readers pull back exactly what they need.
    if ( view.retainResidentState && downloadToHost )
    {
        storage.residentContract.NoteDeviceFieldUpdate(
            owner->FieldGeneration() );
    }
    else if ( view.retainResidentState && view.deviceAuthoritative )
    {
        storage.residentContract.NoteDeviceFieldAdvance(
            owner->FieldGeneration() );
    }
    else
    {
        storage.residentContract.Invalidate();
    }
}

BoundaryBufferView HipFluxBackend::PackBoundary(
    const void * cacheKey,
    const BoundaryFieldExchangePlan & plan,
    BoundaryTransferMode mode )
{
    plan.Validate();
    return this->PackBoundary(
        cacheKey,
        { plan.sourceCellCount, plan.sendCellIds },
        plan.nEquations,
        mode );
}

BoundaryBufferView HipFluxBackend::PackBoundary(
    const void * cacheKey,
    const BoundaryEndpointMap & endpoint,
    int nEquations,
    BoundaryTransferMode mode )
{
    return this->PackBoundary(
        cacheKey, endpoint, nEquations, mode, { -1, -1, 0 } );
}

BoundaryBufferView HipFluxBackend::PackBoundary(
    const void * cacheKey,
    const BoundaryEndpointMap & endpoint,
    int nEquations,
    BoundaryTransferMode mode,
    const HipBoundaryEndpointKey & endpointKey )
{
    endpoint.Validate();
    if ( nEquations <= 0 )
    {
        throw std::invalid_argument(
            "HIP boundary pack equation count must be positive" );
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP boundary pack requires a bound Ns3D state" );
    }
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error(
            "bound HIP boundary pack state is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const int cellCount = storage.residentContract.Cells()
        + storage.residentContract.GhostCells();
    if ( ! storage.residentContract.Resident()
         || cellCount != endpoint.cellCount
         || binding->second->problem.nEquations != nEquations )
    {
        throw std::logic_error(
            "HIP boundary pack requires a compatible resident state" );
    }

    HipBoundaryEndpointStorage & boundary =
        storage.boundarySend[ endpointKey ];
    const std::uint64_t topologyGeneration =
        binding->second->TopologyGeneration();
    const bool endpointChanged =
        boundary.topologyGeneration != topologyGeneration
        || boundary.cellCount != endpoint.cellCount
        || boundary.nEquations != nEquations
        || boundary.role != endpoint.role
        || boundary.cachedCellIds != endpoint.cellIds;
    if ( endpointChanged )
    {
        if ( ! endpoint.cellIds.empty() )
        {
            boundary.cellIds.CopyFromHost(
                endpoint.cellIds.data(), endpoint.cellIds.size() );
        }
        else
        {
            boundary.cellIds.Resize( 0 );
        }
        boundary.cachedCellIds = endpoint.cellIds;
        boundary.topologyGeneration = topologyGeneration;
        boundary.cellCount = endpoint.cellCount;
        boundary.nEquations = nEquations;
        boundary.role = endpoint.role;
    }
    const std::size_t valueCount =
        endpoint.cellIds.size() * static_cast< std::size_t >( nEquations );
    boundary.values.Resize( valueCount );
    if ( valueCount != 0 )
    {
        const int kernelValueCount = static_cast< int >( valueCount );
        const int blockSize = 256;
        hipLaunchKernelGGL(
            PackBoundaryKernel,
            dim3( ( kernelValueCount + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0,
            storage.cellState.Data(), boundary.cellIds.Data(),
            cellCount, static_cast< int >( endpoint.cellIds.size() ),
            nEquations, boundary.values.Data() );
        CheckHip( hipGetLastError(), "boundary pack kernel launch" );
    }

    if ( mode == BoundaryTransferMode::DeviceDirect )
    {
        CheckHip(
            hipDeviceSynchronize(),
            "boundary pack kernel synchronize" );
        return {
            boundary.values.Data(),
            valueCount,
            BoundaryMemorySpace::Device };
    }

    if ( boundary.hostValues.size() != valueCount )
        boundary.hostValues.resize( valueCount );
    if ( valueCount != 0 )
    {
        boundary.values.CopyToHost( boundary.hostValues.data(), valueCount );
    }
    return {
        boundary.hostValues.data(),
        valueCount,
        BoundaryMemorySpace::Host };
}

void HipFluxBackend::UploadBoundaryState(
    const void * cacheKey,
    int nCells,
    int nGhostCells,
    int nEquations,
    const Real * const * components )
{
    if ( cacheKey == nullptr || nCells <= 0 || nGhostCells < 0
         || ( nEquations != 3 && nEquations != 5 )
         || components == nullptr )
    {
        throw std::invalid_argument(
            "invalid HIP boundary state upload" );
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP boundary state upload requires a bound Ns3D state" );
    }
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error(
            "bound HIP boundary upload state is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const int totalCells = nCells + nGhostCells;
    storage.cellState.Resize(
        static_cast< std::size_t >( nEquations ) * totalCells );
    const void * ownerToken = binding->second->OwnerToken();
    const std::uint64_t fieldGeneration =
        binding->second->FieldGeneration();
    const std::uint64_t topologyGeneration =
        binding->second->TopologyGeneration();
    // Packing only reads owned cells, so a pending ghost-only refresh does not
    // force a full upload here.  A changed owned-cell generation or a new
    // binding does -- unless the contract shows the device already owns that
    // newer generation (device-authoritative state update).  Without this the
    // interface exchange would overwrite the resident field with the stale
    // host mirror once per neighbor and per stage.
    const bool generationSatisfied =
        storage.residentContract.Resident()
        && storage.residentContract.FieldGeneration() == fieldGeneration;
    const bool deviceOwnsNewerField =
        storage.residentContract.Resident()
        && storage.residentContract.RequiredSync( fieldGeneration )
            == ResidentStateSync::None;
    bool sameResidentState =
        storage.residentContract.Matches(
            ownerToken, topologyGeneration, nCells, nGhostCells,
            nEquations )
        && storage.residentContract.Resident()
        && ( generationSatisfied || deviceOwnsNewerField );
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        if ( components[ equation ] == nullptr )
        {
            throw std::invalid_argument(
                "HIP boundary state component is null" );
        }
        if ( storage.cachedStateQ[ equation ] != components[ equation ] )
            sameResidentState = false;
    }
    if ( sameResidentState )
    {
        // The HIP RK path already owns the current q state.  Interface
        // exchange must consume that resident allocation instead of copying
        // the complete field once for every neighbor.
        return;
    }
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        AccelRuntime::Instance().Backend().Copy(
            storage.cellState.Data()
                + static_cast< std::size_t >( equation ) * totalCells,
            components[ equation ],
            static_cast< std::size_t >( totalCells ) * sizeof( Real ),
            AccelCopyKind::HostToDevice );
        storage.cachedStateQ[ equation ] = components[ equation ];
    }
    storage.residentContract.Bind(
        ownerToken, topologyGeneration, nCells, nGhostCells, nEquations );
    storage.residentContract.NoteFullFieldSync( fieldGeneration );
    storage.cachedCells = nCells;
    storage.cachedGhostCells = nGhostCells;
    storage.cachedEquations = nEquations;
}

BoundaryBufferView HipFluxBackend::PrepareBoundaryReceive(
    const void * cacheKey,
    const BoundaryFieldExchangePlan & plan,
    BoundaryTransferMode mode )
{
    plan.Validate();
    return this->PrepareBoundaryReceive(
        cacheKey,
        { plan.targetCellCount, plan.receiveCellIds },
        plan.nEquations,
        mode );
}

BoundaryBufferView HipFluxBackend::PrepareBoundaryReceive(
    const void * cacheKey,
    const BoundaryEndpointMap & endpoint,
    int nEquations,
    BoundaryTransferMode mode )
{
    return this->PrepareBoundaryReceive(
        cacheKey, endpoint, nEquations, mode, { -1, -1, 0 } );
}

BoundaryBufferView HipFluxBackend::PrepareBoundaryReceive(
    const void * cacheKey,
    const BoundaryEndpointMap & endpoint,
    int nEquations,
    BoundaryTransferMode mode,
    const HipBoundaryEndpointKey & endpointKey )
{
    endpoint.Validate();
    if ( nEquations <= 0 )
    {
        throw std::invalid_argument(
            "HIP boundary receive equation count must be positive" );
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP boundary receive requires a bound Ns3D state" );
    }
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error(
            "bound HIP boundary receive state is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const int cellCount = storage.residentContract.Cells()
        + storage.residentContract.GhostCells();
    if ( ! storage.residentContract.Resident()
         || cellCount != endpoint.cellCount
         || binding->second->problem.nEquations != nEquations )
    {
        throw std::logic_error(
            "HIP boundary receive requires a compatible resident state" );
    }

    HipBoundaryEndpointStorage & boundary =
        storage.boundaryReceive[ endpointKey ];
    const std::uint64_t topologyGeneration =
        binding->second->TopologyGeneration();
    const bool endpointChanged =
        boundary.topologyGeneration != topologyGeneration
        || boundary.cellCount != endpoint.cellCount
        || boundary.nEquations != nEquations
        || boundary.role != endpoint.role
        || boundary.cachedCellIds != endpoint.cellIds;
    if ( endpointChanged )
    {
        if ( ! endpoint.cellIds.empty() )
        {
            boundary.cellIds.CopyFromHost(
                endpoint.cellIds.data(), endpoint.cellIds.size() );
        }
        else
        {
            boundary.cellIds.Resize( 0 );
        }
        boundary.cachedCellIds = endpoint.cellIds;
        boundary.topologyGeneration = topologyGeneration;
        boundary.cellCount = endpoint.cellCount;
        boundary.nEquations = nEquations;
        boundary.role = endpoint.role;
    }
    const std::size_t valueCount =
        endpoint.cellIds.size() * static_cast< std::size_t >( nEquations );
    boundary.values.Resize( valueCount );
    if ( mode == BoundaryTransferMode::DeviceDirect )
    {
        return {
            boundary.values.Data(),
            valueCount,
            BoundaryMemorySpace::Device };
    }

    if ( boundary.hostValues.size() != valueCount )
        boundary.hostValues.resize( valueCount );
    return {
        boundary.hostValues.data(),
        valueCount,
        BoundaryMemorySpace::Host };
}

void HipFluxBackend::UnpackBoundary(
    const void * cacheKey,
    const BoundaryFieldExchangePlan & plan,
    BoundaryTransferMode mode )
{
    plan.Validate();
    this->UnpackBoundary(
        cacheKey,
        { plan.targetCellCount, plan.receiveCellIds },
        plan.nEquations,
        mode );
}

void HipFluxBackend::UnpackBoundary(
    const void * cacheKey,
    const BoundaryEndpointMap & endpoint,
    int nEquations,
    BoundaryTransferMode mode )
{
    this->UnpackBoundary(
        cacheKey, endpoint, nEquations, mode, { -1, -1, 0 } );
}

void HipFluxBackend::UnpackBoundary(
    const void * cacheKey,
    const BoundaryEndpointMap & endpoint,
    int nEquations,
    BoundaryTransferMode mode,
    const HipBoundaryEndpointKey & endpointKey )
{
    endpoint.Validate();
    if ( nEquations <= 0 )
    {
        throw std::invalid_argument(
            "HIP boundary unpack equation count must be positive" );
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP boundary unpack requires a bound Ns3D state" );
    }
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error(
            "bound HIP boundary unpack state is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const int cellCount = storage.residentContract.Cells()
        + storage.residentContract.GhostCells();
    const std::size_t valueCount =
        endpoint.cellIds.size() * static_cast< std::size_t >( nEquations );
    const auto boundaryIter =
        storage.boundaryReceive.find( endpointKey );
    if ( boundaryIter == storage.boundaryReceive.end() )
    {
        throw std::logic_error(
            "HIP boundary unpack has no persistent receive endpoint" );
    }
    HipBoundaryEndpointStorage & boundary = boundaryIter->second;
    if ( ! storage.residentContract.Resident()
         || cellCount != endpoint.cellCount
         || binding->second->problem.nEquations != nEquations
         || boundary.topologyGeneration
             != binding->second->TopologyGeneration()
         || boundary.cellCount != endpoint.cellCount
         || boundary.nEquations != nEquations
         || boundary.role != endpoint.role
         || boundary.cachedCellIds != endpoint.cellIds
         || boundary.values.Size() != valueCount
         || boundary.cellIds.Size() != endpoint.cellIds.size() )
    {
        throw std::logic_error(
            "HIP boundary unpack requires a prepared receive buffer" );
    }
    if ( mode == BoundaryTransferMode::HostStaged )
    {
        if ( boundary.hostValues.size() != valueCount )
        {
            throw std::logic_error(
                "HIP host-staged boundary receive buffer is not prepared" );
        }
        if ( valueCount != 0 )
        {
            boundary.values.CopyFromHost(
                boundary.hostValues.data(), valueCount );
        }
    }

    if ( valueCount != 0 )
    {
        const int kernelValueCount = static_cast< int >( valueCount );
        const int blockSize = 256;
        hipLaunchKernelGGL(
            UnpackBoundaryKernel,
            dim3( ( kernelValueCount + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0,
            boundary.values.Data(),
            boundary.cellIds.Data(),
            cellCount, static_cast< int >( endpoint.cellIds.size() ),
            nEquations, storage.cellState.Data() );
        CheckHip( hipGetLastError(), "boundary unpack kernel launch" );
        CheckHip(
            hipDeviceSynchronize(),
            "boundary unpack kernel synchronize" );
    }
    // The paired host exchange wrote the same ghost cells, and those values
    // were just received into the state-owned device q.  The ghost epoch is now
    // current for this endpoint; any host mutation that was not covered by a
    // paired exchange leaves the pending flag set and is repaired by the next
    // contract sync.
    storage.residentContract.NoteDeviceGhostUpdate();
}

void HipFluxBackend::NotifyHostFieldStateModified( const void * cacheKey )
{
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() ) return;
    // Owned cells changed outside the device path: bump the owner generation
    // so the next device consumer performs a full upload.
    binding->second->MarkFieldChanged();
}

void HipFluxBackend::NotifyHostGhostStateModified( const void * cacheKey )
{
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() ) return;
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr ) return;
    hipState->Gradient().residentContract.NoteHostGhostMutation();
}

bool HipFluxBackend::SyncResidentGhostStateFromHost(
    const void * cacheKey,
    int nCells,
    int nGhostCells,
    int nEquations,
    const Real * const * components )
{
    if ( cacheKey == nullptr || nCells <= 0 || nGhostCells < 0
         || nEquations <= 0 || components == nullptr )
    {
        throw std::invalid_argument(
            "invalid HIP resident ghost sync request" );
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() ) return false;
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr ) return false;
    HipGradientStorage & storage = hipState->Gradient();
    ResidentStateContract & contract = storage.residentContract;
    if ( ! contract.Matches(
             binding->second->OwnerToken(),
             binding->second->TopologyGeneration(),
             nCells, nGhostCells, nEquations ) )
    {
        return false;
    }
    if ( ! contract.HostGhostPending() ) return true;
    if ( nGhostCells == 0 )
    {
        contract.NoteGhostSync();
        return true;
    }

    const int totalCells = nCells + nGhostCells;
    if ( storage.cellState.Size()
         < static_cast< std::size_t >( nEquations ) * totalCells )
    {
        return false;
    }
    { ScopedStageTimer syncTimer( "halo_ghost_sync" );
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        if ( components[ equation ] == nullptr )
        {
            throw std::invalid_argument(
                "HIP resident ghost sync component is null" );
        }
        AccelRuntime::Instance().Backend().Copy(
            storage.cellState.Data()
                + static_cast< std::size_t >( equation ) * totalCells
                + nCells,
            components[ equation ] + nCells,
            static_cast< std::size_t >( nGhostCells ) * sizeof( Real ),
            AccelCopyKind::HostToDevice );
        storage.cachedStateQ[ equation ] = components[ equation ];
    }
    }
    contract.NoteGhostSync();
    return true;
}

bool HipFluxBackend::HostMirrorStale( const void * cacheKey ) const
{
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() ) return false;
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr ) return false;
    const ResidentStateContract & contract =
        hipState->Gradient().residentContract;
    return contract.Resident() && ! contract.HostMirrorCurrent();
}

bool HipFluxBackend::RefreshHostPrimitiveState(
    const void * cacheKey,
    int nCells,
    int nGhostCells,
    int nEquations,
    Real * const * components )
{
    if ( cacheKey == nullptr || nCells <= 0 || nGhostCells < 0
         || nEquations <= 0 || components == nullptr )
    {
        throw std::invalid_argument(
            "invalid HIP host primitive refresh request" );
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() ) return false;
    Ns3DDeviceState * owner = binding->second;
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        owner->BackendState() );
    if ( hipState == nullptr ) return false;
    HipGradientStorage & storage = hipState->Gradient();
    ResidentStateContract & contract = storage.residentContract;
    if ( ! contract.Matches(
             owner->OwnerToken(), owner->TopologyGeneration(),
             nCells, nGhostCells, nEquations )
         || ! contract.Resident() )
    {
        return false;
    }
    if ( contract.HostMirrorCurrent() ) return false;

    const int totalCells = nCells + nGhostCells;
    if ( storage.cellState.Size()
         < static_cast< std::size_t >( nEquations ) * totalCells )
    {
        return false;
    }
    { ScopedStageTimer refreshTimer( "host_mirror_refresh" );
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        if ( components[ equation ] == nullptr )
        {
            throw std::invalid_argument(
                "HIP host primitive refresh component is null" );
        }
        AccelRuntime::Instance().Backend().Copy(
            components[ equation ],
            storage.cellState.Data()
                + static_cast< std::size_t >( equation ) * totalCells,
            static_cast< std::size_t >( totalCells ) * sizeof( Real ),
            AccelCopyKind::DeviceToHost );
        storage.cachedStateQ[ equation ] = components[ equation ];
    }
    }
    contract.NoteHostMirrorRefresh( owner->FieldGeneration() );
    return true;
}

void HipFluxBackend::MirrorPrimitiveCellsToHost(
    const void * cacheKey,
    int nEquations,
    const int * cellIds,
    int cellCount,
    Real * const * components )
{
    if ( cacheKey == nullptr || nEquations <= 0 || cellIds == nullptr
         || cellCount <= 0 || components == nullptr )
    {
        return;
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() ) return;
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr ) return;
    HipGradientStorage & storage = hipState->Gradient();
    const ResidentStateContract & contract = storage.residentContract;
    if ( ! contract.Resident() || contract.HostMirrorCurrent() ) return;
    const int totalCells = contract.Cells() + contract.GhostCells();
    if ( totalCells <= 0
         || storage.cellState.Size()
             < static_cast< std::size_t >( nEquations ) * totalCells )
    {
        return;
    }
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        if ( components[ equation ] == nullptr ) return;
    }
    const std::size_t packedValues =
        static_cast< std::size_t >( nEquations ) * cellCount;
    storage.stateMirrorPack.Resize( packedValues );
    storage.stateMirrorIds.CopyFromHost( cellIds, cellCount );
    const int blockSize = 256;
    hipLaunchKernelGGL(
        PackStateCellsKernel,
        dim3( ( static_cast< int >( packedValues ) + blockSize - 1 )
              / blockSize ),
        dim3( blockSize ), 0, 0,
        storage.cellState.Data(), storage.stateMirrorIds.Data(),
        cellCount, nEquations, totalCells,
        storage.stateMirrorPack.Data() );
    CheckHip( hipGetLastError(), "state cell pack kernel launch" );
    // One blocking D2H covers every equation; the scatter below is host work.
    Real * staging = LocalPinnedStaging().Get( packedValues );
    AccelRuntime::Instance().Backend().Copy(
        staging, storage.stateMirrorPack.Data(),
        packedValues * sizeof( Real ), AccelCopyKind::DeviceToHost );
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        const Real * source =
            staging + static_cast< std::size_t >( equation ) * cellCount;
        Real * target = components[ equation ];
        for ( int index = 0; index < cellCount; ++ index )
        {
            const int cell = cellIds[ index ];
            if ( cell < 0 || cell >= contract.Cells() )
            {
                throw std::logic_error(
                    "HIP host mirror cell is outside the owned range" );
            }
            target[ cell ] = source[ index ];
        }
    }
}

void HipFluxBackend::MarkResidentHostMirrorStale( const void * cacheKey )
{
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() ) return;
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr ) return;
    hipState->Gradient().residentContract.NoteHostMirrorStale();
}


// Residual diagnostics: one block per contiguous cell chunk, serial inside the
// chunk so the summation order is fixed for a given chunking.  Only scalars
// leave the device; the residual field itself stays resident.
__global__ void ResidualMetricsKernel(
    const Real * residual,
    const Real * q,
    int nCells,
    int totalCells,
    int nEquations,
    int chunkCells,
    Real * chunkSum,
    Real * chunkMax,
    int * chunkIndex,
    Real * chunkMin )
{
    const int chunk = blockIdx.x;
    const int chunks = gridDim.x;
    if ( threadIdx.x != 0 ) return;
    const int begin = chunk * chunkCells;
    const int end = min( begin + chunkCells, nCells );
    if ( begin >= end ) return;
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        Real sum = 0.0;
        Real maximum = 0.0;
        int index = begin;
        const Real * values = residual + equation * nCells;
        for ( int cell = begin; cell < end; ++ cell )
        {
            const Real value = values[ cell ];
            sum += value * value;
            const Real magnitude = fabs( value );
            if ( magnitude > maximum )
            {
                maximum = magnitude;
                index = cell;
            }
        }
        chunkSum[ equation * chunks + chunk ] = sum;
        chunkMax[ equation * chunks + chunk ] = maximum;
        chunkIndex[ equation * chunks + chunk ] = index;
    }
    // Primitive positivity extremes (rho = component 0, p = last component).
    // Primitive state is [equation][owned + ghost cell]; only the owned cells
    // participate in the positivity report.
    Real minDensity = q[ begin ];
    Real minPressure = q[ ( nEquations - 1 ) * totalCells + begin ];
    for ( int cell = begin + 1; cell < end; ++ cell )
    {
        const Real density = q[ cell ];
        const Real pressure = q[ ( nEquations - 1 ) * totalCells + cell ];
        if ( density < minDensity ) minDensity = density;
        if ( pressure < minPressure ) minPressure = pressure;
    }
    chunkMin[ chunk ] = minDensity;
    chunkMin[ chunks + chunk ] = minPressure;
}

bool HipFluxBackend::ComputeResidualMetrics(
    const void * cacheKey,
    int nCells,
    int nEquations,
    Real gamma,
    HipResidualMetrics & metrics )
{
    ( void ) gamma;
    metrics.valid = false;
    if ( cacheKey == nullptr || nCells <= 0 || nEquations <= 0 ) return false;
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() ) return false;
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr ) return false;
    HipGradientStorage & storage = hipState->Gradient();
    if ( ! storage.residualMetricsValid || storage.residualCells != nCells
         || storage.residualEquations != nEquations
         || storage.residual.Size()
             < static_cast< std::size_t >( nCells ) * nEquations )
    {
        return false;
    }
    ScopedStageTimer timer( "residual_metrics" );
    const int chunkCells = 256;
    const int chunks = ( nCells + chunkCells - 1 ) / chunkCells;
    if ( storage.residualMetricChunks != chunks
         || storage.residualMetricEquations != nEquations )
    {
        storage.residualMetricChunkSum.Resize(
            static_cast< std::size_t >( chunks ) * nEquations );
        storage.residualMetricChunkMax.Resize(
            static_cast< std::size_t >( chunks ) * nEquations );
        storage.residualMetricChunkIndex.Resize(
            static_cast< std::size_t >( chunks ) * nEquations );
        storage.residualMetricChunkMin.Resize(
            static_cast< std::size_t >( 2 ) * chunks );
        storage.residualMetricChunks = chunks;
        storage.residualMetricEquations = nEquations;
    }
    hipLaunchKernelGGL(
        ResidualMetricsKernel, dim3( chunks ), dim3( 1 ), 0, 0,
        storage.residual.Data(), storage.cellState.Data(), nCells,
        storage.residentContract.Cells() + storage.residentContract.GhostCells(),
        nEquations, chunkCells, storage.residualMetricChunkSum.Data(),
        storage.residualMetricChunkMax.Data(),
        storage.residualMetricChunkIndex.Data(),
        storage.residualMetricChunkMin.Data() );
    CheckHip( hipGetLastError(), "residual metrics kernel launch" );
    CheckHip( hipDeviceSynchronize(), "residual metrics kernel synchronize" );

    std::vector< Real > chunkSum(
        static_cast< std::size_t >( chunks ) * nEquations );
    std::vector< Real > chunkMax( chunkSum.size() );
    std::vector< int > chunkIndex( chunkSum.size() );
    std::vector< Real > chunkMin( static_cast< std::size_t >( 2 ) * chunks );
    storage.residualMetricChunkSum.CopyToHost( chunkSum.data(), chunkSum.size() );
    storage.residualMetricChunkMax.CopyToHost( chunkMax.data(), chunkMax.size() );
    storage.residualMetricChunkIndex.CopyToHost(
        chunkIndex.data(), chunkIndex.size() );
    storage.residualMetricChunkMin.CopyToHost( chunkMin.data(), chunkMin.size() );

    metrics.sumOfSquares.assign( static_cast< std::size_t >( nEquations ), 0.0 );
    metrics.maxAbsolute.assign( static_cast< std::size_t >( nEquations ), 0.0 );
    metrics.argmaxCell.assign( static_cast< std::size_t >( nEquations ), 0 );
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        Real sum = 0.0;
        Real maximum = 0.0;
        int index = 0;
        for ( int chunk = 0; chunk < chunks; ++ chunk )
        {
            const std::size_t at =
                static_cast< std::size_t >( equation ) * chunks + chunk;
            sum += chunkSum[ at ];
            if ( chunkMax[ at ] > maximum )
            {
                maximum = chunkMax[ at ];
                index = chunkIndex[ at ];
            }
        }
        metrics.sumOfSquares[ equation ] = sum;
        metrics.maxAbsolute[ equation ] = maximum;
        metrics.argmaxCell[ equation ] = index;
    }
    metrics.minDensity = chunkMin[ 0 ];
    metrics.minPressure = chunkMin[ chunks ];
    for ( int chunk = 1; chunk < chunks; ++ chunk )
    {
        if ( chunkMin[ chunk ] < metrics.minDensity )
            metrics.minDensity = chunkMin[ chunk ];
        if ( chunkMin[ chunks + chunk ] < metrics.minPressure )
            metrics.minPressure = chunkMin[ chunks + chunk ];
    }
    metrics.nonPhysicalCells = 0;
    if ( ! ( metrics.minDensity > 0.0 ) || ! ( metrics.minPressure > 0.0 ) )
        metrics.nonPhysicalCells = 1;
    metrics.valid = true;
    return true;
}

void HipFluxBackend::MirrorGradientCellsToHost(
    const void * cacheKey,
    int nEquations,
    const std::vector< int > & cellIds,
    Real * const * dqdx,
    Real * const * dqdy,
    Real * const * dqdz )
{
    if ( cacheKey == nullptr || nEquations <= 0 || dqdx == nullptr
         || dqdy == nullptr || dqdz == nullptr )
    {
        throw std::invalid_argument(
            "invalid HIP gradient mirror request" );
    }
    if ( cellIds.empty() ) return;
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP gradient mirror requires a bound Ns3D state" );
    }
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error(
            "bound HIP gradient mirror state is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const int totalCells = storage.residentContract.Cells()
        + storage.residentContract.GhostCells();
    if ( totalCells <= 0
         || storage.cellGradient.Size()
             < static_cast< std::size_t >( totalCells ) * 3 * nEquations )
    {
        throw std::logic_error(
            "HIP gradient mirror requires a resident gradient" );
    }
    {
        ScopedStageTimer timer( "gradient_mirror" );
        const int cellCount = static_cast< int >( cellIds.size() );
        if ( storage.cachedGradientMirrorCellIds != cellIds
             || storage.cachedGradientMirrorEquations != nEquations )
        {
            storage.gradientMirrorCellIds.CopyFromHost(
                cellIds.data(), cellIds.size() );
            storage.cachedGradientMirrorCellIds = cellIds;
            storage.cachedGradientMirrorEquations = nEquations;
        }
        const int planes = 3 * nEquations;
        const std::size_t valueCount =
            static_cast< std::size_t >( planes ) * cellCount;
        storage.gradientMirrorValues.Resize( valueCount );
        const int blockSize = 256;
        hipLaunchKernelGGL(
            PackBoundaryKernel,
            dim3( ( static_cast< int >( valueCount ) + blockSize - 1 )
                  / blockSize ),
            dim3( blockSize ), 0, 0,
            storage.cellGradient.Data(),
            storage.gradientMirrorCellIds.Data(),
            totalCells, cellCount, planes,
            storage.gradientMirrorValues.Data() );
        CheckHip( hipGetLastError(), "gradient mirror pack kernel launch" );
        Real * pinned = LocalPinnedStaging().Get( valueCount );
        storage.gradientMirrorValues.CopyToHost( pinned, valueCount );
        for ( int plane = 0; plane < planes; ++ plane )
        {
            const int equation = plane / 3;
            Real * destination = plane % 3 == 0
                ? dqdx[ equation ]
                : ( plane % 3 == 1 ? dqdy[ equation ] : dqdz[ equation ] );
            const Real * source =
                pinned + static_cast< std::size_t >( plane ) * cellCount;
            for ( int index = 0; index < cellCount; ++ index )
            {
                destination[ cellIds[ index ] ] = source[ index ];
            }
        }
    }
}

void HipFluxBackend::PackGradientFacesToHost(
    const void * cacheKey, const int * cellIds, int nFaces, int nEquations,
    Real * hostPacked )
{
    if ( cacheKey == nullptr || cellIds == nullptr || nFaces <= 0
         || nEquations <= 0 || hostPacked == nullptr )
    {
        throw std::invalid_argument( "invalid HIP gradient face pack request" );
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error( "HIP gradient face pack requires a bound state" );
    }
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error( "bound HIP state is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const int totalCells = storage.residentContract.Cells()
        + storage.residentContract.GhostCells();
    if ( totalCells <= 0
         || storage.cellGradient.Size()
             < static_cast< std::size_t >( totalCells ) * 3 * nEquations )
    {
        throw std::logic_error( "HIP gradient face pack has no resident gradient" );
    }
    const int valueCount = nFaces * 3 * nEquations;
    storage.gradientPackIds.CopyFromHost( cellIds, nFaces );
    storage.gradientPackValues.Resize( valueCount );
    const int blockSize = 256;
    hipLaunchKernelGGL(
        PackGradientFacesKernel,
        dim3( ( valueCount + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage.cellGradient.Data(), storage.gradientPackIds.Data(),
        nFaces, nEquations, totalCells,
        storage.gradientPackValues.Data() );
    CheckHip( hipGetLastError(), "gradient face pack kernel launch" );
    // CopyToHost is a blocking hipMemcpy on the same stream, so it already
    // waits for the pack kernel.  An explicit device synchronize here only
    // serialized the pipeline once more per interface ghost level and stage.
    storage.gradientPackValues.CopyToHost( hostPacked, valueCount );
}

void HipFluxBackend::UploadGradientFacesToDevice(
    const void * cacheKey, int nFaces, int nEquations, const Real * hostPacked )
{
    if ( cacheKey == nullptr || nFaces <= 0 || nEquations <= 0
         || hostPacked == nullptr )
    {
        throw std::invalid_argument(
            "invalid HIP gradient receive upload request" );
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP gradient receive upload requires a bound state" );
    }
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error( "bound HIP state is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const std::size_t valueCount =
        static_cast< std::size_t >( nFaces ) * 3 * nEquations;
    storage.gradientPackValues.CopyFromHost( hostPacked, valueCount );
}

void HipFluxBackend::UnpackGradientFacesFromDevice(
    const void * cacheKey, const int * ghostCellIds, int nFaces,
    int nEquations )
{
    if ( cacheKey == nullptr || ghostCellIds == nullptr || nFaces <= 0
         || nEquations <= 0 )
    {
        throw std::invalid_argument(
            "invalid HIP gradient receive unpack request" );
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP gradient receive unpack requires a bound state" );
    }
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error( "bound HIP state is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const int totalCells = storage.residentContract.Cells()
        + storage.residentContract.GhostCells();
    if ( totalCells <= 0
         || storage.cellGradient.Size()
             < static_cast< std::size_t >( totalCells ) * 3 * nEquations )
    {
        throw std::logic_error(
            "HIP gradient receive unpack has no resident gradient" );
    }
    const std::size_t needed =
        static_cast< std::size_t >( nFaces ) * 3 * nEquations;
    if ( storage.gradientPackValues.Size() < needed )
    {
        // Before the first exchange round there is no payload; the host path
        // scattered a zero-initialised receive buffer, so mirror that.
        storage.gradientPackValues.Resize( needed );
        CheckHip( hipMemset( storage.gradientPackValues.Data(), 0,
                             needed * sizeof( Real ) ),
                  "gradient receive zero fill" );
    }
    storage.gradientPackIds.CopyFromHost( ghostCellIds, nFaces );
    const int valueCount = nFaces * 3 * nEquations;
    const int blockSize = 256;
    hipLaunchKernelGGL(
        UnpackGradientFacesKernel,
        dim3( ( valueCount + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage.cellGradient.Data(), storage.gradientPackIds.Data(),
        nFaces, nEquations, totalCells,
        storage.gradientPackValues.Data() );
    CheckHip( hipGetLastError(), "gradient receive unpack kernel launch" );
}

void HipFluxBackend::UnpackGradientFacesFromHost(
    const void * cacheKey, const int * ghostCellIds, int nFaces,
    int nEquations, const Real * hostPacked )
{
    if ( cacheKey == nullptr || ghostCellIds == nullptr || nFaces <= 0
         || nEquations <= 0 || hostPacked == nullptr )
    {
        throw std::invalid_argument( "invalid HIP gradient face unpack request" );
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP gradient face unpack requires a bound state" );
    }
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error( "bound HIP state is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const int totalCells = storage.residentContract.Cells()
        + storage.residentContract.GhostCells();
    if ( totalCells <= 0
         || storage.cellGradient.Size()
             < static_cast< std::size_t >( totalCells ) * 3 * nEquations )
    {
        throw std::logic_error(
            "HIP gradient face unpack has no resident gradient" );
    }
    const int valueCount = nFaces * 3 * nEquations;
    storage.gradientPackIds.CopyFromHost( ghostCellIds, nFaces );
    storage.gradientPackValues.CopyFromHost( hostPacked, valueCount );
    const int blockSize = 256;
    hipLaunchKernelGGL(
        UnpackGradientFacesKernel,
        dim3( ( valueCount + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        storage.cellGradient.Data(), storage.gradientPackIds.Data(),
        nFaces, nEquations, totalCells,
        storage.gradientPackValues.Data() );
    CheckHip( hipGetLastError(), "gradient face unpack kernel launch" );
    // No trailing synchronize: the next device operation is stream-ordered
    // after the unpack, and no host reader touches the ghost gradient here.
    // The old synchronize turned every stage into a full pipeline drain.
}

bool HipFluxBackend::VerifyResidentStateAgainstHost(
    const void * cacheKey,
    int nCells,
    int nGhostCells,
    int nEquations,
    const Real * const * components,
    Real tolerance,
    Real * maxAbsoluteDifference,
    long * mismatchedValues )
{
    if ( maxAbsoluteDifference != nullptr ) *maxAbsoluteDifference = 0.0;
    if ( mismatchedValues != nullptr ) *mismatchedValues = 0;
    if ( cacheKey == nullptr || nCells <= 0 || nGhostCells < 0
         || nEquations <= 0 || components == nullptr )
    {
        return false;
    }
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() ) return false;
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        binding->second->BackendState() );
    if ( hipState == nullptr ) return false;
    HipGradientStorage & storage = hipState->Gradient();
    if ( storage.cellState.Size()
         < static_cast< std::size_t >( nCells + nGhostCells ) * nEquations )
    {
        return false;
    }
    VerifyResidentPrimitiveState(
        storage, nCells, nGhostCells, nEquations, components, tolerance,
        maxAbsoluteDifference, mismatchedValues );
    return true;
}

HipFluxBackend & HipFluxBackend::Shared()
{
    if ( ! sharedHipFluxBackend )
    {
        if ( ! AccelRuntime::Instance().IsInitialized() )
        {
            throw std::runtime_error(
                "cannot create shared HIP flux backend before runtime initialization." );
        }
        sharedHipFluxBackend = std::make_unique< HipFluxBackend >();
    }
    return *sharedHipFluxBackend;
}

void HipFluxBackend::ReleaseShared()
{
    sharedHipFluxBackend.reset();
}

void HipFluxBackend::CalcInvFlux(
    const FaceStateView & state,
    FaceFluxView & flux,
    int scheme )
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
        qLeft_.CopyFromHost( state.qLeft, nFaces );
        qRight_.CopyFromHost( state.qRight, nFaces );
        if ( state.xNormal ) xNormal_.CopyFromHost( state.xNormal, nFaces );
        if ( state.faceArea ) area_.CopyFromHost( state.faceArea, nFaces );
        deviceFlux_.Resize( nFaces );

        const int blockSize = 256;
        const int gridSize = ( nFaces + blockSize - 1 ) / blockSize;
        hipLaunchKernelGGL(
            CalcInvFluxScalarKernel,
            dim3( gridSize ),
            dim3( blockSize ),
            0,
            0,
            qLeft_.Data(),
            qRight_.Data(),
            state.xNormal ? xNormal_.Data() : nullptr,
            state.faceArea ? area_.Data() : nullptr,
            deviceFlux_.Data(),
            nFaces );
        CheckHip( hipGetLastError(), "CalcInvFluxScalarKernel launch" );
        CheckHip( hipDeviceSynchronize(), "CalcInvFluxScalarKernel synchronize" );
        deviceFlux_.CopyToHost( flux.values, nFaces );
        deviceFluxFaces_ = nFaces;
        deviceFluxEquations_ = nEq;
    }
    else if ( nEq >= 3 )
    {
        qLeft_.CopyFromHost( state.qLeft, totalValues );
        qRight_.CopyFromHost( state.qRight, totalValues );
        if ( state.xNormal ) xNormal_.CopyFromHost( state.xNormal, nFaces );
        if ( state.yNormal ) yNormal_.CopyFromHost( state.yNormal, nFaces );
        if ( state.zNormal ) zNormal_.CopyFromHost( state.zNormal, nFaces );
        if ( state.meshVelocityNormal )
            meshVelocityNormal_.CopyFromHost( state.meshVelocityNormal, nFaces );
        if ( state.faceArea ) area_.CopyFromHost( state.faceArea, nFaces );
        deviceFlux_.Resize( totalValues );

        const int blockSize = 256;
        const int gridSize = ( nFaces + blockSize - 1 ) / blockSize;
        hipLaunchKernelGGL(
            CalcInvFluxEulerKernel,
            dim3( gridSize ),
            dim3( blockSize ),
            0,
            0,
            qLeft_.Data(),
            qRight_.Data(),
            state.xNormal ? xNormal_.Data() : nullptr,
            state.yNormal ? yNormal_.Data() : nullptr,
            state.zNormal ? zNormal_.Data() : nullptr,
            state.meshVelocityNormal ? meshVelocityNormal_.Data() : nullptr,
            state.faceArea ? area_.Data() : nullptr,
            deviceFlux_.Data(),
            nFaces,
            nEq,
            state.gamma,
            scheme );
        CheckHip( hipGetLastError(), "CalcInvFluxEulerKernel launch" );
        CheckHip( hipDeviceSynchronize(), "CalcInvFluxEulerKernel synchronize" );
        deviceFlux_.CopyToHost( flux.values, totalValues );
        deviceFluxFaces_ = nFaces;
        deviceFluxEquations_ = nEq;
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
        DeviceBuffer< unsigned char > boundaryMask( nFaces );
        DeviceBuffer< Real > deviceResidual( nCells );

        deviceFlux.CopyFromHost( flux.values, nFaces );
        leftCell.CopyFromHost( connectivity.leftCell, nFaces );
        rightCell.CopyFromHost( connectivity.rightCell, nFaces );
        if ( connectivity.boundaryMask )
            boundaryMask.CopyFromHost( connectivity.boundaryMask, nFaces );
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
            connectivity.boundaryMask ? boundaryMask.Data() : nullptr,
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
        DeviceBuffer< unsigned char > boundaryMask( nFaces );
        DeviceBuffer< Real > deviceResidual( totalResidual );

        deviceFlux.CopyFromHost( flux.values, totalFlux );
        leftCell.CopyFromHost( connectivity.leftCell, nFaces );
        rightCell.CopyFromHost( connectivity.rightCell, nFaces );
        if ( connectivity.boundaryMask )
            boundaryMask.CopyFromHost( connectivity.boundaryMask, nFaces );
        deviceResidual.CopyFromHost( residual.values, totalResidual );

        const int blockSize = 256;
        const int gridSize = ( nFaces + blockSize - 1 ) / blockSize;
        hipLaunchKernelGGL(
            ScatterFaceFluxEulerKernel,
            dim3( gridSize ),
            dim3( blockSize ),
            0,
            0,
            deviceFlux.Data(),
            leftCell.Data(),
            rightCell.Data(),
            connectivity.boundaryMask ? boundaryMask.Data() : nullptr,
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

void HipFluxBackend::AddCurrentFaceFlux(
    const FaceConnectivityView & connectivity,
    ResidualView & residual )
{
    ValidateFaceConnectivityView( connectivity );
    if ( deviceFluxFaces_ != connectivity.nFaces
         || deviceFluxEquations_ != residual.nEquations
         || deviceFluxEquations_ < 3
         || residual.nCells <= 0
         || residual.values == nullptr )
    {
        throw std::invalid_argument(
            "HIP cached residual request does not match the current flux." );
    }

    const int nFaces = connectivity.nFaces;
    const int nCells = residual.nCells;
    const int nEq = residual.nEquations;
    const int totalResidual = nCells * nEq;

    const bool connectivityCacheMatches =
        geometryCacheValid_
        && cachedFaces_ == nFaces
        && cachedBoundaryFaces_ == connectivity.nBoundaryFaces
        && cachedLeftCell_ == connectivity.leftCell
        && cachedRightCell_ == connectivity.rightCell
        && cachedBoundaryMask_ == connectivity.boundaryMask;
    if ( ! connectivityCacheMatches )
    {
        ScopedStageTimer timer( "geometry_connectivity_H2D" );
        leftCell_.CopyFromHost( connectivity.leftCell, nFaces );
        rightCell_.CopyFromHost( connectivity.rightCell, nFaces );
        if ( connectivity.boundaryMask )
            boundaryMask_.CopyFromHost( connectivity.boundaryMask, nFaces );
        cachedLeftCell_ = connectivity.leftCell;
        cachedRightCell_ = connectivity.rightCell;
        cachedBoundaryMask_ = connectivity.boundaryMask;
        cachedFaces_ = nFaces;
        cachedBoundaryFaces_ = connectivity.nBoundaryFaces;
        geometryCacheValid_ = true;
    }
    deviceResidual_.Resize( totalResidual );
    { ScopedStageTimer timer( "H2D" );
        deviceResidual_.CopyFromHost( residual.values, totalResidual );
    }

    const int blockSize = 256;
    const int gridSize = ( nFaces + blockSize - 1 ) / blockSize;
    { ScopedStageTimer timer( "HIP kernel" );
    hipLaunchKernelGGL(
        ScatterFaceFluxEulerKernel,
        dim3( gridSize ),
        dim3( blockSize ),
        0,
        0,
        deviceFlux_.Data(),
        leftCell_.Data(),
        rightCell_.Data(),
        connectivity.boundaryMask ? boundaryMask_.Data() : nullptr,
        connectivity.nBoundaryFaces,
        nFaces,
        nCells,
        nEq,
        deviceResidual_.Data() );
    CheckHip( hipGetLastError(), "ScatterFaceFluxEulerKernel launch" );
    CheckHip(
        hipDeviceSynchronize(),
        "ScatterFaceFluxEulerKernel synchronize" );
    }
    { ScopedStageTimer timer( "D2H" );
        deviceResidual_.CopyToHost( residual.values, totalResidual );
    }
}

void HipFluxBackend::CalcAndAddPrimitiveFaceFlux(
    const PrimitiveFaceStateView & state,
    const FaceConnectivityView & connectivity,
    ResidualView & residual,
    int scheme,
    FaceFluxView * hostFlux,
    const void * cacheKey )
{
    const bool componentInput =
        state.primitiveLeft == nullptr && state.primitiveRight == nullptr
        && state.primitiveLeftComponents[ 0 ] != nullptr
        && state.primitiveRightComponents[ 0 ] != nullptr;
    if ( state.nFaces <= 0 || state.nEquations < 3
         || ( !componentInput
              && ( state.primitiveLeft == nullptr
                   || state.primitiveRight == nullptr ) )
         || state.gamma <= 1.0
         || state.nFaces != connectivity.nFaces
         || residual.nCells <= 0
         || residual.nEquations != state.nEquations
         || ( residual.values == nullptr
              && residual.components[ 0 ] == nullptr ) )
    {
        throw std::invalid_argument(
            "invalid fused HIP primitive residual request." );
    }
    if ( hostFlux != nullptr
         && ( hostFlux->nFaces != state.nFaces
              || hostFlux->nEquations != state.nEquations
              || hostFlux->values == nullptr ) )
    {
        throw std::invalid_argument(
            "invalid fused HIP host flux request." );
    }
    ValidateFaceConnectivityView( connectivity );

    const int nFaces = state.nFaces;
    const int nCells = residual.nCells;
    const int nEq = state.nEquations;
    const int totalFlux = nFaces * nEq;
    const int totalResidual = nCells * nEq;

    HipGradientStorage * boundStorage = nullptr;
    if ( cacheKey != nullptr )
    {
        const auto binding = boundStates_.find( cacheKey );
        if ( binding != boundStates_.end() )
        {
            auto * hipState = dynamic_cast< HipNs3DBackendState * >(
                binding->second->BackendState() );
            if ( hipState == nullptr )
            {
                throw std::logic_error(
                    "bound HIP primitive state is not HIP-owned" );
            }
            boundStorage = &hipState->Gradient();
        }
    }

    DeviceBuffer< Real > * leftBuffer = &qLeft_;
    DeviceBuffer< Real > * rightBuffer = &qRight_;
    DeviceBuffer< Real > * xNormalBuffer = &xNormal_;
    DeviceBuffer< Real > * yNormalBuffer = &yNormal_;
    DeviceBuffer< Real > * zNormalBuffer = &zNormal_;
    DeviceBuffer< Real > * meshVelocityBuffer = &meshVelocityNormal_;
    DeviceBuffer< Real > * areaBuffer = &area_;
    DeviceBuffer< int > * leftCellBuffer = &leftCell_;
    DeviceBuffer< int > * rightCellBuffer = &rightCell_;
    DeviceBuffer< unsigned char > * boundaryMaskBuffer = &boundaryMask_;
    DeviceBuffer< Real > * residualBuffer = &deviceResidual_;
    DeviceBuffer< Real > * fluxBuffer = &deviceFlux_;
    if ( boundStorage != nullptr )
    {
        leftBuffer = &boundStorage->qf1;
        rightBuffer = &boundStorage->qf2;
        xNormalBuffer = &boundStorage->xNormal;
        yNormalBuffer = &boundStorage->yNormal;
        zNormalBuffer = &boundStorage->zNormal;
        meshVelocityBuffer = &boundStorage->meshVelocityNormal;
        areaBuffer = &boundStorage->faceArea;
        leftCellBuffer = &boundStorage->leftCell;
        rightCellBuffer = &boundStorage->rightCell;
        boundaryMaskBuffer = &boundStorage->boundaryMask;
        residualBuffer = &boundStorage->residual;
        fluxBuffer = &boundStorage->faceFlux;
    }

    if ( componentInput )
    {
        ScopedStageTimer timer( "H2D" );
        leftBuffer->Resize( totalFlux );
        rightBuffer->Resize( totalFlux );
        for ( int equation = 0; equation < nEq; ++ equation )
        {
            AccelRuntime::Instance().Backend().Copy(
                leftBuffer->Data() + equation * nFaces,
                state.primitiveLeftComponents[ equation ],
                nFaces * sizeof( Real ), AccelCopyKind::HostToDevice );
            AccelRuntime::Instance().Backend().Copy(
                rightBuffer->Data() + equation * nFaces,
                state.primitiveRightComponents[ equation ],
                nFaces * sizeof( Real ), AccelCopyKind::HostToDevice );
        }
    }
    else
    {
        ScopedStageTimer timer( "H2D" );
        leftBuffer->CopyFromHost( state.primitiveLeft, totalFlux );
        rightBuffer->CopyFromHost( state.primitiveRight, totalFlux );
    }

    const void * ownerToken = cacheKey;
    const bool geometryCacheMatches =
        boundStorage != nullptr
        ? ( boundStorage->fluxGeometryCacheValid
            && boundStorage->cachedFluxFaces == nFaces
            && boundStorage->cachedFluxBoundaryFaces
                == connectivity.nBoundaryFaces
            && boundStorage->cachedFluxOwnerToken == ownerToken
            && boundStorage->cachedFluxXNormal == state.xNormal
            && boundStorage->cachedFluxYNormal == state.yNormal
            && boundStorage->cachedFluxZNormal == state.zNormal
            && boundStorage->cachedFluxMeshVelocityNormal
                == state.meshVelocityNormal
            && boundStorage->cachedFluxFaceArea == state.faceArea
            && boundStorage->cachedFluxLeftCell == connectivity.leftCell
            && boundStorage->cachedFluxRightCell == connectivity.rightCell
            && boundStorage->cachedFluxBoundaryMask
                == connectivity.boundaryMask )
        : ( geometryCacheValid_
            && cachedFaces_ == nFaces
            && cachedBoundaryFaces_ == connectivity.nBoundaryFaces
            && cachedOwner_ == cacheKey
            && cachedXNormal_ == state.xNormal
            && cachedYNormal_ == state.yNormal
            && cachedZNormal_ == state.zNormal
            && cachedMeshVelocityNormal_ == state.meshVelocityNormal
            && cachedArea_ == state.faceArea
            && cachedLeftCell_ == connectivity.leftCell
            && cachedRightCell_ == connectivity.rightCell
            && cachedBoundaryMask_ == connectivity.boundaryMask );
    if ( ! geometryCacheMatches )
    {
        ScopedStageTimer timer( "geometry_connectivity_H2D" );
        if ( state.xNormal )
            xNormalBuffer->CopyFromHost( state.xNormal, nFaces );
        if ( state.yNormal )
            yNormalBuffer->CopyFromHost( state.yNormal, nFaces );
        if ( state.zNormal )
            zNormalBuffer->CopyFromHost( state.zNormal, nFaces );
        if ( state.meshVelocityNormal )
        {
            meshVelocityBuffer->CopyFromHost(
                state.meshVelocityNormal, nFaces );
        }
        if ( state.faceArea )
            areaBuffer->CopyFromHost( state.faceArea, nFaces );
        leftCellBuffer->CopyFromHost( connectivity.leftCell, nFaces );
        rightCellBuffer->CopyFromHost( connectivity.rightCell, nFaces );
        if ( connectivity.boundaryMask )
        {
            boundaryMaskBuffer->CopyFromHost(
                connectivity.boundaryMask, nFaces );
        }
        if ( boundStorage != nullptr )
        {
            boundStorage->cachedFluxXNormal = state.xNormal;
            boundStorage->cachedFluxYNormal = state.yNormal;
            boundStorage->cachedFluxZNormal = state.zNormal;
            boundStorage->cachedFluxMeshVelocityNormal =
                state.meshVelocityNormal;
            boundStorage->cachedFluxFaceArea = state.faceArea;
            boundStorage->cachedFluxLeftCell = connectivity.leftCell;
            boundStorage->cachedFluxRightCell = connectivity.rightCell;
            boundStorage->cachedFluxBoundaryMask =
                connectivity.boundaryMask;
            boundStorage->cachedFluxOwnerToken = ownerToken;
            boundStorage->cachedFluxFaces = nFaces;
            boundStorage->cachedFluxBoundaryFaces =
                connectivity.nBoundaryFaces;
            boundStorage->fluxGeometryCacheValid = true;
        }
        else
        {
            cachedXNormal_ = state.xNormal;
            cachedYNormal_ = state.yNormal;
            cachedZNormal_ = state.zNormal;
            cachedMeshVelocityNormal_ = state.meshVelocityNormal;
            cachedArea_ = state.faceArea;
            cachedLeftCell_ = connectivity.leftCell;
            cachedRightCell_ = connectivity.rightCell;
            cachedBoundaryMask_ = connectivity.boundaryMask;
            cachedOwner_ = cacheKey;
            cachedFaces_ = nFaces;
            cachedBoundaryFaces_ = connectivity.nBoundaryFaces;
            geometryCacheValid_ = true;
        }
    }
    residualBuffer->Resize( totalResidual );
    { ScopedStageTimer timer( "H2D" );
    if ( residual.values != nullptr )
    {
        residualBuffer->CopyFromHost( residual.values, totalResidual );
    }
    else
    {
        for ( int equation = 0; equation < nEq; ++ equation )
        {
            AccelRuntime::Instance().Backend().Copy(
                residualBuffer->Data() + equation * nCells,
                residual.components[ equation ],
                nCells * sizeof( Real ), AccelCopyKind::HostToDevice );
        }
    }
    fluxBuffer->Resize( totalFlux );
    }

    HipGradientStorage * cellStorage =
        boundStorage != nullptr ? boundStorage : &fallbackGradientStorage_;
    EnsureCellFaceConnectivity(
        *cellStorage, connectivity.leftCell, connectivity.rightCell,
        nCells, nFaces, connectivity.nBoundaryFaces );

    const int blockSize = 256;
    const int gridSize = ( nFaces + blockSize - 1 ) / blockSize;
    { ScopedStageTimer kernelTimer( "HIP kernel" );
    hipLaunchKernelGGL(
        CalcAndScatterPrimitiveEulerKernel,
        dim3( gridSize ),
        dim3( blockSize ),
        0,
        0,
        leftBuffer->Data(),
        rightBuffer->Data(),
        state.xNormal ? xNormalBuffer->Data() : nullptr,
        state.yNormal ? yNormalBuffer->Data() : nullptr,
        state.zNormal ? zNormalBuffer->Data() : nullptr,
        state.meshVelocityNormal ? meshVelocityBuffer->Data() : nullptr,
        state.faceArea ? areaBuffer->Data() : nullptr,
        leftCellBuffer->Data(),
        rightCellBuffer->Data(),
        connectivity.boundaryMask ? boundaryMaskBuffer->Data() : nullptr,
        connectivity.nBoundaryFaces,
        nFaces,
        nCells,
        nEq,
        state.gamma,
        scheme,
        fluxBuffer->Data(),
        residualBuffer->Data() );
    hipLaunchKernelGGL(
        GatherFaceFluxEulerKernel,
        dim3( ( nCells + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        fluxBuffer->Data(),
        leftCellBuffer->Data(), rightCellBuffer->Data(),
        connectivity.boundaryMask ? boundaryMaskBuffer->Data() : nullptr,
        cellStorage->cellFaceOffsets.Data(),
        cellStorage->cellFaceIds.Data(),
        connectivity.nBoundaryFaces, nFaces, nCells, nEq,
        residualBuffer->Data() );
    CheckHip(
        hipGetLastError(), "CalcAndScatterPrimitiveEulerKernel launch" );
    CheckHip(
        hipDeviceSynchronize(),
        "CalcAndScatterPrimitiveEulerKernel synchronize" );
    }
    { ScopedStageTimer timer( "D2H" );
    if ( residual.values != nullptr )
    {
        residualBuffer->CopyToHost( residual.values, totalResidual );
    }
    else
    {
        for ( int equation = 0; equation < nEq; ++ equation )
        {
            AccelRuntime::Instance().Backend().Copy(
                residual.components[ equation ],
                residualBuffer->Data() + equation * nCells,
                nCells * sizeof( Real ), AccelCopyKind::DeviceToHost );
        }
    }
    }
    if ( hostFlux != nullptr )
    {
        ScopedStageTimer timer( "D2H" );
        fluxBuffer->CopyToHost( hostFlux->values, totalFlux );
        if ( boundStorage == nullptr )
        {
            deviceFluxFaces_ = nFaces;
            deviceFluxEquations_ = nEq;
        }
    }
}

}
