/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "HipFluxBackend.h"

#include "AccelRuntime.h"
// Only PrimitiveFaceStateView is used here, and that struct is execution-space
// independent. Including the CPU adapter for it made the accelerator backend
// depend on a host header it never calls into.
#include "PrimitiveFaceContract.h"
#include "DeviceBuffer.h"
#include "HipProductionPolicy.h"
#include "Ns3DDomain.h"
#include "StageProfiler.h"
#include "ViscousFluxArithmetic.h"

#include <hip/hip_runtime.h>

#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cstddef>
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

// Every kernel batch in the device-resident path used to be followed by a full
// queue drain, including the batches whose next statement only launches another
// kernel and whose result is never read back on the host in between.  Eight
// ranks share one card, so those drains are paid eight times per stage.  Each
// drain here is an error check only: every real host read goes through
// hipMemcpy, which blocks on its own.  ONEFLOW_HIP_DEFER_SYNC=1 keeps the
// launch-error check (hipGetLastError) and drops the drain.  The default stays
// off so the validated timing basis does not move until the gates rerun with
// the flag enabled.
bool HipDeferSyncEnabled()
{
    // Default on since the 2026-11-07 paired run: same node, 3 interleaved
    // repeats each side, 15.57 -> 14.55 ms/step with the trace and the 200-step
    // cross-backend gate green.  =0 restores the per-stage drain.
    static const bool enabled = [] {
        const char * value = std::getenv( "ONEFLOW_HIP_DEFER_SYNC" );
        return value == nullptr || value[ 0 ] != '0';
    }();
    return enabled;
}

// Deferring every drain lets the software queue grow without bound: with 16 or
// 32 ranks on one card the driver aborted inside libgalaxyhip (AQL packet dump,
// signal 6) because the hardware queue filled up.  Bound the number of kernel
// batches that may be outstanding; the default keeps the per-stage savings but
// never lets the queue overflow.
int MaxLaunchesPerDrain()
{
    static const int limit = [] {
        const char * value =
            std::getenv( "ONEFLOW_HIP_MAX_LAUNCHES_PER_DRAIN" );
        if ( value == nullptr || value[ 0 ] == '\0' ) return 32;
        const int parsed = std::atoi( value );
        return parsed > 0 ? parsed : 32;
    }();
    return limit;
}

void HipDrain( const char * operation )
{
    // Per-thread, not per-process.  With several ranks in one process they share
    // a plain static, so one rank's launches push another's counter past the
    // limit and drains land at the wrong point: sometimes early (paying a sync
    // that was not due), sometimes late (letting the hardware queue grow).  The
    // pinned staging buffer in this same file is already thread_local, so the
    // intent was per-thread here too.
    static thread_local int sinceDrain = 0;
    if ( HipDeferSyncEnabled() && ++ sinceDrain < MaxLaunchesPerDrain() ) return;
    sinceDrain = 0;
    CheckHip( hipDeviceSynchronize(), operation );
}

// Scalar status reads are full round trips: the driver has to drain the queue
// to return four bytes.  The hot path did one per time step and one per stage
// update.  With ONEFLOW_HIP_DEFER_STATUS=1 the flag stays on the device and the
// macro step checks it once through CheckResidentStatus().
bool HipDeferStatusEnabled()
{
    return HipStatusCheckDeferred();
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


// Geometry and connectivity as the resident cell kernels see them.  Both the
// mean-flow gradient and the derived temperature gradient consume the identical
// set, so one cache keyed on these pointers serves both instead of a second
// copy of the same face/cell geometry.
void EnsureCellFaceConnectivity(
    HipGradientStorage & storage,
    const int * leftCell,
    const int * rightCell,
    int nCells,
    int nFaces,
    int nBoundaryFaces );

struct ResidentCellGeometry
{
    const void * ownerToken = nullptr;
    std::uint64_t topologyGeneration = 0;
    int nCells = 0;
    int nGhostCells = 0;
    int nFaces = 0;
    int nBoundaryFaces = 0;
    const Real * xFace = nullptr;
    const Real * yFace = nullptr;
    const Real * zFace = nullptr;
    const Real * xNormal = nullptr;
    const Real * yNormal = nullptr;
    const Real * zNormal = nullptr;
    const Real * faceArea = nullptr;
    const Real * xCell = nullptr;
    const Real * yCell = nullptr;
    const Real * zCell = nullptr;
    const Real * cellVolume = nullptr;
    const int * leftCell = nullptr;
    const int * rightCell = nullptr;
};

// Upload face/cell geometry and the cell-to-face CSR only when the resident
// copy is stale for this binding.  The mean-flow and temperature gradients call
// this with the same view, so the second one is a pointer comparison.
bool EnsureResidentCellGeometry(
    HipGradientStorage & storage,
    const ResidentCellGeometry & view,
    const char * stageName )
{
    const bool cacheMatches = storage.geometryCacheValid
        && storage.cachedOwnerToken == view.ownerToken
        && storage.cachedTopologyGeneration == view.topologyGeneration
        && storage.cachedXFace == view.xFace
        && storage.cachedYFace == view.yFace
        && storage.cachedZFace == view.zFace
        && storage.cachedXNormal == view.xNormal
        && storage.cachedYNormal == view.yNormal
        && storage.cachedZNormal == view.zNormal
        && storage.cachedFaceArea == view.faceArea
        && storage.cachedXCell == view.xCell
        && storage.cachedYCell == view.yCell
        && storage.cachedZCell == view.zCell
        && storage.cachedVolume == view.cellVolume
        && storage.cachedLeftCell == view.leftCell
        && storage.cachedRightCell == view.rightCell
        && storage.cachedCells == view.nCells
        && storage.cachedGhostCells == view.nGhostCells
        && storage.cachedFaces == view.nFaces
        && storage.cachedBoundaryFaces == view.nBoundaryFaces;
    if ( ! cacheMatches )
    {
        ScopedStageTimer geomTimer( stageName );
        const int totalCells = view.nCells + view.nGhostCells;
        storage.xFaceCenter.CopyFromHost( view.xFace, view.nFaces );
        storage.yFaceCenter.CopyFromHost( view.yFace, view.nFaces );
        storage.zFaceCenter.CopyFromHost( view.zFace, view.nFaces );
        storage.xNormal.CopyFromHost( view.xNormal, view.nFaces );
        storage.yNormal.CopyFromHost( view.yNormal, view.nFaces );
        storage.zNormal.CopyFromHost( view.zNormal, view.nFaces );
        storage.faceArea.CopyFromHost( view.faceArea, view.nFaces );
        storage.xCellCenter.CopyFromHost( view.xCell, totalCells );
        storage.yCellCenter.CopyFromHost( view.yCell, totalCells );
        storage.zCellCenter.CopyFromHost( view.zCell, totalCells );
        storage.cellVolume.CopyFromHost( view.cellVolume, view.nCells );
        storage.leftCell.CopyFromHost( view.leftCell, view.nFaces );
        storage.rightCell.CopyFromHost( view.rightCell, view.nFaces );
        storage.cachedOwnerToken = view.ownerToken;
        storage.cachedTopologyGeneration = view.topologyGeneration;
        storage.cachedXFace = view.xFace;
        storage.cachedYFace = view.yFace;
        storage.cachedZFace = view.zFace;
        storage.cachedXNormal = view.xNormal;
        storage.cachedYNormal = view.yNormal;
        storage.cachedZNormal = view.zNormal;
        storage.cachedFaceArea = view.faceArea;
        storage.cachedXCell = view.xCell;
        storage.cachedYCell = view.yCell;
        storage.cachedZCell = view.zCell;
        storage.cachedVolume = view.cellVolume;
        storage.cachedLeftCell = view.leftCell;
        storage.cachedRightCell = view.rightCell;
        storage.cachedCells = view.nCells;
        storage.cachedGhostCells = view.nGhostCells;
        storage.cachedFaces = view.nFaces;
        storage.cachedBoundaryFaces = view.nBoundaryFaces;
        storage.geometryCacheValid = true;
    }
    EnsureCellFaceConnectivity(
        storage, view.leftCell, view.rightCell,
        view.nCells, view.nFaces, view.nBoundaryFaces );
    return ! cacheMatches;
}

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

// Roe flux for the five-equation system, mirroring NsInvFlux::Roe and the
// NsInvFluxPolicy helpers term for term: Roe-averaged state, Harten entropy
// fix on the three wave speeds, then the five-equation dissipation.  Host and
// device must agree to round-off, so the operation order follows the host
// implementation instead of being algebraically "simplified".
struct RoePrimitive
{
    Real density;
    Real u;
    Real v;
    Real w;
    Real pressure;
    Real enthalpy;
};

__device__ void RoeFaceFlux(
    const Primitive & left,
    const Primitive & right,
    Real gamma,
    Real nx, Real ny, Real nz, Real vfn,
    Real entropy1, Real entropy2, bool entropyFix,
    Real * flux )
{
    RoePrimitive pL;
    RoePrimitive pR;
    pL.density = left.density;
    pL.u = left.u;
    pL.v = left.v;
    pL.w = left.w;
    pL.pressure = left.pressure;
    pR.density = right.density;
    pR.u = right.u;
    pR.v = right.v;
    pR.w = right.w;
    pR.pressure = right.pressure;

    const Real gamm1 = gamma - 1.0;
    pL.enthalpy = gamma / gamm1 * pL.pressure / pL.density
        + 0.5 * ( pL.u * pL.u + pL.v * pL.v + pL.w * pL.w );
    pR.enthalpy = gamma / gamm1 * pR.pressure / pR.density
        + 0.5 * ( pR.u * pR.u + pR.v * pR.v + pR.w * pR.w );

    // RoeAverageWeights / RoeAverageComponent
    const Real ratio = sqrt( pR.density / pL.density );
    const Real coef = 1.0 / ( 1.0 + ratio );
    const Real densityMid = sqrt( pL.density * pR.density );
    const Real uMid = ( pL.u + pR.u * ratio ) * coef;
    const Real vMid = ( pL.v + pR.v * ratio ) * coef;
    const Real wMid = ( pL.w + pR.w * ratio ) * coef;
    const Real pMid = ( pL.pressure + pR.pressure * ratio ) * coef;
    const Real hMid = ( pL.enthalpy + pR.enthalpy * ratio ) * coef;

    const Real speedSquared = uMid * uMid + vMid * vMid + wMid * wMid;
    const Real vnFlow = nx * uMid + ny * vMid + nz * wMid;
    const Real vnRel = vnFlow - vfn;
    const Real c2 = gamm1 * ( hMid - 0.5 * speedSquared );
    const Real c = sqrt( fabs( c2 ) );

    // RoeAbsoluteEigenvalues + ApplyHartenEntropyFix
    Real absVn = fabs( vnRel );
    Real absVnPlusC = fabs( vnRel + c );
    Real absVnMinusC = fabs( vnRel - c );
    if ( entropyFix )
    {
        const Real maxEigenvalue =
            absVnPlusC > absVnMinusC ? absVnPlusC : absVnMinusC;
        const Real m1 = maxEigenvalue * entropy1;
        const Real m2 = maxEigenvalue * entropy2;
        absVn = m1 > absVn ? m1 : absVn;
        absVnPlusC = m2 > absVnPlusC ? m2 : absVnPlusC;
        absVnMinusC = m2 > absVnMinusC ? m2 : absVnMinusC;
    }

    // RoeDissipationWaveCoefficients
    const Real xi1 = ( 2.0 * absVn - absVnPlusC - absVnMinusC )
        / ( 2.0 * c2 );
    const Real xi2 = ( absVnPlusC - absVnMinusC ) / ( 2.0 * c );

    // Conserved jumps, matching PrimToQ on both sides.
    Real q1[ 5 ] = {};
    Real q2[ 5 ] = {};
    const Real pLc[ 5 ] = { pL.density, pL.u, pL.v, pL.w, pL.pressure };
    const Real pRc[ 5 ] = { pR.density, pR.u, pR.v, pR.w, pR.pressure };
    const Real primLeft[ 5 ] = { pLc[ 0 ], pLc[ 1 ], pLc[ 2 ], pLc[ 3 ], pLc[ 4 ] };
    (void) primLeft;
    const Real eL = pL.pressure / gamm1
        + 0.5 * pL.density * ( pL.u * pL.u + pL.v * pL.v + pL.w * pL.w );
    const Real eR = pR.pressure / gamm1
        + 0.5 * pR.density * ( pR.u * pR.u + pR.v * pR.v + pR.w * pR.w );
    q1[ 0 ] = pL.density;
    q1[ 1 ] = pL.density * pL.u;
    q1[ 2 ] = pL.density * pL.v;
    q1[ 3 ] = pL.density * pL.w;
    q1[ 4 ] = eL;
    q2[ 0 ] = pR.density;
    q2[ 1 ] = pR.density * pR.u;
    q2[ 2 ] = pR.density * pR.v;
    q2[ 3 ] = pR.density * pR.w;
    q2[ 4 ] = eR;
    Real dq[ 5 ] = {};
    for ( int eq = 0; eq < 5; ++ eq ) dq[ eq ] = q2[ eq ] - q1[ eq ];

    // Fluxes
    Real f[ 5 ] = {};
    const Real massL = pL.density * ( nx * pL.u + ny * pL.v + nz * pL.w - vfn );
    const Real massR = pR.density * ( nx * pR.u + ny * pR.v + nz * pR.w - vfn );
    f[ 0 ] = massL + massR;
    f[ 1 ] = massL * pL.u + nx * pL.pressure + massR * pR.u + nx * pR.pressure;
    f[ 2 ] = massL * pL.v + ny * pL.pressure + massR * pR.v + ny * pR.pressure;
    f[ 3 ] = massL * pL.w + nz * pL.pressure + massR * pR.w + nz * pR.pressure;
    f[ 4 ] = massL * pL.enthalpy + vfn * pL.pressure
        + massR * pR.enthalpy + vfn * pR.pressure;

    // RoeContinuityNormalJump, RoeTotalEnthalpyJump, coupling terms
    const Real dc = vnFlow * dq[ 0 ] - nx * dq[ 1 ] - ny * dq[ 2 ] - nz * dq[ 3 ];
    const Real ae = gamma - 1.0;
    Real dh = -ae * ( uMid * dq[ 1 ] + vMid * dq[ 2 ] + wMid * dq[ 3 ] - dq[ 4 ] );
    dh += 0.5 * ae * speedSquared * dq[ 0 ];
    const Real term1 = dh * xi1 + dc * xi2;
    const Real term2 = c2 * dc * xi1 + dh * xi2;

    f[ 0 ] -= ( absVn * dq[ 0 ] - term1 );
    f[ 1 ] -= ( absVn * dq[ 1 ] - uMid * term1 + nx * term2 );
    f[ 2 ] -= ( absVn * dq[ 2 ] - vMid * term1 + ny * term2 );
    f[ 3 ] -= ( absVn * dq[ 3 ] - wMid * term1 + nz * term2 );
    f[ 4 ] -= ( absVn * dq[ 4 ] - hMid * term1 + vnFlow * term2 );

    for ( int eq = 0; eq < 5; ++ eq ) flux[ eq ] = 0.5 * f[ eq ];
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
    Real entropy1,
    Real entropy2,
    int entropyFix,
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
    if ( scheme == 2 )
    {
        Real roeFlux[ 5 ] = {};
        RoeFaceFlux( pL, pR, gamma, nx, ny, nz, vfn,
                     entropy1, entropy2, entropyFix != 0, roeFlux );
        for ( int eq = 0; eq < nEq; ++ eq )
        {
            flux[ eq * nFaces + face ] = faceArea * roeFlux[ eq ];
        }
        return;
    }
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

// Reset a device buffer without a host round trip.  The explicit device chain
// owns its residual accumulation base; uploading a host field that the stage
// did not produce would either cost a full H2D or, if skipped, accumulate on
// top of the previous stage.
__global__ void ZeroValuesKernel( Real * values, int nValues )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index < nValues ) values[ index ] = 0.0;
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

// Temperature from the resident primitive state.
//
// This reproduces NsCalcGamaT's relation term for term:
//
//     T[cell] = p[cell] / ( stateCoefficient * rho[cell] * oamw )
//
// where the host reads rho from q[IDX::IR], p from q[IDX::IP], and oamw is the
// molecular-weight ratio (the local `one` outside the chemistry path).  The
// operation order matches the host expression, so an owned cell agrees to
// round-off rather than merely closely.
//
// Gamma does not appear as a literal.  The host gets it through
// nscom.statecoef = 1 / ( gama_ref * mach_ref^2 ), i.e. gamma enters as a
// configurable field, so the kernel takes the solved coefficient instead of
// substituting 1.4.
__global__ void TemperatureFromPrimitiveKernel(
    const Real * density,
    const Real * pressure,
    int totalCells,
    Real stateCoefficient,
    Real molecularWeightRatio,
    Real * temperature )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= totalCells ) return;
    temperature[ cell ] =
        pressure[ cell ]
        / ( stateCoefficient * density[ cell ] * molecularWeightRatio );
}

// Cell gradient of the temperature field.  Deliberately the same discrete
// operator as GradientCellKernel, because the viscous chain differentiates
// temperature exactly the way it differentiates the mean-flow variables and a
// different stencil would show up as a heat-flux mismatch.  The inner face
// interpolation weights and the normal-area projection are copied unchanged;
// only the scalar field being accumulated differs.
__global__ void TemperatureGradientCellKernel(
    const Real * temperature, Real * gradient,
    const Real * xFace, const Real * yFace, const Real * zFace,
    const Real * xNormal, const Real * yNormal, const Real * zNormal,
    const Real * faceArea, const Real * xCell, const Real * yCell,
    const Real * zCell, const Real * volume,
    const int * leftCell, const int * rightCell,
    const int * cellFaceOffsets, const int * cellFaceIds,
    int nCells, int nBoundaryFaces, int nFaces )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;
    const int totalCells = nCells + nBoundaryFaces;
    Real accumulated[ 3 ] = {};
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
        const Real value = cl * temperature[ left ]
            + cr * temperature[ right ];
        accumulated[ 0 ] += sign * ax * value;
        accumulated[ 1 ] += sign * ay * value;
        accumulated[ 2 ] += sign * az * value;
    }
    const Real inverseVolume = 1.0 / volume[ cell ];
    gradient[ cell ] = accumulated[ 0 ] * inverseVolume;
    gradient[ totalCells + cell ] = accumulated[ 1 ] * inverseVolume;
    gradient[ 2 * totalCells + cell ] = accumulated[ 2 ] * inverseVolume;
}

// Boundary ghost cells have no control volume, so they take the owner cell's
// gradient exactly as GradientCopyBoundaryKernel does for the mean flow.
__global__ void TemperatureGradientCopyBoundaryKernel(
    Real * gradient, const int * leftCell, const int * rightCell,
    int nCells, int nBoundaryFaces )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nBoundaryFaces ) return;
    const int totalCells = nCells + nBoundaryFaces;
    const int left = leftCell[ face ];
    const int ghost = rightCell[ face ];
    gradient[ ghost ] = gradient[ left ];
    gradient[ totalCells + ghost ] = gradient[ totalCells + left ];
    gradient[ 2 * totalCells + ghost ] = gradient[ 2 * totalCells + left ];
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


__device__ Real DeviceBarthFactor(
    Real dqdx, Real dqdy, Real dqdz,
    Real dx, Real dy, Real dz,
    Real nx, Real ny, Real nz,
    Real minv, Real maxv,
    Real orientationSign )
{
    const Real dq = dqdx * dx + dqdy * dy + dqdz * dz;
    const Real ds = sqrt( dx * dx + dy * dy + dz * dz );
    const Real dot = orientationSign * ( nx * dx + ny * dy + nz * dz )
        / ( ds + 1.0e-30 );
    if ( dot < 0.0 ) return 0.0;
    if ( dq > maxv ) return maxv / dq;
    if ( dq < minv ) return minv / dq;
    return 1.0;
}

// Barth limiter factor on the device: neighbour range and the minimum over the
// cell's faces, per equation.  Mirrors the host sequence
//   BuildCellNeighborOffsetBounds -> per-face BarthFunction -> TightenLimiter
// with one difference in form only: the host tightens while walking faces, this
// walks one cell's faces and takes the minimum, which is the same value because
// min() is commutative.
__global__ void BarthLimiterKernel(
    const Real * q,
    const Real * gradient,
    const Real * xFace, const Real * yFace, const Real * zFace,
    const Real * xCell, const Real * yCell, const Real * zCell,
    const Real * xNormal, const Real * yNormal, const Real * zNormal,
    const int * leftCell, const int * rightCell,
    const int * faceOffsets, const int * faceIds,
    const unsigned char * boundaryMask,
    const ReconstructionBoundaryOperation * boundaryOperation,
    int nCells, int totalCells, int nFaces, int nEquations,
    Real * phiField )
{
    // One thread per cell with all equations in registers: the face adjacency is
    // walked once, not once per equation. The per-(cell,equation) form traversed
    // the connectivity list nEquations times and made the stage ~45x slower than
    // the unlimited path.
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;

    Real lo[ 5 ];
    Real hi[ 5 ];
    for ( int eq = 0; eq < nEquations; ++ eq )
    {
        lo[ eq ] = q[ eq * totalCells + cell ];
        hi[ eq ] = lo[ eq ];
    }

    const int begin = faceOffsets[ cell ];
    const int end = faceOffsets[ cell + 1 ];

    for ( int k = begin; k < end; ++ k )
    {
        const int face = faceIds[ k ];
        const bool isLeft = leftCell[ face ] == cell;
        if ( isLeft && boundaryMask != nullptr && boundaryMask[ face ] != 0 )
        {
            const ReconstructionBoundaryOperation op =
                boundaryOperation != nullptr
                ? boundaryOperation[ face ]
                : ReconstructionBoundaryOperation::Preserve;
            if ( op != ReconstructionBoundaryOperation::Preserve ) continue;
        }
        const int neighbour = isLeft ? rightCell[ face ] : leftCell[ face ];
        for ( int eq = 0; eq < nEquations; ++ eq )
        {
            const Real value = q[ eq * totalCells + neighbour ];
            lo[ eq ] = value < lo[ eq ] ? value : lo[ eq ];
            hi[ eq ] = value > hi[ eq ] ? value : hi[ eq ];
        }
    }

    Real phi[ 5 ];
    for ( int eq = 0; eq < nEquations; ++ eq ) phi[ eq ] = 1.0;

    const Real cellX = xCell[ cell ];
    const Real cellY = yCell[ cell ];
    const Real cellZ = zCell[ cell ];
    for ( int k = begin; k < end; ++ k )
    {
        const int face = faceIds[ k ];
        const bool isLeft = leftCell[ face ] == cell;
        const Real dx = xFace[ face ] - cellX;
        const Real dy = yFace[ face ] - cellY;
        const Real dz = zFace[ face ] - cellZ;
        const Real nx = xNormal ? xNormal[ face ] : 1.0;
        const Real ny = yNormal ? yNormal[ face ] : 0.0;
        const Real nz = zNormal ? zNormal[ face ] : 0.0;
        const Real sign = isLeft ? 1.0 : -1.0;
        for ( int eq = 0; eq < nEquations; ++ eq )
        {
            const Real minv = lo[ eq ] - q[ eq * totalCells + cell ];
            const Real maxv = hi[ eq ] - q[ eq * totalCells + cell ];
            const Real factor = DeviceBarthFactor(
                gradient[ eq * 3 * totalCells + cell ],
                gradient[ eq * 3 * totalCells + totalCells + cell ],
                gradient[ eq * 3 * totalCells + 2 * totalCells + cell ],
                dx, dy, dz, nx, ny, nz, minv, maxv, sign );
            phi[ eq ] = factor < phi[ eq ] ? factor : phi[ eq ];
        }
    }

    for ( int eq = 0; eq < nEquations; ++ eq )
        phiField[ eq * totalCells + cell ] = phi[ eq ];
    // Ghost factor: the host initialises the whole field to one and only
    // physical cells are tightened.
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
    const Real * limiter,
    int limiterMode,
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
        // The limiter scales the reconstruction increment, not the cell value.
        // limiterMode 0 keeps the unlimited second-order projection that the
        // device path has always used.
        // LimField::CalcFaceValue blends both sides with min(phi_l, phi_r) and
        // applies that single factor to both reconstructions; matching it here
        // is what keeps the device face values equal to the CPU oracle.
        const Real phiLeft = limiterMode != 0 ? limiter[ qBase + left ] : 1.0;
        const Real phiRight = limiterMode != 0 ? limiter[ qBase + right ] : 1.0;
        const Real phiFace = phiLeft < phiRight ? phiLeft : phiRight;
        leftCandidate[ eq ] = qLeft + phiFace * (
            gradient[ gradientBase + left ] * dxl
            + gradient[ gradientBase + totalCells + left ] * dyl
            + gradient[ gradientBase + 2 * totalCells + left ] * dzl );
        rightCandidate[ eq ] = qRight + phiFace * (
            gradient[ gradientBase + right ] * dxr
            + gradient[ gradientBase + totalCells + right ] * dyr
            + gradient[ gradientBase + 2 * totalCells + right ] * dzr );
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

// Copy a cell-scalar component into a staging buffer.
//
// Defined once instead of repeating the AccelRuntime::Backend().Copy spell-out:
// the resident-state copies in this file are already verbose enough that a
// mixed-up source/destination is hard to see on review.
void CopyRealComponent(
    DeviceBuffer< Real > & destination, std::size_t offset, const Real * source,
    int count, AccelCopyKind kind )
{
    AccelRuntime::Instance().Backend().Copy(
        destination.Data() + offset, source,
        static_cast< std::size_t >( count ) * sizeof( Real ), kind );
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

namespace
{

__global__ void EvaluateDeviceRoeKernel(
    const Real * left,
    const Real * right,
    const Real * normal,
    Real meshVelocityNormal,
    Real gamma,
    Real entropy1,
    Real entropy2,
    int entropyFix,
    Real * flux )
{
    Primitive pL;
    Primitive pR;
    pL.density = left[ 0 ];
    pL.u = left[ 1 ];
    pL.v = left[ 2 ];
    pL.w = left[ 3 ];
    pL.pressure = left[ 4 ];
    pR.density = right[ 0 ];
    pR.u = right[ 1 ];
    pR.v = right[ 2 ];
    pR.w = right[ 3 ];
    pR.pressure = right[ 4 ];
    Real out[ 5 ] = {};
    RoeFaceFlux( pL, pR, gamma, normal[ 0 ], normal[ 1 ], normal[ 2 ],
                 meshVelocityNormal, entropy1, entropy2,
                 entropyFix != 0, out );
    for ( int eq = 0; eq < 5; ++ eq ) flux[ eq ] = out[ eq ];
}

}

namespace
{

// Mirrors BarthFunction + ReconstructFaceIncrement + FaceSideDirectionDot in
// codes/usolver/src/ULimiterPolicy.cpp.  Kept as one function so the device and
// the host cannot drift on the normal-dot sign convention.

__global__ void EvaluateDeviceBarthKernel(
    Real dqdx, Real dqdy, Real dqdz,
    Real dx, Real dy, Real dz,
    Real nx, Real ny, Real nz,
    Real minv, Real maxv,
    Real orientationSign,
    Real * result )
{
    *result = DeviceBarthFactor( dqdx, dqdy, dqdz, dx, dy, dz,
                                 nx, ny, nz, minv, maxv, orientationSign );
}

}

void HipFluxBackend::CopyLimiterStateToHost(
    const void * cacheKey, int nCells, int nFaces, int nEquations,
    Real * phi, Real * faceLeft, Real * faceRight )
{
    const auto binding = boundStates_.find( cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP limiter dump requires a bound state" );
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
    auto copyComponent = []( DeviceBuffer< Real > & buffer,
                             Real * destination, std::size_t componentStride,
                             std::size_t component, std::size_t count )
    {
        AccelRuntime::Instance().Backend().Copy(
            destination,
            buffer.Data() + component * componentStride,
            count * sizeof( Real ),
            AccelCopyKind::DeviceToHost );
    };
    if ( phi != nullptr
         && storage.limiter.Size()
                >= static_cast< std::size_t >( nEquations ) * totalCells )
    {
        for ( int equation = 0; equation < nEquations; ++ equation )
            copyComponent( storage.limiter,
                           phi + static_cast< std::size_t >( equation ) * nCells,
                           totalCells, static_cast< std::size_t >( equation ),
                           static_cast< std::size_t >( nCells ) );
    }
    if ( faceLeft != nullptr
         && storage.qf1.Size()
                >= static_cast< std::size_t >( nEquations ) * nFaces )
    {
        for ( int equation = 0; equation < nEquations; ++ equation )
            copyComponent( storage.qf1,
                           faceLeft
                               + static_cast< std::size_t >( equation ) * nFaces,
                           nFaces, static_cast< std::size_t >( equation ),
                           static_cast< std::size_t >( nFaces ) );
    }
    if ( faceRight != nullptr
         && storage.qf2.Size()
                >= static_cast< std::size_t >( nEquations ) * nFaces )
    {
        for ( int equation = 0; equation < nEquations; ++ equation )
            copyComponent( storage.qf2,
                           faceRight
                               + static_cast< std::size_t >( equation ) * nFaces,
                           nFaces, static_cast< std::size_t >( equation ),
                           static_cast< std::size_t >( nFaces ) );
    }
}

Real HipFluxBackend::EvaluateDeviceBarth(
    Real dqdx, Real dqdy, Real dqdz,
    Real dx, Real dy, Real dz,
    Real nx, Real ny, Real nz,
    Real minv, Real maxv,
    Real orientationSign ) const
{
    Real * deviceResult = nullptr;
    CheckHip( hipMalloc( &deviceResult, sizeof( Real ) ), "barth contract alloc" );
    hipLaunchKernelGGL(
        EvaluateDeviceBarthKernel, dim3( 1 ), dim3( 1 ), 0, 0,
        dqdx, dqdy, dqdz, dx, dy, dz, nx, ny, nz,
        minv, maxv, orientationSign, deviceResult );
    CheckHip( hipGetLastError(), "barth contract kernel launch" );
    CheckHip( hipDeviceSynchronize(), "barth contract kernel synchronize" );
    Real result = 0.0;
    CheckHip( hipMemcpy( &result, deviceResult, sizeof( Real ),
                         hipMemcpyDeviceToHost ), "barth contract copy" );
    CheckHip( hipFree( deviceResult ), "barth contract free" );
    return result;
}

void HipFluxBackend::EvaluateDeviceRoe(
    const DeviceRoeSample & sample, Real flux[ 5 ] ) const
{
    Real * deviceFlux = nullptr;
    CheckHip( hipMalloc( &deviceFlux, 5 * sizeof( Real ) ), "roe contract alloc" );
    Real * deviceLeft = nullptr;
    Real * deviceRight = nullptr;
    Real * deviceNormal = nullptr;
    CheckHip( hipMalloc( &deviceLeft, 5 * sizeof( Real ) ), "roe contract alloc" );
    CheckHip( hipMalloc( &deviceRight, 5 * sizeof( Real ) ), "roe contract alloc" );
    CheckHip( hipMalloc( &deviceNormal, 3 * sizeof( Real ) ), "roe contract alloc" );
    CheckHip( hipMemcpy( deviceLeft, sample.left, 5 * sizeof( Real ),
                         hipMemcpyHostToDevice ), "roe contract copy" );
    CheckHip( hipMemcpy( deviceRight, sample.right, 5 * sizeof( Real ),
                         hipMemcpyHostToDevice ), "roe contract copy" );
    CheckHip( hipMemcpy( deviceNormal, sample.normal, 3 * sizeof( Real ),
                         hipMemcpyHostToDevice ), "roe contract copy" );
    hipLaunchKernelGGL(
        EvaluateDeviceRoeKernel, dim3( 1 ), dim3( 1 ), 0, 0,
        deviceLeft, deviceRight, deviceNormal,
        sample.meshVelocityNormal, sample.gamma,
        sample.entropy1, sample.entropy2, sample.entropyFix ? 1 : 0,
        deviceFlux );
    CheckHip( hipGetLastError(), "roe contract kernel launch" );
    CheckHip( hipDeviceSynchronize(), "roe contract kernel synchronize" );
    CheckHip( hipMemcpy( flux, deviceFlux, 5 * sizeof( Real ),
                         hipMemcpyDeviceToHost ), "roe contract copy" );
    CheckHip( hipFree( deviceLeft ), "roe contract free" );
    CheckHip( hipFree( deviceRight ), "roe contract free" );
    CheckHip( hipFree( deviceNormal ), "roe contract free" );
    CheckHip( hipFree( deviceFlux ), "roe contract free" );
}

void HipFluxBackend::CalcGradient( const CellGradientView & view )
{
    // Per-call host preparation is not free at this scale: the breakdown of a
    // 1000-step run shows the parent category costing more than the kernel it
    // wraps.  Split the entry into prep / contract / geometry so the next round
    // can name the piece instead of guessing.
    HipGradientStorage * storage = &fallbackGradientStorage_;
    Ns3DDeviceState * owner = nullptr;
    const void * ownerToken = view.cacheKey;
    std::uint64_t fieldGeneration = 0;
    std::uint64_t topologyGeneration = 0;
    int totalCells = 0;
    {
        ScopedStageTimer prepTimer( "g_prep" );
        ValidateCellGradientView( view );
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
        totalCells = view.nCells + view.nGhostCells;
        storage->cellState.Resize( view.nEquations * totalCells );
        ownerToken = owner != nullptr ? owner->OwnerToken() : view.cacheKey;
        fieldGeneration = owner != nullptr ? owner->FieldGeneration() : 0;
        topologyGeneration = owner != nullptr ? owner->TopologyGeneration() : 0;
    }
    const int totalGradient = view.nEquations * 3 * totalCells;
    ResidentStateSync sync = ResidentStateSync::None;
    {
        ScopedStageTimer contractTimer( "g_contract" );
        sync = RefreshResidentPrimitiveState(
            * storage, owner, ownerToken, topologyGeneration, fieldGeneration,
            view.nCells, view.nGhostCells, view.nEquations, view.q );
    }
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
    ResidentCellGeometry geometry;
    geometry.ownerToken = ownerToken;
    geometry.topologyGeneration = topologyGeneration;
    geometry.nCells = view.nCells;
    geometry.nGhostCells = view.nGhostCells;
    geometry.nFaces = view.nFaces;
    geometry.nBoundaryFaces = view.nBoundaryFaces;
    geometry.xFace = view.xFace;
    geometry.yFace = view.yFace;
    geometry.zFace = view.zFace;
    geometry.xNormal = view.xNormal;
    geometry.yNormal = view.yNormal;
    geometry.zNormal = view.zNormal;
    geometry.faceArea = view.faceArea;
    geometry.xCell = view.xCell;
    geometry.yCell = view.yCell;
    geometry.zCell = view.zCell;
    geometry.cellVolume = view.cellVolume;
    geometry.leftCell = view.leftCell;
    geometry.rightCell = view.rightCell;
    EnsureResidentCellGeometry( *storage, geometry, "g_geom" );
    storage->cellGradient.Resize( totalGradient );
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
    HipDrain( "gradient kernel synchronize" );
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

void HipFluxBackend::CalcCellTemperature( const CellTemperatureView & view )
{
    ValidateCellTemperatureView( view );

    // Resolve the resident storage the same way CalcGradient does: the bound
    // state owns it so temperature survives across stages, and the fallback
    // keeps the unbound diagnostic path working.
    HipGradientStorage * storage = &fallbackGradientStorage_;
    Ns3DDeviceState * owner = nullptr;
    const void * ownerToken = view.cacheKey;
    std::uint64_t fieldGeneration = 0;
    std::uint64_t topologyGeneration = 0;
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
        ownerToken = owner->OwnerToken();
        fieldGeneration = owner->FieldGeneration();
        topologyGeneration = owner->TopologyGeneration();
    }

    const int totalCells = view.nCells + view.nGhostCells;
    storage->cellState.Resize(
        static_cast< std::size_t >( view.nEquations ) * totalCells );

    // Temperature must be built from the newest primitive state, so consume the
    // same resident-state contract the gradient uses.  This also makes a stale
    // host mirror fail loudly rather than silently producing a temperature
    // field from an old q.
    RefreshResidentPrimitiveState(
        * storage, owner, ownerToken, topologyGeneration, fieldGeneration,
        view.nCells, view.nGhostCells, view.nEquations, view.q );

    // The density and pressure planes are already resident, so the kernel reads
    // them from the device copy instead of uploading them again.
    const Real * density = storage->cellState.Data();
    const Real * pressure =
        storage->cellState.Data() + 4 * static_cast< std::size_t >( totalCells );

    storage->temperature.Resize( totalCells );
    storage->temperatureGradient.Resize(
        static_cast< std::size_t >( 3 ) * totalCells );

    ResidentCellGeometry geometry;
    geometry.ownerToken = ownerToken;
    geometry.topologyGeneration = topologyGeneration;
    geometry.nCells = view.nCells;
    geometry.nGhostCells = view.nGhostCells;
    geometry.nFaces = view.nFaces;
    geometry.nBoundaryFaces = view.nBoundaryFaces;
    geometry.xFace = view.xFace;
    geometry.yFace = view.yFace;
    geometry.zFace = view.zFace;
    geometry.xNormal = view.xNormal;
    geometry.yNormal = view.yNormal;
    geometry.zNormal = view.zNormal;
    geometry.faceArea = view.faceArea;
    geometry.xCell = view.xCell;
    geometry.yCell = view.yCell;
    geometry.zCell = view.zCell;
    geometry.cellVolume = view.cellVolume;
    geometry.leftCell = view.leftCell;
    geometry.rightCell = view.rightCell;
    // Reuses the geometry and cell-face CSR the mean-flow gradient already
    // uploaded for the same binding.
    EnsureResidentCellGeometry( * storage, geometry, "t_geom" );

    const int blockSize = 256;
    {
        ScopedStageTimer timer( "temperature_kernel" );
        hipLaunchKernelGGL(
            TemperatureFromPrimitiveKernel,
            dim3( ( totalCells + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0,
            density, pressure, totalCells,
            view.stateCoefficient, view.molecularWeightRatio,
            storage->temperature.Data() );
        hipLaunchKernelGGL(
            TemperatureGradientCellKernel,
            dim3( ( view.nCells + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0,
            storage->temperature.Data(), storage->temperatureGradient.Data(),
            storage->xFaceCenter.Data(), storage->yFaceCenter.Data(),
            storage->zFaceCenter.Data(), storage->xNormal.Data(),
            storage->yNormal.Data(), storage->zNormal.Data(),
            storage->faceArea.Data(), storage->xCellCenter.Data(),
            storage->yCellCenter.Data(), storage->zCellCenter.Data(),
            storage->cellVolume.Data(),
            storage->leftCell.Data(), storage->rightCell.Data(),
            storage->cellFaceOffsets.Data(), storage->cellFaceIds.Data(),
            view.nCells, view.nBoundaryFaces, view.nFaces );
        if ( view.nBoundaryFaces > 0 )
        {
            hipLaunchKernelGGL(
                TemperatureGradientCopyBoundaryKernel,
                dim3( ( view.nBoundaryFaces + blockSize - 1 ) / blockSize ),
                dim3( blockSize ), 0, 0,
                storage->temperatureGradient.Data(),
                storage->leftCell.Data(), storage->rightCell.Data(),
                view.nCells, view.nBoundaryFaces );
        }
        CheckHip( hipGetLastError(), "temperature kernel launch" );
        HipDrain( "temperature kernel synchronize" );
    }

    // Mirror the derived fields into the state-owned host storage so a host
    // reader (and the contract test) can compare them without reaching into
    // device buffers.  The device gradient is [component][cell] while the state
    // stores one plane per direction, so one contiguous copy is split on the
    // host.  A caller that also asked for the diagnostic mirror uses the same
    // staging buffer rather than repeating the transfer.
    //
    // Nothing is transferred when neither consumer asked for it.  The previous
    // version pulled the three gradient planes unconditionally, which is
    // 3 * totalCells * 8 bytes - about 7 MB for a 300k-cell zone - on every
    // call even when the caller wanted the field only on the device.  No host
    // code reads owner->temperature* today, so that transfer was pure waste.
    const bool mirrorWanted = ( owner != nullptr ) || view.downloadToHost;
    const std::size_t totalValues = static_cast< std::size_t >( 3 ) * totalCells;
    const Real * staging = nullptr;
    if ( mirrorWanted )
    {
        staging = LocalPinnedStaging().Get( totalValues );
        AccelRuntime::Instance().Backend().Copy(
            const_cast< Real * >( staging ),
            storage->temperatureGradient.Data(),
            totalValues * sizeof( Real ), AccelCopyKind::DeviceToHost );
    }
    if ( owner != nullptr )
    {
        ScopedStageTimer timer( "t_download" );
        // The state sizes these in Upload(); a reservation-only state or a
        // grid whose ghost extent grew would otherwise be too short to receive
        // the mirror, so make the destination match the device extent first.
        const std::size_t mirrorSize = static_cast< std::size_t >( totalCells );
        if ( owner->temperature.size() != mirrorSize )
        {
            owner->temperature.assign( mirrorSize, 0.0 );
            owner->temperatureGradientX.assign( mirrorSize, 0.0 );
            owner->temperatureGradientY.assign( mirrorSize, 0.0 );
            owner->temperatureGradientZ.assign( mirrorSize, 0.0 );
        }
        storage->temperature.CopyToHost(
            owner->temperature.data(), totalCells );
        std::copy(
            staging, staging + totalCells,
            owner->temperatureGradientX.data() );
        std::copy(
            staging + totalCells, staging + 2 * totalCells,
            owner->temperatureGradientY.data() );
        std::copy(
            staging + 2 * totalCells, staging + 3 * totalCells,
            owner->temperatureGradientZ.data() );
    }

    if ( ! view.downloadToHost ) return;
    storage->temperature.CopyToHost( view.temperature, totalCells );
    std::copy( staging, staging + totalCells, view.dTdx );
    std::copy( staging + totalCells, staging + 2 * totalCells, view.dTdy );
    std::copy( staging + 2 * totalCells, staging + 3 * totalCells, view.dTdz );
}

// ---------------------------------------------------------------------------
// Viscous mean-flow flux (laminar + turbulent stress and heat conduction)
//
// The host chain this mirrors, in execution order:
//
//   UNsVisFlux::CalcFlux            per-face loop, per RK stage
//   UNsVisFlux::PrepareFaceValue    geometry, side gradients, averaged viscosity
//   UNsVisFlux::StdMethod           VIS_STD (ivischeme == 0) face-gradient build
//   UNsVisFlux::CalcFaceVisFlux     heat flux -> stress -> NS viscous flux
//   UNsVisFlux::UpdateFaceVisFlux   coeff = -invReynolds * faceArea
//
// Every arithmetic step keeps the host operation order, because the contract is
// round-off agreement with the host face flux, not merely "the same formula".
//
// Packed index convention for cell fields, totalCells = nCells + nGhostCells:
//   cellState[ eq * totalCells + cell ]                     rho,u,v,w,p
//   cellGradient[ (eq*3 + component) * totalCells + cell ]  component 0=x,1=y,2=z
//   temperatureGradient[ component * totalCells + cell ]
//   temperature[ cell ]
// ---------------------------------------------------------------------------

namespace
{

// Sutherland viscosity: the nondimensional sqrt( t^3 ) form.
//
// The expression lives in ViscousFluxArithmetic.h so this kernel and the host
// reference compile the same function rather than two transcriptions.  The
// header is plain inline arithmetic, which the device toolchain accepts from a
// __device__ body; the host policy headers are not.
__device__ Real DeviceSutherlandViscosity(
    Real temperatureNondim, Real sutherlandTsNondim )
{
    return ONEFLOW::SutherlandViscosityNondimDevice(
        temperatureNondim, sutherlandTsNondim );
}

// One face of UNsVisFlux::CalcFaceVisFlux under the VIS_STD face preparation.
//
// Host correspondence, term by term:
//   PrepareFaceValue   gcom normals/area, side gradients, visl1/2, vist1/2,
//                      the two averages, the face state q1/q2
//   AverGrad           VisGrad::AverGrad, arithmetic mean of the two sides
//   CalcGradCoef       VisGradGeom::CalcGradCoef (unit cell-centre displacement)
//   AverFaceValue      VisGrad::AverFaceValue (arithmetic mean, both fields)
//   CorrectFaceGrad    CorrectFaceGradient (UGradPolicy.cpp:86), both fields
//   CalcNormalGrad     VisGrad::CalcNormalGrad (both fields)
//   StokesStressTensor / ProjectStressMomentumFlux / ViscousEnergyFlux
//   EffectiveThermalConductivity / ConductiveHeatFluxNormal
//   UpdateFaceVisFlux  coeff = -invReynolds * faceArea
//
// Two orders that matter and are easy to get wrong:
//   * AverGrad runs before CorrectFaceGrad, and CalcNormalGrad runs after it,
//     so the normal gradient is taken from the corrected gradient.  Swapping
//     either pair changes the heat flux.
//   * nscom.vis is visl + vist of the *averaged* coefficients, so each field is
//     averaged between the two cells and the results are added.
//
// The temperature field is needed even though no temperature average enters the
// flux: CorrectFaceGrad corrects the temperature gradient using the two cells'
// temperature values, exactly as VisGrad::CorrectFaceGrad does for visT.
//
// qNormal is returned unscaled; the caller applies SaveHeatFlux's
// -oreynolds * qNormal for solid faces.
__device__ void ViscousFaceFlux(
    const Real * cellState,
    const Real * cellGradient,
    const Real * temperature,
    const Real * temperatureGradient,
    const Real * laminarViscosity,
    const Real * turbulentViscosity,
    int totalCells,
    Real xfn, Real yfn, Real zfn, Real faceArea,
    Real xCellL, Real yCellL, Real zCellL,
    Real xCellR, Real yCellR, Real zCellR,
    int left, int right,
    Real inversePrandtlLaminar, Real inversePrandtlTurbulent,
    Real specificHeatCp, Real inverseReynolds,
    Real flux[ 5 ], Real & qNormal )
{
    // --- AverGrad -----------------------------------------------------------
    // The stress tensor reads du/dx..dw/dz, i.e. the three momentum equations
    // (1, 2, 3).  The density and pressure gradients are not read by any term
    // of the host chain either.
    Real dqdx[ 3 ] = {};
    Real dqdy[ 3 ] = {};
    Real dqdz[ 3 ] = {};
    for ( int m = 0; m < 3; ++ m )
    {
        const int eq = 1 + m;
        const std::size_t base =
            static_cast< std::size_t >( eq ) * 3 * totalCells;
        dqdx[ m ] = 0.5 * ( cellGradient[ base + left ]
                          + cellGradient[ base + right ] );
        dqdy[ m ] = 0.5 * ( cellGradient[ base + totalCells + left ]
                          + cellGradient[ base + totalCells + right ] );
        dqdz[ m ] = 0.5 * ( cellGradient[ base + 2 * totalCells + left ]
                          + cellGradient[ base + 2 * totalCells + right ] );
    }

    // --- CalcGradCoef -------------------------------------------------------
    Real dx = xCellR - xCellL;
    Real dy = yCellR - yCellL;
    Real dz = zCellR - zCellL;
    const Real ods = ViscousUnitizeDisplacement( dx, dy, dz );

    // --- AverFaceValue ------------------------------------------------------
    // The energy flux uses the face-averaged velocity, not either side's.
    const Real uFace = 0.5 * ( cellState[ 1 * totalCells + left ]
                             + cellState[ 1 * totalCells + right ] );
    const Real vFace = 0.5 * ( cellState[ 2 * totalCells + left ]
                             + cellState[ 2 * totalCells + right ] );
    const Real wFace = 0.5 * ( cellState[ 3 * totalCells + left ]
                             + cellState[ 3 * totalCells + right ] );

    Real dTdx = 0.5 * ( temperatureGradient[ left ]
                      + temperatureGradient[ right ] );
    Real dTdy = 0.5 * ( temperatureGradient[ totalCells + left ]
                      + temperatureGradient[ totalCells + right ] );
    Real dTdz = 0.5 * ( temperatureGradient[ 2 * totalCells + left ]
                      + temperatureGradient[ 2 * totalCells + right ] );

    // --- CorrectFaceGrad ----------------------------------------------------
    // CorrectFaceGradient (codes/usolver/src/UGradPolicy.cpp:86): project the
    // averaged gradient onto the cell-centre direction, replace that component
    // with the finite difference across the face, and add the correction back
    // along the same unit direction.
    for ( int m = 0; m < 3; ++ m )
    {
        const int eq = 1 + m;
        CorrectViscousFaceGradient(
            cellState[ eq * totalCells + left ],
            cellState[ eq * totalCells + right ],
            dqdx[ m ], dqdy[ m ], dqdz[ m ],
            dx, dy, dz, ods );
    }
    CorrectViscousFaceGradient(
        temperature[ left ], temperature[ right ],
        dTdx, dTdy, dTdz, dx, dy, dz, ods );

    // --- CalcNormalGrad -----------------------------------------------------
    const Real dTdn = xfn * dTdx + yfn * dTdy + zfn * dTdz;

    // --- PrepareFaceValue: viscosity averages -------------------------------
    const Real visl = 0.5 * ( laminarViscosity[ left ]
                            + laminarViscosity[ right ] );
    const Real vist = turbulentViscosity != nullptr
        ? 0.5 * ( turbulentViscosity[ left ] + turbulentViscosity[ right ] )
        : 0.0;
    const Real vis = visl + vist;

    // --- StokesStressTensor / ProjectStressMomentumFlux ----------------------
    Real txx = 0.0;
    Real tyy = 0.0;
    Real tzz = 0.0;
    Real txy = 0.0;
    Real txz = 0.0;
    Real tyz = 0.0;
    ViscousStokesStress(
        vis,
        dqdx[ 0 ], dqdy[ 0 ], dqdz[ 0 ],
        dqdx[ 1 ], dqdy[ 1 ], dqdz[ 1 ],
        dqdx[ 2 ], dqdy[ 2 ], dqdz[ 2 ],
        txx, tyy, tzz, txy, txz, tyz );
    Real fluxMx = 0.0;
    Real fluxMy = 0.0;
    Real fluxMz = 0.0;
    ViscousProjectStressMomentumFlux(
        xfn, yfn, zfn,
        txx, tyy, tzz, txy, txz, tyz,
        fluxMx, fluxMy, fluxMz );

    // --- AddHeatFlux ---------------------------------------------------------
    // The host zeroes qx/qy/qz (AddChemHeatFlux is empty for chemModel != 1),
    // adds FaceNormalDot( n, 0, 0, 0 ) == 0, then adds the conduction term, so
    // qNormal reduces to kcp * dTdn.  Both dropped terms are identically zero
    // in this path, not unmodelled physics.
    const Real kcp = ViscousEffectiveThermalConductivity(
        visl, vist, inversePrandtlLaminar, inversePrandtlTurbulent, specificHeatCp );
    qNormal = ViscousConductiveHeatFluxNormal( kcp, dTdn );

    // --- UpdateFaceVisFlux ---------------------------------------------------
    // The host sets fvis[IR] = 0 and scales every component, so the energy term
    // consumes the pre-coefficient momentum flux.
    const Real coeff = ViscousFaceFluxCoefficientDevice(
        inverseReynolds, faceArea );
    flux[ 0 ] = coeff * 0.0;
    flux[ 1 ] = coeff * fluxMx;
    flux[ 2 ] = coeff * fluxMy;
    flux[ 3 ] = coeff * fluxMz;
    flux[ 4 ] = coeff * ViscousEnergyFluxDevice(
        uFace, vFace, wFace, fluxMx, fluxMy, fluxMz, qNormal );
}

// Face normal orientation.
//
// The host's GCom::SetGeometry has an rc == cId swap that reverses the normal,
// but the viscous chain never calls it: UNsVisFlux::PrepareFaceValue assigns
// gcom.xfn straight from ug.xfn[fId] and both cell ids from ug.lcf/ug.rcf.  The
// kernel therefore uses the stored normal unchanged, which is what keeps the
// stress and heat projections on the host sign.
//
// Whole-face-array kernel.  The per-face body is the same function the
// single-face contract entry calls, so the two cannot drift.
__global__ void ViscousFluxKernel(
    const Real * cellState,
    const Real * cellGradient,
    const Real * temperature,
    const Real * temperatureGradient,
    const Real * laminarViscosity,
    const Real * turbulentViscosity,
    int totalCells,
    const Real * xNormal,
    const Real * yNormal,
    const Real * zNormal,
    const Real * faceArea,
    const Real * xCell,
    const Real * yCell,
    const Real * zCell,
    const int * leftCell,
    const int * rightCell,
    int nEquations,
    int nFaces,
    Real inversePrandtlLaminar,
    Real inversePrandtlTurbulent,
    Real specificHeatCp,
    Real inverseReynolds,
    Real * faceFlux,
    int faceFluxStride,
    Real * wallHeatFlux,
    const int * wallSlot,
    int nWallFaces )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nFaces ) return;

    const int left = leftCell[ face ];
    const int right = rightCell[ face ];
    Real flux[ 5 ] = {};
    Real qNormal = 0.0;
    ViscousFaceFlux(
        cellState, cellGradient, temperature, temperatureGradient,
        laminarViscosity, turbulentViscosity, totalCells,
        xNormal[ face ], yNormal[ face ], zNormal[ face ], faceArea[ face ],
        xCell[ left ], yCell[ left ], zCell[ left ],
        xCell[ right ], yCell[ right ], zCell[ right ],
        left, right,
        inversePrandtlLaminar, inversePrandtlTurbulent,
        specificHeatCp, inverseReynolds,
        flux, qNormal );

    for ( int eq = 0; eq < nEquations; ++ eq )
    {
        faceFlux[ eq * faceFluxStride + face ] = flux[ eq ];
    }

    if ( wallHeatFlux != nullptr && wallSlot != nullptr )
    {
        const int slot = wallSlot[ face ];
        // UNsVisFlux::SaveHeatFlux: non_dim_heatflux = -oreynolds * qNormal.
        if ( slot >= 0 && slot < nWallFaces )
        {
            wallHeatFlux[ slot ] = ViscousWallHeatFlux(
                inverseReynolds, qNormal );
        }
    }
}

__global__ void EvaluateDeviceSutherlandKernel(
    Real temperatureNondim, Real sutherlandTsNondim, Real * result )
{
    *result = DeviceSutherlandViscosity(
        temperatureNondim, sutherlandTsNondim );
}

}

void HipFluxBackend::CalcViscousFaceFlux( const ViscousFluxView & view )
{
    // The FluxBackend entry has no boundary classification, so it cannot
    // produce wall heat flux.  Fail closed rather than leave the caller's
    // buffer untouched, which would hand the solver stale values.
    if ( view.wallHeatFlux != nullptr || view.nWallFaces != 0 )
    {
        throw std::invalid_argument(
            "HIP viscous flux: wall heat flux needs the solid face ids, so use "
            "the CalcViscousFaceFlux( view, wallFaceIds ) overload" );
    }
    CalcViscousFaceFlux( view, nullptr );
}

void HipFluxBackend::CalcViscousFaceFlux(
    const ViscousFluxView & view, const int * wallFaceIds )
{
    ValidateViscousFluxView( view );
    if ( ( view.wallHeatFlux != nullptr ) != ( wallFaceIds != nullptr ) )
    {
        throw std::invalid_argument(
            "HIP viscous flux: wall heat flux needs both an output buffer and "
            "the solid face ids" );
    }
    if ( view.nEquations != 5 )
    {
        // The stress tensor and the energy flux are the five-equation NS form.
        // A three-equation state would silently index past its own planes.
        throw std::invalid_argument(
            "HIP viscous flux requires the five-equation primitive state" );
    }
    if ( view.nWallFaces < 0 || view.nWallFaces > view.nFaces )
    {
        throw std::invalid_argument(
            "HIP viscous flux: wall face count is out of range" );
    }
    if ( view.faceFluxStride < view.nEquations )
    {
        throw std::invalid_argument(
            "HIP viscous flux: output stride is shorter than the equation "
            "count" );
    }
    // ValidateViscousFluxView deliberately stops at the fields every execution
    // space needs; the cell-centre coordinates are only read by the face
    // construction here (corrected gradient direction), and a null one would
    // read garbage rather than fail.
    if ( view.xCell == nullptr || view.yCell == nullptr
         || view.zCell == nullptr )
    {
        throw std::invalid_argument(
            "HIP viscous flux requires cell-centre coordinates" );
    }
    if ( view.q[ 1 ] == nullptr || view.q[ 2 ] == nullptr
         || view.q[ 3 ] == nullptr )
    {
        throw std::invalid_argument(
            "HIP viscous flux requires the full five-component primitive "
            "state" );
    }
    for ( int eq = 0; eq < view.nEquations; ++ eq )
    {
        if ( view.dqdx[ eq ] == nullptr || view.dqdy[ eq ] == nullptr
             || view.dqdz[ eq ] == nullptr )
        {
            throw std::invalid_argument(
                "HIP viscous flux requires a complete primitive gradient" );
        }
    }

    const int nFaces = view.nFaces;
    const int nEq = view.nEquations;
    const int totalCells = view.nCells + view.nGhostCells;

    // Resident path: when the caller has bound a state for this cache key, the
    // device already owns the primitive state, its gradient, the temperature
    // and its gradient, and all of the geometry.  Uploading them again is the
    // whole cost of this call - 26 cell-sized planes for a 300k-cell zone, so
    // roughly 60 MB and 5-9 ms per call at measured copy rates - which is more
    // than the host per-face loop it replaces.  Use the resident storage and
    // skip the upload entirely.
    //
    // The fallback below stays for callers that legitimately have no bound
    // state (the contract tests drive the kernel directly with hand-built
    // arrays), and because an unbound call must still work rather than read
    // whatever an earlier call left behind.
    HipGradientStorage * resident = nullptr;
    if ( view.cacheKey != nullptr )
    {
        const auto binding = boundStates_.find( view.cacheKey );
        if ( binding != boundStates_.end() )
        {
            auto * hipState = dynamic_cast< HipNs3DBackendState * >(
                binding->second->BackendState() );
            if ( hipState != nullptr )
            {
                HipGradientStorage & candidate = hipState->Gradient();
                const std::size_t stateValues =
                    static_cast< std::size_t >( nEq ) * totalCells;
                const std::size_t gradientValues =
                    static_cast< std::size_t >( nEq ) * 3 * totalCells;
                const bool geometryReady =
                    candidate.xNormal.Size() >= static_cast< std::size_t >( nFaces )
                    && candidate.faceArea.Size()
                        >= static_cast< std::size_t >( nFaces )
                    && candidate.leftCell.Size()
                        >= static_cast< std::size_t >( nFaces )
                    && candidate.rightCell.Size()
                        >= static_cast< std::size_t >( nFaces )
                    && candidate.xCellCenter.Size()
                        >= static_cast< std::size_t >( totalCells );
                const bool stateReady =
                    candidate.cellState.Size() >= stateValues
                    && candidate.cellGradient.Size() >= gradientValues
                    && candidate.temperature.Size()
                        >= static_cast< std::size_t >( totalCells )
                    && candidate.temperatureGradient.Size()
                        >= static_cast< std::size_t >( 3 ) * totalCells;
                if ( geometryReady && stateReady )
                {
                    resident = &candidate;
                }
            }
        }
    }

    // The primitive gradient is uploaded on both paths rather than taken from
    // the resident storage.
    //
    // storage->cellGradient is produced by GradientCellKernel from the resident
    // state, which is a different computation from the dqdx/dqdy/dqdz the host
    // chain hands the viscous flux: the solver's UNsGrad owns those, with its own
    // limiter and face correction.  Reusing the resident gradient to skip this
    // upload changed the answer.  Isolating the hooks one at a time measured it:
    // the viscous term alone deviated 6.7e-3 relative from the host arm, while
    // the four SA hooks were bit-identical.
    viscousCellGradient_.Resize(
        static_cast< std::size_t >( nEq ) * 3 * totalCells );
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const std::size_t base =
            static_cast< std::size_t >( eq ) * 3 * totalCells;
        CopyRealComponent( viscousCellGradient_, base, view.dqdx[ eq ],
                           totalCells, AccelCopyKind::HostToDevice );
        CopyRealComponent( viscousCellGradient_, base + totalCells,
                           view.dqdy[ eq ], totalCells,
                           AccelCopyKind::HostToDevice );
        CopyRealComponent( viscousCellGradient_, base + 2 * totalCells,
                           view.dqdz[ eq ], totalCells,
                           AccelCopyKind::HostToDevice );
    }

    // State, gradient, temperature and its gradient all come from the caller on
    // both paths.
    //
    // These are not the same values the device holds.  Upload() puts a snapshot
    // of the state there, but the host chain advances state, temperature and
    // gradients at different points in the step, and the resident copies are
    // whatever the last device update left.  Using them measured 6.7e-3
    // relative error against the host arm on the viscous hook alone, with the
    // four SA hooks bit-identical - so the drift was here, not in a formula.
    // Only the geometry is shared, because topology does not change per step.
    viscousCellState_.Resize(
        static_cast< std::size_t >( nEq ) * totalCells );
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        CopyRealComponent( viscousCellState_,
                           static_cast< std::size_t >( eq ) * totalCells,
                           view.q[ eq ], totalCells, AccelCopyKind::HostToDevice );
    }
    viscousTemperature_.CopyFromHost( view.temperature, totalCells );
    viscousTemperatureGradient_.Resize(
        static_cast< std::size_t >( 3 ) * totalCells );
    CopyRealComponent( viscousTemperatureGradient_, 0, view.dTdx,
                       totalCells, AccelCopyKind::HostToDevice );
    CopyRealComponent( viscousTemperatureGradient_, totalCells, view.dTdy,
                       totalCells, AccelCopyKind::HostToDevice );
    CopyRealComponent( viscousTemperatureGradient_, 2 * totalCells, view.dTdz,
                       totalCells, AccelCopyKind::HostToDevice );

    if ( resident == nullptr )
    {
    { ScopedStageTimer timer( "v_upload" );
    viscousFaceNormals_.Resize( static_cast< std::size_t >( 3 ) * nFaces );
    const Real * normals[ 3 ] = { view.xNormal, view.yNormal, view.zNormal };
    std::vector< Real > zeroNormal;
    for ( int component = 0; component < 3; ++ component )
    {
        // The validator only requires xNormal, so an absent y/z normal is legal
        // and means zero.  Upload it explicitly: DeviceBuffer::Resize is
        // grow-only and does not zero a fresh allocation.
        if ( normals[ component ] == nullptr )
        {
            zeroNormal.assign( nFaces, 0.0 );
        }
        CopyRealComponent(
            viscousFaceNormals_,
            static_cast< std::size_t >( component ) * nFaces,
            normals[ component ] != nullptr ? normals[ component ]
                                            : zeroNormal.data(),
            nFaces, AccelCopyKind::HostToDevice );
    }
    viscousFaceArea_.CopyFromHost( view.faceArea, nFaces );
    viscousCellCenter_.Resize( static_cast< std::size_t >( 3 ) * totalCells );
    CopyRealComponent( viscousCellCenter_, 0, view.xCell, totalCells,
                       AccelCopyKind::HostToDevice );
    CopyRealComponent( viscousCellCenter_, totalCells, view.yCell, totalCells,
                       AccelCopyKind::HostToDevice );
    CopyRealComponent( viscousCellCenter_, 2 * totalCells, view.zCell,
                       totalCells, AccelCopyKind::HostToDevice );
    viscousFaceCells_.Resize( static_cast< std::size_t >( 2 ) * nFaces );
    AccelRuntime::Instance().Backend().Copy(
        viscousFaceCells_.Data(), view.leftCell, nFaces * sizeof( int ),
        AccelCopyKind::HostToDevice );
    AccelRuntime::Instance().Backend().Copy(
        viscousFaceCells_.Data() + nFaces, view.rightCell,
        nFaces * sizeof( int ), AccelCopyKind::HostToDevice );
    }
    }

    // Source pointers.  Only the geometry may come from the resident storage:
    // Upload() put exactly those values there and they are static per topology.
    // The state, the gradient and the temperature all come from the caller,
    // because each is updated by a different part of the host chain at a
    // different point in the step, and the resident copies are whatever the last
    // device update left behind rather than what this call was handed.
    const Real * deviceState = viscousCellState_.Data();
    const Real * deviceGradient = viscousCellGradient_.Data();
    const Real * deviceTemperature = viscousTemperature_.Data();
    const Real * deviceTemperatureGradient = viscousTemperatureGradient_.Data();
    // Normals and cell centres are three separate allocations in the resident
    // storage, each holding exactly one plane.  The fallback staging buffers
    // pack them contiguously, so the kernel is handed base + nFaces and
    // base + 2*nFaces.  Doing that to the resident pointers read past the end of
    // the x plane into whatever the allocator placed next; it only looked right
    // when the buffers happened to be adjacent.
    const Real * deviceNormalX = resident != nullptr
        ? resident->xNormal.Data() : viscousFaceNormals_.Data();
    const Real * deviceNormalY = resident != nullptr
        ? resident->yNormal.Data() : viscousFaceNormals_.Data() + nFaces;
    const Real * deviceNormalZ = resident != nullptr
        ? resident->zNormal.Data() : viscousFaceNormals_.Data() + 2 * nFaces;
    const Real * deviceFaceArea = resident != nullptr
        ? resident->faceArea.Data() : viscousFaceArea_.Data();
    const Real * deviceCellX = resident != nullptr
        ? resident->xCellCenter.Data() : viscousCellCenter_.Data();
    const Real * deviceCellY = resident != nullptr
        ? resident->yCellCenter.Data() : viscousCellCenter_.Data() + totalCells;
    const Real * deviceCellZ = resident != nullptr
        ? resident->zCellCenter.Data()
        : viscousCellCenter_.Data() + 2 * totalCells;
    const int * deviceLeftCell = resident != nullptr
        ? resident->leftCell.Data() : viscousFaceCells_.Data();
    const int * deviceRightCell = resident != nullptr
        ? resident->rightCell.Data() : viscousFaceCells_.Data() + nFaces;

    // Viscosity is a per-call input, not resident state: it is Sutherland on the
    // laminar side and the SA solve's output on the turbulent side, so it has no
    // slot in HipGradientStorage.  These two planes stay uploaded in either
    // path.  They are 2 of the 26 planes the resident path avoids, which is why
    // this is still a large win rather than a rewrite.
    viscousLaminar_.CopyFromHost( view.laminarViscosity, totalCells );
    if ( view.turbulentViscosity != nullptr )
    {
        viscousTurbulent_.CopyFromHost( view.turbulentViscosity, totalCells );
    }
    else
    {
        viscousTurbulent_.Reset();
    }
    const Real * deviceLaminar = viscousLaminar_.Data();
    const Real * deviceTurbulent = view.turbulentViscosity != nullptr
        ? viscousTurbulent_.Data() : nullptr;

    // Wall classification: face -> slot in the caller's wall list.  Built once
    // per call so the kernel does not scan the wall list per face.
    std::vector< int > wallSlot;
    DeviceBuffer< int > deviceWallSlot;
    if ( wallFaceIds != nullptr && view.nWallFaces > 0 )
    {
        wallSlot.assign( nFaces, -1 );
        for ( int slot = 0; slot < view.nWallFaces; ++ slot )
        {
            const int face = wallFaceIds[ slot ];
            if ( face < 0 || face >= nFaces )
            {
                throw std::invalid_argument(
                    "HIP viscous flux: wall face id is out of range" );
            }
            if ( wallSlot[ face ] != -1 )
            {
                throw std::invalid_argument(
                    "HIP viscous flux: duplicate wall face id" );
            }
            wallSlot[ face ] = slot;
        }
        deviceWallSlot.CopyFromHost( wallSlot.data(), nFaces );
    }
    // The kernel writes the wall heat flux on the device: the view's buffer is
    // a host pointer and cannot be dereferenced by a kernel (VMFault).
    if ( view.wallHeatFlux != nullptr && view.nWallFaces > 0 )
    {
        viscousWallHeatFlux_.Resize(
            static_cast< std::size_t >( view.nWallFaces ) );
    }

    // The kernel writes equation-major with the caller's stride, i.e. the
    // highest index is ( nEquations - 1 ) * stride + ( nFaces - 1 ).  Sizing by
    // stride * nFaces is only correct when stride == nEquations; for a wider
    // stride it under-allocates and the last planes land outside the buffer.
    viscousFaceFlux_.Resize(
        static_cast< std::size_t >( nEq - 1 ) * view.faceFluxStride
        + static_cast< std::size_t >( nFaces ) );
    {
        ScopedStageTimer timer( "v_kernel" );
        const int blockSize = 256;
        const int gridSize = ( nFaces + blockSize - 1 ) / blockSize;
        hipLaunchKernelGGL(
            ViscousFluxKernel,
            dim3( gridSize ), dim3( blockSize ), 0, 0,
            deviceState,
            deviceGradient,
            deviceTemperature,
            deviceTemperatureGradient,
            deviceLaminar,
            view.turbulentViscosity != nullptr
                ? deviceTurbulent : nullptr,
            totalCells,
            deviceNormalX,
            deviceNormalY,
            deviceNormalZ,
            deviceFaceArea,
            deviceCellX,
            deviceCellY,
            deviceCellZ,
            deviceLeftCell,
            deviceRightCell,
            nEq, nFaces,
            view.inversePrandtlLaminar, view.inversePrandtlTurbulent,
            view.specificHeatCp, view.inverseReynolds,
            viscousFaceFlux_.Data(), view.faceFluxStride,
            view.wallHeatFlux != nullptr && view.nWallFaces > 0
                ? viscousWallHeatFlux_.Data() : nullptr,
            deviceWallSlot.Size() > 0 ? deviceWallSlot.Data() : nullptr,
            view.nWallFaces );
        CheckHip( hipGetLastError(), "ViscousFluxKernel launch" );
        HipDrain( "ViscousFluxKernel synchronize" );
    }

    { ScopedStageTimer timer( "v_download" );
    // Copy only the nFaces columns of each equation plane: the caller's stride
    // may exceed nEquations, and a whole-plane copy would overrun its buffer.
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        AccelRuntime::Instance().Backend().Copy(
            view.faceFlux
                + static_cast< std::size_t >( eq ) * view.faceFluxStride,
            viscousFaceFlux_.Data()
                + static_cast< std::size_t >( eq ) * view.faceFluxStride,
            nFaces * sizeof( Real ), AccelCopyKind::DeviceToHost );
    }
    if ( view.wallHeatFlux != nullptr && view.nWallFaces > 0 )
    {
        viscousWallHeatFlux_.CopyToHost(
            view.wallHeatFlux, view.nWallFaces );
    }
    }
}

Real HipFluxBackend::EvaluateDeviceSutherland(
    Real temperatureNondim, Real sutherlandTsNondim ) const
{
    Real * deviceResult = nullptr;
    CheckHip( hipMalloc( &deviceResult, sizeof( Real ) ),
              "sutherland contract alloc" );
    hipLaunchKernelGGL(
        EvaluateDeviceSutherlandKernel, dim3( 1 ), dim3( 1 ), 0, 0,
        temperatureNondim, sutherlandTsNondim, deviceResult );
    CheckHip( hipGetLastError(), "sutherland contract kernel launch" );
    CheckHip( hipDeviceSynchronize(),
              "sutherland contract kernel synchronize" );
    Real result = 0.0;
    CheckHip( hipMemcpy( &result, deviceResult, sizeof( Real ),
                         hipMemcpyDeviceToHost ), "sutherland contract copy" );
    CheckHip( hipFree( deviceResult ), "sutherland contract free" );
    return result;
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
    if ( ! HipDeferStatusEnabled() )
    {
        storage->timeStepStatus.CopyToHost( & status, 1 );
        if ( status != 0 )
        {
            throw std::runtime_error(
                "HIP local timestep produced an invalid spectral radius" );
        }
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
    HipDrain( "local timestep kernel synchronize" );

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
    bool inflowChanged = storage->cachedInflowCount != view.nEquations;
    for ( int equation = 0;
          ! inflowChanged && equation < view.nEquations; ++ equation )
    {
        inflowChanged = storage->cachedInflowValues[ equation ]
            != view.inflow[ equation ];
    }
    if ( inflowChanged )
    {
        storage->inflow.CopyFromHost( view.inflow, view.nEquations );
        for ( int equation = 0; equation < view.nEquations; ++ equation )
            storage->cachedInflowValues[ equation ] = view.inflow[ equation ];
        storage->cachedInflowCount = view.nEquations;
    }
    storage->boundaryFaceState.Resize(
        view.nEquations * view.nBoundaryFaces );
    storage->wallStateValid = view.deviceAuthoritativeHostMirror;
    storage->wallStateFaces = view.nBoundaryFaces;
    storage->wallStateEquations = view.nEquations;
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
    // The boundary stage just produced bcQ on the device.  With a
    // device-authoritative mirror nobody reads the host copy: the next
    // reconstruction uploads it straight back.  Skip both ends of that round
    // trip (this D2H and the matching H2D) instead of shuffling it per stage.
    const bool bcQStaysOnDevice =
        view.deviceAuthoritativeHostMirror && storage->wallStateValid;
    const bool mirrorBcQ =
        view.bcQ[ 0 ] != nullptr && ! bcQStaysOnDevice;
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
    HipGradientStorage * storage = &fallbackGradientStorage_;
    std::uint64_t topologyGeneration = view.topologyGeneration;
    {
        ScopedStageTimer reconPrepTimer( "r_prep" );
        ValidateCellFaceReconstructionView( view );
        const void * cacheKey = view.cacheKey != nullptr
            ? view.cacheKey : view.ownerToken;
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
            topologyGeneration = owner->TopologyGeneration();
        }
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
    }

    // The limiter factor is still produced on the host and uploaded; the
    // device reconstruction consumes exactly that field.
    const bool limiterActive =
        view.limiterMode != ReconstructionLimiterMode::Disabled;
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
    // Set inside the solid branch below and consumed at the kernel launch.
    bool deviceOwnsWallState = false;
    if ( hasSolid )
    {
        // Solid-override faces are a prefix of the face array. Callers pack
        // bcQ as one contiguous equation-major buffer of
        // nEquations*nBoundaryFaces values so this can be one transfer.
        // Indexing that buffer with a full-face stride reads the wrong
        // equation for every face after the first.
        for ( int equation = 1; equation < view.nEquations; ++ equation )
        {
            if ( view.bcQ[ equation ]
                 != view.bcQ[ 0 ]
                    + static_cast< std::ptrdiff_t >( equation )
                        * view.nBoundaryFaces )
            {
                throw std::invalid_argument(
                    "HIP reconstruction bcQ must be contiguous and "
                    "equation-major with stride nBoundaryFaces" );
            }
        }
        ScopedStageTimer bcqTimer( "recon_bcq" );
        const std::size_t bcQValues =
            static_cast< std::size_t >( view.nEquations )
            * view.nBoundaryFaces;
        // The device boundary stage leaves the newest wall state in
        // boundaryFaceState, in the same equation-major [equation][face]
        // layout this kernel reads.  Point at it and skip the host round trip
        // (a D2H mirror plus an H2D re-upload per stage); otherwise fall back
        // to uploading the caller's packed buffer.
        deviceOwnsWallState =
            storage->wallStateValid
            && storage->wallStateFaces == view.nBoundaryFaces
            && storage->wallStateEquations == view.nEquations
            && storage->boundaryFaceState.Size() == bcQValues;
        if ( ! deviceOwnsWallState )
        {
            storage->bcQ.Resize( bcQValues );
            AccelRuntime::Instance().Backend().Copy(
                storage->bcQ.Data(), view.bcQ[ 0 ],
                bcQValues * sizeof( Real ),
                AccelCopyKind::HostToDevice );
        }
    }
    else
    {
        storage->bcQ.Reset();
        storage->wallStateValid = false;
    }
    // Compute the limiter factor on the device from the resident gradient.
    // The host copy is not usable here: with a device-authoritative q the host
    // gradient mirror is deliberately stale, and a factor built from it biased
    // the solution by 8.5e-5 over 200 steps.
    if ( limiterActive )
    {
        // The factor kernel walks each cell's face list, so the cell->face CSR
        // must exist before it runs; the flux stage used to build it later.
        EnsureCellFaceConnectivity(
            *storage, view.leftCell, view.rightCell,
            view.nCells, view.nFaces, view.nBoundaryFaces );
        ScopedStageTimer limiterTimer( "recon_limiter_kernel" );
        const std::size_t limiterValues =
            static_cast< std::size_t >( view.nEquations ) * totalCells;
        storage->limiter.Resize( limiterValues );
        const int limiterBlock = 128;
        hipLaunchKernelGGL(
            BarthLimiterKernel,
            dim3( ( view.nCells + limiterBlock - 1 ) / limiterBlock ),
            dim3( limiterBlock ), 0, 0,
            storage->cellState.Data(), storage->cellGradient.Data(),
            storage->xFaceCenter.Data(), storage->yFaceCenter.Data(),
            storage->zFaceCenter.Data(), storage->xCellCenter.Data(),
            storage->yCellCenter.Data(), storage->zCellCenter.Data(),
            storage->xNormal.Data(), storage->yNormal.Data(),
            storage->zNormal.Data(),
            storage->leftCell.Data(), storage->rightCell.Data(),
            storage->cellFaceOffsets.Data(), storage->cellFaceIds.Data(),
            storage->boundaryMaskValid ? storage->boundaryMask.Data() : nullptr,
            storage->boundaryOperation.Data(),
            view.nCells, totalCells, view.nFaces, view.nEquations,
            storage->limiter.Data() );
        CheckHip( hipGetLastError(), "Barth limiter kernel launch" );
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
        hasSolid ? ( deviceOwnsWallState
                         ? storage->boundaryFaceState.Data()
                         : storage->bcQ.Data() )
                 : nullptr,
        view.nBoundaryFaces,
        limiterActive ? storage->limiter.Data() : nullptr,
        limiterActive ? 1 : 0,
        view.nCells, view.nBoundaryFaces, view.nFaces, view.nEquations,
        view.physicality == ReconstructionPhysicalityPolicy::None ? 0 : 1,
        view.densityComponent, view.pressureComponent,
        storage->qf1.Data(), storage->qf2.Data() );
    CheckHip( hipGetLastError(), "reconstruction kernel launch" );
    HipDrain( "reconstruction kernel synchronize" );

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
    bool downloadResidual,
    ResidualBase residualBase,
    Real entropy1,
    Real entropy2,
    bool entropyFix )
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
    if ( residualBase == ResidualBase::DeviceZero )
    {
        // The device chain owns the base: start from zero on the device.  The
        // host residual field still receives the final values when the caller
        // asks for a download below.
        ScopedStageTimer zeroTimer( "residual_zero" );
        const int zeroBlock = 256;
        hipLaunchKernelGGL(
            ZeroValuesKernel,
            dim3( ( totalResidual + zeroBlock - 1 ) / zeroBlock ),
            dim3( zeroBlock ), 0, 0,
            storage->residual.Data(), totalResidual );
        CheckHip( hipGetLastError(), "residual zero kernel launch" );
    }
    else if ( residual.values != nullptr )
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
        entropy1, entropy2, entropyFix ? 1 : 0,
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
    HipDrain( "reconstructed primitive flux kernel synchronize" );
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
    HipDrain( "state residual scale kernel synchronize" );
    storage.residualMetricsValid = true;
}

void HipFluxBackend::RequestResidentStatusCheck()
{
    // The read itself is a blocking device->host transfer, i.e. a full queue
    // drain.  Doing it at the end of the macro step means the driver has to
    // empty a queue that nothing else needs to wait for.  Record the request
    // and let the host-staged interface pack (which already drains for its own
    // copy) consume it, so the check keeps its per-step cadence without paying
    // a second drain.
    statusCheckPending_ = true;
}

void HipFluxBackend::FlushRequestedResidentStatus()
{
    if ( ! statusCheckPending_ ) return;
    statusCheckPending_ = false;
    CheckResidentStatus();
}

void HipFluxBackend::CheckResidentStatus()
{
    // Only bound states carry the flags; the fallback storage is used by the
    // contract tests, which check their own status synchronously.
    for ( auto & entry : boundStates_ )
    {
        auto * hipState = dynamic_cast< HipNs3DBackendState * >(
            entry.second->BackendState() );
        if ( hipState == nullptr ) continue;
        HipGradientStorage & storage = hipState->Gradient();
        // A zone that never ran the resident time step (or never updated the
        // state) has an empty flag buffer; reading it would trip the copy
        // guard rather than report a status.
        if ( storage.timeStepStatus.Size() < 1
             && storage.stateUpdateStatus.Size() < 1 )
        {
            continue;
        }
        int zero = 0;
        if ( storage.timeStepStatus.Size() >= 1 )
        {
            int timeStepStatus = 0;
            storage.timeStepStatus.CopyToHost( & timeStepStatus, 1 );
            if ( timeStepStatus != 0 )
            {
                throw std::runtime_error(
                    "HIP local timestep produced an invalid spectral radius" );
            }
            storage.timeStepStatus.CopyFromHost( & zero, 1 );
        }
        if ( storage.stateUpdateStatus.Size() >= 1 )
        {
            int updateStatus = 0;
            storage.stateUpdateStatus.CopyToHost( & updateStatus, 1 );
            if ( updateStatus != 0 )
            {
                throw std::runtime_error(
                    "HIP primitive state update produced a non-physical state" );
            }
            storage.stateUpdateStatus.CopyFromHost( & zero, 1 );
        }
    }
}

void HipFluxBackend::ScaleAndUpdatePrimitiveState(
    const CellStateUpdateView & view,
    bool downloadToHost )
{
    // Issuing the scale and the update from one entry removes a host entry and
    // its queue drain per stage without touching either formula: the two
    // kernels keep their own launch geometry and the residual field is written
    // back exactly as before.
    ValidateCellStateUpdateView( view );
    const auto binding = boundStates_.find( view.cacheKey );
    if ( binding == boundStates_.end() )
    {
        throw std::logic_error(
            "HIP fused scale/update requires a bound Ns3D state" );
    }
    Ns3DDeviceState * owner = binding->second;
    auto * hipState = dynamic_cast< HipNs3DBackendState * >(
        owner->BackendState() );
    if ( hipState == nullptr )
    {
        throw std::logic_error(
            "bound HIP state update is not HIP-owned" );
    }
    HipGradientStorage & storage = hipState->Gradient();
    const int totalCells = view.nCells + view.nGhostCells;
    if ( ! storage.residualValid
         || storage.residualCells != view.nCells
         || storage.residualEquations != view.nEquations
         || storage.cachedCells != view.nCells
         || storage.cachedGhostCells != view.nGhostCells
         || storage.cachedVolume != view.cellVolume
         || storage.cellVolume.Size()
             < static_cast< std::size_t >( view.nCells )
         || storage.cellState.Size()
             < static_cast< std::size_t >( totalCells * view.nEquations ) )
    {
        throw std::logic_error(
            "HIP fused scale/update requires resident state and residual" );
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

    const int blockSize = 256;
    {
        ScopedStageTimer timer( "fused_scale" );
        const int total = view.nCells * view.nEquations;
        hipLaunchKernelGGL(
            ScaleResidualKernel,
            dim3( ( total + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0,
            storage.residual.Data(), storage.timeStep.Data(),
            storage.cellVolume.Data(), view.nCells, view.nEquations,
            view.rkCoefficient );
        CheckHip( hipGetLastError(), "state residual scale kernel launch" );
    }
    storage.residualMetricsValid = true;

    int status = 0;
    const bool deferStatus = HipDeferStatusEnabled();
    // Deferred status checking still needs a live flag on the first call.
    // Leave an existing flag untouched between checks so a later bad stage
    // cannot be erased before CheckResidentStatus() observes it.
    if ( ! deferStatus || storage.stateUpdateStatus.Size() < 1 )
    {
        storage.stateUpdateStatus.CopyFromHost( & status, 1 );
    }
    {
        ScopedStageTimer timer( "state_update" );
        hipLaunchKernelGGL(
            UpdatePrimitiveStateKernel,
            dim3( ( view.nCells + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0,
            storage.cellState.Data(), storage.residual.Data(),
            view.nCells, totalCells, view.gamma,
            storage.stateUpdateStatus.Data() );
        CheckHip( hipGetLastError(), "primitive state update kernel launch" );
    }
    if ( ! deferStatus )
    {
        storage.stateUpdateStatus.CopyToHost( & status, 1 );
        if ( status != 0 )
        {
            throw std::runtime_error(
                "HIP primitive state update produced a non-physical state" );
        }
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
    owner->MarkFieldChanged();
    storage.residualValid = false;
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

    if ( storage.residualScalePending )
    {
        CellStateUpdateView fused = view;
        fused.rkCoefficient = storage.residualScaleCoefficient;
        fused.preferResidentTimeStep = true;
        storage.residualScalePending = false;
        ScaleAndUpdatePrimitiveState( fused, downloadToHost );
        return;
    }
    int status = 0;
    const bool deferStatus = HipDeferStatusEnabled();
    // Deferred status checking still needs a live flag on the first call.
    // Leave an existing flag untouched between checks so a later bad stage
    // cannot be erased before CheckResidentStatus() observes it.
    if ( ! deferStatus || storage.stateUpdateStatus.Size() < 1 )
    {
        storage.stateUpdateStatus.CopyFromHost( & status, 1 );
    }
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
    if ( ! deferStatus )
    {
        storage.stateUpdateStatus.CopyToHost( & status, 1 );
        if ( status != 0 )
        {
            throw std::runtime_error(
                "HIP primitive state update produced a non-physical state" );
        }
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
        HipDrain( "boundary pack kernel synchronize" );
        // Same one-shot consume as the host-staged path below: the pending
        // status read is guarded by its own flag, so only the first endpoint
        // after a macro step pays for it.
        FlushRequestedResidentStatus();
        return {
            boundary.values.Data(),
            valueCount,
            BoundaryMemorySpace::Device };
    }

    if ( boundary.hostValues.size() != valueCount )
        boundary.hostValues.resize( valueCount );
    if ( valueCount != 0 )
    {
        // Pageable host memory stages through a driver buffer and measured at
        // roughly 200 MB/s on this platform, which made this per-endpoint copy
        // dominate the interface round.  Transfer into the shared pinned
        // staging block and memcpy the short distance on the host instead; the
        // endpoint buffer keeps its own ownership and lifetime.
        Real * staging = LocalPinnedStaging().Get( valueCount );
        boundary.values.CopyToHost( staging, valueCount );
        std::copy(
            staging, staging + valueCount, boundary.hostValues.begin() );
    }
    // The pack copy above drained the queue, so a status read here costs a
    // four-byte transfer on an idle device instead of a full drain.
    FlushRequestedResidentStatus();
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
            Real * staging = LocalPinnedStaging().Get( valueCount );
            std::copy(
                boundary.hostValues.begin(), boundary.hostValues.end(),
                staging );
            boundary.values.CopyFromHost( staging, valueCount );
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
        HipDrain( "boundary unpack kernel synchronize" );
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

bool HipFluxBackend::PackGradientFacesToHost(
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
        // No resident gradient for this zone yet: the caller falls back to the
        // host mirror. Happens on the first step and for a zone that this rank
        // has not advanced yet (multiple zones per rank).
        return false;
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
    return true;
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
    if ( storage.cachedGradientPackIdsHost != ghostCellIds
         || storage.cachedGradientPackIdCount != nFaces )
    {
        storage.gradientPackIds.CopyFromHost( ghostCellIds, nFaces );
        storage.cachedGradientPackIdsHost = ghostCellIds;
        storage.cachedGradientPackIdCount = nFaces;
    }
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
    const void * cacheKey,
    Real entropy1,
    Real entropy2,
    bool entropyFix )
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
        entropy1, entropy2, entropyFix ? 1 : 0,
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
