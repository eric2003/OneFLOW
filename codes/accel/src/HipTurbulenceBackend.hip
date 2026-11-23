/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
\*---------------------------------------------------------------------------*/

// HIP implementation of the Spalart-Allmaras turbulence seam.
//
// The kernels are thin: every formula they evaluate comes from
// TurbSaArithmetic.h, which is also what a host loop compiles, so the device
// and the host cannot drift into two transcriptions.  This file contributes
// only the launch geometry, the bounds checks, and the two array layouts the
// host chain uses (cell-major scalar planes and face-major scalar planes).
//
// Deliberately a separate translation unit and a separate backend class from
// HipFluxBackend: SA is a one-equation scalar transport with a source, not a
// sixth mean-flow equation, and the five-equation kernels index fixed
// five-element arrays by equation.

#include "HipTurbulenceBackend.h"

#include "AccelRuntime.h"
#include "DeviceBuffer.h"
#include "TurbSaArithmetic.h"

#include <hip/hip_runtime.h>

#include <memory>
#include <stdexcept>
#include <string>

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

int LaunchBlockSize()
{
    return 256;
}

// ---------------------------------------------------------------------------
// Eddy viscosity: TurbSrcFlux::CalcCellVist1Equ:282.
// ---------------------------------------------------------------------------
__global__ void SaEddyViscosityKernel(
    const Real * workingVariable,
    const Real * density,
    const Real * laminarViscosity,
    Real * turbulentViscosity,
    int nCells,
    SaModelConstants constants )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;

    TurbSaArithmetic::EddyViscosityInput input;
    input.density = density[ cell ];
    input.workingVariable = workingVariable[ cell ];
    input.laminarViscosity = laminarViscosity[ cell ];
    input.constants = constants;
    const TurbSaArithmetic::EddyViscosityResult result =
        TurbSaArithmetic::ComputeEddyViscosity( input );
    turbulentViscosity[ cell ] = result.eddyViscosity;
}

// ---------------------------------------------------------------------------
// Source term: TurbCom::CalcSrcSa:747 + UTurbSrcFlux::Update1Equ:242.
//
// The host Update1Equ accumulates res_sa*cvol and spec_sa*cvol into the
// existing fields, so the default path here is `+=`.  `overwrite` exists for a
// device-authoritative stage that owns the accumulation base, exactly like
// ResidualBase::DeviceZero on the mean-flow side.
// ---------------------------------------------------------------------------
__global__ void SaSourceKernel(
    const Real * workingVariable,
    const Real * workingVariableGradientX,
    const Real * workingVariableGradientY,
    const Real * workingVariableGradientZ,
    const Real * density,
    const Real * gradientDudx, const Real * gradientDudy,
    const Real * gradientDudz, const Real * gradientDvdx,
    const Real * gradientDvdy, const Real * gradientDvdz,
    const Real * gradientDwdx, const Real * gradientDwdy,
    const Real * gradientDwdz,
    const Real * wallDistance,
    const Real * lengthScale,
    const Real * laminarViscosity,
    const Real * cellVolume,
    Real * residual,
    Real * implicitSpectrum,
    int nCells,
    SaModelConstants constants,
    int overwrite )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;

    TurbSaArithmetic::SourceInput input;
    input.density = density[ cell ];
    input.workingVariable = workingVariable[ cell ];
    input.laminarViscosity = laminarViscosity[ cell ];
    // Host order is dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
    // see UTurbSrcFlux::PrepareCellValue1Equ:314-322.
    input.velocityGradient[ 0 ] = gradientDudx[ cell ];
    input.velocityGradient[ 1 ] = gradientDudy[ cell ];
    input.velocityGradient[ 2 ] = gradientDudz[ cell ];
    input.velocityGradient[ 3 ] = gradientDvdx[ cell ];
    input.velocityGradient[ 4 ] = gradientDvdy[ cell ];
    input.velocityGradient[ 5 ] = gradientDvdz[ cell ];
    input.velocityGradient[ 6 ] = gradientDwdx[ cell ];
    input.velocityGradient[ 7 ] = gradientDwdy[ cell ];
    input.velocityGradient[ 8 ] = gradientDwdz[ cell ];
    input.workingVariableGradient[ 0 ] = workingVariableGradientX[ cell ];
    input.workingVariableGradient[ 1 ] = workingVariableGradientY[ cell ];
    input.workingVariableGradient[ 2 ] = workingVariableGradientZ[ cell ];
    // CalcLengthScaleOfWallDist copies the wall distance into len_scale for
    // the plain model, so the wall distance is the fallback.
    input.lengthScale = lengthScale != nullptr
        ? lengthScale[ cell ] : wallDistance[ cell ];
    input.cellVolume = cellVolume[ cell ];
    input.constants = constants;

    const TurbSaArithmetic::SourceResult result =
        TurbSaArithmetic::ComputeSource( input );

    const Real scale = cellVolume[ cell ];
    if ( overwrite )
    {
        residual[ cell ] = result.residual * scale;
        implicitSpectrum[ cell ] = result.implicitSpectrum * scale;
    }
    else
    {
        residual[ cell ] += result.residual * scale;
        implicitSpectrum[ cell ] += result.implicitSpectrum * scale;
    }
}

// ---------------------------------------------------------------------------
// Gradient: Green-Gauss, CalcGradGGCellWeight (UGrad.cpp:88).
//
// The host makes three sequential passes: accumulate per face, scale by the
// owner volume, then copy each boundary ghost from its owner.  A per-face
// device launch would have to atomicAdd into the two cells it touches, and the
// summation order would no longer be the host's.  Because the residual gate is
// compared at 1e-17, the summation order is part of the contract, so this
// thread owns one cell and walks the face array in the host's face order.  That
// is O(nCells * nFaces) per stage; it is a correctness-first port, and the
// cell-to-face CSR the mean-flow gradient already builds is the obvious next
// optimization once the arithmetic is pinned.
// ---------------------------------------------------------------------------
__global__ void SaGradientKernel(
    const Real * workingVariable,
    const Real * xFace,
    const Real * yFace,
    const Real * zFace,
    const Real * xNormal,
    const Real * yNormal,
    const Real * zNormal,
    const Real * faceArea,
    const Real * xCell,
    const Real * yCell,
    const Real * zCell,
    const Real * cellVolume,
    const int * leftCell,
    const int * rightCell,
    const int * faceOffsets,
    const int * faceIds,
    int nCells,
    int nBoundaryFaces,
    Real * gradientX,
    Real * gradientY,
    Real * gradientZ )
{
    // One thread per cell walking that cell's own incident faces.  The earlier
    // version walked the whole face array from every thread, which is
    // O(nCells * nFaces) - for the 48x48x128 flagship that is 2.6e11 inner
    // iterations per stage, so it could never be used.  The cell-to-face CSR is
    // the same adjacency the mean-flow gradient already builds.
    //
    // Summation order is preserved exactly: the CSR lists a cell's faces in
    // ascending face id (EnsureCellFaceConnectivity fills it in one face-order
    // pass), which is the order the host accumulates them in UGrad.cpp:88.
    // That is what keeps this a bit-for-bit port rather than a re-summation.
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;

    Real accumulatedX = 0.0;
    Real accumulatedY = 0.0;
    Real accumulatedZ = 0.0;

    for ( int slot = faceOffsets[ cell ]; slot < faceOffsets[ cell + 1 ]; ++ slot )
    {
        const int face = faceIds[ slot ];
        const int left = leftCell[ face ];
        const int right = rightCell[ face ];
        const bool isLeft = left == cell;
        // A boundary face contributes to its owner only (UGrad.cpp:121's
        // `continue` skips the right cell), so a boundary ghost never
        // accumulates.
        if ( ! isLeft && face < nBoundaryFaces ) continue;

        // InverseDistanceFaceWeights, UGradPolicy.cpp:11.
        const Real dxl = xFace[ face ] - xCell[ left ];
        const Real dyl = yFace[ face ] - yCell[ left ];
        const Real dzl = zFace[ face ] - zCell[ left ];
        const Real dxr = xFace[ face ] - xCell[ right ];
        const Real dyr = yFace[ face ] - yCell[ right ];
        const Real dzr = zFace[ face ] - zCell[ right ];
        Real weightLeft = 0.0;
        Real weightRight = 0.0;
        TurbSaArithmetic::InverseDistanceWeights(
            dxl, dyl, dzl, dxr, dyr, dzr, weightLeft, weightRight );

        // WeightedFaceAverage, UGradPolicy.cpp:25.
        const Real value = weightLeft * workingVariable[ left ]
            + weightRight * workingVariable[ right ];

        // AreaWeightedFaceNormal, UGradPolicy.cpp:36.
        const Real snx = xNormal[ face ] * faceArea[ face ];
        const Real sny = yNormal[ face ] * faceArea[ face ];
        const Real snz = zNormal[ face ] * faceArea[ face ];

        if ( isLeft )
        {
            TurbSaArithmetic::AccumulateFaceGradient(
                snx, sny, snz, value, 1.0,
                accumulatedX, accumulatedY, accumulatedZ );
        }
        else
        {
            // Only an interior face reaches here for the right cell, so the
            // negative contribution is unconditional (UGrad.cpp:122).
            TurbSaArithmetic::AccumulateFaceGradient(
                snx, sny, snz, value, -1.0,
                accumulatedX, accumulatedY, accumulatedZ );
        }
    }

    // ScaleGradientByVolume, UGradPolicy.cpp:47.
    TurbSaArithmetic::ScaleGradientByVolume(
        cellVolume[ cell ], accumulatedX, accumulatedY, accumulatedZ );

    gradientX[ cell ] = accumulatedX;
    gradientY[ cell ] = accumulatedY;
    gradientZ[ cell ] = accumulatedZ;
}

// Boundary ghosts inherit the left cell gradient (UGrad.cpp:133).
__global__ void SaGradientBoundaryKernel(
    const int * leftCell,
    const int * rightCell,
    int nCells,
    int totalCells,
    int nBoundaryFaces,
    Real * gradientX,
    Real * gradientY,
    Real * gradientZ )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nBoundaryFaces ) return;
    const int left = leftCell[ face ];
    const int ghost = rightCell[ face ];
    if ( left < 0 || left >= nCells ) return;
    if ( ghost < nCells || ghost >= totalCells ) return;
    gradientX[ ghost ] = gradientX[ left ];
    gradientY[ ghost ] = gradientY[ left ];
    gradientZ[ ghost ] = gradientZ[ left ];
}

// ---------------------------------------------------------------------------
// Face flux: UTurbInvFlux + UTurbVisFlux.
// ---------------------------------------------------------------------------
__global__ void SaFaceFluxKernel(
    const Real * xNormal,
    const Real * yNormal,
    const Real * zNormal,
    const Real * meshVelocityNormal,
    const Real * faceArea,
    const int * leftCell,
    const int * rightCell,
    const Real * faceWorkingVariableLeft,
    const Real * faceWorkingVariableRight,
    const Real * faceDensityLeft,
    const Real * faceDensityRight,
    const Real * faceVelocityLeftU,
    const Real * faceVelocityLeftV,
    const Real * faceVelocityLeftW,
    const Real * faceVelocityRightU,
    const Real * faceVelocityRightV,
    const Real * faceVelocityRightW,
    const Real * xFace,
    const Real * yFace,
    const Real * zFace,
    const Real * xCell,
    const Real * yCell,
    const Real * zCell,
    const Real * cellWorkingVariable,
    const Real * gradientX,
    const Real * gradientY,
    const Real * gradientZ,
    const Real * cellLaminarViscosity,
    const Real * cellDensity,
    const unsigned char * boundaryMask,
    const unsigned char * interfaceMask,
    int nFaces,
    int nBoundaryFaces,
    Real densityWeight,
    Real invSigma,
    Real invReynolds,
    Real * inviscidFlux,
    Real * diffusiveFlux )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nFaces ) return;

    const bool boundary = boundaryMask != nullptr
        ? boundaryMask[ face ] != 0
        : face < nBoundaryFaces;
    const bool isInterface = interfaceMask != nullptr
        && interfaceMask[ face ] != 0;
    const int left = leftCell[ face ];
    const int right = rightCell[ face ];

    if ( inviscidFlux != nullptr )
    {
        // TurbInv::RoeFlux + UTurbInvFlux::UpdateFaceInvFlux.
        TurbSaArithmetic::InviscidFaceInput inviscid;
        inviscid.workingVariableLeft = faceWorkingVariableLeft[ face ];
        inviscid.workingVariableRight = faceWorkingVariableRight[ face ];
        inviscid.densityLeft = faceDensityLeft[ face ];
        inviscid.densityRight = faceDensityRight[ face ];
        inviscid.normal[ 0 ] = xNormal[ face ];
        inviscid.normal[ 1 ] = yNormal[ face ];
        inviscid.normal[ 2 ] = zNormal[ face ];
        inviscid.meshVelocityNormal = meshVelocityNormal[ face ];
        inviscid.faceArea = faceArea[ face ];
        inviscid.densityWeight = densityWeight;
        inviscid.velocityLeft[ 0 ] = faceVelocityLeftU[ face ];
        inviscid.velocityLeft[ 1 ] = faceVelocityLeftV[ face ];
        inviscid.velocityLeft[ 2 ] = faceVelocityLeftW[ face ];
        inviscid.velocityRight[ 0 ] = faceVelocityRightU[ face ];
        inviscid.velocityRight[ 1 ] = faceVelocityRightV[ face ];
        inviscid.velocityRight[ 2 ] = faceVelocityRightW[ face ];
        inviscidFlux[ face ] =
            TurbSaArithmetic::ComputeInviscidFaceFlux( inviscid );
    }

    if ( diffusiveFlux == nullptr ) return;

    // UTurbVisFlux::CalcFaceVisFlux1Equ.  The cell-centred working variable and
    // the cell gradients, not the reconstructed face values.
    TurbSaArithmetic::FaceNormalGradientInput gradient;
    gradient.workingVariableLeft = cellWorkingVariable[ left ];
    gradient.workingVariableRight = cellWorkingVariable[ right ];
    gradient.gradientLeft[ 0 ] = gradientX[ left ];
    gradient.gradientLeft[ 1 ] = gradientY[ left ];
    gradient.gradientLeft[ 2 ] = gradientZ[ left ];
    gradient.gradientRight[ 0 ] = gradientX[ right ];
    gradient.gradientRight[ 1 ] = gradientY[ right ];
    gradient.gradientRight[ 2 ] = gradientZ[ right ];
    gradient.xFace = xFace[ face ];
    gradient.yFace = yFace[ face ];
    gradient.zFace = zFace[ face ];
    gradient.xCellLeft = xCell[ left ];
    gradient.yCellLeft = yCell[ left ];
    gradient.zCellLeft = zCell[ left ];
    gradient.xCellRight = xCell[ right ];
    gradient.yCellRight = yCell[ right ];
    gradient.zCellRight = zCell[ right ];
    gradient.normal[ 0 ] = xNormal[ face ];
    gradient.normal[ 1 ] = yNormal[ face ];
    gradient.normal[ 2 ] = zNormal[ face ];
    gradient.boundary = boundary;
    gradient.isInterface = isInterface;

    const TurbSaArithmetic::FaceNormalGradientResult normalGradient =
        TurbSaArithmetic::ComputeFaceNormalGradient( gradient );

    TurbSaArithmetic::DiffusiveFaceInput diffusion;
    diffusion.workingVariableLeft = cellWorkingVariable[ left ];
    diffusion.workingVariableRight = cellWorkingVariable[ right ];
    diffusion.laminarViscosityLeft = cellLaminarViscosity[ left ];
    diffusion.laminarViscosityRight = cellLaminarViscosity[ right ];
    diffusion.densityLeft = cellDensity[ left ];
    diffusion.densityRight = cellDensity[ right ];
    diffusion.normalGradient = normalGradient.normalGradient;
    diffusion.faceArea = faceArea[ face ];
    diffusion.invSigma = invSigma;
    diffusion.invReynolds = invReynolds;

    diffusiveFlux[ face ] =
        TurbSaArithmetic::ComputeDiffusiveFaceFlux( diffusion );
}

// ---------------------------------------------------------------------------
// Update: TurbUpdate::CalcFlowField1Equ:64.
// ---------------------------------------------------------------------------
__global__ void SaUpdateKernel(
    const Real * oldValue,
    const Real * delta,
    Real * workingVariable,
    int nCells,
    Real relaxationCoefficient )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;
    workingVariable[ cell ] = TurbSaArithmetic::UnderRelaxedUpdate(
        oldValue[ cell ], delta[ cell ], relaxationCoefficient );
}

// ---------------------------------------------------------------------------
// Implicit spectrum: TurbSpectrum::CalcFaceSpectrum1Equ + the unsteady
// diagonal of UTurbSpectrum::CalcUnsteadySpectrum.
//
// The host walks the face array and accumulates each face's radius into both
// of its cells and its matrix entries into the face's own slots.  Those are
// four independent scatter targets per face, so the device uses one atomic add
// per target per face.  That reorders the summation relative to the host, which
// is acceptable for these operator coefficients and unacceptable for the
// residual; the residual path keeps its host order.
// ---------------------------------------------------------------------------
__global__ void SaSpectrumDiagonalKernel(
    const Real * timestep,
    const Real * cellVolume,
    Real turbCflRatio,
    Real unsteadySp2,
    int nCells,
    Real * implicitSpectrum )
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if ( cell >= nCells ) return;
    // CalcUnsteadySpectrum, idualtime == 0 branch: sp2 / (dt * ratio) * vol.
    implicitSpectrum[ cell ] =
        ( unsteadySp2 / ( timestep[ cell ] * turbCflRatio ) ) * cellVolume[ cell ];
}

__global__ void SaFaceSpectrumKernel(
    const Real * xNormal,
    const Real * yNormal,
    const Real * zNormal,
    const Real * meshVelocityNormal,
    const Real * faceArea,
    const int * leftCell,
    const int * rightCell,
    const Real * density,
    const Real * velocityU,
    const Real * velocityV,
    const Real * velocityW,
    const Real * workingVariable,
    const Real * laminarViscosity,
    const Real * cellVolume,
    int nFaces,
    Real invReynolds,
    Real invSigma,
    Real * implicitSpectrum,
    Real * matrixLeft,
    Real * matrixRight )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nFaces ) return;

    const int left = leftCell[ face ];
    const int right = rightCell[ face ];

    TurbSaArithmetic::SaFaceSpectrumInput in;
    in.faceArea = faceArea[ face ];

    // RelativeFaceNormalVelocity, NsInvFluxPolicy.cpp, at each side.
    in.vnLeft = xNormal[ face ] * velocityU[ left ]
        + yNormal[ face ] * velocityV[ left ]
        + zNormal[ face ] * velocityW[ left ]
        - meshVelocityNormal[ face ];
    in.vnRight = xNormal[ face ] * velocityU[ right ]
        + yNormal[ face ] * velocityV[ right ]
        + zNormal[ face ] * velocityW[ right ]
        - meshVelocityNormal[ face ];

    in.densityLeft = density[ left ];
    in.densityRight = density[ right ];
    in.laminarViscosityLeft = laminarViscosity[ left ];
    in.laminarViscosityRight = laminarViscosity[ right ];
    in.workingVariableLeft = workingVariable[ left ];
    in.workingVariableRight = workingVariable[ right ];
    in.cellVolumeLeft = cellVolume[ left ];
    in.cellVolumeRight = cellVolume[ right ];
    in.invReynolds = invReynolds;
    in.invSigma = invSigma;

    const TurbSaArithmetic::SaFaceSpectrumResult result =
        TurbSaArithmetic::ComputeFaceSpectrum( in );

    // UpdateSpectrumRadius (UTurbSpectrum.cpp:147).
    atomicAdd( & implicitSpectrum[ left ], result.radius );
    atomicAdd( & implicitSpectrum[ right ], result.radius );
    atomicAdd( & matrixLeft[ face ], result.matrixLeft );
    atomicAdd( & matrixRight[ face ], result.matrixRight );
}

// ---------------------------------------------------------------------------
// Turbulence boundary conditions: UTurbBcSolver::CalcBc.
//
// One thread per boundary face.  Faces owned by an external exchange
// (interface/periodic/overset) are skipped exactly as the host's updateFlag
// does, and only the ghost working variable is written.
// ---------------------------------------------------------------------------
__global__ void SaBoundaryKernel(
    const int * leftCell,
    const int * rightCell,
    const int * boundaryType,
    const Real * xNormal,
    const Real * yNormal,
    const Real * zNormal,
    const Real * meshVelocityNormal,
    int faceOuterNormal,
    Real * workingVariable,
    const Real * density,
    const Real * pressure,
    const Real * velocityU,
    const Real * velocityV,
    const Real * velocityW,
    Real inflowWorkingVariable,
    Real inflowDensity,
    Real inflowPressure,
    Real inflowVelocityU,
    Real inflowVelocityV,
    Real inflowVelocityW,
    Real gammaRef,
    Real gamma,
    int nBoundaryFaces )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face >= nBoundaryFaces ) return;

    TurbSaArithmetic::SaBoundaryState in;
    in.boundaryType = boundaryType[ face ];
    const int left = leftCell[ face ];
    in.interior = workingVariable[ left ];
    in.inflow = inflowWorkingVariable;

    in.density = density[ left ];
    in.pressure = pressure[ left ];
    in.velocity[ 0 ] = velocityU[ left ];
    in.velocity[ 1 ] = velocityV[ left ];
    in.velocity[ 2 ] = velocityW[ left ];
    in.normal[ 0 ] = xNormal[ face ];
    in.normal[ 1 ] = yNormal[ face ];
    in.normal[ 2 ] = zNormal[ face ];
    in.meshVelocityNormal = meshVelocityNormal[ face ];
    in.faceOuterNormal = static_cast< Real >( faceOuterNormal );

    in.inflowDensity = inflowDensity;
    in.inflowPressure = inflowPressure;
    in.inflowVelocity[ 0 ] = inflowVelocityU;
    in.inflowVelocity[ 1 ] = inflowVelocityV;
    in.inflowVelocity[ 2 ] = inflowVelocityW;
    in.gammaRef = gammaRef;
    in.gamma = gamma;

    workingVariable[ rightCell[ face ] ] =
        TurbSaArithmetic::ComputeBoundaryWorkingVariable( in );
}

// ---------------------------------------------------------------------------
// LU-SGS scalar sweeps: UTurbLusgs::LowerSweep / UpperSweep.
//
// One launch per dependency level.  Every face a level's cells touch points at
// a neighbour on a strictly lower level (in the sweep direction), so the value
// read here is the fully updated one the host would have seen - the sequential
// semantics are preserved, not approximated.  Inside a level the cells are
// mutually independent, so the level's cells fill the launch.
// ---------------------------------------------------------------------------
__global__ void SaSweepLevelKernel(
    const int * cellFaceOffsets,
    const int * cellFaceIds,
    const int * levelOffset,
    const int * levelCells,
    int level,
    const int * leftCell,
    const int * rightCell,
    const Real * rhs,
    const Real * implicitSpectrum,
    const Real * matrixLeft,
    const Real * matrixRight,
    Real * delta,
    int lowerSweep )
{
    const int begin = levelOffset[ level ];
    const int count = levelOffset[ level + 1 ] - begin;
    const int slotInLevel = blockIdx.x * blockDim.x + threadIdx.x;
    if ( slotInLevel >= count ) return;

    const int cell = levelCells[ begin + slotInLevel ];

    // ZeroFluxIncrement + Solve*Cell + AddFluxIncrement: the host accumulates
    // matrix*dqj over the cell's faces into rhs0 before the change formula.
    Real rhsIncrement = 0.0;

    for ( int slot = cellFaceOffsets[ cell ];
          slot < cellFaceOffsets[ cell + 1 ]; ++ slot )
    {
        const int face = cellFaceIds[ slot ];
        const int left = leftCell[ face ];
        const bool isLeft = left == cell;

        // UTurbLusgs::CanNotLowerSolve / CanNotUpperSolve: a face is part of the
        // sweep's triangular half only when the cell is on the expected side.
        //   lower: skip when the *other* cell has the larger id
        //   upper: skip when the *other* cell has the smaller id
        const int neighbour = isLeft ? rightCell[ face ] : left;
        if ( lowerSweep != 0 )
        {
            if ( neighbour > cell ) continue;
        }
        else
        {
            if ( neighbour < cell ) continue;
        }

        // UTurbLusgs::PrepareData picks the matrix by the face's orientation:
        // matrix_l when the cell is the face's right cell (the host reverses
        // lc/rc and sets swapflag), matrix_r when it is the left cell.
        //
        // The matrix already carries the face area (TurbSpectrum folds `farea`
        // into every contribution), so the raw neighbour delta times the matrix
        // is the whole term - TurbLusgs::GetStandardFluxIncrement:124 does
        // exactly this.  Face geometry and mean-flow state are therefore not
        // kernel inputs at all.
        const Real matrix = isLeft ? matrixRight[ face ] : matrixLeft[ face ];
        rhsIncrement += TurbSaArithmetic::LusgsMatrixTimesDelta(
            matrix, delta[ neighbour ] );
    }

    const Real radius = implicitSpectrum[ cell ];
    if ( lowerSweep != 0 )
    {
        delta[ cell ] = TurbSaArithmetic::LusgsLowerDelta(
            rhs[ cell ], rhsIncrement, radius );
    }
    else
    {
        delta[ cell ] = TurbSaArithmetic::LusgsUpperDelta(
            delta[ cell ], rhsIncrement, radius );
    }
}

// ---------------------------------------------------------------------------
// Single-sample evaluators for the contract test.  Each writes one scalar, so
// the test can compare device and host on the exact same inputs without
// building a mesh.  This mirrors HipFluxBackend::EvaluateDeviceRoe.
// ---------------------------------------------------------------------------
__global__ void SaSourceSampleKernel(
    Real density,
    Real workingVariable,
    Real laminarViscosity,
    const Real * velocityGradient,
    const Real * workingVariableGradient,
    Real lengthScale,
    Real cellVolume,
    SaModelConstants constants,
    Real * residual,
    Real * implicitSpectrum )
{
    TurbSaArithmetic::SourceInput input;
    input.density = density;
    input.workingVariable = workingVariable;
    input.laminarViscosity = laminarViscosity;
    for ( int component = 0; component < 9; ++ component )
    {
        input.velocityGradient[ component ] = velocityGradient[ component ];
    }
    for ( int component = 0; component < 3; ++ component )
    {
        input.workingVariableGradient[ component ] =
            workingVariableGradient[ component ];
    }
    input.lengthScale = lengthScale;
    input.cellVolume = cellVolume;
    input.constants = constants;

    const TurbSaArithmetic::SourceResult result =
        TurbSaArithmetic::ComputeSource( input );
    *residual = result.residual;
    *implicitSpectrum = result.implicitSpectrum;
}

__global__ void SaEddyViscositySampleKernel(
    Real density,
    Real workingVariable,
    Real laminarViscosity,
    SaModelConstants constants,
    Real * eddyViscosity )
{
    TurbSaArithmetic::EddyViscosityInput input;
    input.density = density;
    input.workingVariable = workingVariable;
    input.laminarViscosity = laminarViscosity;
    input.constants = constants;
    *eddyViscosity = TurbSaArithmetic::ComputeEddyViscosity( input )
        .eddyViscosity;
}

__global__ void SaInviscidFaceSampleKernel(
    Real workingVariableLeft,
    Real workingVariableRight,
    Real densityLeft,
    Real densityRight,
    const Real * velocityLeft,
    const Real * velocityRight,
    const Real * normal,
    Real meshVelocityNormal,
    Real faceArea,
    Real densityWeight,
    Real * flux )
{
    TurbSaArithmetic::InviscidFaceInput input;
    input.workingVariableLeft = workingVariableLeft;
    input.workingVariableRight = workingVariableRight;
    input.densityLeft = densityLeft;
    input.densityRight = densityRight;
    for ( int component = 0; component < 3; ++ component )
    {
        input.velocityLeft[ component ] = velocityLeft[ component ];
        input.velocityRight[ component ] = velocityRight[ component ];
        input.normal[ component ] = normal[ component ];
    }
    input.meshVelocityNormal = meshVelocityNormal;
    input.faceArea = faceArea;
    input.densityWeight = densityWeight;
    *flux = TurbSaArithmetic::ComputeInviscidFaceFlux( input );
}

__global__ void SaNormalGradientSampleKernel(
    Real workingVariableLeft,
    Real workingVariableRight,
    const Real * gradientLeft,
    const Real * gradientRight,
    const Real * faceCenter,
    const Real * cellLeft,
    const Real * cellRight,
    const Real * normal,
    int boundary,
    int isInterface,
    Real * normalGradient )
{
    TurbSaArithmetic::FaceNormalGradientInput input;
    input.workingVariableLeft = workingVariableLeft;
    input.workingVariableRight = workingVariableRight;
    for ( int component = 0; component < 3; ++ component )
    {
        input.gradientLeft[ component ] = gradientLeft[ component ];
        input.gradientRight[ component ] = gradientRight[ component ];
        input.normal[ component ] = normal[ component ];
    }
    input.xFace = faceCenter[ 0 ];
    input.yFace = faceCenter[ 1 ];
    input.zFace = faceCenter[ 2 ];
    input.xCellLeft = cellLeft[ 0 ];
    input.yCellLeft = cellLeft[ 1 ];
    input.zCellLeft = cellLeft[ 2 ];
    input.xCellRight = cellRight[ 0 ];
    input.yCellRight = cellRight[ 1 ];
    input.zCellRight = cellRight[ 2 ];
    input.boundary = boundary != 0;
    input.isInterface = isInterface != 0;
    *normalGradient =
        TurbSaArithmetic::ComputeFaceNormalGradient( input ).normalGradient;
}

// Copies a small host array to a fresh device allocation for the sample
// kernels.  The sample entry points are test-only, so a plain hipMalloc /
// hipFree pair is clearer than adding them to the resident buffer set.
class SampleScratch
{
public:
    explicit SampleScratch( int count )
    {
        CheckHip( hipMalloc( &pointer_, count * sizeof( Real ) ),
                  "SA sample scratch alloc" );
    }

    ~SampleScratch()
    {
        if ( pointer_ != nullptr ) hipFree( pointer_ );
    }

    SampleScratch( const SampleScratch & ) = delete;
    SampleScratch & operator=( const SampleScratch & ) = delete;

    Real * Pointer() { return pointer_; }

private:
    Real * pointer_ = nullptr;
};

void CopyToDevice( Real * destination, const Real * source, int count )
{
    CheckHip( hipMemcpy( destination, source, count * sizeof( Real ),
                         hipMemcpyHostToDevice ),
              "SA sample copy" );
}

void CopyFromDevice( Real * destination, const Real * source, int count )
{
    CheckHip( hipMemcpy( destination, source, count * sizeof( Real ),
                         hipMemcpyDeviceToHost ),
              "SA sample copy" );
}

// Copy one solver plane onto the device for the duration of a call.
//
// Every view in this file carries host addresses: the SA fields are MRFields the
// solver owns.  A kernel handed a host pointer does not fail at launch - it
// aborts inside libgalaxyhip with "Invalid address access" the first time a
// wavefront touches it, which looks like a device-side bug rather than a
// missing copy.  That is exactly how the four SA entries failed when they were
// first wired in.
using PlaneBuffer = std::unique_ptr< DeviceBuffer< Real > >;

PlaneBuffer UploadPlane( const Real * host, int nCells )
{
    if ( host == nullptr )
    {
        throw std::invalid_argument( "SA plane pointer is null" );
    }
    PlaneBuffer buffer( new DeviceBuffer< Real >() );
    buffer->CopyFromHost( host, static_cast< std::size_t >( nCells ) );
    return buffer;
}

// Integer counterpart, for the connectivity planes (left/right cell, bc type,
// CSR, level schedule).  PlaneBuffer is Real-typed, so an int array handed to
// it silently picks the wrong overload.
using IntPlaneBuffer = std::unique_ptr< DeviceBuffer< int > >;

IntPlaneBuffer UploadIntPlane( const int * host, int count )
{
    if ( host == nullptr )
    {
        throw std::invalid_argument( "SA integer plane pointer is null" );
    }
    IntPlaneBuffer buffer( new DeviceBuffer< int >() );
    buffer->CopyFromHost( host, static_cast< std::size_t >( count ) );
    return buffer;
}

std::unique_ptr< HipTurbulenceBackend > sharedHipTurbulenceBackend;

}

// ---------------------------------------------------------------------------
// Resident-array entries.
//
// These consume host or device pointers supplied by the caller.  Ownership and
// lifetime stay with the solver, exactly as the mean-flow views do; nothing
// here allocates resident storage, because the SA fields already exist in the
// host chain and the device copy of them is a follow-up slice.
// ---------------------------------------------------------------------------

bool HipTurbulenceBackend::SupportsEddyViscosity() const
{
    return true;
}

void HipTurbulenceBackend::CalcEddyViscosity(
    const TurbSaEddyViscosityView & view )
{
    ValidateTurbSaEddyViscosityView( view );
    const int nCells = view.nCells;
    const std::size_t count = static_cast< std::size_t >( nCells );
    // The view's pointers are host addresses from the solver's MRFields.  Copy
    // the three inputs up and compute into a device output, then bring the
    // result back: a kernel cannot dereference host memory, and handing it one
    // aborts inside the driver rather than failing at launch.
    DeviceBuffer< Real > deviceWorkingVariable;
    DeviceBuffer< Real > deviceDensity;
    DeviceBuffer< Real > deviceLaminar;
    DeviceBuffer< Real > deviceOut;
    deviceWorkingVariable.CopyFromHost( view.workingVariable, count );
    deviceDensity.CopyFromHost( view.density, count );
    deviceLaminar.CopyFromHost( view.laminarViscosity, count );
    deviceOut.Resize( count );

    const int blockSize = LaunchBlockSize();
    hipLaunchKernelGGL(
        SaEddyViscosityKernel,
        dim3( ( nCells + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        deviceWorkingVariable.Data(), deviceDensity.Data(), deviceLaminar.Data(),
        deviceOut.Data(), nCells, view.constants );
    CheckHip( hipGetLastError(), "SA eddy viscosity kernel launch" );
    CheckHip( hipDeviceSynchronize(), "SA eddy viscosity kernel" );
    deviceOut.CopyToHost( view.turbulentViscosity, count );
}

bool HipTurbulenceBackend::SupportsSource() const
{
    return true;
}

void HipTurbulenceBackend::CalcSource( const TurbSaSourceView & view )
{
    ValidateTurbSaSourceView( view );
    const TurbSaStateView & state = * view.state;
    const int nCells = view.nCells;
    const std::size_t count = static_cast< std::size_t >( nCells );

    // Every pointer in the view is a host address from the solver's MRFields, so
    // each plane is copied up before the kernel sees it.  A kernel that is
    // handed host memory does not fail at launch: it aborts inside the driver
    // on first touch, which is the failure this replaced.
    auto upload = [ nCells ]( const Real * host )
    {
        auto buffer = std::make_unique< DeviceBuffer< Real > >();
        buffer->CopyFromHost( host, static_cast< std::size_t >( nCells ) );
        return buffer;
    };

    auto workingVariable = upload( state.workingVariable );
    auto gradientX = upload( state.workingVariableGradientX );
    auto gradientY = upload( state.workingVariableGradientY );
    auto gradientZ = upload( state.workingVariableGradientZ );
    auto density = upload( state.density );
    std::unique_ptr< DeviceBuffer< Real > > velocityGradient[ 9 ];
    for ( int plane = 0; plane < 9; ++ plane )
    {
        velocityGradient[ plane ] = upload( state.velocityGradient[ plane ] );
    }
    auto wallDistance = upload( state.wallDistance );
    // A null length scale means "use the wall distance", which is what
    // PrepareCellValue1Equ:364 does for the non-DES model.
    auto lengthScale = state.lengthScale != nullptr
        ? upload( state.lengthScale ) : nullptr;
    auto laminarViscosity = upload( state.laminarViscosity );
    auto cellVolume = upload( state.cellVolume );

    auto deviceResidual = std::make_unique< DeviceBuffer< Real > >();
    auto deviceSpectrum = std::make_unique< DeviceBuffer< Real > >();
    if ( view.overwrite )
    {
        deviceResidual->Resize( count );
        deviceSpectrum->Resize( count );
        CheckHip( hipMemset( deviceResidual->Data(), 0, count * sizeof( Real ) ),
                  "SA source residual zero" );
        CheckHip( hipMemset( deviceSpectrum->Data(), 0, count * sizeof( Real ) ),
                  "SA source spectrum zero" );
    }
    else
    {
        // Accumulate: the kernel adds to what is already there, so the existing
        // values have to be present on the device.  Zeroing would drop the
        // inviscid base that Update1Equ is accumulating onto.
        deviceResidual->CopyFromHost( view.residual, count );
        deviceSpectrum->CopyFromHost( view.implicitSpectrum, count );
    }

    const int blockSize = LaunchBlockSize();
    hipLaunchKernelGGL(
        SaSourceKernel,
        dim3( ( nCells + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        workingVariable->Data(),
        gradientX->Data(),
        gradientY->Data(),
        gradientZ->Data(),
        density->Data(),
        velocityGradient[ 0 ]->Data(), velocityGradient[ 1 ]->Data(),
        velocityGradient[ 2 ]->Data(), velocityGradient[ 3 ]->Data(),
        velocityGradient[ 4 ]->Data(), velocityGradient[ 5 ]->Data(),
        velocityGradient[ 6 ]->Data(), velocityGradient[ 7 ]->Data(),
        velocityGradient[ 8 ]->Data(),
        wallDistance->Data(),
        lengthScale != nullptr ? lengthScale->Data() : nullptr,
        laminarViscosity->Data(),
        cellVolume->Data(),
        deviceResidual->Data(),
        deviceSpectrum->Data(),
        nCells,
        state.constants,
        view.overwrite ? 1 : 0 );
    CheckHip( hipGetLastError(), "SA source kernel launch" );
    CheckHip( hipDeviceSynchronize(), "SA source kernel" );
    // Bring the accumulated sums back.  Without this the solver's residual and
    // spectrum fields never see the device contribution.
    deviceResidual->CopyToHost( view.residual, count );
    deviceSpectrum->CopyToHost( view.implicitSpectrum, count );
}

bool HipTurbulenceBackend::SupportsGradient() const
{
    return true;
}

// Cell-to-face CSR over owned cells.
//
// Built once per partition and reused.  The single pass over the face array
// appends each face to its left cell and, for an interior face, to its right
// cell, so every cell's list comes out in ascending face id.  That ordering is
// the whole point: the host gradient accumulates a cell's faces in face-id
// order (UGrad.cpp:88 walks `fId` ascending), so summing in this order
// reproduces the host's floating-point result instead of merely a close one.
struct HipTurbulenceBackend::CellFaceAdjacency
{
    DeviceBuffer< int > offsets;
    DeviceBuffer< int > ids;
    const int * cachedLeft = nullptr;
    const int * cachedRight = nullptr;
    int cachedCells = -1;
    int cachedFaces = -1;
    int cachedBoundaryFaces = -1;
    bool valid = false;
};

HipTurbulenceBackend::CellFaceAdjacency &
HipTurbulenceBackend::EnsureCellFaceAdjacency(
    const int * leftCell, const int * rightCell,
    int nCells, int nFaces, int nBoundaryFaces )
{
    if ( cellFaceAdjacency_ == nullptr )
    {
        cellFaceAdjacency_.reset( new CellFaceAdjacency() );
    }
    CellFaceAdjacency & adjacency = * cellFaceAdjacency_;

    if ( adjacency.valid
         && adjacency.cachedLeft == leftCell
         && adjacency.cachedRight == rightCell
         && adjacency.cachedCells == nCells
         && adjacency.cachedFaces == nFaces
         && adjacency.cachedBoundaryFaces == nBoundaryFaces )
    {
        return adjacency;
    }

    if ( nCells <= 0 || nFaces <= 0 || leftCell == nullptr
         || rightCell == nullptr )
    {
        throw std::invalid_argument(
            "SA cell-face adjacency requires a non-empty face partition" );
    }

    std::vector< int > offsets( static_cast< std::size_t >( nCells ) + 1, 0 );
    for ( int face = 0; face < nFaces; ++ face )
    {
        const int left = leftCell[ face ];
        if ( left >= 0 && left < nCells ) ++ offsets[ left + 1 ];
        // A boundary face's right cell is a ghost, so it never enters the CSR.
        if ( face >= nBoundaryFaces )
        {
            const int right = rightCell[ face ];
            if ( right >= 0 && right < nCells ) ++ offsets[ right + 1 ];
        }
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
        if ( face >= nBoundaryFaces )
        {
            const int right = rightCell[ face ];
            if ( right >= 0 && right < nCells ) ids[ cursor[ right ] ++ ] = face;
        }
    }

    adjacency.offsets.CopyFromHost(
        offsets.data(), offsets.size() );
    adjacency.ids.CopyFromHost(
        ids.data(), ids.size() );
    adjacency.cachedLeft = leftCell;
    adjacency.cachedRight = rightCell;
    adjacency.cachedCells = nCells;
    adjacency.cachedFaces = nFaces;
    adjacency.cachedBoundaryFaces = nBoundaryFaces;
    adjacency.valid = true;
    return adjacency;
}

void HipTurbulenceBackend::CalcGradient( const TurbSaGradientView & view )
{
    ValidateTurbSaGradientView( view );

    CellFaceAdjacency & adjacency = EnsureCellFaceAdjacency(
        view.leftCell, view.rightCell,
        view.nCells, view.nFaces, view.nBoundaryFaces );

    const int nCells = view.nCells;
    const int nFaces = view.nFaces;
    const int totalCells = nCells + view.nGhostCells;

    // Every view pointer above is a host address from the solver's MRFields.
    // The kernel can only read device memory, so each plane is copied up and
    // the three gradient planes are copied back: the source term and the
    // diffusive flux read them on the host immediately after this call.
    //
    // The cell-indexed planes use the total extent because the CSR reaches the
    // ghost range through rightCell for a boundary face.
    auto workingVariable = UploadPlane( view.workingVariable, totalCells );
    auto xFace = UploadPlane( view.xFace, nFaces );
    auto yFace = UploadPlane( view.yFace, nFaces );
    auto zFace = UploadPlane( view.zFace, nFaces );
    auto xNormal = UploadPlane( view.xNormal, nFaces );
    auto yNormal = UploadPlane( view.yNormal, nFaces );
    auto zNormal = UploadPlane( view.zNormal, nFaces );
    auto faceArea = UploadPlane( view.faceArea, nFaces );
    auto xCell = UploadPlane( view.xCell, totalCells );
    auto yCell = UploadPlane( view.yCell, totalCells );
    auto zCell = UploadPlane( view.zCell, totalCells );
    auto cellVolume = UploadPlane( view.cellVolume, totalCells );
    auto leftCell = UploadIntPlane( view.leftCell, nFaces );
    auto rightCell = UploadIntPlane( view.rightCell, nFaces );

    // The ghost planes are written by the boundary inheriting kernel, so the
    // device output spans the ghosts too and is copied back whole.
    DeviceBuffer< Real > gradientX(
        static_cast< std::size_t >( totalCells ) );
    DeviceBuffer< Real > gradientY(
        static_cast< std::size_t >( totalCells ) );
    DeviceBuffer< Real > gradientZ(
        static_cast< std::size_t >( totalCells ) );

    const int blockSize = LaunchBlockSize();
    hipLaunchKernelGGL(
        SaGradientKernel,
        dim3( ( nCells + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        workingVariable->Data(),
        xFace->Data(), yFace->Data(), zFace->Data(),
        xNormal->Data(), yNormal->Data(), zNormal->Data(),
        faceArea->Data(),
        xCell->Data(), yCell->Data(), zCell->Data(), cellVolume->Data(),
        leftCell->Data(), rightCell->Data(),
        adjacency.offsets.Data(), adjacency.ids.Data(),
        nCells, view.nBoundaryFaces,
        gradientX.Data(), gradientY.Data(), gradientZ.Data() );
    CheckHip( hipGetLastError(), "SA gradient kernel launch" );

    // Boundary ghosts inherit the owner gradient (UGrad.cpp:133).
    if ( view.nBoundaryFaces > 0 )
    {
        hipLaunchKernelGGL(
            SaGradientBoundaryKernel,
            dim3( ( view.nBoundaryFaces + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0,
            leftCell->Data(), rightCell->Data(), nCells,
            totalCells, view.nBoundaryFaces,
            gradientX.Data(), gradientY.Data(), gradientZ.Data() );
        CheckHip( hipGetLastError(), "SA gradient boundary kernel launch" );
    }
    CheckHip( hipDeviceSynchronize(), "SA gradient kernel" );

    const std::size_t count = static_cast< std::size_t >( totalCells );
    gradientX.CopyToHost( view.gradientX, count );
    gradientY.CopyToHost( view.gradientY, count );
    gradientZ.CopyToHost( view.gradientZ, count );
}

bool HipTurbulenceBackend::SupportsFaceFlux() const
{
    return true;
}

void HipTurbulenceBackend::CalcFaceFlux( const TurbSaFaceView & view )
{
    ValidateTurbSaFaceView( view );

    const int nFaces = view.nFaces;
    const int nCells = view.nCells + view.nGhostCells;
    if ( nCells <= 0 )
    {
        throw std::invalid_argument(
            "SA face flux requires the cell extent of its per-cell planes" );
    }

    // Face-length planes.  The reconstructed face values are only present when
    // an inviscid flux was requested; the diffusive-only caller leaves them
    // null, so upload each only when it exists.
    auto fXNormal = UploadPlane( view.xNormal, nFaces );
    auto fYNormal = UploadPlane( view.yNormal, nFaces );
    auto fZNormal = UploadPlane( view.zNormal, nFaces );
    auto fVelocityNormal = UploadPlane( view.meshVelocityNormal, nFaces );
    auto fArea = UploadPlane( view.faceArea, nFaces );
    std::unique_ptr< DeviceBuffer< Real > > fA;
    std::unique_ptr< DeviceBuffer< Real > > fB;
    std::unique_ptr< DeviceBuffer< Real > > fRhoA;
    std::unique_ptr< DeviceBuffer< Real > > fRhoB;
    std::unique_ptr< DeviceBuffer< Real > > fUa, fVa, fWa;
    std::unique_ptr< DeviceBuffer< Real > > fUb, fVb, fWb;
    if ( view.faceWorkingVariableLeft != nullptr )
    {
        fA = UploadPlane( view.faceWorkingVariableLeft, nFaces );
        fB = UploadPlane( view.faceWorkingVariableRight, nFaces );
        fRhoA = UploadPlane( view.faceDensityLeft, nFaces );
        fRhoB = UploadPlane( view.faceDensityRight, nFaces );
        fUa = UploadPlane( view.faceVelocityLeft[ 0 ], nFaces );
        fVa = UploadPlane( view.faceVelocityLeft[ 1 ], nFaces );
        fWa = UploadPlane( view.faceVelocityLeft[ 2 ], nFaces );
        fUb = UploadPlane( view.faceVelocityRight[ 0 ], nFaces );
        fVb = UploadPlane( view.faceVelocityRight[ 1 ], nFaces );
        fWb = UploadPlane( view.faceVelocityRight[ 2 ], nFaces );
    }
    auto fXFace = UploadPlane( view.xFace, nFaces );
    auto fYFace = UploadPlane( view.yFace, nFaces );
    auto fZFace = UploadPlane( view.zFace, nFaces );

    // Cell-length planes, indexed through leftCell/rightCell into the ghosts.
    auto cX = UploadPlane( view.xCell, nCells );
    auto cY = UploadPlane( view.yCell, nCells );
    auto cZ = UploadPlane( view.zCell, nCells );
    auto cWorking = UploadPlane( view.cellWorkingVariable, nCells );
    auto cGradX = UploadPlane( view.workingVariableGradientX, nCells );
    auto cGradY = UploadPlane( view.workingVariableGradientY, nCells );
    auto cGradZ = UploadPlane( view.workingVariableGradientZ, nCells );
    auto cLaminar = UploadPlane( view.cellLaminarViscosity, nCells );
    auto cDensity = UploadPlane( view.cellDensity, nCells );

    // Integer face->cell maps.
    DeviceBuffer< int > dLeft;
    DeviceBuffer< int > dRight;
    dLeft.CopyFromHost( view.leftCell, static_cast< std::size_t >( nFaces ) );
    dRight.CopyFromHost( view.rightCell, static_cast< std::size_t >( nFaces ) );

    // Optional masks.  Absent means "no mask", which the kernel already treats
    // as boundary-by-index and non-interface.
    DeviceBuffer< unsigned char > dBoundary;
    DeviceBuffer< unsigned char > dInterface;
    if ( view.boundaryMask != nullptr )
    {
        dBoundary.CopyFromHost( view.boundaryMask,
                                static_cast< std::size_t >( nFaces ) );
    }
    if ( view.interfaceMask != nullptr )
    {
        dInterface.CopyFromHost( view.interfaceMask,
                                 static_cast< std::size_t >( nFaces ) );
    }

    DeviceBuffer< Real > dInviscid;
    DeviceBuffer< Real > dDiffusive;
    if ( view.inviscidFlux != nullptr ) dInviscid.Resize( nFaces );
    if ( view.diffusiveFlux != nullptr ) dDiffusive.Resize( nFaces );

    const int blockSize = LaunchBlockSize();
    hipLaunchKernelGGL(
        SaFaceFluxKernel,
        dim3( ( nFaces + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        fXNormal->Data(), fYNormal->Data(), fZNormal->Data(),
        fVelocityNormal->Data(), fArea->Data(),
        dLeft.Data(), dRight.Data(),
        fA ? fA->Data() : nullptr, fB ? fB->Data() : nullptr,
        fRhoA ? fRhoA->Data() : nullptr, fRhoB ? fRhoB->Data() : nullptr,
        fUa ? fUa->Data() : nullptr, fVa ? fVa->Data() : nullptr,
        fWa ? fWa->Data() : nullptr,
        fUb ? fUb->Data() : nullptr, fVb ? fVb->Data() : nullptr,
        fWb ? fWb->Data() : nullptr,
        fXFace->Data(), fYFace->Data(), fZFace->Data(),
        cX->Data(), cY->Data(), cZ->Data(),
        cWorking->Data(),
        cGradX->Data(), cGradY->Data(), cGradZ->Data(),
        cLaminar->Data(), cDensity->Data(),
        view.boundaryMask != nullptr ? dBoundary.Data() : nullptr,
        view.interfaceMask != nullptr ? dInterface.Data() : nullptr,
        nFaces, view.nBoundaryFaces,
        view.densityWeight, view.invSigma, view.invReynolds,
        view.inviscidFlux != nullptr ? dInviscid.Data() : nullptr,
        view.diffusiveFlux != nullptr ? dDiffusive.Data() : nullptr );
    CheckHip( hipGetLastError(), "SA face flux kernel launch" );
    CheckHip( hipDeviceSynchronize(), "SA face flux kernel" );

    if ( view.inviscidFlux != nullptr )
    {
        dInviscid.CopyToHost( view.inviscidFlux,
                              static_cast< std::size_t >( nFaces ) );
    }
    if ( view.diffusiveFlux != nullptr )
    {
        dDiffusive.CopyToHost( view.diffusiveFlux,
                               static_cast< std::size_t >( nFaces ) );
    }
}

bool HipTurbulenceBackend::SupportsUpdate() const
{
    return true;
}

void HipTurbulenceBackend::Update( const TurbSaUpdateView & view )
{
    ValidateTurbSaUpdateView( view );
    const std::size_t count = static_cast< std::size_t >( view.nCells );

    // The view hands over host MRField planes.  The kernel can only touch device
    // memory, so oldValue and delta are copied up and the working variable is
    // copied back; the host reads it immediately after this call.
    auto oldValue = UploadPlane( view.oldValue, view.nCells );
    auto delta = UploadPlane( view.delta, view.nCells );

    DeviceBuffer< Real > workingVariable( count );

    const int blockSize = LaunchBlockSize();
    hipLaunchKernelGGL(
        SaUpdateKernel,
        dim3( ( view.nCells + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        oldValue->Data(), delta->Data(), workingVariable.Data(), view.nCells,
        view.relaxationCoefficient );
    CheckHip( hipGetLastError(), "SA update kernel launch" );
    CheckHip( hipDeviceSynchronize(), "SA update kernel" );

    workingVariable.CopyToHost( view.workingVariable, count );
}

bool HipTurbulenceBackend::SupportsSpectrum() const
{
    return true;
}

void HipTurbulenceBackend::CalcSpectrum( const TurbSaSpectrumView & view )
{
    ValidateTurbSaSpectrumView( view );

    const int totalCells = view.nCells + view.nGhostCells;
    const int nFaces = view.nFaces;

    // Per-cell planes are indexed through leftCell/rightCell, which reach into
    // the ghost range, so the cell extent is the total, not the owned count.
    auto density = UploadPlane( view.density, totalCells );
    auto velocityU = UploadPlane( view.velocity[ 0 ], totalCells );
    auto velocityV = UploadPlane( view.velocity[ 1 ], totalCells );
    auto velocityW = UploadPlane( view.velocity[ 2 ], totalCells );
    auto workingVariable = UploadPlane( view.workingVariable, totalCells );
    auto laminarViscosity = UploadPlane( view.laminarViscosity, totalCells );
    auto cellVolume = UploadPlane( view.cellVolume, totalCells );

    auto xNormal = UploadPlane( view.xNormal, nFaces );
    auto yNormal = UploadPlane( view.yNormal, nFaces );
    auto zNormal = UploadPlane( view.zNormal, nFaces );
    auto meshVelocityNormal = UploadPlane( view.meshVelocityNormal, nFaces );
    auto faceArea = UploadPlane( view.faceArea, nFaces );
    auto leftCell = UploadIntPlane( view.leftCell, nFaces );
    auto rightCell = UploadIntPlane( view.rightCell, nFaces );

    // ZeroSpectrum (UTurbSpectrum.cpp:79) before accumulating.
    DeviceBuffer< Real > spectrum( static_cast< std::size_t >( totalCells ) );
    CheckHip( hipMemset( spectrum.Data(), 0, spectrum.Size() * sizeof( Real ) ),
              "SA spectrum zero" );
    DeviceBuffer< Real > matrixLeft( static_cast< std::size_t >( nFaces ) );
    DeviceBuffer< Real > matrixRight( static_cast< std::size_t >( nFaces ) );
    CheckHip( hipMemset( matrixLeft.Data(), 0,
                         matrixLeft.Size() * sizeof( Real ) ),
              "SA spectrum matrix zero" );
    CheckHip( hipMemset( matrixRight.Data(), 0,
                         matrixRight.Size() * sizeof( Real ) ),
              "SA spectrum matrix zero" );

    const int blockSize = LaunchBlockSize();
    if ( view.unsteadyDiagonal )
    {
        PlaneBuffer timestep( new DeviceBuffer< Real >() );
        timestep->CopyFromHost( view.timestep,
                               static_cast< std::size_t >( totalCells ) );
        hipLaunchKernelGGL(
            SaSpectrumDiagonalKernel,
            dim3( ( view.nCells + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0,
            timestep->Data(), cellVolume->Data(),
            view.turbCflRatio, view.unsteadySp2, view.nCells,
            spectrum.Data() );
        CheckHip( hipGetLastError(), "SA spectrum diagonal kernel launch" );
    }

    hipLaunchKernelGGL(
        SaFaceSpectrumKernel,
        dim3( ( nFaces + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        xNormal->Data(), yNormal->Data(), zNormal->Data(),
        meshVelocityNormal->Data(), faceArea->Data(),
        leftCell->Data(), rightCell->Data(),
        density->Data(), velocityU->Data(), velocityV->Data(),
        velocityW->Data(), workingVariable->Data(),
        laminarViscosity->Data(), cellVolume->Data(),
        nFaces, view.invReynolds, view.invSigma,
        spectrum.Data(), matrixLeft.Data(), matrixRight.Data() );
    CheckHip( hipGetLastError(), "SA face spectrum kernel launch" );
    CheckHip( hipDeviceSynchronize(), "SA spectrum kernel" );

    spectrum.CopyToHost( view.implicitSpectrum,
                         static_cast< std::size_t >( totalCells ) );
    matrixLeft.CopyToHost( view.matrixLeft,
                           static_cast< std::size_t >( nFaces ) );
    matrixRight.CopyToHost( view.matrixRight,
                            static_cast< std::size_t >( nFaces ) );
}

bool HipTurbulenceBackend::SupportsBoundary() const
{
    return true;
}

void HipTurbulenceBackend::CalcBoundary( const TurbSaBoundaryView & view )
{
    ValidateTurbSaBoundaryView( view );

    const int totalCells = view.nCells + view.nGhostCells;
    const int nFaces = view.nBoundaryFaces;

    // The working variable is read and written in place.
    DeviceBuffer< Real > workingVariable(
        static_cast< std::size_t >( totalCells ) );
    workingVariable.CopyFromHost(
        view.workingVariable, static_cast< std::size_t >( totalCells ) );

    auto density = UploadPlane( view.density, totalCells );
    auto pressure = UploadPlane( view.pressure, totalCells );
    auto velocityU = UploadPlane( view.velocity[ 0 ], totalCells );
    auto velocityV = UploadPlane( view.velocity[ 1 ], totalCells );
    auto velocityW = UploadPlane( view.velocity[ 2 ], totalCells );

    auto xNormal = UploadPlane( view.xNormal, nFaces );
    auto yNormal = UploadPlane( view.yNormal, nFaces );
    auto zNormal = UploadPlane( view.zNormal, nFaces );
    auto meshVelocityNormal = UploadPlane( view.meshVelocityNormal, nFaces );
    auto leftCell = UploadIntPlane( view.leftCell, nFaces );
    auto rightCell = UploadIntPlane( view.rightCell, nFaces );
    auto boundaryType = UploadIntPlane( view.boundaryType, nFaces );

    const int blockSize = LaunchBlockSize();
    hipLaunchKernelGGL(
        SaBoundaryKernel,
        dim3( ( nFaces + blockSize - 1 ) / blockSize ),
        dim3( blockSize ), 0, 0,
        leftCell->Data(), rightCell->Data(), boundaryType->Data(),
        xNormal->Data(), yNormal->Data(), zNormal->Data(),
        meshVelocityNormal->Data(), view.faceOuterNormal,
        workingVariable.Data(), density->Data(), pressure->Data(),
        velocityU->Data(), velocityV->Data(), velocityW->Data(),
        view.inflowWorkingVariable, view.inflowDensity,
        view.inflowPressure, view.inflowVelocity[ 0 ],
        view.inflowVelocity[ 1 ], view.inflowVelocity[ 2 ],
        view.gammaRef, view.gamma, nFaces );
    CheckHip( hipGetLastError(), "SA boundary kernel launch" );
    CheckHip( hipDeviceSynchronize(), "SA boundary kernel" );

    workingVariable.CopyToHost(
        view.workingVariable, static_cast< std::size_t >( totalCells ) );
}

bool HipTurbulenceBackend::SupportsSweep() const
{
    return true;
}

// Resident sweep inputs.
//
// Only what the sweep actually reads is held: the partition (CSR, level
// schedule, face-cell map) is cached by identity, and the four per-step arrays
// are uploaded on every call.  Face geometry and mean-flow state are not here
// at all - the matrix already carries the face area, so the kernel never needs
// a normal, an area or a velocity.
struct HipTurbulenceBackend::SweepStorage
{
    DeviceBuffer< int > cellFaceOffsets;
    DeviceBuffer< int > cellFaceIds;
    DeviceBuffer< int > levelOffset;
    DeviceBuffer< int > levelCells;
    DeviceBuffer< int > leftCell;
    DeviceBuffer< int > rightCell;
    DeviceBuffer< Real > rhs;
    DeviceBuffer< Real > spectrum;
    DeviceBuffer< Real > matrixLeft;
    DeviceBuffer< Real > matrixRight;
    DeviceBuffer< Real > delta;

    const int * cachedOffsets = nullptr;
    const int * cachedIds = nullptr;
    const int * cachedLevelOffset = nullptr;
    const int * cachedLevelCells = nullptr;
    const int * cachedLeftCell = nullptr;
    const int * cachedRightCell = nullptr;
    int cachedCells = -1;
    int cachedFaces = -1;
    int cachedLevels = -1;
};

void HipTurbulenceBackend::CalcSweep( const TurbSaSweepView & view )
{
    ValidateTurbSaSweepView( view );

    if ( sweepStorage_ == nullptr )
    {
        sweepStorage_.reset( new SweepStorage() );
    }
    SweepStorage & storage = * sweepStorage_;

    const int nCells = view.nCells;
    const int nFaces = view.nFaces;
    const int nLevels = view.nLevels;
    const std::size_t cellCount = static_cast< std::size_t >( nCells );
    // delta spans the ghost range: the upper sweep multiplies a boundary
    // face's ghost delta by matrix_r, exactly as the host does.
    const std::size_t deltaCount =
        static_cast< std::size_t >( nCells + view.nGhostCells );

    auto uploadInt = []( DeviceBuffer< int > & buffer, const int * host,
                         std::size_t count, const int *& key,
                         const int * hostKey )
    {
        if ( key == hostKey && buffer.Size() >= count ) return;
        buffer.CopyFromHost( host, count );
        key = hostKey;
    };

    const bool partitionChanged =
        storage.cachedCells != nCells || storage.cachedFaces != nFaces
        || storage.cachedLevels != nLevels;

    if ( partitionChanged )
    {
        uploadInt( storage.cellFaceOffsets, view.cellFaceOffsets,
                   cellCount + 1, storage.cachedOffsets, view.cellFaceOffsets );
        uploadInt( storage.cellFaceIds, view.cellFaceIds,
                   static_cast< std::size_t >( view.cellFaceOffsets[ nCells ] ),
                   storage.cachedIds, view.cellFaceIds );
        uploadInt( storage.levelOffset, view.levelOffset,
                   static_cast< std::size_t >( nLevels ) + 1,
                   storage.cachedLevelOffset, view.levelOffset );
        uploadInt( storage.levelCells, view.levelCells,
                   static_cast< std::size_t >( view.levelOffset[ nLevels ] ),
                   storage.cachedLevelCells, view.levelCells );
        uploadInt( storage.leftCell, view.leftCell,
                   static_cast< std::size_t >( nFaces ),
                   storage.cachedLeftCell, view.leftCell );
        uploadInt( storage.rightCell, view.rightCell,
                   static_cast< std::size_t >( nFaces ),
                   storage.cachedRightCell, view.rightCell );
        storage.cachedCells = nCells;
        storage.cachedFaces = nFaces;
        storage.cachedLevels = nLevels;
    }

    // Per-call state.  The residual, spectrum and both matrices are rewritten
    // by the time loop even though their host pointers are stable, and delta is
    // read-modify-written by the sweep, so all four are uploaded every call.
    storage.rhs.CopyFromHost( view.rhs, cellCount );
    storage.spectrum.CopyFromHost( view.implicitSpectrum, cellCount );
    storage.matrixLeft.CopyFromHost(
        view.matrixLeft, static_cast< std::size_t >( nFaces ) );
    storage.matrixRight.CopyFromHost(
        view.matrixRight, static_cast< std::size_t >( nFaces ) );
    storage.delta.CopyFromHost( view.delta, deltaCount );

    const int blockSize = LaunchBlockSize();
    const int lower = view.lowerSweep ? 1 : 0;
    for ( int level = 0; level < nLevels; ++ level )
    {
        const int count = view.levelOffset[ level + 1 ]
            - view.levelOffset[ level ];
        if ( count <= 0 ) continue;
        hipLaunchKernelGGL(
            SaSweepLevelKernel,
            dim3( ( count + blockSize - 1 ) / blockSize ),
            dim3( blockSize ), 0, 0,
            storage.cellFaceOffsets.Data(), storage.cellFaceIds.Data(),
            storage.levelOffset.Data(), storage.levelCells.Data(), level,
            storage.leftCell.Data(), storage.rightCell.Data(),
            storage.rhs.Data(), storage.spectrum.Data(),
            storage.matrixLeft.Data(), storage.matrixRight.Data(),
            storage.delta.Data(), lower );
        CheckHip( hipGetLastError(), "SA sweep level kernel launch" );
    }
    CheckHip( hipDeviceSynchronize(), "SA sweep kernels" );

    storage.delta.CopyToHost( view.delta, deltaCount );
}

// ---------------------------------------------------------------------------
// Test-only single-sample entries.
// ---------------------------------------------------------------------------

void HipTurbulenceBackend::EvaluateSourceSample(
    const SourceSample & sample, Real & residual, Real & implicitSpectrum ) const
{
    SampleScratch velocityGradient( 9 );
    SampleScratch workingVariableGradient( 3 );
    Real * deviceResidual = nullptr;
    Real * deviceSpectrum = nullptr;
    CheckHip( hipMalloc( &deviceResidual, sizeof( Real ) ),
              "SA source sample alloc" );
    CheckHip( hipMalloc( &deviceSpectrum, sizeof( Real ) ),
              "SA source sample alloc" );

    CopyToDevice( velocityGradient.Pointer(), sample.velocityGradient, 9 );
    CopyToDevice( workingVariableGradient.Pointer(),
                  sample.workingVariableGradient, 3 );

    hipLaunchKernelGGL(
        SaSourceSampleKernel, dim3( 1 ), dim3( 1 ), 0, 0,
        sample.density, sample.workingVariable, sample.laminarViscosity,
        velocityGradient.Pointer(), workingVariableGradient.Pointer(),
        sample.lengthScale, sample.cellVolume, sample.constants,
        deviceResidual, deviceSpectrum );
    CheckHip( hipGetLastError(), "SA source sample kernel launch" );
    CheckHip( hipDeviceSynchronize(), "SA source sample kernel" );

    CopyFromDevice( &residual, deviceResidual, 1 );
    CopyFromDevice( &implicitSpectrum, deviceSpectrum, 1 );
    CheckHip( hipFree( deviceResidual ), "SA source sample free" );
    CheckHip( hipFree( deviceSpectrum ), "SA source sample free" );
}

void HipTurbulenceBackend::EvaluateEddyViscositySample(
    const EddyViscositySample & sample, Real & eddyViscosity ) const
{
    Real * deviceResult = nullptr;
    CheckHip( hipMalloc( &deviceResult, sizeof( Real ) ),
              "SA eddy viscosity sample alloc" );
    hipLaunchKernelGGL(
        SaEddyViscositySampleKernel, dim3( 1 ), dim3( 1 ), 0, 0,
        sample.density, sample.workingVariable, sample.laminarViscosity,
        sample.constants, deviceResult );
    CheckHip( hipGetLastError(), "SA eddy viscosity sample kernel launch" );
    CheckHip( hipDeviceSynchronize(), "SA eddy viscosity sample kernel" );
    CopyFromDevice( &eddyViscosity, deviceResult, 1 );
    CheckHip( hipFree( deviceResult ), "SA eddy viscosity sample free" );
}

void HipTurbulenceBackend::EvaluateInviscidFaceSample(
    const InviscidFaceSample & sample, Real & flux ) const
{
    SampleScratch velocityLeft( 3 );
    SampleScratch velocityRight( 3 );
    SampleScratch normal( 3 );
    Real * deviceResult = nullptr;
    CheckHip( hipMalloc( &deviceResult, sizeof( Real ) ),
              "SA inviscid sample alloc" );
    CopyToDevice( velocityLeft.Pointer(), sample.velocityLeft, 3 );
    CopyToDevice( velocityRight.Pointer(), sample.velocityRight, 3 );
    CopyToDevice( normal.Pointer(), sample.normal, 3 );
    hipLaunchKernelGGL(
        SaInviscidFaceSampleKernel, dim3( 1 ), dim3( 1 ), 0, 0,
        sample.workingVariableLeft, sample.workingVariableRight,
        sample.densityLeft, sample.densityRight,
        velocityLeft.Pointer(), velocityRight.Pointer(), normal.Pointer(),
        sample.meshVelocityNormal, sample.faceArea, sample.densityWeight,
        deviceResult );
    CheckHip( hipGetLastError(), "SA inviscid sample kernel launch" );
    CheckHip( hipDeviceSynchronize(), "SA inviscid sample kernel" );
    CopyFromDevice( &flux, deviceResult, 1 );
    CheckHip( hipFree( deviceResult ), "SA inviscid sample free" );
}

void HipTurbulenceBackend::EvaluateNormalGradientSample(
    const NormalGradientSample & sample, Real & normalGradient ) const
{
    SampleScratch gradientLeft( 3 );
    SampleScratch gradientRight( 3 );
    SampleScratch faceCenter( 3 );
    SampleScratch cellLeft( 3 );
    SampleScratch cellRight( 3 );
    SampleScratch normal( 3 );
    Real * deviceResult = nullptr;
    CheckHip( hipMalloc( &deviceResult, sizeof( Real ) ),
              "SA normal gradient sample alloc" );
    CopyToDevice( gradientLeft.Pointer(), sample.gradientLeft, 3 );
    CopyToDevice( gradientRight.Pointer(), sample.gradientRight, 3 );
    CopyToDevice( faceCenter.Pointer(), sample.faceCenter, 3 );
    CopyToDevice( cellLeft.Pointer(), sample.cellLeft, 3 );
    CopyToDevice( cellRight.Pointer(), sample.cellRight, 3 );
    CopyToDevice( normal.Pointer(), sample.normal, 3 );
    hipLaunchKernelGGL(
        SaNormalGradientSampleKernel, dim3( 1 ), dim3( 1 ), 0, 0,
        sample.workingVariableLeft, sample.workingVariableRight,
        gradientLeft.Pointer(), gradientRight.Pointer(),
        faceCenter.Pointer(), cellLeft.Pointer(), cellRight.Pointer(),
        normal.Pointer(), sample.boundary ? 1 : 0,
        sample.isInterface ? 1 : 0, deviceResult );
    CheckHip( hipGetLastError(), "SA normal gradient sample kernel launch" );
    CheckHip( hipDeviceSynchronize(), "SA normal gradient sample kernel" );
    CopyFromDevice( &normalGradient, deviceResult, 1 );
    CheckHip( hipFree( deviceResult ), "SA normal gradient sample free" );
}

HipTurbulenceBackend & HipTurbulenceBackend::Shared()
{
    if ( ! sharedHipTurbulenceBackend )
    {
        sharedHipTurbulenceBackend.reset( new HipTurbulenceBackend() );
    }
    return * sharedHipTurbulenceBackend;
}

// Out of line so both PIMPL members are complete at the point unique_ptr
// instantiates its deleter.
HipTurbulenceBackend::HipTurbulenceBackend() = default;

HipTurbulenceBackend::~HipTurbulenceBackend() = default;

void HipTurbulenceBackend::ReleaseShared()
{
    sharedHipTurbulenceBackend.reset();
}

}
