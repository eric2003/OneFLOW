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

#pragma once

#include "HXTypeBasic.h"

#include <cstdint>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

enum class FieldLayout
{
    EquationMajor,
    EntityMajor
};

enum class FieldRepresentation
{
    Primitive,
    Conserved,
    Residual
};

enum class FaceAreaPolicy
{
    BackendMultiplies,
    CallerMultiplies
};

enum class MemorySpace
{
    Host,
    Device
};

enum class ReconstructionLimiterMode
{
    Disabled,
    // The caller supplies the per-cell factor (host-computed field pointers).
    Cell,
    // The backend computes the factor on the device from its resident
    // gradient; the host pointers are then not required. Used for Barth,
    // because the host gradient mirror is stale under a device-authoritative
    // state and a host-side factor would be built from old values.
    DeviceCell
};

// The solver translates its BC types into these operations before dispatch.
// Preserve keeps the reconstructed values (INTERFACE/PERIODIC), Average uses
// the left/right cell average (ordinary boundaries), and SolidOverride copies
// the face-indexed bcQ value to both sides.
enum class ReconstructionBoundaryOperation : unsigned char
{
    Preserve,
    Average,
    SolidOverride
};

// Ghost ownership is separate from face reconstruction. ExternalHalo means
// interface/periodic data is supplied by an exchange path; a local boundary
// kernel must not overwrite that ghost. LocalBoundaryModel means the
// solver-owned boundary model computes the ghost and face boundary state.
enum class BoundaryGhostUpdatePolicy : unsigned char
{
    ExternalHalo,
    LocalBoundaryModel
};

enum class ReconstructionPhysicalityPolicy
{
    None,
    PositiveDensityPressure
};

// Solver-side capability metadata. Ownership remains with the solver; halo
// exchange remains outside the numerical kernel.
struct SolverDomainCapabilities
{
    int nEquations = 0;
    int nGhostCells = 0;
    int nHaloLayers = 0;
    bool hasFaceGeometry = false;
    bool hasFaceConnectivity = false;
    bool supportsStateUpload = false;
    bool supportsResidualAdd = false;

    bool SupportsEulerState() const
    {
        return nEquations == 3 || nEquations == 5;
    }
};

struct SolverFieldView
{
    int nEntities = 0;
    int nComponents = 0;
    Real * values = nullptr;
    FieldLayout layout = FieldLayout::EquationMajor;
    FieldRepresentation representation = FieldRepresentation::Conserved;
};

struct SolverConstFieldView
{
    int nEntities = 0;
    int nComponents = 0;
    const Real * values = nullptr;
    FieldLayout layout = FieldLayout::EquationMajor;
    FieldRepresentation representation = FieldRepresentation::Conserved;
};

struct FaceGeometryView
{
    int nFaces = 0;
    const Real * xNormal = nullptr;
    const Real * yNormal = nullptr;
    const Real * zNormal = nullptr;
    const Real * meshVelocityNormal = nullptr;
    const Real * faceArea = nullptr;
    FaceAreaPolicy areaPolicy = FaceAreaPolicy::BackendMultiplies;
};

// These views are deliberately backend-neutral. The owning solver remains
// responsible for lifetime and layout; a future HIP/CUDA/Kokkos adapter only
// receives pointers and extents, not MRField internals.
//
// Data layout convention (multi-equation):
//   Equation-major: data[eq * nFaces + face]  (or data[eq * nCells + cell])
//   This matches the main solver's MRField and the port's CI(c,i,nx) convention.
//
// For nEquations == 1 (scalar convection), qLeft[face] is the scalar at that face.
// For nEquations >= 3 (Euler/NS), qLeft stores conserved variables
//   [rho, rho*u, rho*v, rho*w, rho*E] for each face, equation-major.
struct FaceStateView
{
    int nFaces = 0;
    int nEquations = 0;
    const Real * qLeft = nullptr;
    const Real * qRight = nullptr;
    const Real * xNormal = nullptr;
    const Real * yNormal = nullptr;
    const Real * zNormal = nullptr;
    const Real * meshVelocityNormal = nullptr;
    const Real * faceArea = nullptr;
    Real gamma = 1.4;  // ratio of specific heats (used when nEquations >= 3)
    FieldLayout layout = FieldLayout::EquationMajor;
    FieldRepresentation representation = FieldRepresentation::Conserved;
    FaceAreaPolicy areaPolicy = FaceAreaPolicy::BackendMultiplies;
};

inline void ValidateFaceStateView( const FaceStateView & state )
{
    if ( state.nFaces <= 0 || state.nEquations <= 0
         || state.qLeft == nullptr || state.qRight == nullptr
         || state.gamma <= 1.0 )
    {
        throw std::invalid_argument( "invalid solver face state view" );
    }
    if ( state.layout != FieldLayout::EquationMajor
         || state.representation != FieldRepresentation::Conserved
         || state.areaPolicy != FaceAreaPolicy::BackendMultiplies )
    {
        throw std::invalid_argument(
            "unsupported solver face state contract" );
    }
}

struct FaceFluxView
{
    int nFaces = 0;
    int nEquations = 0;
    Real * values = nullptr;
};

struct FaceConnectivityView
{
    int nFaces = 0;
    int nBoundaryFaces = 0;
    const int * leftCell = nullptr;
    const int * rightCell = nullptr;
    // Optional explicit boundary mask. Without it, boundary faces occupy
    // [0, nBoundaryFaces) for backward compatibility.
    const unsigned char * boundaryMask = nullptr;
};

inline void ValidateFaceConnectivityView(
    const FaceConnectivityView & connectivity )
{
    if ( connectivity.nFaces <= 0
         || connectivity.nBoundaryFaces < 0
         || connectivity.nBoundaryFaces > connectivity.nFaces
         || connectivity.leftCell == nullptr
         || connectivity.rightCell == nullptr )
    {
        throw std::invalid_argument(
            "invalid solver face connectivity view" );
    }
}

struct CellGradientView
{
    int nCells = 0;
    int nGhostCells = 0;
    int nFaces = 0;
    int nBoundaryFaces = 0;
    int nEquations = 0;
    const Real * q[ 5 ] = {};
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
    Real * dqdx[ 5 ] = {};
    Real * dqdy[ 5 ] = {};
    Real * dqdz[ 5 ] = {};
    const void * cacheKey = nullptr;
    // Keep the device gradient resident for production HIP reconstruction.
    // The default preserves host-oracle and diagnostic behavior.
    bool downloadToHost = true;
};

// Boundary update is separate from reconstruction. The solver provides
// host-owned primitive state and boundary metadata; an accelerator backend may
// reuse a compatible resident state, evaluate supported local boundary models,
// and copy only the ghost/interface results needed by the CPU oracle.
struct BoundaryStateView
{
    int nCells = 0;
    int nGhostCells = 0;
    int nBoundaryFaces = 0;
    int nEquations = 0;
    Real * q[ 5 ] = {};
    Real * bcQ[ 5 ] = {};
    const int * leftCell = nullptr;
    const int * rightCell = nullptr;
    const int * boundaryType = nullptr;
    const Real * xNormal = nullptr;
    const Real * yNormal = nullptr;
    const Real * zNormal = nullptr;
    const Real * meshVelocityNormal = nullptr;
    const Real * inflow = nullptr;
    Real gamma = 1.4;
    int faceOuterNormal = 1;
    const void * cacheKey = nullptr;
    // Reuse the bound device state when it is current and compatible.  The
    // backend falls back to the host primitive arrays during initialization or
    // after any host-only update.
    bool preferResidentState = false;
    // The device state owns the newest q and the host mirror is stale by
    // design.  The backend then keeps the ghost cells on the device instead of
    // copying them back to the host every RK stage.
    bool deviceAuthoritativeHostMirror = false;
};

inline void ValidateBoundaryStateView( const BoundaryStateView & view )
{
    if ( view.nCells <= 0 || view.nGhostCells != view.nBoundaryFaces
         || view.nBoundaryFaces < 0 || view.nEquations != 5
         || view.leftCell == nullptr || view.rightCell == nullptr
         || view.boundaryType == nullptr || view.xNormal == nullptr
         || view.yNormal == nullptr || view.zNormal == nullptr
         || view.meshVelocityNormal == nullptr || view.inflow == nullptr
         || view.gamma <= 1.0 || view.cacheKey == nullptr )
    {
        throw std::invalid_argument( "invalid solver boundary state view" );
    }
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        if ( view.q[ equation ] == nullptr )
        {
            throw std::invalid_argument(
                "boundary state view has a missing primitive component" );
        }
    }
}

inline void ValidateCellGradientView( const CellGradientView & view )
{
    if ( view.nCells <= 0 || view.nGhostCells < 0
         || view.nGhostCells != view.nBoundaryFaces || view.nFaces <= 0
         || view.nBoundaryFaces < 0 || view.nBoundaryFaces > view.nFaces
         || ( view.nEquations != 3 && view.nEquations != 5 )
         || view.xFace == nullptr || view.yFace == nullptr
         || view.zFace == nullptr || view.xNormal == nullptr
         || view.yNormal == nullptr || view.zNormal == nullptr
         || view.faceArea == nullptr || view.xCell == nullptr
         || view.yCell == nullptr || view.zCell == nullptr
         || view.cellVolume == nullptr || view.leftCell == nullptr
         || view.rightCell == nullptr )
    {
        throw std::invalid_argument( "invalid solver cell gradient view" );
    }
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        if ( view.q[ equation ] == nullptr || view.dqdx[ equation ] == nullptr
             || view.dqdy[ equation ] == nullptr
             || view.dqdz[ equation ] == nullptr )
        {
            throw std::invalid_argument(
                "cell gradient view has a missing equation component" );
        }
    }
}

// Local inviscid timestep contract for the explicit 3D Euler/RK path.  The
// backend may reuse resident primitive state and cached face geometry, compute
// the face spectral radius and cell timestep, and keep the timestep resident
// for the following residual scaling stage.
struct CellTimeStepView
{
    int nCells = 0;
    int nGhostCells = 0;
    int nFaces = 0;
    int nBoundaryFaces = 0;
    int nEquations = 0;
    const Real * q[ 5 ] = {};
    const Real * xNormal = nullptr;
    const Real * yNormal = nullptr;
    const Real * zNormal = nullptr;
    const Real * meshVelocityNormal = nullptr;
    const Real * faceArea = nullptr;
    const Real * cellVolume = nullptr;
    const int * leftCell = nullptr;
    const int * rightCell = nullptr;
    Real * inviscidSpectrum = nullptr;
    Real * timeStep = nullptr;
    Real gamma = 1.4;
    Real cfl = 0.0;
    Real maxTimeRatio = 0.0;
    FieldLayout layout = FieldLayout::EquationMajor;
    FieldRepresentation representation = FieldRepresentation::Primitive;
    const void * cacheKey = nullptr;
    bool preferResidentState = false;
    bool downloadToHost = true;
};

inline void ValidateCellTimeStepView( const CellTimeStepView & view )
{
    if ( view.nCells <= 0 || view.nGhostCells < 0
         || view.nGhostCells != view.nBoundaryFaces || view.nFaces <= 0
         || view.nBoundaryFaces < 0 || view.nBoundaryFaces > view.nFaces
         || view.nEquations != 5 || view.xNormal == nullptr
         || view.yNormal == nullptr || view.zNormal == nullptr
         || view.meshVelocityNormal == nullptr || view.faceArea == nullptr
         || view.cellVolume == nullptr || view.leftCell == nullptr
         || view.rightCell == nullptr || view.timeStep == nullptr
         || view.gamma <= 1.0 || view.cfl <= 0.0
         || view.layout != FieldLayout::EquationMajor
         || view.representation != FieldRepresentation::Primitive
         || view.cacheKey == nullptr )
    {
        throw std::invalid_argument(
            "invalid solver cell timestep view" );
    }
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        if ( view.q[ equation ] == nullptr )
        {
            throw std::invalid_argument(
                "cell timestep view has a missing primitive component" );
        }
    }
}

// Backend-neutral reconstruction contract. Pointers may refer to host or
// device memory as selected by memorySpace; ownership remains with the
// per-solver/zone/grid state. All component arrays are equation-major by
// construction: component[equation][entity].
struct CellFaceReconstructionView
{
    int nCells = 0;
    int nGhostCells = 0;
    int nFaces = 0;
    int nBoundaryFaces = 0;
    int nEquations = 0;
    const Real * q[ 5 ] = {};
    const Real * dqdx[ 5 ] = {};
    const Real * dqdy[ 5 ] = {};
    const Real * dqdz[ 5 ] = {};
    const Real * limiter[ 5 ] = {};
    const Real * xFace = nullptr;
    const Real * yFace = nullptr;
    const Real * zFace = nullptr;
    const Real * xCell = nullptr;
    const Real * yCell = nullptr;
    const Real * zCell = nullptr;
    const int * leftCell = nullptr;
    const int * rightCell = nullptr;
    const unsigned char * boundaryMask = nullptr;
    const ReconstructionBoundaryOperation * boundaryOperation = nullptr;
    const BoundaryGhostUpdatePolicy * boundaryGhostPolicy = nullptr;
    const Real * bcQ[ 5 ] = {};
    Real * qLeft[ 5 ] = {};
    Real * qRight[ 5 ] = {};
    FieldLayout layout = FieldLayout::EquationMajor;
    FieldRepresentation representation = FieldRepresentation::Primitive;
    MemorySpace memorySpace = MemorySpace::Host;
    ReconstructionLimiterMode limiterMode =
        ReconstructionLimiterMode::Disabled;
    ReconstructionPhysicalityPolicy physicality =
        ReconstructionPhysicalityPolicy::PositiveDensityPressure;
    int densityComponent = 0;
    int pressureComponent = 4;
    bool hasSolidBoundary = false;
    // A cross-zone gradient halo is owned by the interface exchange, which
    // runs on the host for the host-staged path.  When the solver has already
    // applied that exchange to the host gradient mirror, the backend must
    // refresh the device ghost gradient before it reconstructs interface
    // faces; otherwise the ghost side of a Preserve face would use the local
    // owner gradient instead of the neighbor's.
    bool syncGradientGhostsFromHost = false;
    const void * ownerToken = nullptr;
    std::uint64_t topologyGeneration = 0;
    std::uint64_t fieldGeneration = 0;
    // Backend binding key (for example the local grid identity).
    const void * cacheKey = nullptr;
};

inline void ValidateCellFaceReconstructionView(
    const CellFaceReconstructionView & view )
{
    if ( view.nCells <= 0 || view.nGhostCells < 0
         || view.nGhostCells != view.nBoundaryFaces || view.nFaces <= 0
         || view.nBoundaryFaces < 0 || view.nBoundaryFaces > view.nFaces
         || ( view.nEquations != 3 && view.nEquations != 5 )
         || view.xFace == nullptr || view.yFace == nullptr
         || view.zFace == nullptr || view.xCell == nullptr
         || view.yCell == nullptr || view.zCell == nullptr
         || view.leftCell == nullptr || view.rightCell == nullptr
         || ( view.nBoundaryFaces > 0
              && ( view.boundaryOperation == nullptr
                   || view.boundaryGhostPolicy == nullptr ) )
         || view.layout != FieldLayout::EquationMajor
         || view.representation != FieldRepresentation::Primitive
         || view.ownerToken == nullptr || view.topologyGeneration == 0
         || view.fieldGeneration == 0 )
    {
        throw std::invalid_argument(
            "invalid solver cell/face reconstruction view" );
    }
    if ( view.physicality
             == ReconstructionPhysicalityPolicy::PositiveDensityPressure
         && ( view.densityComponent < 0
              || view.densityComponent >= view.nEquations
              || view.pressureComponent < 0
              || view.pressureComponent >= view.nEquations ) )
    {
        throw std::invalid_argument(
            "invalid reconstruction physicality components" );
    }
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        if ( view.q[ equation ] == nullptr
             || view.dqdx[ equation ] == nullptr
             || view.dqdy[ equation ] == nullptr
             || view.dqdz[ equation ] == nullptr
             || view.qLeft[ equation ] == nullptr
             || view.qRight[ equation ] == nullptr
             || ( view.limiterMode == ReconstructionLimiterMode::Cell
                  && view.limiter[ equation ] == nullptr )
             || ( view.hasSolidBoundary && view.bcQ[ equation ] == nullptr ) )
        {
            throw std::invalid_argument(
                "cell/face reconstruction view has a missing component" );
        }
    }

    // Face connectivity is part of the reconstruction contract, not merely
    // metadata for the later flux stage.  Boundary faces may point at their
    // corresponding ghost cells; internal faces must point at real cells.
    // Host views can be checked here.  Device views are validated by the
    // backend that owns the device topology; dereferencing device pointers
    // from this host-side contract validator would itself be invalid.
    if ( view.memorySpace == MemorySpace::Host )
    {
        int explicitBoundaryFaces = 0;
        for ( int face = 0; face < view.nFaces; ++ face )
        {
            const bool boundary = view.boundaryMask != nullptr
                ? view.boundaryMask[ face ] != 0
                : face < view.nBoundaryFaces;
            if ( boundary ) ++ explicitBoundaryFaces;

            if ( view.leftCell[ face ] < 0
                 || view.leftCell[ face ] >= view.nCells )
            {
                throw std::invalid_argument(
                    "reconstruction left-cell index is out of range" );
            }

            const int rightCell = view.rightCell[ face ];
            if ( boundary )
            {
                if ( rightCell < view.nCells
                     || rightCell >= view.nCells + view.nGhostCells )
                {
                    throw std::invalid_argument(
                        "reconstruction boundary right-cell is not a ghost" );
                }

                const BoundaryGhostUpdatePolicy ghostPolicy =
                    view.boundaryGhostPolicy[ face ];
                const ReconstructionBoundaryOperation operation =
                    view.boundaryOperation[ face ];
                if ( ghostPolicy == BoundaryGhostUpdatePolicy::ExternalHalo
                     && operation != ReconstructionBoundaryOperation::Preserve )
                {
                    throw std::invalid_argument(
                        "external-halo boundary must preserve reconstruction" );
                }
                if ( ghostPolicy
                         == BoundaryGhostUpdatePolicy::LocalBoundaryModel
                     && operation == ReconstructionBoundaryOperation::Preserve )
                {
                    throw std::invalid_argument(
                        "local boundary model cannot preserve reconstruction" );
                }
            }
            else if ( rightCell < 0 || rightCell >= view.nCells )
            {
                throw std::invalid_argument(
                    "reconstruction internal right-cell index is out of range" );
            }
        }

        if ( view.boundaryMask != nullptr
             && explicitBoundaryFaces != view.nBoundaryFaces )
        {
            throw std::invalid_argument(
                "reconstruction boundary mask disagrees with boundary count" );
        }
    }
}

struct ResidualView
{
    int nCells = 0;
    int nEquations = 0;
    Real * values = nullptr;
    Real * components[ 5 ] = {};
};

// Where the accumulating residual starts before the face scatter.  The fused
// flux add reads a base, then each face adds its contribution.
//   HostUpload  - the host residual field is the base (legacy host RK, host
//                 fallback, multigrid defect, any path that still runs
//                 LOAD_RESIDUALS).
//   DeviceZero  - the host field was not produced for this stage; the device
//                 residual is zeroed instead of uploading it.  Only valid
//                 when the caller owns the same base semantics on the device.
enum class ResidualBase
{
    HostUpload,
    DeviceZero
};

// Contract for the explicit Euler/RK state-update slice.  The HIP backend
// consumes the residual left resident by the fused flux path, scales it with
// the host-computed local timestep and cell volume, updates its resident
// primitive state, and optionally mirrors the internal-cell primitive values
// back to the solver-owned MRField.
struct CellStateUpdateView
{
    int nCells = 0;
    int nGhostCells = 0;
    int nEquations = 0;
    const Real * timeStep = nullptr;
    const Real * cellVolume = nullptr;
    Real * primitive[ 5 ] = {};
    Real gamma = 1.4;
    Real rkCoefficient = 0.0;
    FieldLayout layout = FieldLayout::EquationMajor;
    FieldRepresentation representation = FieldRepresentation::Primitive;
    const void * cacheKey = nullptr;
    // Keep the updated device state available for the following local boundary
    // update and RK stage.  Host mirroring remains controlled by the backend
    // call so legacy task consumers can stay authoritative during migration.
    bool retainResidentState = false;
    // Reuse the device timestep produced by the preceding timestep stage.
    // The default preserves the host-computed/uploaded path.
    bool preferResidentTimeStep = false;
    // The backend may keep the updated field on the device without mirroring
    // it back to the host.  Requires retainResidentState; host readers are
    // responsible for pulling the cells they need.
    bool deviceAuthoritative = false;
};

inline void ValidateCellStateUpdateView(
    const CellStateUpdateView & view )
{
    if ( view.nCells <= 0 || view.nGhostCells < 0
         || view.nEquations != 5 || view.timeStep == nullptr
         || view.cellVolume == nullptr || view.gamma <= 1.0
         || view.rkCoefficient <= 0.0
         || view.layout != FieldLayout::EquationMajor
         || view.representation != FieldRepresentation::Primitive
         || view.cacheKey == nullptr )
    {
        throw std::invalid_argument(
            "invalid solver cell state update view" );
    }
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        if ( view.primitive[ equation ] == nullptr )
        {
            throw std::invalid_argument(
                "cell state update view has a missing primitive component" );
        }
    }
    if ( view.deviceAuthoritative && ! view.retainResidentState )
    {
        throw std::invalid_argument(
            "device-authoritative state update requires a retained resident "
            "state" );
    }
}

inline void ValidateSolverFieldView( const SolverConstFieldView & field )
{
    if ( field.nEntities <= 0 || field.nComponents <= 0
         || field.values == nullptr )
    {
        throw std::invalid_argument( "invalid solver field view" );
    }
}

// Viscous and heat-conduction slice of the mean-flow residual.
//
// This is the backend-neutral contract for what the host UNsVisFlux chain
// computes per RK stage: a temperature gradient, then per-face stress, energy
// and conduction fluxes.  It exists as a view rather than as methods on a
// concrete backend because the production HIP path currently bypasses every
// abstraction and calls HipFluxBackend::Shared() directly; adding viscous flux
// in that style would hard-wire CUDA and Kokkos out of the feature entirely.
//
// Fields deliberately mirror the host inputs one for one, including the ones
// the device does not have yet: the resident temperature field and its
// gradient are new storage the backend must allocate, because
// CellGradientView is a fixed five-slot view owned by a single bound state.
//
// Everything here is plain pointers and scalars.  No HIP/CUDA/Kokkos type
// appears, so an implementation can be written for any execution space.
struct ViscousFluxView
{
    int nCells = 0;
    int nGhostCells = 0;
    int nFaces = 0;
    int nBoundaryFaces = 0;
    int nEquations = 0;

    // Geometry.  The device already carries equivalents; the view keeps them
    // so an implementation is not required to reach into backend storage.
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

    // Primitive state and its gradient, as consumed by the stress and energy
    // flux terms.  Equation-major, length totalCells = nCells + nGhostCells.
    const Real * q[ 5 ] = {};
    const Real * dqdx[ 5 ] = {};
    const Real * dqdy[ 5 ] = {};
    const Real * dqdz[ 5 ] = {};

    // Temperature and its gradient.  The device does not build either today:
    // UpdatePrimitiveStateKernel derives pressure from total energy and never
    // constructs tempr, so a viscous backend has to own this storage.
    const Real * temperature = nullptr;
    const Real * dTdx = nullptr;
    const Real * dTdy = nullptr;
    const Real * dTdz = nullptr;

    // Transport coefficients.  laminarViscosity is Sutherland on the device
    // call site; turbulentViscosity arrives from the SA solve when a
    // turbulence model is active.
    const Real * laminarViscosity = nullptr;
    const Real * turbulentViscosity = nullptr;

    // Thermodynamic and non-dimensional constants the host chain reads from
    // ambient globals.  Passing them explicitly keeps the view self-contained
    // and is the shape the E1 de-globalisation work wants.
    Real gamma = 1.4;
    Real gasConstant = 0.0;
    Real prandtlLaminar = 0.72;
    Real prandtlTurbulent = 0.9;
    Real specificHeatCp = 0.0;
    Real inverseReynolds = 0.0;
    Real wallTemperature = 0.0;

    // Face tangents for the skew correction.  Optional: an implementation may
    // reconstruct them from the normal and cell centres, but supplying them
    // matches the host path exactly.
    const Real * t1x = nullptr;
    const Real * t1y = nullptr;
    const Real * t1z = nullptr;
    const Real * t2x = nullptr;
    const Real * t2y = nullptr;
    const Real * t2z = nullptr;

    // Output: per-face viscous flux in the same equation-major face layout the
    // inviscid fused-flux path already produces, so the existing add-face-flux
    // scatter can consume it unchanged.
    Real * faceFlux = nullptr;
    int faceFluxStride = 0;

    // Optional wall heat flux, one value per solid boundary face.
    Real * wallHeatFlux = nullptr;
    int nWallFaces = 0;

    const void * cacheKey = nullptr;
};

// A viscous configuration is only meaningful when every field the host chain
// reads is present.  Rejecting early keeps a mis-wired backend from producing a
// residual that silently omits a term, which is exactly the failure mode that
// makes an inviscid-only device result diverge from the host oracle.
inline void ValidateViscousFluxView( const ViscousFluxView & view )
{
    if ( view.nCells <= 0 || view.nFaces <= 0 || view.nEquations <= 0 )
    {
        throw std::invalid_argument(
            "viscous flux view requires cells, faces and equations" );
    }
    if ( view.q[ 0 ] == nullptr || view.dqdx[ 0 ] == nullptr
         || view.dqdy[ 0 ] == nullptr || view.dqdz[ 0 ] == nullptr )
    {
        throw std::invalid_argument(
            "viscous flux view requires primitive state and its gradient" );
    }
    if ( view.temperature == nullptr || view.dTdx == nullptr
         || view.dTdy == nullptr || view.dTdz == nullptr )
    {
        throw std::invalid_argument(
            "viscous flux view requires temperature and its gradient" );
    }
    if ( view.laminarViscosity == nullptr )
    {
        throw std::invalid_argument(
            "viscous flux view requires a laminar viscosity field" );
    }
    if ( view.leftCell == nullptr || view.rightCell == nullptr
         || view.faceArea == nullptr || view.xNormal == nullptr )
    {
        throw std::invalid_argument(
            "viscous flux view requires face geometry and connectivity" );
    }
    if ( view.faceFlux == nullptr || view.faceFluxStride < view.nEquations )
    {
        throw std::invalid_argument(
            "viscous flux view requires an output buffer of at least nEquations "
            "stride" );
    }
}

inline void ValidateFaceGeometryView( const FaceGeometryView & geometry )
{
    if ( geometry.nFaces <= 0 || geometry.faceArea == nullptr )
    {
        throw std::invalid_argument( "invalid solver face geometry view" );
    }
}

EndNameSpace
