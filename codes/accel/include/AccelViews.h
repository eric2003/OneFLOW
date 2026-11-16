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

// Derived thermodynamic field for the viscous and heat-conduction port.
//
// Temperature is not a sixth conserved variable and must not be folded into
// CellGradientView, which is a fixed five-slot view sized nEquations*3 and
// owned by one bound state.  The host builds temperature from the primitive
// state once per boundary stage (NsCalcGamaT) and the viscous chain then
// differentiates it for heat conduction, so the device needs the same derived
// field plus its cell gradient as separate storage.
//
// The relation is exactly what NsCalcGamaT writes for every cell it owns,
// inner and ghost:
//
//     T[cell] = p[cell] / ( stateCoefficient * rho[cell]
//                           * molecularWeightRatio )
//
// stateCoefficient is the solver's nscom.statecoef, initialised to
// 1 / ( gama_ref * mach_ref^2 ).  Gamma therefore reaches this relation through
// the coefficient rather than as a compile-time constant: the caller must pass
// the value the solver actually initialised, because gama_ref is an input
// field.  molecularWeightRatio is the host routine's `oamw` (the local `one`
// for the non-chemistry path).
struct CellTemperatureView
{
    int nCells = 0;
    int nGhostCells = 0;
    int nFaces = 0;
    int nBoundaryFaces = 0;
    int nEquations = 5;

    // Primitive state, equation-major, length totalCells = nCells+nGhostCells.
    // Only the density (component 0) and pressure (component 4) planes are
    // read; the remaining planes document the shared primitive layout.
    const Real * q[ 5 ] = {};

    // Face geometry and cell-face connectivity.  Identical to CellGradientView
    // so the temperature gradient reuses the same resident geometry cache and
    // cell-to-face CSR instead of building a second one.
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

    // Thermodynamic coefficient of the relation above.
    Real stateCoefficient = 0.0;
    Real molecularWeightRatio = 1.0;

    // Outputs, length totalCells each.  Keep the gradient device-resident for
    // the viscous flux unless downloadToHost requests the diagnostic mirror.
    Real * temperature = nullptr;
    Real * dTdx = nullptr;
    Real * dTdy = nullptr;
    Real * dTdz = nullptr;

    const void * cacheKey = nullptr;
    bool downloadToHost = true;
};

inline void ValidateCellTemperatureView( const CellTemperatureView & view )
{
    if ( view.nCells <= 0 || view.nGhostCells < 0
         || view.nGhostCells != view.nBoundaryFaces || view.nFaces <= 0
         || view.nBoundaryFaces < 0 || view.nBoundaryFaces > view.nFaces
         || view.nEquations != 5
         || view.xFace == nullptr || view.yFace == nullptr
         || view.zFace == nullptr || view.xNormal == nullptr
         || view.yNormal == nullptr || view.zNormal == nullptr
         || view.faceArea == nullptr || view.xCell == nullptr
         || view.yCell == nullptr || view.zCell == nullptr
         || view.cellVolume == nullptr || view.leftCell == nullptr
         || view.rightCell == nullptr
         || view.stateCoefficient <= 0.0
         || view.molecularWeightRatio <= 0.0
         || view.temperature == nullptr || view.dTdx == nullptr
         || view.dTdy == nullptr || view.dTdz == nullptr
         || view.cacheKey == nullptr )
    {
        throw std::invalid_argument(
            "invalid solver cell temperature view" );
    }
    // All five primitive planes are required even though the relation only
    // reads density and pressure: a backend that reuses its resident state
    // binds the whole state, and a missing plane there would fail later with a
    // much less specific error.
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        if ( view.q[ equation ] == nullptr )
        {
            throw std::invalid_argument(
                "cell temperature view has a missing primitive component" );
        }
    }
}

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

// ---------------------------------------------------------------------------
// Spalart-Allmaras one-equation turbulence transport.
//
// SA is deliberately NOT a sixth slot in the Euler/NS views.  Two properties
// of the mean-flow contract do not hold for it:
//
//   * several fixed five-slot arrays are indexed by equation up to nEquations,
//     so a sixth slot writes past them;
//   * UpdatePrimitiveStateKernel derives pressure from total energy, which is
//     only meaningful for the conserved mean-flow set.  The SA working
//     variable is a transported scalar with a source term, so treating it as a
//     conserved sixth equation is numerically wrong, not merely awkward.
//
// The host chain is UTurbSrcFlux / UTurbInvFlux / UTurbVisFlux over the
// `turb` solver's own fields, so these views mirror those fields one for one.
// ---------------------------------------------------------------------------

// SaModelConstants is the derived constant set a view and a kernel must agree
// on, so it is declared here rather than in TurbSaArithmetic.h.  Keeping it in
// this shared header is what lets an execution space carry the constants
// without including the arithmetic: AccelViews.h is reachable from every
// backend, so a dependency the other way would drag the SA kernels into CUDA's
// and Kokkos's include graph for a plain POD.
//
// Every value is the derived constant the host carries in `turbcom` after
// TurbCom::InitConst (codes/turb/src/TurbCom.cpp:119-357), not a raw model
// coefficient that this header would re-derive.  That way the device provably
// uses whatever the host initialised, including `cw1` and `cw1k`, which depend
// on karm and sigma.
struct SaModelConstants
{
    Real cb1 = 0.1355;              // cb1
    Real cb2OverSigma = 0.0;        // cb2s = osigma * cb2
    Real cv1Cubed = 0.0;            // cv13 = cv1^3, cv1 = 7.1
    Real cw1 = 0.0;                 // cb1*okarm2 + (1 + cb2)/sigma
    Real cw1K = 0.0;                // cw1k = cw1 * karm2
    Real cw2 = 0.3;
    Real cw3ToSix = 64.0;           // cw36 = cw3^6, cw3 = 2.0
    Real ct3 = 1.2;
    Real ct4 = 0.5;
    Real karm2 = 0.0;               // karm^2, karm = 0.41
    Real okarm2 = 0.0;              // 1/karm2
    Real sac2 = 0.7;
    Real sac3 = 0.9;
    Real fwStar = 0.424;
    Real productionRatio = 1.0;     // rprod
    Real invReynolds = 0.0;         // oreynolds = 1/reynolds
    Real invSigma = 1.5;            // osigma = 1/sigma, sigma = 2/3
    Real maxViscosityRatio = 0.0;   // max_vis_ratio
    bool ft2Enabled = false;        // ft2_flag != 0
    // TurbCom::CalcSaProd branch selector (iprod_sa).  The flagship case sets
    // iprod_sa = 2, which matches none of the branches and leaves work4 at the
    // vorticity magnitude.
    int productionOption = 2;
};

// Fields the SA residual chain reads for one zone.  Every array is cell-major
// with length nCells + nGhostCells, matching the host MRField layout for a
// single-equation field; the velocity-gradient planes are the nine mean-flow
// gradients in host order dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy,
// dwdz.
struct TurbSaStateView
{
    int nCells = 0;
    int nGhostCells = 0;

    // Turbulence solution field and its gradient (uturbf.q / dqdx..dqdz).
    const Real * workingVariable = nullptr;
    const Real * workingVariableGradientX = nullptr;
    const Real * workingVariableGradientY = nullptr;
    const Real * workingVariableGradientZ = nullptr;

    // Mean-flow coupling (uturbf.meanFlow).
    const Real * density = nullptr;
    const Real * velocity[ 3 ] = {};
    // Index convention, and it is NOT the transpose: 0..2 are the x-velocity's
    // three derivatives, 3..5 the y-velocity's, 6..8 the z-velocity's.
    //
    //   [0]=dudx [1]=dudy [2]=dudz
    //   [3]=dvdx [4]=dvdy [5]=dvdz
    //   [6]=dwdx [7]=dwdy [8]=dwdz
    //
    // This matches the host order in UTurbSrcFlux::PrepareCellValue1Equ:341-349.
    // It is easy to get backwards here because the source arrays are
    // `mr[ velocityComponent ][ direction ]`: taking `mr[IU][0], mr[IV][0],
    // mr[IW][0]` for the first three slots yields the transpose, which compiles,
    // passes a view-constructed contract test, and produces wrong residuals
    // without failing anything.
    const Real * velocityGradient[ 9 ] = {};

    // Wall distance (grid CellMesh().dist), one value per owned cell.  It is
    // precomputed at setup (WallDist.cpp) and never updated per step.
    const Real * wallDistance = nullptr;
    // Optional SA length scale (uturbf.len_scale).  When null, wallDistance is
    // used, which is what CalcLengthScaleOfWallDist writes for the plain model.
    const Real * lengthScale = nullptr;

    // Transport coefficients: laminar from Sutherland, eddy an output here.
    const Real * laminarViscosity = nullptr;
    Real * turbulentViscosity = nullptr;

    const Real * cellVolume = nullptr;

    SaModelConstants constants;
    const void * cacheKey = nullptr;
};

inline void ValidateTurbSaStateView( const TurbSaStateView & view )
{
    if ( view.nCells <= 0 || view.nGhostCells < 0
         || view.workingVariable == nullptr
         || view.workingVariableGradientX == nullptr
         || view.workingVariableGradientY == nullptr
         || view.workingVariableGradientZ == nullptr
         || view.density == nullptr || view.wallDistance == nullptr
         || view.laminarViscosity == nullptr
         || view.turbulentViscosity == nullptr
         || view.cellVolume == nullptr || view.density == nullptr
         || view.constants.cv1Cubed <= 0.0 || view.constants.karm2 <= 0.0
         || view.constants.invReynolds <= 0.0 )
    {
        throw std::invalid_argument(
            "invalid Spalart-Allmaras state view" );
    }
    for ( int component = 0; component < 3; ++ component )
    {
        if ( view.velocity[ component ] == nullptr )
        {
            throw std::invalid_argument(
                "SA state view has a missing velocity component" );
        }
    }
    for ( int component = 0; component < 9; ++ component )
    {
        if ( view.velocityGradient[ component ] == nullptr )
        {
            throw std::invalid_argument(
                "SA state view has a missing velocity-gradient component" );
        }
    }
}

// Per-cell source term: the inputs TurbCom::CalcSrcSa reads and the two values
// UTurbSrcFlux::Update1Equ accumulates into the residual and the implicit
// spectrum, already multiplied by the cell volume.
struct TurbSaSourceView
{
    int nCells = 0;
    const TurbSaStateView * state = nullptr;

    // Output, length nCells.  residual receives res_sa * cellVolume and
    // implicitSpectrum receives spec_sa * cellVolume, matching
    // UTurbSrcFlux::Update1Equ:242.
    Real * residual = nullptr;
    Real * implicitSpectrum = nullptr;

    // false accumulates into the existing values, true overwrites them.  The
    // host always accumulates, so the default matches Update1Equ.
    bool overwrite = false;

    const void * cacheKey = nullptr;
};

inline void ValidateTurbSaSourceView( const TurbSaSourceView & view )
{
    if ( view.nCells <= 0 || view.state == nullptr
         || view.residual == nullptr || view.implicitSpectrum == nullptr )
    {
        throw std::invalid_argument( "invalid SA source view" );
    }
    if ( view.state->nCells < view.nCells )
    {
        throw std::invalid_argument(
            "SA source view covers more cells than its state" );
    }
}

// Per-cell eddy viscosity.  Same inputs as TurbSaStateView carries, minus the
// gradient the closure does not read; it is a separate view because the host
// calls CalcVist at least twice per step while the source runs once, and a
// backend may want to keep the two entries and their caching independent.
struct TurbSaEddyViscosityView
{
    int nCells = 0;
    const Real * workingVariable = nullptr;
    const Real * density = nullptr;
    const Real * laminarViscosity = nullptr;
    Real * turbulentViscosity = nullptr;
    SaModelConstants constants;
    const void * cacheKey = nullptr;
};

inline void ValidateTurbSaEddyViscosityView(
    const TurbSaEddyViscosityView & view )
{
    if ( view.nCells <= 0 || view.workingVariable == nullptr
         || view.density == nullptr || view.laminarViscosity == nullptr
         || view.turbulentViscosity == nullptr
         || view.constants.cv1Cubed <= 0.0
         || view.constants.maxViscosityRatio <= 0.0 )
    {
        throw std::invalid_argument(
            "invalid Spalart-Allmaras eddy-viscosity view" );
    }
}

// Green-Gauss cell gradient of the working variable, the turbulence solver's
// counterpart of CellGradientView.
struct TurbSaGradientView
{
    int nCells = 0;
    int nGhostCells = 0;
    int nFaces = 0;
    int nBoundaryFaces = 0;
    const Real * workingVariable = nullptr;
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
    Real * gradientX = nullptr;
    Real * gradientY = nullptr;
    Real * gradientZ = nullptr;
    const void * cacheKey = nullptr;
};

inline void ValidateTurbSaGradientView( const TurbSaGradientView & view )
{
    if ( view.nCells <= 0 || view.nGhostCells < 0 || view.nFaces <= 0
         || view.nBoundaryFaces < 0 || view.nBoundaryFaces > view.nFaces
         || view.workingVariable == nullptr || view.xFace == nullptr
         || view.yFace == nullptr || view.zFace == nullptr
         || view.xNormal == nullptr || view.yNormal == nullptr
         || view.zNormal == nullptr || view.faceArea == nullptr
         || view.xCell == nullptr || view.yCell == nullptr
         || view.zCell == nullptr || view.cellVolume == nullptr
         || view.leftCell == nullptr || view.rightCell == nullptr
         || view.gradientX == nullptr || view.gradientY == nullptr
         || view.gradientZ == nullptr )
    {
        throw std::invalid_argument(
            "invalid Spalart-Allmaras gradient view" );
    }
}

// Explicit update of the turbulence solution field.
struct TurbSaUpdateView
{
    int nCells = 0;
    const Real * oldValue = nullptr;
    const Real * delta = nullptr;
    Real * workingVariable = nullptr;
    Real relaxationCoefficient = 1.0;
    const void * cacheKey = nullptr;
};

inline void ValidateTurbSaUpdateView( const TurbSaUpdateView & view )
{
    if ( view.nCells <= 0 || view.oldValue == nullptr
         || view.delta == nullptr || view.workingVariable == nullptr
         || view.relaxationCoefficient <= 0.0
         || view.relaxationCoefficient > 1.0 )
    {
        throw std::invalid_argument(
            "invalid Spalart-Allmaras update view" );
    }
}

// Per-face scalar flux.  leftCell/rightCell follow the same convention as
// FaceConnectivityView: boundary faces occupy [0, nBoundaryFaces) by default,
// and their right cell is a ghost.
struct TurbSaFaceView
{
    int nFaces = 0;
    int nBoundaryFaces = 0;
    const Real * xNormal = nullptr;
    const Real * yNormal = nullptr;
    const Real * zNormal = nullptr;
    const Real * meshVelocityNormal = nullptr;
    const Real * faceArea = nullptr;
    const int * leftCell = nullptr;
    const int * rightCell = nullptr;
    const unsigned char * boundaryMask = nullptr;

    // Face geometry needed by the scalar diffusion normal gradient.
    const Real * xFace = nullptr;
    const Real * yFace = nullptr;
    const Real * zFace = nullptr;
    const Real * xCell = nullptr;
    const Real * yCell = nullptr;
    const Real * zCell = nullptr;

    // Reconstructed face values: the turbulence limiter output.  These feed the
    // inviscid flux only, because UTurbInvFlux::PrepareFaceValue:126 reads
    // qf1/qf2 and the NS limiter's face densities and velocities.
    const Real * faceWorkingVariableLeft = nullptr;
    const Real * faceWorkingVariableRight = nullptr;
    const Real * faceDensityLeft = nullptr;
    const Real * faceDensityRight = nullptr;
    const Real * faceVelocityLeft[ 3 ] = {};
    const Real * faceVelocityRight[ 3 ] = {};

    // Cell-centred values, cell-major over nCells + nGhostCells.  The diffusive
    // flux reads these at each face's left/right cell, not the reconstructed
    // ones: UTurbVisFlux::PrepareFaceValue copies `(*uturbf.q)[IKE][lc]` into
    // visTurb.q1 (:303-311), and for the one-equation model IKE == ISA.
    // Its laminar viscosity and density are also per cell (`visl[0][lc]`,
    // meanFlow.q[IR][lc]).  Mixing the two sources would change the face
    // diffusivity.
    const Real * cellWorkingVariable = nullptr;
    const Real * cellLaminarViscosity = nullptr;
    const Real * cellDensity = nullptr;

    // Cell gradient of the working variable, cell-major.  Required when a
    // diffusive flux is requested: the VIS_TEST normal gradient reads both
    // cell gradients and both cell centres.
    const Real * workingVariableGradientX = nullptr;
    const Real * workingVariableGradientY = nullptr;
    const Real * workingVariableGradientZ = nullptr;

    // rho_coef: 1 for SA, 0 for the two-equation pair (TurbInv::Init:47).
    Real densityWeight = 1.0;

    // Diffusion constants.  The backend forms the per-side kinematic laminar
    // viscosity as visl * (1/rho) exactly as CalcFaceVisFlux1Equ:339 does.
    Real invSigma = 1.5;
    Real invReynolds = 0.0;
    // iturb_visflux: the host maps 2 to the VIS_TEST normal-gradient variant,
    // which is what the flagship uses.  A backend that cannot reproduce the
    // requested variant must refuse rather than substitute another one.
    int viscousFluxVariant = 2;

    // Optional explicit boundary mask for "this face is an interface".  The
    // host's CalcC1C2 (VisGrad.cpp:206) forces c1 = 1 for a boundary face whose
    // bc type is not INTERFACE; without this mask every boundary face is
    // treated as non-interface, which is correct for solid/farfield/inflow but
    // not for a multizone interface face.
    const unsigned char * interfaceMask = nullptr;

    // Output, length nFaces.  Kept separate so a stage can add the inviscid and
    // diffusive contributions independently, as UTurbInvFlux and UTurbVisFlux
    // do on the host.
    Real * inviscidFlux = nullptr;
    Real * diffusiveFlux = nullptr;

    const void * cacheKey = nullptr;
};

inline void ValidateTurbSaFaceView( const TurbSaFaceView & view )
{
    if ( view.nFaces <= 0 || view.nBoundaryFaces < 0
         || view.nBoundaryFaces > view.nFaces
         || view.xNormal == nullptr || view.yNormal == nullptr
         || view.zNormal == nullptr || view.meshVelocityNormal == nullptr
         || view.faceArea == nullptr || view.leftCell == nullptr
         || view.rightCell == nullptr )
    {
        throw std::invalid_argument( "invalid SA face view" );
    }
    if ( view.inviscidFlux == nullptr && view.diffusiveFlux == nullptr )
    {
        throw std::invalid_argument( "SA face view has no output buffer" );
    }
    if ( view.inviscidFlux != nullptr )
    {
        if ( view.faceWorkingVariableLeft == nullptr
             || view.faceWorkingVariableRight == nullptr
             || view.faceDensityLeft == nullptr
             || view.faceDensityRight == nullptr )
        {
            throw std::invalid_argument(
                "SA inviscid flux requires reconstructed face values and face "
                "densities" );
        }
        for ( int component = 0; component < 3; ++ component )
        {
            if ( view.faceVelocityLeft[ component ] == nullptr
                 || view.faceVelocityRight[ component ] == nullptr )
            {
                throw std::invalid_argument(
                    "SA inviscid flux has a missing face velocity component" );
            }
        }
    }
    if ( view.diffusiveFlux != nullptr )
    {
        if ( view.xFace == nullptr || view.yFace == nullptr
             || view.zFace == nullptr || view.xCell == nullptr
             || view.yCell == nullptr || view.zCell == nullptr
             || view.cellWorkingVariable == nullptr
             || view.cellLaminarViscosity == nullptr
             || view.cellDensity == nullptr
             || view.workingVariableGradientX == nullptr
             || view.workingVariableGradientY == nullptr
             || view.workingVariableGradientZ == nullptr )
        {
            throw std::invalid_argument(
                "SA diffusive flux requires face geometry, the cell working "
                "variable, its cell gradient, the cell laminar viscosity and "
                "the cell density" );
        }
        if ( view.invReynolds <= 0.0 )
        {
            throw std::invalid_argument(
                "SA diffusive flux requires a positive inverse Reynolds number" );
        }
        if ( view.viscousFluxVariant != 2 )
        {
            throw std::invalid_argument(
                "only the VIS_TEST (iturb_visflux = 2) normal-gradient variant "
                "is implemented; refusing to substitute another one" );
        }
    }
}

EndNameSpace
