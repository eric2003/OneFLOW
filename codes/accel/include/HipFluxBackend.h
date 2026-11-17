/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\---------------------------------------------------------------------------*/

#pragma once

#include "FluxBackend.h"
#include "BoundaryExchangeContract.h"
#include "DeviceBuffer.h"
#include "EulerDomain.h"
#include "ResidentStateContract.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <unordered_map>
#include <vector>

namespace ONEFLOW
{

struct PrimitiveFaceStateView;

// Per-zone residual diagnostics computed from the state-owned device residual.
// The solver consumes the same quantities the host residual task used to
// derive by looping the MRField: sum of squares, max |residual| with its cell
// index, and the primitive positivity extremes.  Only these values cross the
// PCIe boundary; the residual field itself stays on the device.
struct HipResidualMetrics
{
    std::vector< Real > sumOfSquares;
    std::vector< Real > maxAbsolute;
    std::vector< int > argmaxCell;
    Real minDensity = 0.0;
    Real minPressure = 0.0;
    int nonPhysicalCells = 0;
    bool valid = false;
};

struct HipBoundaryEndpointKey
{
    int localZone = -1;
    int neighborZone = -1;
    int ghostLevel = 0;

    bool operator==( const HipBoundaryEndpointKey & other ) const noexcept
    {
        return localZone == other.localZone
            && neighborZone == other.neighborZone
            && ghostLevel == other.ghostLevel;
    }
};

struct HipBoundaryEndpointKeyHash
{
    std::size_t operator()(
        const HipBoundaryEndpointKey & key ) const noexcept;
};

struct HipBoundaryEndpointStorage
{
    DeviceBuffer< int > cellIds;
    DeviceBuffer< Real > values;
    std::vector< int > cachedCellIds;
    std::vector< Real > hostValues;
    std::uint64_t topologyGeneration = 0;
    int cellCount = 0;
    int nEquations = 0;
    BoundaryEndpointRole role = BoundaryEndpointRole::Unspecified;
};

// HIP-only storage owned either by a solver/grid Ns3DBackendState or by the
// compatibility fallback in HipFluxBackend. Keeping this type in the HIP
// header avoids leaking DeviceBuffer through backend-neutral contracts.
struct HipGradientStorage
{
    DeviceBuffer< Real > cellState;
    DeviceBuffer< Real > cellGradient;
    DeviceBuffer< Real > xFaceCenter;
    DeviceBuffer< Real > yFaceCenter;
    DeviceBuffer< Real > zFaceCenter;
    DeviceBuffer< Real > xNormal;
    DeviceBuffer< Real > yNormal;
    DeviceBuffer< Real > zNormal;
    DeviceBuffer< Real > meshVelocityNormal;
    DeviceBuffer< Real > faceArea;
    DeviceBuffer< Real > xCellCenter;
    DeviceBuffer< Real > yCellCenter;
    DeviceBuffer< Real > zCellCenter;
    DeviceBuffer< Real > cellVolume;
    DeviceBuffer< int > leftCell;
    DeviceBuffer< int > rightCell;
    // Cell-to-face CSR built once per topology.  The cell-based gradient and
    // residual kernels consume this list; it removes the face-parallel
    // atomicAdd scatter from the production path.
    DeviceBuffer< int > cellFaceOffsets;
    DeviceBuffer< int > cellFaceIds;
    const int * cachedCellFaceLeft = nullptr;
    const int * cachedCellFaceRight = nullptr;
    int cachedCellFaceCells = -1;
    int cachedCellFaceFaces = -1;
    int cachedCellFaceBoundaryFaces = -1;
    bool cellFaceValid = false;
    const void * cachedOwnerToken = nullptr;
    std::uint64_t cachedTopologyGeneration = 0;
    const Real * cachedXFace = nullptr;
    const Real * cachedYFace = nullptr;
    const Real * cachedZFace = nullptr;
    const Real * cachedXNormal = nullptr;
    const Real * cachedYNormal = nullptr;
    const Real * cachedZNormal = nullptr;
    const Real * cachedFaceArea = nullptr;
    const Real * cachedXCell = nullptr;
    const Real * cachedYCell = nullptr;
    const Real * cachedZCell = nullptr;
    const Real * cachedVolume = nullptr;
    const int * cachedLeftCell = nullptr;
    const int * cachedRightCell = nullptr;
    const Real * cachedMeshVelocityNormal = nullptr;
    int cachedCells = 0;
    int cachedGhostCells = -1;
    int cachedFaces = 0;
    int cachedBoundaryFaces = -1;
    bool geometryCacheValid = false;
    DeviceBuffer< Real > qf1;
    DeviceBuffer< Real > qf2;
    // Face flux is part of the state-owned reconstructed production chain.
    // The backend-level deviceFlux_ remains a compatibility buffer for the
    // legacy host-staged path only.
    DeviceBuffer< Real > faceFlux;
    DeviceBuffer< Real > residual;
    DeviceBuffer< Real > inviscidSpectrum;
    DeviceBuffer< Real > timeStep;
    DeviceBuffer< Real > minTimeStep;
    DeviceBuffer< int > timeStepStatus;
    DeviceBuffer< int > stateUpdateStatus;
    DeviceBuffer< unsigned char > boundaryMask;
    DeviceBuffer< ReconstructionBoundaryOperation > boundaryOperation;
    DeviceBuffer< Real > bcQ;
    DeviceBuffer< Real > boundaryFaceState;
    // Packed [equation][ghost cell] slice used for a single D2H per stage.
    DeviceBuffer< Real > boundaryGhostPack;
    // Packed [equation][owned cell] slice for on-demand host readers.
    DeviceBuffer< Real > stateMirrorPack;
    DeviceBuffer< int > stateMirrorIds;
    DeviceBuffer< int > boundaryType;
    DeviceBuffer< Real > inflow;
    DeviceBuffer< Real > boundaryXNormal;
    DeviceBuffer< Real > boundaryYNormal;
    DeviceBuffer< Real > boundaryZNormal;
    DeviceBuffer< Real > boundaryMeshVelocityNormal;
    DeviceBuffer< int > boundaryLeftCell;
    DeviceBuffer< int > boundaryRightCell;
    // Boundary geometry/topology is immutable for a bound grid across RK
    // stages. Keep its host-view identity so UpdateBoundaryState does not
    // re-upload the same metadata on every stage. Dynamic inflow values are
    // intentionally not part of this cache.
    const void * cachedBoundaryOwnerToken = nullptr;
    std::uint64_t cachedBoundaryGeometryGeneration = 0;
    const Real * cachedBoundaryXNormal = nullptr;
    const Real * cachedBoundaryYNormal = nullptr;
    const Real * cachedBoundaryZNormal = nullptr;
    const Real * cachedBoundaryMeshVelocityNormal = nullptr;
    const int * cachedBoundaryLeftCell = nullptr;
    const int * cachedBoundaryRightCell = nullptr;
    const int * cachedBoundaryType = nullptr;
    bool boundaryGeometryCacheValid = false;
    std::unordered_map<
        HipBoundaryEndpointKey,
        HipBoundaryEndpointStorage,
        HipBoundaryEndpointKeyHash > boundarySend;
    std::unordered_map<
        HipBoundaryEndpointKey,
        HipBoundaryEndpointStorage,
        HipBoundaryEndpointKeyHash > boundaryReceive;
    const ReconstructionBoundaryOperation * cachedBoundaryOperationHost =
        nullptr;
    std::size_t cachedBoundaryOperationCount = 0;
    std::uint64_t cachedBoundaryTopologyGeneration = 0;
    int cachedEquations = 0;
    int residualCells = 0;
    int residualEquations = 0;
    const Real * cachedTimeStepXNormal = nullptr;
    const Real * cachedTimeStepYNormal = nullptr;
    const Real * cachedTimeStepZNormal = nullptr;
    const Real * cachedTimeStepMeshVelocityNormal = nullptr;
    const Real * cachedTimeStepFaceArea = nullptr;
    const Real * cachedTimeStepVolume = nullptr;
    const int * cachedTimeStepLeftCell = nullptr;
    const int * cachedTimeStepRightCell = nullptr;
    const Real * cachedTimeStepHost = nullptr;
    int cachedTimeStepCells = 0;
    int cachedTimeStepGhostCells = -1;
    int cachedTimeStepFaces = 0;
    int cachedTimeStepBoundaryFaces = -1;
    bool timeStepGeometryCacheValid = false;
    bool residentTimeStepReady = false;
    // Generation of the reconstruction metadata/field view.
    std::uint64_t cachedFieldGeneration = 0;
    // Single source of truth for "is the resident q current, and for which
    // owner/generation".  Ghost-only host mutations are tracked separately so
    // an interface round trip does not force a full field upload.
    ResidentStateContract residentContract;
    bool boundaryMaskValid = false;
    bool residualValid = false;
    // True between CALC_LHS (device residual scaling) and the next flux write.
    // The state update consumes residualValid, but the diagnostics still need
    // the scaled residual, so the two flags are tracked separately.
    bool residualMetricsValid = false;
    const Real * cachedStateQ[ 5 ] = {};
    // Cells whose gradient is mirrored to the host for the host-side interface
    // exchange.  Kept resident so the mirror does not re-upload ids per stage.
    DeviceBuffer< int > gradientMirrorCellIds;
    DeviceBuffer< int > gradientPackIds;
    DeviceBuffer< Real > gradientPackValues;
    std::vector< int > cachedGradientMirrorCellIds;
    DeviceBuffer< Real > gradientMirrorValues;
    int cachedGradientMirrorEquations = 0;
    // Residual diagnostic reduction buffers (small, reused every outer step).
    DeviceBuffer< Real > residualMetricChunkSum;
    DeviceBuffer< Real > residualMetricChunkMax;
    DeviceBuffer< int > residualMetricChunkIndex;
    DeviceBuffer< Real > residualMetricChunkMin;
    int residualMetricChunks = 0;
    int residualMetricEquations = 0;
    // Cache metadata for the bound-state host-staged primitive path. The
    // unbound direct/smoke path keeps using HipFluxBackend's fallback cache.
    const void * cachedFluxOwnerToken = nullptr;
    const Real * cachedFluxXNormal = nullptr;
    const Real * cachedFluxYNormal = nullptr;
    const Real * cachedFluxZNormal = nullptr;
    const Real * cachedFluxMeshVelocityNormal = nullptr;
    const Real * cachedFluxFaceArea = nullptr;
    const int * cachedFluxLeftCell = nullptr;
    const int * cachedFluxRightCell = nullptr;
    const unsigned char * cachedFluxBoundaryMask = nullptr;
    int cachedFluxFaces = 0;
    int cachedFluxBoundaryFaces = -1;
    bool fluxGeometryCacheValid = false;
};

std::unique_ptr< Ns3DBackendState > CreateHipNs3DBackendState(
    int deviceId );

class HipFluxBackend final : public FluxBackend
{
public:
    static HipFluxBackend & Shared();
    static void ReleaseShared();

    void CalcInvFlux(
        const FaceStateView & state,
        FaceFluxView & flux,
        int scheme ) override;

    void AddFaceFlux(
        const FaceFluxView & flux,
        const FaceConnectivityView & connectivity,
        ResidualView & residual ) override;

    // Add the flux produced by the immediately preceding CalcInvFlux call
    // without uploading the host flux array again. This is the production
    // seam used by the 3D main solver HIP batch path.
    void AddCurrentFaceFlux(
        const FaceConnectivityView & connectivity,
        ResidualView & residual );

    // Fused production path for primitive face data. The kernel converts
    // primitive states, evaluates the numerical flux, and scatters directly
    // into the device residual. hostFlux may be null in NoTrace mode.
    void CalcAndAddPrimitiveFaceFlux(
        const PrimitiveFaceStateView & state,
        const FaceConnectivityView & connectivity,
        ResidualView & residual,
        int scheme,
        FaceFluxView * hostFlux = nullptr,
        const void * cacheKey = nullptr );

    void CalcGradient( const CellGradientView & view );

    // Compute the local inviscid spectral radius and timestep from resident
    // primitive state. The resulting timestep can remain device-resident for
    // ScaleCurrentResidual.
    void CalcLocalInviscidTimeStep( const CellTimeStepView & view );

    // Update supported local inviscid boundary models on the resident state.
    void UpdateBoundaryState( const BoundaryStateView & view );

    // Reconstruct face primitive values into state-owned device buffers. The
    // host view supplies metadata and optional diagnostic destinations; the
    // numerical path consumes q/gradient/geometry already resident on device.
    void ReconstructFaceValues(
        const CellFaceReconstructionView & view,
        bool downloadToHost = false );

    // Consume the device qf1/qf2 produced by ReconstructFaceValues and scatter
    // the resulting flux directly into the device residual. hostFlux is only
    // populated for explicit trace/diagnostic requests.
    void CalcAndAddReconstructedFaceFlux(
        const FaceConnectivityView & connectivity,
        ResidualView & residual,
        int scheme,
        Real gamma,
        FaceFluxView * hostFlux = nullptr,
        const void * cacheKey = nullptr,
        bool downloadResidual = true );

    // Continue the reconstructed-flux path without a residual round trip.
    // ScaleCurrentResidual implements CALC_LHS. UpdatePrimitiveState implements
    // the five-equation primitive -> conserved -> primitive flow update and
    // mirrors the internal-cell primitive state to host when requested.
    void ScaleCurrentResidual( const CellStateUpdateView & view );
    void UpdatePrimitiveState(
        const CellStateUpdateView & view,
        bool downloadToHost = true );

    BoundaryBufferView PackBoundary(
        const void * cacheKey,
        const BoundaryFieldExchangePlan & plan,
        BoundaryTransferMode mode );
    BoundaryBufferView PackBoundary(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode );
    BoundaryBufferView PackBoundary(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode,
        const HipBoundaryEndpointKey & endpointKey );
    void UploadBoundaryState(
        const void * cacheKey,
        int nCells,
        int nGhostCells,
        int nEquations,
        const Real * const * components );
    BoundaryBufferView PrepareBoundaryReceive(
        const void * cacheKey,
        const BoundaryFieldExchangePlan & plan,
        BoundaryTransferMode mode );
    BoundaryBufferView PrepareBoundaryReceive(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode );
    BoundaryBufferView PrepareBoundaryReceive(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode,
        const HipBoundaryEndpointKey & endpointKey );
    void UnpackBoundary(
        const void * cacheKey,
        const BoundaryFieldExchangePlan & plan,
        BoundaryTransferMode mode );
    void UnpackBoundary(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode );
    void UnpackBoundary(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode,
        const HipBoundaryEndpointKey & endpointKey );

    void BindState( const void * cacheKey, Ns3DDeviceState & state );
    void UnbindState( const void * cacheKey ) noexcept;

    // --- Resident state ownership contract -------------------------------
    // The solver calls these when it writes the host mirror outside the
    // device path.  Owned-cell writes invalidate the whole resident field;
    // ghost writes only invalidate the ghost region.
    void NotifyHostFieldStateModified( const void * cacheKey );
    void NotifyHostGhostStateModified( const void * cacheKey );

    // Refresh the ghost region of the state-owned device q from the host
    // mirror when (and only when) the contract says it is stale.  This is the
    // host-staged halo binding: the exchanged host ghost remains the transfer
    // source of truth, while the device state owns the resident copy.
    bool SyncResidentGhostStateFromHost(
        const void * cacheKey,
        int nCells,
        int nGhostCells,
        int nEquations,
        const Real * const * components );

    // Mirror the resident gradient at the given cells into the solver-owned
    // host components.  A device reconstruction only needs the host mirror at
    // the boundary-adjacent owned cells (the interface send values and the
    // face-indexed boundary gradient), not the whole field: the interface
    // exchange, including its one-round lag, still runs through
    // Grad::SwapBcGrad on that mirror.
    void MirrorGradientCellsToHost(
        const void * cacheKey,
        int nEquations,
        const std::vector< int > & cellIds,
        Real * const * dqdx,
        Real * const * dqdy,
        Real * const * dqdz );

    // Pack the resident gradient at the interface send cells into a
    // [component][equation][face] host buffer.  The host then copies it into
    // the interface send storage in exchange order, so the endpoint cell
    // gather runs on the device instead of the host.
    void PackGradientFacesToHost(
        const void * cacheKey,
        const int * cellIds,
        int nFaces,
        int nEquations,
        Real * hostPacked );

    // Read the interface receive storage into a packed host buffer and scatter
    // it directly into the resident ghost gradient.  This replaces the host
    // ghost scatter plus the separate ghost-gradient H2D sync.
    void UnpackGradientFacesFromHost(
        const void * cacheKey,
        const int * ghostCellIds,
        int nFaces,
        int nEquations,
        const Real * hostPacked );

    // Upload the interface receive payload for the device gradient into the
    // resident pack buffer, once per exchange round.  Every RK stage then
    // unpacks it on the device instead of repeating the same host->device copy.
    void UploadGradientFacesToDevice(
        const void * cacheKey,
        int nFaces,
        int nEquations,
        const Real * hostPacked );

    // Scatter the uploaded receive payload into the resident ghost gradient.
    // The payload buffer is the one UploadGradientFacesToDevice filled.
    void UnpackGradientFacesFromDevice(
        const void * cacheKey,
        const int * ghostCellIds,
        int nFaces,
        int nEquations );

    // Reduce the resident (already LHS-scaled) device residual into the
    // per-zone diagnostic scalars the solver's residual task reports.  Returns
    // false when the grid has no resident residual, so callers keep the host
    // loop as the fallback.
    bool ComputeResidualMetrics(
        const void * cacheKey,
        int nCells,
        int nEquations,
        Real gamma,
        HipResidualMetrics & metrics );

    // Diagnostic contract check: compare the complete resident primitive
    // state (owned + ghost cells) against the host mirror.  Returns false when
    // the grid has no bound resident state.
    bool VerifyResidentStateAgainstHost(
        const void * cacheKey,
        int nCells,
        int nGhostCells,
        int nEquations,
        const Real * const * components,
        Real tolerance,
        Real * maxAbsoluteDifference,
        long * mismatchedValues );

    // True when the resident device state is newer than the host mirror, i.e.
    // the device is authoritative for q and host readers must pull the values
    // back before using them.
    bool HostMirrorStale( const void * cacheKey ) const;

    // Pull the complete resident field (owned + ghost cells) back into the
    // host q mirror.  No-op when the mirror is already current.  Returns true
    // when values were transferred.
    bool RefreshHostPrimitiveState(
        const void * cacheKey,
        int nCells,
        int nGhostCells,
        int nEquations,
        Real * const * components );

    // Pull the resident field at the given owned cells back into the host
    // mirror.  Used by the per-step host readers that only touch
    // boundary-adjacent cells (aerodynamic force, pressure coefficient).
    void MirrorPrimitiveCellsToHost(
        const void * cacheKey,
        int nEquations,
        const int * cellIds,
        int cellCount,
        Real * const * components );

    // The device wrote cells the host mirror does not reflect (an interface
    // unpack whose paired host exchange was skipped).
    void MarkResidentHostMirrorStale( const void * cacheKey );

private:
    DeviceBuffer< Real > qLeft_;
    DeviceBuffer< Real > qRight_;
    DeviceBuffer< Real > xNormal_;
    DeviceBuffer< Real > yNormal_;
    DeviceBuffer< Real > zNormal_;
    DeviceBuffer< Real > meshVelocityNormal_;
    DeviceBuffer< Real > area_;
    DeviceBuffer< Real > deviceFlux_;
    DeviceBuffer< int > leftCell_;
    DeviceBuffer< int > rightCell_;
    DeviceBuffer< unsigned char > boundaryMask_;
    DeviceBuffer< Real > deviceResidual_;
    HipGradientStorage fallbackGradientStorage_;
    std::unordered_map< const void *, Ns3DDeviceState * > boundStates_;
    int deviceFluxFaces_ = 0;
    int deviceFluxEquations_ = 0;
    const Real * cachedXNormal_ = nullptr;
    const Real * cachedYNormal_ = nullptr;
    const Real * cachedZNormal_ = nullptr;
    const Real * cachedMeshVelocityNormal_ = nullptr;
    const Real * cachedArea_ = nullptr;
    const int * cachedLeftCell_ = nullptr;
    const int * cachedRightCell_ = nullptr;
    const unsigned char * cachedBoundaryMask_ = nullptr;
    const void * cachedOwner_ = nullptr;
    int cachedFaces_ = 0;
    int cachedBoundaryFaces_ = -1;
    bool geometryCacheValid_ = false;
};

}
